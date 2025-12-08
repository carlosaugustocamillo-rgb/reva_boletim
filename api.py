
from fastapi import FastAPI, BackgroundTasks, Request, Response
from fastapi.responses import HTMLResponse, StreamingResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import json
import uuid

from boletim_service import rodar_boletim
from campanha_service import criar_campanha_tematica
from revamais_service import criar_campanha_revamais
from simple_agent import run_agent, AgentInput
from pydantic import BaseModel
import firebase_admin
from firebase_admin import credentials, storage

class CampanhaInput(BaseModel):
    tema: str = ""

# Force rebuild for Python 3.11
app = FastAPI(title="RevaCast Boletim Service")

# Configuração de CORS PRIMEIRO (será "envelopada" pelo middleware manual abaixo)
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://localhost:8000",
        "https://revalidatie.com.br",
        "https://www.revalidatie.com.br",
        "https://reva-boletim-production.up.railway.app",
        "https://stackblitz.com",
        "https://cors-proxy-production.stackblitz.workers.dev",
    ],
    # allow_origin_regex removed to avoid conflicts since manual middleware handles it
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
 
@app.middleware("http")
async def cors_override_middleware(request: Request, call_next):
    origin = request.headers.get("origin")
    
    # Lista de domínios permitidos dinamicamente
    allowed_domains = ["webcontainer-api.io", "stackblitz.com", "stackblitz.io", "stackblitz.workers.dev", "localhost"]
    
    # Verifica se a origem bate com algum dos domínios permitidos
    is_allowed = False
    if origin:
        for domain in allowed_domains:
            if domain in origin:
                is_allowed = True
                break
    else:
        # Se não tem origin (server-to-server ou local), as vezes é bom liberar ou tratar como seguro dependendo do caso.
        # Mas para browser, geralmente tem origin. Se for null (redirects), liberamos.
        if origin is None or origin == "null":
             is_allowed = True
             origin = "*" # Fallback

    # Se for requisição OPTIONS vinda de origem permitida, responde direto
    if request.method == "OPTIONS" and is_allowed:
        response = Response(status_code=200)
        response.headers["Access-Control-Allow-Origin"] = origin
        response.headers["Access-Control-Allow-Credentials"] = "true"
        response.headers["Access-Control-Allow-Methods"] = "*"
        response.headers["Access-Control-Allow-Headers"] = "*"
        return response

    response = await call_next(request)
    
    # Adiciona headers nas rotas normais também
    if is_allowed:
        response.headers["Access-Control-Allow-Origin"] = origin
        response.headers["Access-Control-Allow-Credentials"] = "true"
        response.headers["Access-Control-Allow-Methods"] = "*"
        response.headers["Access-Control-Allow-Headers"] = "*"
        
    return response

@app.get("/ping")
def ping():
    return {"status": "ok", "version": "with_cors_fix_v2", "cors": "enabled"}

# --- Armazenamento em ARQUIVO para tarefas em background (Persistência) ---
import os
import json

TASKS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "tasks")
os.makedirs(TASKS_DIR, exist_ok=True)

def save_task(task_id, data):
    """Salva o estado da tarefa em um arquivo JSON."""
    try:
        filepath = os.path.join(TASKS_DIR, f"{task_id}.json")
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"Erro ao salvar tarefa {task_id}: {e}")

def load_task(task_id):
    """Carrega o estado da tarefa do arquivo JSON."""
    try:
        filepath = os.path.join(TASKS_DIR, f"{task_id}.json")
        if not os.path.exists(filepath):
            return None
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        print(f"Erro ao carregar tarefa {task_id}: {e}")
        return None

def processar_boletim_background(task_id: str, opcoes: dict):
    """Função wrapper que roda o boletim e salva logs em arquivo."""
    # Estado inicial
    task_state = {"status": "running", "logs": [], "result": None}
    save_task(task_id, task_state)
    
    try:
        for log_msg in rodar_boletim(opcoes):
            # Verifica se foi solicitado cancelamento
            current_state = load_task(task_id)
            if current_state and current_state.get("status") == "canceling":
                task_state["status"] = "canceled"
                task_state["logs"].append("🛑 Processo cancelado pelo usuário.")
                save_task(task_id, task_state)
                return # Interrompe imediatamente

            if isinstance(log_msg, dict):
                task_state["result"] = log_msg
            else:
                task_state["logs"].append(str(log_msg))
            
            # Salva periodicamente
            save_task(task_id, task_state)
        
        task_state["status"] = "completed"
        task_state["logs"].append("✅ Processo finalizado com sucesso.")
        save_task(task_id, task_state)
        
        task_state["status"] = "error"
        task_state["logs"].append(error_msg)
        save_task(task_id, task_state)
        print(error_msg)

def processar_revamais_background(task_id: str, tema: str):
    """Função wrapper para Reva+ em background."""
    task_state = {"status": "running", "logs": [], "result": None}
    save_task(task_id, task_state)
    
    def log_callback(msg):
        task_state["logs"].append(msg)
        save_task(task_id, task_state)
        
    def check_cancel():
        current = load_task(task_id)
        if current and current.get("status") == "canceling":
            task_state["status"] = "canceled"
            task_state["logs"].append("🛑 Cancelado pelo usuário.")
            save_task(task_id, task_state)
            return True
        return False
        
    try:
        from revamais_service import criar_campanha_revamais
        resultado = criar_campanha_revamais(tema, log_callback=log_callback, check_cancel=check_cancel)
        
        # Se retornou sucesso/erro (dict)
        if isinstance(resultado, dict):
            if resultado.get("status") == "error":
                task_state["status"] = "error"
            else:
                task_state["status"] = "completed"
            task_state["result"] = resultado
            task_state["logs"].append("✅ Processo Reva+ finalizado.")
        else:
            task_state["status"] = "completed"
            
        save_task(task_id, task_state)
        
    except Exception as e:
        if str(e) == "CANCELADO_PELO_USUARIO":
             return # Já salvou status canceled no check_cancel
             
        import traceback
        error_msg = f"❌ Erro fatal: {str(e)}\n{traceback.format_exc()}"
        task_state["status"] = "error"
        task_state["logs"].append(error_msg)
        save_task(task_id, task_state)


@app.post("/cancelar-boletim/{task_id}")
def cancelar_boletim(task_id: str):
    task = load_task(task_id)
    if not task:
        return {"status": "not_found", "message": "Tarefa não encontrada."}
    
    if task["status"] in ["running", "queued"]:
        task["status"] = "canceling"
        task["logs"].append("⚠️ Solicitando cancelamento...")
        save_task(task_id, task)
        return {"status": "canceling", "message": "Solicitação de cancelamento enviada."}
    
    return {"status": task["status"], "message": "Tarefa já finalizada ou cancelada."}

@app.post("/atualizar-feed")
def trigger_feed_update():
    """Força a atualização do Feed RSS com base nos arquivos do Firebase."""
    try:
        from force_feed_update import force_update_feed
        force_update_feed()
        return {"status": "success", "message": "Feed RSS atualizado com sucesso!"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

@app.get("/listar-episodios")
def listar_episodios():
    """Lista os últimos 5 episódios do Firebase com link de download."""
    try:
        # Garante inicialização (caso não tenha rodado via boletim_service)
        if not firebase_admin._apps:
            cred_path = os.environ.get("FIREBASE_CREDENTIALS_JSON", "firebase_credentials.json")
            if os.path.exists(cred_path):
                cred = credentials.Certificate(cred_path)
                firebase_admin.initialize_app(cred, {
                    'storageBucket': os.environ.get("FIREBASE_BUCKET_NAME")
                })
        
        bucket = storage.bucket()
        # Debug: Listar tudo para ver se o prefixo está certo
        # blobs = list(bucket.list_blobs(prefix="episodios/"))
        # print(f"DEBUG: Encontrados {len(blobs)} blobs em episodios/")
        
        blobs = bucket.list_blobs(prefix="episodios/")
        lista = []
        for blob in blobs:
            if blob.name.endswith(".mp3"):
                # Garante que é público
                try:
                    blob.make_public()
                except: pass
                
                lista.append({
                    "nome": blob.name.split("/")[-1],
                    "url": blob.public_url,
                    "data": blob.time_created.isoformat(),
                    "tamanho": f"{blob.size / 1024 / 1024:.2f} MB"
                })
        
        # Ordena por data (mais recente primeiro)
        lista.sort(key=lambda x: x['data'], reverse=True)
        return lista[:5] # Retorna só os 5 últimos
    except Exception as e:
        return {"error": str(e)}

@app.get("/")
def root():
    return {"message": "Serviço do boletim científico está no ar. Acesse /painel para interface gráfica."}

@app.get("/health")
def health_check():
    return {"status": "ok", "python_version": "unknown"}

@app.get("/painel", response_class=HTMLResponse)
def painel():
    with open("dashboard.html", "r", encoding="utf-8") as f:
        return f.read()

# --- Novos Endpoints Background ---

@app.post("/iniciar-boletim")
def iniciar_boletim(
    background_tasks: BackgroundTasks,
    resumos: bool = True,
    roteiro: bool = True,
    audio: bool = True,
    mailchimp: bool = True,
    firebase: bool = True
):
    task_id = str(uuid.uuid4())
    opcoes = {
        'resumos': resumos,
        'roteiro': roteiro,
        'audio': audio,
        'mailchimp': mailchimp,
        'firebase': firebase
    }
    
    # Cria o arquivo inicial
    save_task(task_id, {"status": "queued", "logs": ["⏳ Iniciando..."], "result": None})
    
    # Inicia a tarefa em background
    background_tasks.add_task(processar_boletim_background, task_id, opcoes)
    
    return {
        "task_id": task_id,
        "status": "started",
        "message": "Boletim iniciado em segundo plano. Verifique o status com o ID fornecido."
    }

@app.post("/iniciar-revamais")
def iniciar_revamais(
    background_tasks: BackgroundTasks,
    input_data: CampanhaInput
):
    task_id = str(uuid.uuid4())
    
    # Cria estado inicial
    save_task(task_id, {"status": "queued", "logs": ["⏳ Iniciando Reva+ ..."], "result": None})
    
    background_tasks.add_task(processar_revamais_background, task_id, input_data.tema)
    
    return {
        "task_id": task_id,
        "status": "started",
        "message": "Reva+ iniciado em background."
    }
    
@app.post("/cancelar-tarefa/{task_id}")
def cancelar_tarefa_generica(task_id: str):
    return cancelar_boletim(task_id) # Reutiliza a mesma lógica

@app.get("/status-boletim/{task_id}")
def get_status_boletim(task_id: str):
    task = load_task(task_id)
    if not task:
        return {"status": "not_found", "message": "Tarefa não encontrada."}
    return task

@app.post("/criar-revamais")
def criar_revamais_endpoint(input_data: CampanhaInput):
    """
    Cria uma edição do Reva + (newsletter temática).
    """
    print(f"🔎 DEBUG: Request recebido no endpoint /criar-revamais. Tema: {input_data.tema}")
    try:
        from revamais_service import criar_campanha_revamais
        resultado = criar_campanha_revamais(input_data.tema)
        
        # Força Header CORS manual (Cinto e Suspensórios)
        response = JSONResponse(content={"success": True, "resultado": resultado})
        return response
    except Exception as e:
        print(f"❌ DEBUG: Erro no endpoint: {e}")
        response = JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )
        return response

# --- Endpoint Antigo (Mantido para compatibilidade) ---
@app.get("/rodar-boletim-stream")
def rodar_boletim_stream(
    resumos: bool = True,
    roteiro: bool = True,
    audio: bool = True,
    mailchimp: bool = True,
    firebase: bool = True
):
    """
    Endpoint de streaming que envia atualizações de progresso via Server-Sent Events (SSE).
    Aceita query params para controlar quais etapas executar.
    """
    opcoes = {
        'resumos': resumos,
        'roteiro': roteiro,
        'audio': audio,
        'mailchimp': mailchimp,
        'firebase': firebase
    }

    def iter_boletim():
        try:
            for passo in rodar_boletim(opcoes):
                if isinstance(passo, dict):
                    # Resultado final
                    yield f"data: {json.dumps({'status': 'done', 'result': passo})}\n\n"
                else:
                    # Mensagem de progresso
                    yield f"data: {json.dumps({'status': 'progress', 'message': passo})}\n\n"
        except Exception as e:
            yield f"data: {json.dumps({'status': 'error', 'message': str(e)})}\n\n"

    return StreamingResponse(iter_boletim(), media_type="text/event-stream")


@app.post("/rodar-boletim")
def rodar_boletim_endpoint():
    """
    Endpoint legado (mantido para compatibilidade, mas agora consome o gerador até o fim).
    """
    try:
        # Consome o gerador até o último item (que é o resultado)
        resultado = None
        for item in rodar_boletim():
            resultado = item
        
        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )


@app.post("/criar-campanha")
def criar_campanha_endpoint(input_data: CampanhaInput):
    """
    Cria uma campanha temática completa: Deep Research -> Imagem -> Texto -> Mailchimp.
    """
    try:
        resultado = criar_campanha_tematica(input_data.tema)
        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )


@app.post("/criar-revamais")
def criar_revamais_endpoint(input_data: CampanhaInput):
    """
    Cria uma edição do boletim 'Reva +' (foco em pacientes).
    """
    try:
        resultado = criar_campanha_revamais(input_data.tema)
        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )


@app.post("/agente-revacast")
async def agente_revacast(input_data: AgentInput):
    """
    Endpoint novo, que conversa com o AGENTE Customizado (sem Agent Kit).
    
    - Recebe a mensagem do usuário e histórico.
    - O agente decide se chama a tool rodar_boletim.
    - Retorna a resposta textual e o histórico atualizado.
    """
    try:
        result = run_agent(input_data)
        return result
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"error": str(e)}
        )


@app.get("/debug-import")
def debug_import():
    try:
        import boletim_service
        return {"status": "ok", "message": "Importação bem sucedida!"}
    except Exception as e:
        import traceback
        return {
            "status": "error",
            "message": str(e),
            "traceback": traceback.format_exc()
        }


@app.get("/teste-audio-curto")
def teste_audio_curto():
    """
    Gera um áudio de teste curto (2 falas) e faz upload para o Firebase
    para validar o fluxo completo sem gastar muitos créditos.
    """
    try:
        import os
        from elevenlabs.client import ElevenLabs
        from elevenlabs import VoiceSettings
        from pydub import AudioSegment
        from firebase_service import upload_file
        
        # Configuração
        API_KEY = os.environ.get("ELEVENLABS_API_KEY")
        VOICE_HOST = os.environ.get("ELEVEN_VOICE_ID_HOST", "p5oveq8dCbyBIAaD6gzR")
        VOICE_COHOST = os.environ.get("ELEVEN_VOICE_ID_COHOST", "x3mAOLD9WzlmrFCwA1S3")
        
        if not API_KEY:
            return {"status": "error", "message": "Sem chave ElevenLabs configurada."}
            
        client = ElevenLabs(api_key=API_KEY)
        
        # Roteiro curto
        roteiro = [
            {"speaker": "HOST", "text": "Olá, este é um teste rápido do RevaCast para validar o sistema."},
            {"speaker": "COHOST", "text": "Exato! Estamos testando a velocidade e o upload para o Firebase."}
        ]
        
        audios_temp = []
        base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        os.makedirs(base_dir, exist_ok=True)
        
        # Gera áudios
        for i, fala in enumerate(roteiro):
            voice_id = VOICE_HOST if fala['speaker'] == "HOST" else VOICE_COHOST
            model = "eleven_turbo_v2_5" if fala['speaker'] == "HOST" else "eleven_multilingual_v2"
            
            audio_gen = client.text_to_speech.convert(
                voice_id=voice_id,
                text=fala['text'],
                model_id=model,
                voice_settings=VoiceSettings(stability=0.4, similarity_boost=0.8, style=0.6, use_speaker_boost=True)
            )
            
            temp_path = os.path.join(base_dir, f"temp_test_{i}.mp3")
            with open(temp_path, "wb") as f:
                for chunk in audio_gen: f.write(chunk)
            
            # Aplica velocidade
            seg = AudioSegment.from_file(temp_path)
            speed = 1.0 if fala['speaker'] == "HOST" else 1.30
            
            # Só aplica speedup se for diferente de 1.0 para evitar processamento desnecessário/ruído
            if speed != 1.0:
                seg_fast = seg.speedup(playback_speed=speed)
                seg_fast.export(temp_path, format="mp3")
                audios_temp.append(seg_fast) # Append the modified segment
            else:
                # Se for 1.0, mantém o arquivo original (ElevenLabs direto)
                # No caso de 1.0, seg_fast não é criado, então usamos 'seg'
                audios_temp.append(seg)
            
        # Junta tudo
        final_audio = AudioSegment.empty()
        
        # Adiciona Intro se existir
        intro_path = os.path.join(base_dir, "intro_guto.mp3")
        if os.path.exists(intro_path):
            try:
                intro = AudioSegment.from_file(intro_path)
                final_audio += intro + AudioSegment.silent(duration=1000)
            except Exception as e:
                print(f"Erro ao carregar intro no teste: {e}")
        
        for a in audios_temp:
            final_audio += a + AudioSegment.silent(duration=500)
            
        final_path = os.path.join(base_dir, "teste_fluxo_completo.mp3")
        final_audio.export(final_path, format="mp3")
        
        # Upload
        url = upload_file(final_path, "testes/teste_fluxo_completo.mp3")
        
        return {
            "status": "success",
            "message": "Áudio gerado e enviado com sucesso!",
            "url": url,
            "host_voice": VOICE_HOST,
            "cohost_voice": VOICE_COHOST
        }
        
    except Exception as e:
        import traceback
        return {
            "status": "error",
            "message": str(e),
            "traceback": traceback.format_exc()
        }


@app.get("/rescue-audio")
def rescue_audio():
    """
    Tenta encontrar arquivos MP3 na pasta data/ e fazer upload manual para o Firebase.
    Útil se o processo falhou logo após gerar o áudio mas antes de subir.
    """
    try:
        import os
        import glob
        from firebase_service import upload_file
        
        base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        mp3_files = glob.glob(os.path.join(base_dir, "*.mp3"))
        
        results = []
        
        if not mp3_files:
            return {"status": "warning", "message": "Nenhum arquivo MP3 encontrado na pasta data/."}
            
        for mp3_path in mp3_files:
            filename = os.path.basename(mp3_path)
            # Evita subir a intro se ela estiver lá
            if "intro" in filename: continue
            
            destination = f"resgate/{filename}"
            url = upload_file(mp3_path, destination)
            results.append({"file": filename, "url": url})
            
        return {
            "status": "success",
            "message": f"{len(results)} arquivos resgatados.",
            "files": results
        }
        
    except Exception as e:
        import traceback
        return {
            "status": "error",
            "message": str(e),
            "traceback": traceback.format_exc()
        }
