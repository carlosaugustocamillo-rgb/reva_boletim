
from fastapi import FastAPI, BackgroundTasks
from fastapi.responses import HTMLResponse, StreamingResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import json
import uuid

from boletim_service import rodar_boletim
from campanha_service import criar_campanha_tematica
from simple_agent import run_agent, AgentInput
from pydantic import BaseModel

class CampanhaInput(BaseModel):
    tema: str

# Force rebuild for Python 3.11
app = FastAPI(title="RevaCast Boletim Service")

# Configuração de CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Armazenamento em memória para tarefas em background ---
tasks_db = {}

def processar_boletim_background(task_id: str, opcoes: dict):
    """Função wrapper que roda o boletim e salva logs na memória."""
    tasks_db[task_id] = {"status": "running", "logs": [], "result": None}
    try:
        for log_msg in rodar_boletim(opcoes):
            if isinstance(log_msg, dict):
                tasks_db[task_id]["result"] = log_msg
            else:
                tasks_db[task_id]["logs"].append(str(log_msg))
        
        tasks_db[task_id]["status"] = "completed"
        tasks_db[task_id]["logs"].append("✅ Processo finalizado com sucesso.")
        
    except Exception as e:
        import traceback
        error_msg = f"❌ Erro fatal: {str(e)}\n{traceback.format_exc()}"
        tasks_db[task_id]["status"] = "error"
        tasks_db[task_id]["logs"].append(error_msg)
        print(error_msg)

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
    
    # Inicia a tarefa em background
    background_tasks.add_task(processar_boletim_background, task_id, opcoes)
    
    return {
        "task_id": task_id,
        "status": "started",
        "message": "Boletim iniciado em segundo plano. Verifique o status com o ID fornecido."
    }

@app.get("/status-boletim/{task_id}")
def get_status_boletim(task_id: str):
    task = tasks_db.get(task_id)
    if not task:
        return {"status": "not_found", "message": "Tarefa não encontrada."}
    return task

# --- Endpoint Antigo (Mantido para compatibilidade, mas não recomendado para conexões instáveis) ---
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
