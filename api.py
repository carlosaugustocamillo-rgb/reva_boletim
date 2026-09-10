
from fastapi import FastAPI, BackgroundTasks, Request, Response
from fastapi.responses import HTMLResponse, StreamingResponse, JSONResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware
import base64
import json
import uuid

from boletim_service import rodar_boletim
from campanha_service import criar_campanha_tematica
from revamais_service import criar_campanha_revamais
from simple_agent import run_agent, AgentInput
from pydantic import BaseModel, Field
import firebase_admin
from firebase_admin import credentials, storage
from datetime import date, datetime
import re
from elevenlabs_utils import diagnosticar_erro_elevenlabs
from connected_papers import build_manual_report, parse_bibtex
from pubmed_related import PubMedRelatedClient, RelatedConfig, prefilter

def simple_slugify(text):
    text = text.lower().strip()
    text = re.sub(r'[^\w\s-]', '', text)
    text = re.sub(r'[\s_-]+', '-', text)
    return text


def model_to_dict(model):
    if hasattr(model, "model_dump"):
        return model.model_dump()
    return model.dict()

class CampanhaInput(BaseModel):
    tema: str = ""


class ReferenciaSelecionadaInput(BaseModel):
    pmid: str | None = None
    texto: str
    link: str
    resumo: str = ""
    jif: float = 0.0
    journal: str = ""
    fonte: str = ""
    source_label: str = ""
    doi: str = ""
    tipo_estudo: str = ""
    achado_principal: str = ""
    consensus_claim: str = ""

class NewsPayload(BaseModel):
    titulo: str
    html_content: str
    data_publicacao: str
    imagem_capa: str = None
    resumo: str = None
    autor: str = "Revalidatie"


class ConnectedPapersImportInput(BaseModel):
    content: str
    filename: str = "connected_papers.bib"


class ManualPodcastContextInput(BaseModel):
    main_articles: list[dict] = Field(default_factory=list)
    selected_by_anchor: dict[str, list[dict]] = Field(default_factory=dict)


class PubMedSimilarInput(BaseModel):
    article: dict = Field(default_factory=dict)


# Force rebuild for Python 3.11
app = FastAPI(title="Revahub")

# Configuração de CORS Permissiva (Resolve problemas com WebContainers e Ambientes de Dev)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Permite qualquer origem
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
@app.get("/ping")
def ping():
    return {"status": "ok", "version": "with_cors_fix_v2", "cors": "enabled"}

# --- Armazenamento em ARQUIVO para tarefas em background (Persistência) ---
import os
import json

# Define diretório global de tarefas
TASK_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "tasks")
os.makedirs(TASK_DIR, exist_ok=True)

def save_task(task_id, data):
    """Salva o estado da tarefa em um arquivo JSON."""
    try:
        filepath = os.path.join(TASK_DIR, f"{task_id}.json")
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"Erro ao salvar tarefa {task_id}: {e}")

def load_task(task_id):
    """Carrega o estado da tarefa do arquivo JSON."""
    try:
        filepath = os.path.join(TASK_DIR, f"{task_id}.json")
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
        
    except Exception as e:
        import traceback
        error_msg = f"❌ Erro fatal: {str(e)}\n{traceback.format_exc()}"
        task_state["status"] = "error"
        task_state["logs"].append(error_msg)
        save_task(task_id, task_state)
        print(error_msg)

        print(error_msg)

# Helper factories to cleaner code
def make_log_callback(task_id):
    def log(msg):
        # Carrega estado atual para não sobrescrever logs anteriores
        state = load_task(task_id) or {"logs": [], "status": "running"}
        state["logs"].append(msg)
        save_task(task_id, state)
    return log

def make_check_cancel(task_id):
    def check():
        state = load_task(task_id)
        if state and state.get("status") == "canceling":
            # Atualiza para canceled
            state["status"] = "canceled"
            state["logs"].append("🛑 Cancelado pelo usuário.")
            save_task(task_id, state)
            return True
        return False
    return check

def update_task_status(task_id, status, log_msg=None):
    state = load_task(task_id) or {"logs": []}
    state["status"] = status
    if log_msg:
        state["logs"].append(log_msg)
    save_task(task_id, state)

def processar_revamais_background(task_id: str, opcoes: dict):
    """Função wrapper que roda o Reva+ e salva logs em arquivo. Agora com suporte a flags."""
    log_callback = make_log_callback(task_id)
    check_cancel = make_check_cancel(task_id)
    
    # Se recebeu string antiga (legacy), converte para dict
    if isinstance(opcoes, str):
        opcoes = {"tema_usuario": opcoes}
        
    tema = opcoes.get("tema_usuario")
    
    # Atualiza status inicial
    update_task_status(task_id, "running", f"🚀 Iniciando geração do Reva+ para: {tema}")
        
    try:
        from revamais_service import criar_campanha_revamais
        # Chama o serviço com as opções desempacotadas
        resultado = criar_campanha_revamais(
            tema_usuario=tema,
            calendar_index=opcoes.get("calendar_index"),
            gerar_midia=opcoes.get("gerar_midia", True),
            gerar_instagram=opcoes.get("gerar_instagram", True),
            enviar_email=opcoes.get("enviar_email", True),
            referencias_selecionadas=opcoes.get("referencias_selecionadas", []),
            relatorio_consensus=opcoes.get("relatorio_consensus", ""),
            log_callback=log_callback,
            check_cancel=check_cancel
        )

        # Se retornou sucesso/erro (dict)
        if isinstance(resultado, dict):
            if resultado.get("status") == "error":
                task_state = load_task(task_id) or {"logs": []} # reload state
                task_state["status"] = "error"
                task_state["result"] = resultado
                save_task(task_id, task_state)
            else:
                task_state = load_task(task_id) or {"logs": []}
                task_state["status"] = "completed"
                task_state["result"] = resultado
                task_state["logs"].append("✅ Processo Reva+ finalizado.")
                save_task(task_id, task_state)
        else:
            task_state = load_task(task_id) or {"logs": []}
            task_state["status"] = "completed"
            save_task(task_id, task_state)
        
    except Exception as e:
        if str(e) == "CANCELADO_PELO_USUARIO":
             return # Já salvou status canceled no check_cancel
             
        import traceback
        error_msg = f"❌ Erro fatal: {str(e)}\n{traceback.format_exc()}"
        task_state = load_task(task_id) or {"logs": []}
        task_state["status"] = "error"
        task_state["logs"].append(error_msg)
        save_task(task_id, task_state)


def ensure_storage_ready():
    if firebase_admin._apps:
        return

    cred_path = os.environ.get("FIREBASE_CREDENTIALS_JSON", "firebase_credentials.json")
    if os.path.exists(cred_path):
        cred = credentials.Certificate(cred_path)
        firebase_admin.initialize_app(cred, {
            'storageBucket': os.environ.get("FIREBASE_BUCKET_NAME")
        })


def serialize_storage_blob(blob, tipo):
    try:
        blob.make_public()
    except Exception:
        pass

    return {
        "tipo": tipo,
        "nome": blob.name.split("/")[-1],
        "url": blob.public_url,
        "data": blob.time_created.isoformat() if getattr(blob, "time_created", None) else "",
        "tamanho": f"{blob.size / 1024 / 1024:.2f} MB" if getattr(blob, "size", None) else "",
    }


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
        ensure_storage_ready()
        
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


@app.get("/listar-instagram-revamais")
def listar_instagram_revamais():
    """Lista os últimos reels, roteiros e carrosséis gerados no Reva+."""
    try:
        ensure_storage_ready()

        bucket = storage.bucket()
        reels = []
        roteiros = []
        carrosseis = []

        for blob in bucket.list_blobs(prefix="instagram/"):
            nome = blob.name.split("/")[-1]
            if nome.startswith("instagram_reel_") and nome.endswith(".md"):
                reels.append(serialize_storage_blob(blob, "reel"))
            elif nome.startswith("roteiro_carrossel_") and nome.endswith(".md"):
                roteiros.append(serialize_storage_blob(blob, "roteiro"))
            elif nome.startswith("instagram_carrossel_") and nome.endswith(".zip"):
                carrosseis.append(serialize_storage_blob(blob, "carrossel"))

        reels.sort(key=lambda item: item["data"], reverse=True)
        roteiros.sort(key=lambda item: item["data"], reverse=True)
        carrosseis.sort(key=lambda item: item["data"], reverse=True)

        return {
            "reels": reels[:5],
            "roteiros": roteiros[:5],
            "carrosseis": carrosseis[:5],
        }
    except Exception as e:
        return {"error": str(e)}

@app.get("/")
def root():
    return {"message": "Serviço do boletim científico está no ar. Acesse /painel para interface gráfica."}

@app.get("/health")
def health_check():
    return {"status": "ok", "python_version": "unknown"}


@app.api_route("/favicon.ico", methods=["GET", "HEAD"])
def favicon():
    return Response(status_code=204)


@app.api_route("/.well-known/appspecific/com.chrome.devtools.json", methods=["GET", "HEAD"])
def chrome_devtools_metadata():
    return JSONResponse(content={})

@app.get("/painel", response_class=HTMLResponse)
def painel():
    with open("dashboard.html", "r", encoding="utf-8") as f:
        return f.read()

# --- Novos Endpoints Background ---

def processar_medico_background(task_id: str, dry_run: bool):
    """Função wrapper que roda o boletim MEDICO e salva logs."""
    log_callback = make_log_callback(task_id)
    
    update_task_status(task_id, "running", "🚀 Iniciando Boletim MEDICO...")
    
    try:
        from medico_boletim_service import rodar_boletim_medico
        rodar_boletim_medico(dry_run=dry_run, log_callback=log_callback)
        
        update_task_status(task_id, "completed", "✅ Boletim MEDICO finalizado.")
        
    except Exception as e:
        import traceback
        error_msg = f"❌ Erro fatal: {str(e)}\n{traceback.format_exc()}"
        state = load_task(task_id) or {"logs": []}
        state["status"] = "error"
        state["logs"].append(error_msg)
        save_task(task_id, state)

@app.post("/iniciar-boletim-medico")
def iniciar_boletim_medico(
    background_tasks: BackgroundTasks,
    dry_run: bool = False
):
    """
    Inicia a geração do Boletim Médico Semanal.
    Use dry_run=True para teste (sem Mailchimp).
    """
    task_id = str(uuid.uuid4())
    save_task(task_id, {"status": "queued", "logs": ["⏳ Iniciando..."], "result": None})
    
    background_tasks.add_task(processar_medico_background, task_id, dry_run)
    
    return {
        "task_id": task_id,
        "status": "started",
        "message": "Boletim Médico iniciado. Use /status-boletim/{task_id} para acompanhar."
    }

@app.post("/iniciar-boletim")
def iniciar_boletim(
    background_tasks: BackgroundTasks,
    resumos: bool = True,
    roteiro: bool = True,
    revisao_roteiro: bool = True,
    brief_spotify: bool = True,
    audio: bool = True,
    mailchimp: bool = True,
    firebase: bool = True,
    referencias_pubmed: bool = True,
    payload: dict | None = None,
):
    body = payload or {}
    query_options = {
        'resumos': resumos,
        'roteiro': roteiro,
        'revisao_roteiro': revisao_roteiro,
        'brief_spotify': brief_spotify,
        'audio': audio,
        'mailchimp': mailchimp,
        'firebase': firebase,
        'referencias_pubmed': referencias_pubmed,
    }
    for key in query_options:
        value = body.get(key)
        if value is not None:
            query_options[key] = value
    task_id = str(uuid.uuid4())
    opcoes = query_options
    opcoes['somente_curadoria'] = bool(body.get('somente_curadoria'))
    opcoes['artigos_podcast_aprovados'] = body.get('artigos_podcast_aprovados') or []
    opcoes['contexto_pubmed_manual'] = body.get('contexto_pubmed_manual')
    
    # Cria o arquivo inicial
    save_task(task_id, {"status": "queued", "logs": ["⏳ Iniciando..."], "result": None})
    
    # Inicia a tarefa em background
    background_tasks.add_task(processar_boletim_background, task_id, opcoes)
    
    return {
        "task_id": task_id,
        "status": "started",
        "message": "Boletim iniciado em segundo plano. Verifique o status com o ID fornecido."
    }


@app.post("/importar-connected-papers")
def importar_connected_papers(payload: ConnectedPapersImportInput):
    """Analisa um BibTeX exportado manualmente do Connected Papers."""
    try:
        result = parse_bibtex(payload.content)
        result["arquivo_nome"] = payload.filename
        return result
    except ValueError as error:
        return JSONResponse(status_code=400, content={"error": str(error)})


@app.post("/preparar-contexto-manual")
def preparar_contexto_manual(payload: ManualPodcastContextInput):
    """Valida a seleção manual e monta o contrato usado pelo roteirista."""
    if not payload.main_articles:
        return JSONResponse(status_code=400, content={"error": "Selecione ao menos um artigo principal."})
    return build_manual_report(payload.main_articles, payload.selected_by_anchor)


@app.post("/buscar-similares-pubmed")
def buscar_similares_pubmed(payload: PubMedSimilarInput):
    """Busca similares no PubMed sob demanda, para substituir um grafo indisponível."""
    article = payload.article or {}
    pmid = str(article.get("pmid", "")).strip()
    if not re.fullmatch(r"\d+", pmid):
        return JSONResponse(status_code=400, content={"error": "O artigo precisa ter um PMID válido."})

    cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "cache", "pubmed_related")
    pubmed = PubMedRelatedClient(
        cache_dir,
        email=os.environ.get("ENTREZ_EMAIL", ""),
        api_key=os.environ.get("NCBI_API_KEY"),
        config=RelatedConfig.from_env(),
    )
    try:
        candidates, cached = pubmed.candidates(pmid)
        eligible, rejected = prefilter(candidates, article, set(), date.today())
        return {
            "versao": 1,
            "fonte": "PubMed Similar Articles (fallback manual)",
            "ancora": article,
            "candidatos": eligible,
            "excluidos": rejected,
            "cache": cached,
            "total_candidatos": len(eligible),
            "mensagem": (
                "Nenhum similar anterior com resumo disponível foi encontrado."
                if not eligible else "Selecione manualmente os artigos que deseja usar como contexto."
            ),
        }
    except Exception as error:
        return JSONResponse(status_code=502, content={"error": f"PubMed indisponível: {type(error).__name__}."})
    finally:
        pubmed.session.close()

class RevaMaisInput(BaseModel):
    tema: str = ""
    calendar_index: int | None = None
    gerar_midia: bool = True
    gerar_instagram: bool = True
    enviar_email: bool = True
    referencias_selecionadas: list[ReferenciaSelecionadaInput] = Field(default_factory=list)
    relatorio_consensus: str = ""


class RevaMaisPrepareInput(BaseModel):
    tema: str = ""
    calendar_index: int | None = None
    quantidade_referencias: int = 8


class RevaMaisConsensusInput(BaseModel):
    tema: str = ""
    calendar_index: int | None = None
    relatorio: str
    quantidade_referencias: int | None = None


class RevaMaisConsensusPdfInput(BaseModel):
    tema: str = ""
    calendar_index: int | None = None
    filename: str = "consensus.pdf"
    pdf_base64: str
    quantidade_referencias: int | None = None

@app.post("/iniciar-revamais")
def iniciar_revamais(
    background_tasks: BackgroundTasks,
    input_data: RevaMaisInput
):
    try:
        task_id = str(uuid.uuid4())
        
        # Cria estado inicial
        save_task(task_id, {"status": "queued", "logs": ["⏳ Iniciando Reva+ ..."], "result": None})
        
        # Prepara opções
        opcoes = {
            "tema_usuario": input_data.tema,
            "calendar_index": input_data.calendar_index,
            "gerar_midia": input_data.gerar_midia,
            "gerar_instagram": input_data.gerar_instagram,
            "enviar_email": input_data.enviar_email,
            "referencias_selecionadas": [model_to_dict(ref) for ref in input_data.referencias_selecionadas],
            "relatorio_consensus": input_data.relatorio_consensus,
        }
        
        background_tasks.add_task(processar_revamais_background, task_id, opcoes)
        
        return {
            "task_id": task_id,
            "status": "started",
            "message": "Reva+ iniciado em background."
        }
    except Exception as e:
        import traceback
        print(f"❌ CRITICAL ERROR IN /iniciar-revamais: {e}\n{traceback.format_exc()}")
        return JSONResponse(
            status_code=500,
            content={
                "status": "error", 
                "message": "Critical Server Error during initialization",
                "error": str(e), 
                "traceback": traceback.format_exc()
            }
        )

@app.post("/cancelar-tarefa/{task_id}")
def cancelar_tarefa_generica(task_id: str):
    return cancelar_boletim(task_id) # Reutiliza a mesma lógica


@app.get("/calendario-revamais")
def calendario_revamais_endpoint():
    try:
        from revamais_service import listar_calendario_revamais

        resultado = listar_calendario_revamais()
        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )


@app.post("/preparar-revamais")
def preparar_revamais_endpoint(input_data: RevaMaisPrepareInput):
    """
    Resolve o tema e retorna os artigos candidatos para curadoria manual
    antes da geração final do Reva+.
    """
    try:
        from revamais_service import preparar_referencias_revamais

        resultado = preparar_referencias_revamais(
            tema_usuario=input_data.tema,
            calendar_index=input_data.calendar_index,
            quantidade_referencias=input_data.quantidade_referencias,
        )

        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )


@app.post("/importar-consensus-revamais")
def importar_consensus_revamais_endpoint(input_data: RevaMaisConsensusInput):
    """
    Importa um relatório colado/exportado do Consensus.app e retorna referências
    no mesmo formato usado pela curadoria manual do Reva+.
    """
    try:
        from revamais_service import preparar_referencias_consensus_revamais

        resultado = preparar_referencias_consensus_revamais(
            tema_usuario=input_data.tema,
            calendar_index=input_data.calendar_index,
            relatorio_texto=input_data.relatorio,
            quantidade_referencias=input_data.quantidade_referencias,
        )

        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )


@app.post("/importar-consensus-revamais-pdf")
def importar_consensus_revamais_pdf_endpoint(input_data: RevaMaisConsensusPdfInput):
    """
    Recebe um PDF do Consensus em base64, extrai o texto e devolve o mesmo
    payload usado pela curadoria do Reva+.
    """
    try:
        from revamais_service import (
            extrair_texto_pdf_consensus,
            preparar_referencias_consensus_revamais,
        )

        pdf_base64 = input_data.pdf_base64.strip()
        if "," in pdf_base64:
            pdf_base64 = pdf_base64.split(",", 1)[1]

        pdf_bytes = base64.b64decode(pdf_base64)
        relatorio_texto = extrair_texto_pdf_consensus(pdf_bytes)

        resultado = preparar_referencias_consensus_revamais(
            tema_usuario=input_data.tema,
            calendar_index=input_data.calendar_index,
            relatorio_texto=relatorio_texto,
            quantidade_referencias=input_data.quantidade_referencias,
        )
        resultado["arquivo_importado"] = input_data.filename

        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)},
        )

@app.get("/status-boletim/{task_id}")
async def status_boletim(task_id: str):
    status_file = os.path.join(TASK_DIR, f"{task_id}.json")
    if not os.path.exists(status_file):
        # Tenta achar log antigo ou retorna 404
        return {"status": "not_found"}
        
    with open(status_file, "r") as f:
        data = json.load(f)
    return data

@app.get("/ultimo-roteiro")
async def get_ultimo_roteiro():
    """Retorna o JSON do último roteiro gerado para inspeção."""
    roteiro_dir = os.path.join("data", "roteiros")
    if not os.path.exists(roteiro_dir):
        return {"error": "Pasta de roteiros não encontrada."}
    
    arquivos = sorted([f for f in os.listdir(roteiro_dir) if f.endswith(".json")], reverse=True)
    if not arquivos:
        return {"error": "Nenhum roteiro encontrado."}
    
    ultimo_arquivo = os.path.join(roteiro_dir, arquivos[0])
    with open(ultimo_arquivo, "r", encoding="utf-8") as f:
        conteudo = json.load(f)
        
    return {"arquivo": arquivos[0], "conteudo": conteudo}


@app.get("/baixar-brief/{data_ref}")
async def baixar_brief_spotify(data_ref: str):
    """
    Download do brief gerado para Spotify no formato brief_spotify_YYYY-MM-DD.txt
    """
    if not re.match(r"^\d{4}-\d{2}-\d{2}$", data_ref):
        return JSONResponse(
            status_code=400,
            content={"error": "Formato de data inválido. Use YYYY-MM-DD."}
        )

    base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    brief_path = os.path.join(base_dir, f"brief_spotify_{data_ref}.txt")

    if not os.path.exists(brief_path):
        return JSONResponse(
            status_code=404,
            content={"error": f"Brief não encontrado para {data_ref}."}
        )

    return FileResponse(
        path=brief_path,
        filename=f"brief_spotify_{data_ref}.txt",
        media_type="text/plain; charset=utf-8"
    )

@app.get("/baixar-referencias-podcast/{data_ref}")
async def baixar_referencias_podcast(data_ref: str):
    """Relatório de fontes e triagem, separado do boletim enviado por e-mail."""
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", data_ref):
        return JSONResponse(status_code=400, content={"error": "Use YYYY-MM-DD."})
    try:
        datetime.strptime(data_ref, "%Y-%m-%d")
    except ValueError:
        return JSONResponse(status_code=400, content={"error": "Data inválida."})
    filename = f"contexto_pubmed_{data_ref}.json"
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "referencias", filename)
    if not os.path.isfile(path):
        return JSONResponse(status_code=404, content={"error": "Referências não encontradas para esta edição."})
    return FileResponse(path=path, filename=filename, media_type="application/json")


@app.get("/testar-firebase")
async def start_firebase_test():
    """Rota para diagnosticar conexão com Firebase (Storage e Firestore)."""
    results = {"storage": "pending", "firestore": "pending", "logs": []}
    
    # 1. Teste Firestore
    try:
        from firebase_service import save_firestore_document
        test_data = {"teste": True, "data": datetime.now().isoformat(), "msg": "Verificação de escrita"}
        sucesso = save_firestore_document("debug_testes", "teste_conexao", test_data)
        if sucesso:
            results["firestore"] = "✅ SUCESSO! Coleção 'debug_testes' criada."
        else:
            results["firestore"] = "❌ FALHA. Verifique logs do servidor."
    except Exception as e:
        results["firestore"] = f"❌ ERRO EXCEÇÃO: {str(e)}"

    # 2. Teste Storage
    try:
        from firebase_service import upload_file
        # Cria arquivo temp
        temp_file = "teste_storage.txt"
        with open(temp_file, "w") as f: f.write("Teste de upload do RevaCast.")
        
        url = upload_file(temp_file, "test-uploads/teste_conexao.txt")
        if url:
             results["storage"] = f"✅ SUCESSO! Arquivo salvo em test-uploads/. URL: {url}"
        else:
             results["storage"] = "❌ FALHA no Upload (URL vazia)."
        
        if os.path.exists(temp_file): os.remove(temp_file)
    except Exception as e:
        results["storage"] = f"❌ ERRO EXCEÇÃO: {str(e)}"
        
    return results
    return results

@app.post("/remixar-audio")
def endpoint_remixar_audio():
    """Tenta recuperar e remixar áudios do dia direto do Firebase."""
    try:
        from boletim_service import remixar_audio_from_firebase
        resultado = remixar_audio_from_firebase()
        return resultado
    except Exception as e:
        return {"success": False, "error": str(e)}

@app.post("/criar-revamais")
def criar_revamais_endpoint(input_data: RevaMaisInput):
    """
    Cria uma edição do Reva + (newsletter temática).
    """
    print(f"🔎 DEBUG: Request recebido no endpoint /criar-revamais. Tema: {input_data.tema}")
    try:
        from revamais_service import criar_campanha_revamais
        resultado = criar_campanha_revamais(
            tema_usuario=input_data.tema,
            calendar_index=input_data.calendar_index,
            gerar_midia=input_data.gerar_midia,
            gerar_instagram=input_data.gerar_instagram,
            enviar_email=input_data.enviar_email,
            referencias_selecionadas=[model_to_dict(ref) for ref in input_data.referencias_selecionadas],
            relatorio_consensus=input_data.relatorio_consensus,
        )
        
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

# --- Novo Endpoint para Publicar no Site ---
@app.post("/publicar-noticia")
async def publicar_noticia_website(payload: NewsPayload):
    """
    Recebe o conteúdo aprovado pelo usuário e salva na coleção 'news' do Firestore.
    Isso permitirá que o site exiba a notícia automaticamente.
    """
    try:
        from firebase_service import save_firestore_document
        
        # Gera ID amigável (Slug)
        slug = simple_slugify(payload.titulo)
        # Adiciona timestamp curto para garantir unicidade
        doc_id = f"{slug}-{datetime.now().strftime('%Y%m%d')}"
        
        doc_data = {
            "title": payload.titulo,
            "content": payload.html_content, # HTML interno
            "publishedAt": payload.data_publicacao, # String DD/MM/YYYY
            "createdAt": datetime.now().isoformat(),
            "coverImage": payload.imagem_capa,
            "summary": payload.resumo or "",
            "author": payload.autor,
            "status": "published",
            "slug": doc_id
        }
        
        # Salva no Firestore
        success = save_firestore_document("news", doc_id, doc_data)
        
        if success:
            return {"status": "success", "message": "Notícia publicada com sucesso!", "id": doc_id, "slug": doc_id}
        else:
            return JSONResponse(status_code=500, content={"status": "error", "message": "Falha ao salvar no Firestore."})
            
    except Exception as e:
        return JSONResponse(status_code=500, content={"status": "error", "message": str(e)})

# --- WhatsApp Automation Endpoints ---

class WhatsAppSendPayload(BaseModel):
    draft_id: str
    to_number: str = None  # If None, sends to default/channel if configured, or error
    override_text: str = None # Allows editing before sending

class WhatsAppSchedulePayload(BaseModel):
    draft_id: str
    schedule_time: str # ISO format

@app.get("/whatsapp/drafts")
def list_whatsapp_drafts(status: str = None):
    """Lists generated WhatsApp message drafts."""
    try:
        from whatsapp_service import list_drafts
        drafts = list_drafts(status)
        return drafts
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})

@app.post("/whatsapp/send")
def send_whatsapp_message(payload: WhatsAppSendPayload):
    """Sends a WhatsApp message from a draft."""
    try:
        from whatsapp_service import send_message, update_draft_status, list_drafts
        
        # 1. Get draft
        # Note: list_drafts is not efficient for getting one, but for now it works. 
        # Ideal: get_draft(id) in service.
        # Let's assume the frontend passes the text or we fetch it.
        # For MVP, let's fetch all and filter (bad performance but safe for <100 drafts).
        # Better: Add get_draft to service. But to save steps, I'll use the one I have or generic firestore get.
        
        from firebase_service import read_json_from_storage # Wrong tool.
        # Let's use generic firestore get if needed or trust the payload.
        
        # Actually, let's trust the payload override_text if present, else we need to fetch.
        text_to_send = payload.override_text
        
        # If no text provided, we MUST fetch existing.
        if not text_to_send:
            # Quick hack: use the list drafts
             all_drafts = list_drafts()
             draft = next((d for d in all_drafts if d['id'] == payload.draft_id), None)
             if not draft:
                 return JSONResponse(status_code=404, content={"error": "Draft not found"})
             text_to_send = draft.get('generated_text')

        if not payload.to_number:
            # TODO: Add Channel ID or Admin Number in .env
            # For now, require it or fail.
             return JSONResponse(status_code=400, content={"error": "target phone number (to_number) is required for now."})

        # 2. Send
        result = send_message(payload.to_number, text_to_send)
        
        if "error" in result:
             return JSONResponse(status_code=500, content=result)
             
        # 3. Update status
        update_draft_status(payload.draft_id, "SENT")
        
        return {"status": "success", "api_response": result}
        
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})

@app.post("/whatsapp/schedule")
def schedule_whatsapp_message(payload: WhatsAppSchedulePayload):
    """Schedules a draft (just updates metadata for now)."""
    try:
        from whatsapp_service import update_draft_status
        success = update_draft_status(payload.draft_id, "SCHEDULED", payload.schedule_time)
        return {"success": success}
    except Exception as e:
         return JSONResponse(status_code=500, content={"error": str(e)})


# --- Endpoint Antigo (Mantido para compatibilidade) ---
@app.get("/rodar-boletim-stream")
def rodar_boletim_stream(
    resumos: bool = True,
    roteiro: bool = True,
    revisao_roteiro: bool = True,
    brief_spotify: bool = True,
    audio: bool = True,
    mailchimp: bool = True,
    firebase: bool = True,
    referencias_pubmed: bool = True,
):
    """
    Endpoint de streaming que envia atualizações de progresso via Server-Sent Events (SSE).
    Aceita query params para controlar quais etapas executar.
    """
    opcoes = {
        'resumos': resumos,
        'roteiro': roteiro,
        'revisao_roteiro': revisao_roteiro,
        'brief_spotify': brief_spotify,
        'audio': audio,
        'mailchimp': mailchimp,
        'firebase': firebase
    }
    opcoes['referencias_pubmed'] = referencias_pubmed

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
    Testa o mesmo Text to Dialogue/Eleven v3 usado pelo RevaCast Weekly
    e faz upload para o Firebase com baixo consumo de caracteres.
    """
    try:
        import os
        from boletim_service import gerar_dialogo_com_eleven_v3
        from firebase_service import upload_file

        api_key = os.environ.get("ELEVENLABS_API_KEY")
        VOICE_HOST = os.environ.get("ELEVEN_VOICE_ID_HOST", "p5oveq8dCbyBIAaD6gzR")
        VOICE_COHOST = os.environ.get("ELEVEN_VOICE_ID_COHOST", "tnSpp4vdxKPjI9w0GnoV")

        if not api_key:
            return {"status": "error", "message": "Sem chave ElevenLabs configurada."}

        roteiro = [
            {"speaker": "HOST", "text": "Olá, este é um teste rápido do RevaCast para validar o sistema."},
            {"speaker": "COHOST", "text": "Exato! Estamos testando a velocidade e o upload para o Firebase."}
        ]

        base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
        os.makedirs(base_dir, exist_ok=True)
        final_path = os.path.join(base_dir, "teste_fluxo_completo.mp3")
        gerar_dialogo_com_eleven_v3(roteiro, final_path)
        url = upload_file(final_path, "testes/teste_fluxo_completo.mp3")

        return {
            "status": "success",
            "message": "Áudio gerado e enviado com sucesso!",
            "url": url,
            "model": os.environ.get("ELEVEN_AUDIO_DIALOGUE_MODEL", "eleven_v3"),
            "host_voice": VOICE_HOST,
            "cohost_voice": VOICE_COHOST
        }

    except Exception as e:
        diagnostico = diagnosticar_erro_elevenlabs(e)
        return JSONResponse(status_code=424, content={
            "status": "error",
            **diagnostico,
        })


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
