
from fastapi import FastAPI
from fastapi.responses import HTMLResponse, StreamingResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import json

from boletim_service import rodar_boletim
from campanha_service import criar_campanha_tematica
from simple_agent import run_agent, AgentInput
from pydantic import BaseModel

class CampanhaInput(BaseModel):
    tema: str

app = FastAPI(title="RevaCast Boletim Service")

# Configuração de CORS para permitir requisições do frontend (localhost:3000, etc)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Em produção, substitua "*" pela URL do seu site
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
def root():
    return {"message": "Serviço do boletim científico está no ar. Acesse /painel para interface gráfica."}


@app.get("/painel", response_class=HTMLResponse)
def painel():
    with open("dashboard.html", "r", encoding="utf-8") as f:
        return f.read()


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
