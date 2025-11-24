from fastapi import FastAPI
from fastapi.responses import JSONResponse

from boletim_service import rodar_boletim
from simple_agent import run_agent, AgentInput

app = FastAPI(title="RevaCast Boletim Service")


@app.get("/")
def root():
    return {"message": "Serviço do boletim científico está no ar."}


@app.post("/rodar-boletim")
def rodar_boletim_endpoint():
    """
    Endpoint 'antigo' que dispara diretamente a função local boletim_service.rodar_boletim().
    Continua funcionando como antes.
    """
    try:
        resultado = rodar_boletim()
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
