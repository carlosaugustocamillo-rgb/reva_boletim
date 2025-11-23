# api.py

from fastapi import FastAPI
from fastapi.responses import JSONResponse

from boletim_service import rodar_boletim
from revacast_workflow import WorkflowInput, run_workflow

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
async def agente_revacast(input: WorkflowInput):
    """
    Endpoint novo, que conversa com o AGENTE.

    - Recebe a mensagem do usuário em input.input_as_text
    - Roda o workflow do agente
    - O agente decide quando chamar a tool rodar_boletim(tipo_execucao=...)
    - A tool chama o endpoint /rodar-boletim do serviço de boletim (via HTTP)
    """
    return await run_workflow(input)
