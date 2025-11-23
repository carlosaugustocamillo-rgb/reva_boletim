# api.py
#
# API simples com FastAPI para expor a função rodar_boletim()
# do arquivo boletim_service.py

from fastapi import FastAPI
from fastapi.responses import JSONResponse
from boletim_service import rodar_boletim

app = FastAPI(title="RevaCast Boletim Service")

@app.get("/")
def root():
    return {"message": "Serviço do boletim científico está no ar."}

@app.post("/rodar-boletim")
def rodar_boletim_endpoint():
    try:
        resultado = rodar_boletim()
        return JSONResponse(content={"success": True, "resultado": resultado})
    except Exception as e:
        return JSONResponse(
            status_code=500,
            content={"success": False, "erro": str(e)}
        )
