import requests
import json
import sys

def chat_with_agent(message):
    url = "http://127.0.0.1:8000/agente-revacast"
    
    payload = {
        "message": message,
        "conversation_history": []
    }
    
    try:
        print(f"Você: {message}")
        print("Agente: (pensando...)")
        
        response = requests.post(url, json=payload)
        response.raise_for_status()
        
        data = response.json()
        agent_response = data.get("response", "Sem resposta")
        
        print(f"Agente: {agent_response}")
        return data.get("history", [])
        
    except requests.exceptions.ConnectionError:
        print("\n❌ Erro: Não foi possível conectar ao servidor.")
        print("Certifique-se de que o servidor está rodando com o comando:")
        print("uvicorn main:app --reload")
    except Exception as e:
        print(f"\n❌ Erro: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        msg = " ".join(sys.argv[1:])
    else:
        msg = "Gere o boletim semanal"
        
    chat_with_agent(msg)
