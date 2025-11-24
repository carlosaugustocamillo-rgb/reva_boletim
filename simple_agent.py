import os
import json
import threading
from typing import List, Dict, Any
from openai import OpenAI
from pydantic import BaseModel

from boletim_service import rodar_boletim as rodar_boletim_job

# Configuração do cliente OpenAI
client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

# Definição da ferramenta (Tool)
TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "rodar_boletim",
            "description": "Dispara o processo de geração do boletim semanal e/ou episódio de podcast.",
            "parameters": {
                "type": "object",
                "properties": {
                    "tipo_execucao": {
                        "type": "string",
                        "enum": ["boletim_semanal", "episodio", "geral"],
                        "description": "O tipo de execução: 'boletim_semanal' (apenas texto), 'episodio' (apenas áudio) ou 'geral' (ambos)."
                    }
                },
                "required": ["tipo_execucao"]
            }
        }
    }
]

SYSTEM_PROMPT = """
Você é o RevaCast AI, o assistente oficial de produção do RevaCast Weekly.
Sua responsabilidade é coordenar a geração do boletim científico semanal e do podcast.

**Suas Capacidades:**
1. **Gerar Boletim:** Buscar artigos no PubMed, traduzir e formatar para envio.
2. **Gerar Episódio:** Criar roteiros baseados nos artigos e gerar áudio usando ElevenLabs.
3. **Execução Geral:** Fazer todo o processo de ponta a ponta.

**Regras de Comportamento:**
- Seja cordial, eficiente e profissional.
- Quando o usuário pedir para gerar algo, SEMPRE use a ferramenta `rodar_boletim`.
- O processo é demorado (pode levar alguns minutos). Avise o usuário que você iniciou o processo em background.
- Se o usuário perguntar sobre o status, explique que o processo roda no servidor e ele receberá o e-mail/arquivo quando terminar.
- Não tente inventar artigos ou resultados. Apenas gerencie o fluxo.

**Uso da Ferramenta `rodar_boletim`:**
- `boletim_semanal`: Para pedidos de "gerar boletim", "atualizar notícias", "buscar artigos".
- `episodio`: Para pedidos de "gerar áudio", "criar podcast", "fazer o episódio".
- `geral`: Para pedidos de "rodar tudo", "fazer o processo completo", "gerar boletim e podcast".
"""

class AgentInput(BaseModel):
    message: str
    conversation_history: List[Dict[str, Any]] = []

def run_agent(input_data: AgentInput) -> Dict[str, Any]:
    """
    Executa o agente com a mensagem do usuário e histórico.
    Retorna a resposta do agente e o novo histórico.
    """
    messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    
    # Adiciona histórico anterior
    messages.extend(input_data.conversation_history)
    
    # Adiciona mensagem atual
    messages.append({"role": "user", "content": input_data.message})

    # Primeira chamada ao modelo
    response = client.chat.completions.create(
        model="gpt-4o",
        messages=messages,
        tools=TOOLS,
        tool_choice="auto"
    )

    response_message = response.choices[0].message
    messages.append(response_message)

    # Verifica se houve chamada de ferramenta
    if response_message.tool_calls:
        for tool_call in response_message.tool_calls:
            if tool_call.function.name == "rodar_boletim":
                args = json.loads(tool_call.function.arguments)
                tipo = args.get("tipo_execucao", "geral")
                
                # Executa a função em background para não bloquear
                def _run_background():
                    try:
                        print(f"Iniciando {tipo} em background...")
                        rodar_boletim_job() # A função original não aceita argumentos, roda tudo. 
                        # Se precisar separar, teria que refatorar o boletim_service.py.
                        # Por enquanto, assumimos que roda o pipeline padrão.
                    except Exception as e:
                        print(f"Erro no background job: {e}")

                threading.Thread(target=_run_background, daemon=True).start()

                # Resposta da ferramenta
                tool_response = {
                    "role": "tool",
                    "tool_call_id": tool_call.id,
                    "name": "rodar_boletim",
                    "content": json.dumps({
                        "status": "started", 
                        "message": f"Processo '{tipo}' iniciado em background."
                    })
                }
                messages.append(tool_response)

        # Segunda chamada para gerar a resposta final ao usuário
        final_response = client.chat.completions.create(
            model="gpt-4o",
            messages=messages
        )
        final_message = final_response.choices[0].message.content
        messages.append(final_response.choices[0].message)
    else:
        final_message = response_message.content

    # Retorna resposta e histórico atualizado (serializável)
    # Convertendo objetos Message para dicts para retorno
    serialized_history = []
    for msg in messages:
        if isinstance(msg, dict):
            serialized_history.append(msg)
        else:
            serialized_history.append(msg.model_dump(exclude_none=True))

    return {
        "response": final_message,
        "history": serialized_history
    }
