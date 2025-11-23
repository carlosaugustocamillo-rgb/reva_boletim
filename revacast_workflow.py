import threading
from typing import Dict

from agents import (
    Agent,
    ModelSettings,
    RunConfig,
    Runner,
    TResponseInputItem,
    function_tool,
    trace,
)
from pydantic import BaseModel
from boletim_service import rodar_boletim as rodar_boletim_job

# Tool definitions
@function_tool
def rodar_boletim(tipo_execucao: str) -> Dict[str, str]:
    """
    Dispara o backend do REVABOLETIM para gerar boletim/episódio.

    Args:
        tipo_execucao: "boletim_semanal", "episodio" ou "geral".
    """
    # Dispara o job pesado em background para não bloquear a chamada da tool.
    def _run_boletim_background():
        try:
            rodar_boletim_job()
        except Exception as exc:  # pragma: no cover - apenas log de proteção
            print(f"Erro ao rodar boletim em background ({tipo_execucao}): {exc}")

    threading.Thread(target=_run_boletim_background, daemon=True).start()

    return {
        "status": "started",
        "tipo_execucao": tipo_execucao,
        "message": "O processo foi iniciado em background.",
    }


revacast_agent = Agent(
    name="Revacast Weekly",
    instructions="""Você é o agente oficial do REVACAST Weekly.
Sua função é acionar o backend sempre que o usuário pedir para gerar o boletim semanal,
gerar episódio, atualizar boletim, criar podcast ou expressões semelhantes.

Você NUNCA deve tentar realizar processamento interno.
Sempre que o usuário fizer um pedido relacionado ao boletim, você deve acionar a ferramenta `rodar_boletim`.

Quando usar a ferramenta `rodar_boletim`, você deve sempre preencher o campo `tipo_execucao` assim:
- Use \"boletim_semanal\" quando o usuário pedir para gerar ou atualizar o boletim semanal.
- Use \"episodio\" quando o usuário pedir para criar o episódio de podcast com base no boletim.
- Use \"geral\" quando o usuário pedir para rodar tudo (boletim + episódio).

Após chamar a ferramenta, aguarde a resposta e informe ao usuário o resultado.
""",
    model="gpt-4.1-mini",
    tools=[rodar_boletim],
    model_settings=ModelSettings(
        temperature=1,
        top_p=1,
        parallel_tool_calls=True,
        max_tokens=2048,
        store=True,
    )
)


class WorkflowInput(BaseModel):
    input_as_text: str


# Main code entrypoint
async def run_workflow(workflow_input: WorkflowInput) -> Dict[str, str]:
    with trace("New workflow"):
        conversation_history: list[TResponseInputItem] = [
            {
                "role": "user",
                "content": [
                    {
                        "type": "input_text",
                        "text": workflow_input.input_as_text,
                    }
                ]
            }
        ]

        revacast_result_temp = await Runner.run(
            revacast_agent,
            input=[*conversation_history],
            run_config=RunConfig(
                trace_metadata={
                    "__trace_source__": "agent-builder",
                    "workflow_id": "wf_69234ded373c8190b1ece2f8f6fcec9307f7235d790729ee",
                }
            )
        )

        conversation_history.extend(
            [item.to_input_item() for item in revacast_result_temp.new_items]
        )

        revacast_result = {
            "output_text": revacast_result_temp.final_output_as(str)
        }

        return revacast_result
