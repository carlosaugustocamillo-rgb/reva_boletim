"""Utilitarios pequenos para transformar erros da ElevenLabs em diagnosticos uteis."""

from __future__ import annotations

from typing import Any


def _error_detail(error: Any) -> dict:
    body = getattr(error, "body", None)
    if isinstance(body, dict):
        detail = body.get("detail", body)
        return detail if isinstance(detail, dict) else {}
    return {}


def diagnosticar_erro_elevenlabs(error: Any) -> dict:
    """Retorna um erro seguro e estruturado, sem expor chave ou traceback."""
    detail = _error_detail(error)
    raw = " ".join(
        str(value)
        for value in (
            detail.get("code"),
            detail.get("status"),
            detail.get("message"),
            str(error),
        )
        if value
    ).lower()

    request_id = detail.get("request_id")
    status_code = getattr(error, "status_code", None)

    if any(token in raw for token in (
        "subscription_required",
        "paid_plan_required",
        "ivc_not_permitted",
        "only_for_creator",
        "professional voices require",
        "free users cannot use library voices",
    )):
        code = "elevenlabs_plan_required"
        message = (
            "O plano atual da ElevenLabs nao permite uma ou mais vozes configuradas. "
            "Ative o plano Creator (ou superior) ou troque por vozes permitidas no plano atual."
        )
    elif any(token in raw for token in ("quota_exceeded", "insufficient_credits", "usage_limit")):
        code = "elevenlabs_quota_exceeded"
        message = "A cota/creditos da ElevenLabs acabou. Verifique o uso e o limite da assinatura."
    elif any(token in raw for token in ("invalid_api_key", "unauthorized", "authentication_error")):
        code = "elevenlabs_authentication_failed"
        message = "A chave da ElevenLabs e invalida ou nao possui a permissao necessaria."
    elif any(token in raw for token in ("voice_not_found", "voice does not exist")):
        code = "elevenlabs_voice_unavailable"
        message = "Uma das vozes configuradas nao existe ou nao esta acessivel para esta chave."
    else:
        code = "elevenlabs_generation_failed"
        original = detail.get("message")
        message = f"A ElevenLabs recusou a geracao: {original}" if original else "Falha ao gerar audio na ElevenLabs."

    return {
        "code": code,
        "message": message,
        "provider_status": status_code,
        "request_id": request_id,
    }


def formatar_erro_elevenlabs(error: Any) -> str:
    diagnostico = diagnosticar_erro_elevenlabs(error)
    request_suffix = f" (request_id: {diagnostico['request_id']})" if diagnostico["request_id"] else ""
    return f"{diagnostico['message']}{request_suffix}"
