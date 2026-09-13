"""Episode-level podcast writing. No newsletter, retrieval, publishing or TTS side effects."""
import hashlib
import json
import os
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path

VERSION = "editorial-1"
DEFAULT_MODEL = "gpt-6-astra"
PROMPTS = Path(__file__).parent / "prompts" / "podcast"


def obj(**properties):
    return {"type": "object", "properties": properties, "required": list(properties), "additionalProperties": False}


def array(items):
    return {"type": "array", "items": items}


STR = {"type": "string"}
SUPPORT = obj(source_id=STR, quote=STR)
OBSERVATION = obj(
    kind={"type": "string", "enum": ["strength", "documented_limit", "not_reported", "design_scope"]},
    explanation=STR, supports=array(SUPPORT),
)
PLAN_SCHEMA = obj(
    central_question=STR,
    studies=array(obj(pmid=STR, question=STR, finding=STR, supports=array(SUPPORT),
                      appraisal=array(OBSERVATION), caution=STR, connection_to_next=STR)),
)
TURN = obj(speaker={"type": "string", "enum": ["HOST", "COHOST"]}, text=STR, source_ids=array(STR))
SCRIPT_SCHEMA = obj(
    title=STR, opening=array(TURN),
    studies=array(obj(pmid=STR, dialogue=array(TURN), spoken_caution=STR)),
    closing=array(TURN),
)
AUDIT_SCHEMA = obj(issues=array(obj(
    severity={"type": "string", "enum": ["blocking", "note"]}, location=STR, reason=STR,
)))


class EditorialError(ValueError):
    pass


def _clean(value):
    return " ".join(str(value or "").split())


def evidence_packet(articles, context):
    """Only source metadata/abstracts, never similarity decisions as scientific evidence."""
    if not 1 <= len(articles) <= 6:
        raise EditorialError("Selecione de um a seis estudos principais.")
    context_by_id = {str(s.get("pmid_ancora")): s for s in (context or {}).get("estudos", [])}
    packet = []
    seen = set()
    for article in articles:
        pmid = str(article.get("pmid", ""))
        if not pmid.isdigit() or pmid in seen:
            raise EditorialError("Estudo principal sem PMID válido ou duplicado.")
        seen.add(pmid)
        group = context_by_id.get(pmid, {})
        sources = []
        refs = group.get("referencias", [])
        if len(refs) > 12:
            raise EditorialError("Limite editorial: até 12 referências por estudo principal.")
        for index, item in enumerate([article, *refs]):
            original = _clean(item.get("resumo_original"))
            # A translation is not a full text; keep its provenance explicit.
            abstract = original or (_clean(item.get("resumo_traduzido")) if index == 0 else "")
            if len(abstract) > 16000:
                raise EditorialError("Resumo excede o limite editorial; não será truncado silenciosamente.")
            sources.append({
                "source_id": f"{pmid}:main" if index == 0 else f"{pmid}:ref{index}",
                "role": "main" if index == 0 else "context",
                "pmid": str(item.get("pmid") or ""), "doi": _clean(item.get("doi")),
                "title": _clean(item.get("titulo")),
                "authors": item.get("autores") or [item.get("primeiro_autor") or "Autor não informado"],
                "publication_date": _clean(item.get("data_publicacao")),
                "study_types": item.get("tipos") or [], "abstract": abstract,
                "material": "original_abstract" if original else "translated_abstract" if abstract else "metadata_only",
            })
        if not sources[0]["abstract"]:
            raise EditorialError(f"PMID {pmid} sem resumo: falta material para o roteiro.")
        packet.append({"pmid": pmid, "selection_mode": group.get("modo", "automatic"), "sources": sources})
    return packet


def _request(client, model, stage, schema, payload, usage):
    prompt = (PROMPTS / f"{stage}.txt").read_text(encoding="utf-8")
    # Keep the existing Chat Completions API. Reasoning models do not accept temperature.
    response = client.with_options(timeout=180.0, max_retries=1).chat.completions.create(
        model=model, messages=[{"role": "developer", "content": prompt},
                               {"role": "user", "content": json.dumps(payload, ensure_ascii=False)}],
        response_format={"type": "json_schema", "json_schema": {"name": f"podcast_{stage}", "strict": True, "schema": schema}},
        max_completion_tokens=12000,
    )
    choice = response.choices[0] if response.choices else None
    if not choice or choice.finish_reason != "stop" or getattr(choice.message, "refusal", None):
        raise EditorialError(f"Etapa {stage} recusada ou incompleta; nenhum texto parcial será narrado.")
    usage.append({"stage": stage, "model": getattr(response, "model", model),
                  "input_tokens": getattr(response.usage, "prompt_tokens", 0),
                  "output_tokens": getattr(response.usage, "completion_tokens", 0)})
    try:
        return json.loads(choice.message.content)
    except (TypeError, json.JSONDecodeError) as error:
        raise EditorialError(f"JSON inválido na etapa {stage}.") from error


def _ordered_studies(value, packet):
    expected = [s["pmid"] for s in packet]
    if [s.get("pmid") for s in value.get("studies", [])] != expected:
        raise EditorialError("O modelo omitiu, duplicou ou mudou a ordem dos estudos aprovados.")


def validate_plan(plan, packet):
    _ordered_studies(plan, packet)
    for study, evidence in zip(plan["studies"], packet):
        sources = {s["source_id"]: s for s in evidence["sources"]}
        def check_support(supports, required=True):
            if required and not supports:
                raise EditorialError("Afirmação sem trecho de apoio no resumo.")
            for support in supports:
                source = sources.get(support.get("source_id"))
                quote = _clean(support.get("quote"))
                if not source or not quote or quote not in source["abstract"]:
                    raise EditorialError("Trecho de apoio ausente da fonte indicada.")
        check_support(study["supports"])
        if f'{study["pmid"]}:main' not in {s["source_id"] for s in study["supports"]}:
            raise EditorialError("O achado principal não está vinculado ao estudo principal.")
        if not study["appraisal"] or not _clean(study["caution"]):
            raise EditorialError("Falta avaliação crítica ou ressalva do estudo.")
        for observation in study["appraisal"]:
            check_support(observation["supports"], observation["kind"] != "not_reported")


def validate_script(script, packet):
    _ordered_studies(script, packet)
    all_ids = {s["source_id"] for p in packet for s in p["sources"]}
    total_words = 0
    def check_turns(turns, allowed):
        nonlocal total_words
        if not turns:
            raise EditorialError("Bloco de roteiro vazio.")
        for turn in turns:
            if turn.get("speaker") not in {"HOST", "COHOST"} or not _clean(turn.get("text")):
                raise EditorialError("Fala ou locutor inválido.")
            if not set(turn["source_ids"]) <= allowed:
                raise EditorialError("Fala atribuída a fonte desconhecida ou de outro estudo.")
            total_words += len(turn["text"].split())
    check_turns(script["opening"], all_ids)
    check_turns(script["closing"], all_ids)
    for study, evidence in zip(script["studies"], packet):
        check_turns(study["dialogue"], {s["source_id"] for s in evidence["sources"]})
        if {t["speaker"] for t in study["dialogue"]} != {"HOST", "COHOST"}:
            raise EditorialError("Os dois apresentadores precisam participar de cada discussão.")
        caution = _clean(study["spoken_caution"])
        if not caution or not any(caution in _clean(t["text"]) for t in study["dialogue"]):
            raise EditorialError("A ressalva científica não aparece na fala do estudo.")
    if total_words > 350 + 420 * len(packet):
        raise EditorialError("Roteiro excede o limite de duração; revise antes do áudio.")


def dialogues(draft):
    script = draft["script"]
    blocks = [script["opening"], *(s["dialogue"] for s in script["studies"]), script["closing"]]
    return [[{"speaker": t["speaker"], "text": t["text"]} for t in block] for block in blocks]


def transcript(draft):
    names = {"HOST": "Ivo", "COHOST": "Manu"}
    return "\n\n".join(f'{names[t["speaker"]]}: {t["text"]}' for block in dialogues(draft) for t in block)


def fingerprint(draft):
    protected = {k: draft[k] for k in ("version", "model", "evidence", "plan", "script", "audit")}
    return hashlib.sha256(json.dumps(protected, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def generate_episode_steps(articles, context, client):
    packet = evidence_packet(articles, context)
    model = os.environ.get("PODCAST_SCRIPT_MODEL", DEFAULT_MODEL).strip()
    usage = []
    yield f"🧠 Planejando o episódio completo com {model} e avaliação crítica baseada nos resumos..."
    plan = _request(client, model, "plan", PLAN_SCHEMA, {"evidence": packet}, usage)
    validate_plan(plan, packet)
    yield "✍️ Escrevendo a conversa completa, incluindo abertura, ressalvas e encerramento..."
    script = _request(client, model, "write", SCRIPT_SCHEMA,
                      {"evidence": packet, "plan": plan, "target_words": 180 + 280 * len(packet)}, usage)
    validate_script(script, packet)
    yield "🔎 Conferindo afirmações, números, comparações e ressalvas contra as fontes..."
    audit = _request(client, model, "audit", AUDIT_SCHEMA, {"evidence": packet, "plan": plan, "script": script}, usage)
    draft = {"id": uuid.uuid4().hex, "version": VERSION, "model": model,
             "created_at": datetime.now(timezone.utc).isoformat(), "status": "pending_review",
             "evidence": packet, "plan": plan, "script": script, "audit": audit, "usage": usage}
    draft["sha256"] = fingerprint(draft)
    return draft


def generate_episode(articles, context, client, log=lambda message: None):
    steps = generate_episode_steps(articles, context, client)
    while True:
        try:
            log(next(steps))
        except StopIteration as completed:
            return completed.value


def draft_path(base_dir, draft_id):
    if not re.fullmatch(r"[a-f0-9]{32}", str(draft_id)):
        raise EditorialError("Identificador de roteiro inválido.")
    return Path(base_dir) / "editorial" / f"{draft_id}.json"


def save_draft(base_dir, draft):
    path = draft_path(base_dir, draft["id"])
    path.parent.mkdir(parents=True, exist_ok=True)
    # Atomic replacement, never expose partially written approval/version data.
    temp = path.with_suffix(f".{uuid.uuid4().hex}.tmp")
    try:
        temp.write_text(json.dumps(draft, ensure_ascii=False, indent=2), encoding="utf-8")
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()
    path.with_suffix(".txt").write_text(transcript(draft), encoding="utf-8")


def load_draft(base_dir, draft_id):
    draft = json.loads(draft_path(base_dir, draft_id).read_text(encoding="utf-8"))
    if draft["sha256"] != fingerprint(draft):
        raise EditorialError("Roteiro alterado após a revisão; gere uma nova versão.")
    return draft


def approve_draft(base_dir, draft_id, sha256):
    draft = load_draft(base_dir, draft_id)
    if sha256 != draft["sha256"]:
        raise EditorialError("Esta versão não é a exibida na revisão. Recarregue o roteiro.")
    if any(issue["severity"] == "blocking" for issue in draft["audit"]["issues"]):
        raise EditorialError("Há pendências científicas bloqueantes. Gere um novo roteiro antes de aprovar.")
    draft.update(status="approved", approved_at=datetime.now(timezone.utc).isoformat())
    save_draft(base_dir, draft)
    return draft


def approved_audio(base_dir, draft_id, sha256):
    draft = load_draft(base_dir, draft_id)
    if draft["status"] != "approved" or draft["sha256"] != sha256:
        raise EditorialError("Revise e aprove esta versão do roteiro antes de gerar áudio.")
    if any(i["severity"] == "blocking" for i in draft["audit"]["issues"]):
        raise EditorialError("Roteiro com pendências científicas bloqueantes.")
    return draft


def review_payload(draft):
    return {"id": draft["id"], "sha256": draft["sha256"], "status": draft["status"],
            "model": draft["model"], "title": draft["script"]["title"], "text": transcript(draft),
            "plan": draft["plan"], "audit": draft["audit"], "usage": draft["usage"],
            "evidence": draft["evidence"],
            "can_approve": not any(i["severity"] == "blocking" for i in draft["audit"]["issues"]),
            "download_url": f'/podcast-roteiro/{draft["id"]}/texto'}
