"""Adapt Reva+ Consensus extraction to podcast sources without inventing abstracts."""
import hashlib


def build_consensus_import(report, parsed):
    report = str(report or "").strip()
    if not 80 <= len(report) <= 300_000:
        raise ValueError("Forneça um relatório Consensus entre 80 e 300 mil caracteres.")
    digest = hashlib.sha256(report.encode()).hexdigest()[:16]
    # Keep the complete report available even when reference extraction is empty.
    candidates = [{
        "source_id": f"consensus-report:{digest}",
        "titulo": "Relatório Consensus — síntese complementar",
        "fonte": "Consensus", "material": "consensus_report",
        "evidence_content": report,
    }]
    for index, ref in enumerate(parsed.get("referencias") or []):
        candidates.append({
            **ref,
            "source_id": f"consensus:{digest}:{index}",
            "titulo": ref.get("titulo") or ref.get("texto") or f"Referência {index + 1}",
            "resumo_original": "",
            "material": ref.get("material") or "consensus_extract",
            "evidence_content": ref.get("evidence_content") or ref.get("consensus_claim") or "",
            "fonte": "Consensus",
        })
    return {"fonte": "Consensus", "candidatos": candidates,
            "parser_warning": parsed.get("erro") or "", "total_candidatos": len(candidates)}
