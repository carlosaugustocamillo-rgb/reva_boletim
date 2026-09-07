"""Parser e normalização de exportações BibTeX do Connected Papers.

O Connected Papers exporta o artigo de origem como a primeira entrada e os
artigos relacionados nas entradas seguintes. O arquivo não traz um score
numérico confiável; por isso a ordem é preservada, mas nunca é tratada como
evidência de qualidade ou concordância clínica.
"""

from __future__ import annotations

import html
import re
from typing import Any


_FIELD_RE = re.compile(r"([A-Za-z][A-Za-z0-9_-]*)\s*=\s*", re.MULTILINE)
_NULL_VALUES = {"", "null", "none", "undefined", "not defined within the journal database."}


def _clean_value(value: str) -> str:
    value = html.unescape(str(value or "")).strip()
    if value.casefold() in _NULL_VALUES:
        return ""
    return value


def _read_value(text: str, start: int) -> tuple[str, int]:
    """Lê valor BibTeX com chaves/aspas aninhadas e retorna fim exclusivo."""
    while start < len(text) and text[start].isspace():
        start += 1
    if start >= len(text):
        return "", start
    opener = text[start]
    if opener == "{":
        depth = 1
        idx = start + 1
        chunks: list[str] = []
        while idx < len(text) and depth:
            char = text[idx]
            if char == "\\" and idx + 1 < len(text):
                chunks.extend((char, text[idx + 1]))
                idx += 2
                continue
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    return "".join(chunks), idx + 1
            chunks.append(char)
            idx += 1
        raise ValueError("Valor BibTeX com chaves não fechadas")
    if opener == '"':
        idx = start + 1
        chunks: list[str] = []
        while idx < len(text):
            char = text[idx]
            if char == "\\" and idx + 1 < len(text):
                chunks.extend((char, text[idx + 1]))
                idx += 2
                continue
            if char == '"':
                return "".join(chunks), idx + 1
            chunks.append(char)
            idx += 1
        raise ValueError("Valor BibTeX com aspas não fechadas")

    idx = start
    while idx < len(text) and text[idx] not in ",\n}":
        idx += 1
    return text[start:idx].strip(), idx


def _parse_entry(entry: str) -> dict[str, str]:
    match = re.match(r"@([A-Za-z]+)\s*\{\s*([^,]+),", entry, re.DOTALL)
    if not match:
        raise ValueError("Entrada BibTeX inválida")
    kind = match.group(1).casefold()
    source_id = match.group(2).strip()
    fields: dict[str, str] = {"tipo": kind, "source_id": source_id}
    body_start = match.end()
    for field in _FIELD_RE.finditer(entry, body_start):
        if field.start() < body_start:
            continue
        name = field.group(1).casefold()
        value, end = _read_value(entry, field.end())
        fields[name] = _clean_value(value)
        if end <= field.end():
            break
    return fields


def _entries(text: str) -> list[str]:
    entries: list[str] = []
    for match in re.finditer(r"@[A-Za-z]+\s*\{", text):
        start = match.start()
        depth = 0
        quoted = False
        escaped = False
        for idx in range(match.end() - 1, len(text)):
            char = text[idx]
            if escaped:
                escaped = False
                continue
            if char == "\\":
                escaped = True
                continue
            if char == '"':
                quoted = not quoted
                continue
            if quoted:
                continue
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    entries.append(text[start:idx + 1])
                    break
        else:
            raise ValueError("Entrada BibTeX sem fechamento")
    return entries


def _authors(value: str) -> list[str]:
    return [part.strip() for part in re.split(r"\s+and\s+", value or "", flags=re.IGNORECASE) if part.strip()]


def _normalize_entry(fields: dict[str, str], position: int) -> dict[str, Any]:
    pmid = re.sub(r"\D", "", fields.get("pmid", ""))
    doi = fields.get("doi", "").strip()
    if doi.casefold() in _NULL_VALUES:
        doi = ""
    url = fields.get("url", "").strip()
    year_match = re.search(r"\b(?:19|20)\d{2}\b", fields.get("year", ""))
    year = year_match.group(0) if year_match else ""
    return {
        "source_id": fields.get("source_id", ""),
        "posicao_connected_papers": position,
        "titulo": fields.get("title", "").strip(),
        "ano": year,
        "data_publicacao": year,
        "url_semantic_scholar": url,
        "url": url,
        "resumo_original": fields.get("abstract", "").strip(),
        "autores": _authors(fields.get("author", "")),
        "journal": fields.get("journal", "").strip(),
        "volume": fields.get("volume", "").strip(),
        "pages": fields.get("pages", "").strip(),
        "doi": doi,
        "pmid": pmid,
        "fonte": "Connected Papers",
    }


def parse_bibtex(text: str) -> dict[str, Any]:
    """Converte uma exportação BibTeX em âncora e candidatos normalizados."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Arquivo BibTeX vazio")
    parsed = [_parse_entry(entry) for entry in _entries(text)]
    if not parsed:
        raise ValueError("Nenhuma entrada BibTeX encontrada")
    articles = [_normalize_entry(item, index) for index, item in enumerate(parsed)]
    anchor = articles[0]
    candidates: list[dict[str, Any]] = []
    seen: set[str] = set()
    for item in articles[1:]:
        key = item.get("pmid") or item.get("doi") or item.get("source_id") or item.get("titulo")
        if not key or key in seen:
            continue
        seen.add(key)
        candidates.append(item)
    return {
        "versao": 1,
        "fonte": "Connected Papers",
        "arquivo_tipo": "bibtex",
        "ancora": anchor,
        "candidatos": candidates,
        "total_entradas": len(articles),
    }


def article_identity(article: dict[str, Any]) -> str:
    return str(article.get("pmid") or article.get("doi") or article.get("source_id") or article.get("titulo") or "").strip()


def build_manual_report(
    main_articles: list[dict[str, Any]],
    selected_by_anchor: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    """Monta o mesmo contrato do relatório PubMed para referências manuais."""
    studies = []
    for anchor in main_articles:
        anchor_id = article_identity(anchor)
        references = []
        for reference in selected_by_anchor.get(anchor_id, []):
            item = dict(reference)
            item["triagem"] = {
                "incluir": True,
                "populacao": "revisao_manual",
                "intervencao": "revisao_manual",
                "desfecho": "revisao_manual",
                "comparador": "nao_informado",
                "desenho": "nao_informado",
                "justificativa": item.get("justificativa_manual", "Selecionada manualmente no Connected Papers."),
                "trecho_ancora": "",
                "trecho_candidato": "",
                "diferenca_importante": item.get("limitacao_manual", "Revisão editorial manual necessária."),
                "origem": "curadoria_manual",
            }
            references.append(item)
        studies.append({
            "pmid_ancora": str(anchor.get("pmid", "")),
            "titulo_ancora": anchor.get("titulo", ""),
            "referencias": references,
            "status": "ok" if references else "sem_referencias",
            "modo": "curadoria_manual",
        })
    return {
        "versao": 1,
        "fonte": "Connected Papers (curadoria manual)",
        "habilitado": bool(selected_by_anchor),
        "modo": "manual",
        "status": "concluido",
        "artigos_principais": main_articles,
        "estudos": studies,
        "tokens_triagem": {"prompt_tokens": 0, "completion_tokens": 0},
    }
