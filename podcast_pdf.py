"""Persistência e extração segura de artigos PDF usados somente no podcast."""
from __future__ import annotations

import hashlib
import json
import os
import re
import uuid
from pathlib import Path

from pypdf import PdfReader

MAX_PDF_BYTES = 25 * 1024 * 1024
MAX_PDF_CHARS = 300_000


class PodcastPdfError(ValueError):
    pass


def _root(base_dir: str | Path) -> Path:
    path = Path(base_dir) / "podcast_pdfs"
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_pdf(base_dir: str | Path, data: bytes, filename: str, source: dict) -> dict:
    if not data or len(data) > MAX_PDF_BYTES:
        raise PodcastPdfError("O PDF deve ter entre 1 byte e 25 MB.")
    if not data.startswith(b"%PDF-"):
        raise PodcastPdfError("O arquivo enviado não é um PDF válido.")

    document_id = uuid.uuid4().hex
    pdf_path = _root(base_dir) / f"{document_id}.pdf"
    pdf_path.write_bytes(data)
    try:
        reader = PdfReader(str(pdf_path))
        if reader.is_encrypted:
            raise PodcastPdfError("PDF protegido por senha não pode ser analisado.")
        pages = []
        total_chars = 0
        for number, page in enumerate(reader.pages, start=1):
            text = " ".join((page.extract_text() or "").split())
            total_chars += len(text)
            if total_chars > MAX_PDF_CHARS:
                raise PodcastPdfError(
                    "O texto do PDF excede 300 mil caracteres; envie uma versão menor do artigo."
                )
            pages.append({"page": number, "text": text})
        extracted_pages = [page for page in pages if page["text"]]
        if not extracted_pages:
            raise PodcastPdfError(
                "Não foi possível extrair texto. Use um PDF pesquisável (com OCR), não apenas imagens."
            )
        metadata = {
            "id": document_id,
            "filename": Path(filename or "artigo.pdf").name,
            "sha256": hashlib.sha256(data).hexdigest(),
            "page_count": len(reader.pages),
            "extracted_page_count": len(extracted_pages),
            "character_count": total_chars,
            "source": source,
            "pages": extracted_pages,
        }
        (_root(base_dir) / f"{document_id}.json").write_text(
            json.dumps(metadata, ensure_ascii=False), encoding="utf-8"
        )
        return {key: metadata[key] for key in metadata if key != "pages"}
    except PodcastPdfError:
        pdf_path.unlink(missing_ok=True)
        raise
    except Exception as error:
        pdf_path.unlink(missing_ok=True)
        raise PodcastPdfError("Não foi possível abrir ou interpretar este PDF.") from error


def load_pdf(base_dir: str | Path, document_id: str) -> dict:
    if not re.fullmatch(r"[a-f0-9]{32}", str(document_id or "")):
        raise PodcastPdfError("Identificador de PDF inválido.")
    path = _root(base_dir) / f"{document_id}.json"
    if not path.is_file():
        raise PodcastPdfError("PDF não encontrado no armazenamento do podcast.")
    return json.loads(path.read_text(encoding="utf-8"))


def source_text(document: dict) -> str:
    return "\n\n".join(
        f"[[PÁGINA {page['page']}]]\n{page['text']}" for page in document.get("pages", [])
    )
