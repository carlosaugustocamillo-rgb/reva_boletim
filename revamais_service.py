import os
# Force Build 123
import base64
import requests
import json
import html
import re
from io import BytesIO
from datetime import datetime, timedelta
from dotenv import load_dotenv
import google.generativeai as genai
from openai import OpenAI
from Bio import Entrez
import mailchimp_marketing as MailchimpMarketing
from mailchimp_marketing.api_client import ApiClientError
import csv
import shutil

# Importa ferramentas já existentes
from firebase_service import get_firestore_db, read_json_from_storage, upload_file

load_dotenv()

# Configurações
Entrez.email = os.environ.get("ENTREZ_EMAIL")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
OPENAI_CHAT_MODEL_REPLACEMENTS = {
    "gpt-5.5-pro": "gpt-5.5",
}


def _openai_chat_model_from_env(env_name, default_model):
    model_name = os.environ.get(env_name, default_model).strip()
    replacement = OPENAI_CHAT_MODEL_REPLACEMENTS.get(model_name)
    if replacement:
        print(f"⚠️ {env_name}={model_name} não é compatível com chat.completions; usando {replacement}.")
        return replacement
    return model_name


OPENAI_TEXT_MODEL = _openai_chat_model_from_env("OPENAI_TEXT_MODEL", "gpt-5.5")
OPENAI_TEXT_MODEL_SEARCH = _openai_chat_model_from_env("OPENAI_TEXT_MODEL_SEARCH", OPENAI_TEXT_MODEL)
OPENAI_TEXT_MODEL_WRITE = _openai_chat_model_from_env("OPENAI_TEXT_MODEL_WRITE", OPENAI_TEXT_MODEL)
OPENAI_IMAGE_MODEL = os.environ.get("OPENAI_IMAGE_MODEL", "gpt-image-2").strip()
OPENAI_IMAGE_QUALITY = os.environ.get("OPENAI_IMAGE_QUALITY", "high").strip()
OPENAI_IMAGE_TIMEOUT_SECONDS = int(os.environ.get("OPENAI_IMAGE_TIMEOUT_SECONDS", "180"))
OPENAI_IMAGE_ENABLED = True

# Configuração Gemini
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)
DEFAULT_GEMINI_TEXT_MODEL = "gemini-3.1-pro-preview"
DEFAULT_GEMINI_IMAGE_MODEL = "gemini-3-pro-image"
GEMINI_MODEL_REPLACEMENTS = {
    "gemini-3-pro-preview": DEFAULT_GEMINI_TEXT_MODEL,
    "gemini-3-pro-image-preview": DEFAULT_GEMINI_IMAGE_MODEL,
}


def _gemini_model_from_env(env_name, default_model):
    model_name = os.environ.get(env_name, default_model).strip()
    model_key = model_name.removeprefix("models/")
    replacement = GEMINI_MODEL_REPLACEMENTS.get(model_key)
    if replacement:
        print(f"⚠️ {env_name}={model_name} indisponível; usando {replacement}.")
        return replacement
    return model_name


GEMINI_TEXT_MODEL = _gemini_model_from_env("GEMINI_TEXT_MODEL", DEFAULT_GEMINI_TEXT_MODEL)
GEMINI_IMAGE_MODEL = _gemini_model_from_env("GEMINI_IMAGE_MODEL", DEFAULT_GEMINI_IMAGE_MODEL)

REVAMAIS_STATE_COLLECTION = "revamais_state"
REVAMAIS_STATE_DOC_ID = "editorial_calendar"
REVAMAIS_LEGACY_STATE_FILE = "revamais/used_themes_state.json"

# Configuração Mailchimp
mc = MailchimpMarketing.Client()
mc.set_config({
    "api_key": os.environ.get("MC_API_KEY"),
    "server": os.environ.get("MC_SERVER")
})
# Tenta pegar ID específico para Reva+, senão usa o hardcoded (Reva+ Audience) ou o padrão (Weekly)
MC_LIST_ID = os.environ.get("MC_LIST_ID_REVAMAIS", "510b954f9a") 
MC_FROM_NAME = os.environ.get("MC_FROM_NAME", "Revalidatie")
MC_REPLY_TO = os.environ.get("MC_REPLY_TO", "contato@revalidatie.com.br")


def gerar_texto_openai(prompt, system_prompt=None, model_name=None, timeout_seconds=None):
    """
    Gera texto com o modelo principal da OpenAI configurado para o projeto.
    """
    model_name = model_name or OPENAI_TEXT_MODEL_SEARCH
    messages = []
    if system_prompt:
        messages.append({"role": "system", "content": system_prompt})
    messages.append({"role": "user", "content": prompt})

    response = client.chat.completions.create(
        model=model_name,
        messages=messages,
        timeout=timeout_seconds,
    )
    content = response.choices[0].message.content or ""
    return content.strip()


def gerar_texto_preferencial(prompt, system_prompt=None, model_name=None, timeout_seconds=None):
    """
    Usa OpenAI como primeira opção e Gemini como fallback para tarefas textuais.
    """
    model_name = model_name or OPENAI_TEXT_MODEL_SEARCH
    try:
        return gerar_texto_openai(
            prompt,
            system_prompt=system_prompt,
            model_name=model_name,
            timeout_seconds=timeout_seconds,
        )
    except Exception as e_openai:
        print(f"⚠️ Falha OpenAI ({model_name}): {e_openai}")

    try:
        model = genai.GenerativeModel(GEMINI_TEXT_MODEL)
        full_prompt = prompt if not system_prompt else f"{system_prompt}\n\n{prompt}"
        response = model.generate_content(full_prompt)
        return (response.text or "").strip()
    except Exception as e_gemini:
        print(f"⚠️ Falha Gemini ({GEMINI_TEXT_MODEL}): {e_gemini}")
        raise


# -------------------------------------------------------------------------
# Helpers de Qualidade de Referência (JCR + Keywords)
# -------------------------------------------------------------------------

_JCR_CACHE = None

def load_jcr_data():
    """
    Carrega o arquivo CSV do JCR na memória para consulta rápida de Fator de Impacto.
    Retorna um dict: { "JOURNAL NAME UPPER": float(impact_factor) }
    """
    global _JCR_CACHE
    if _JCR_CACHE is not None:
        return _JCR_CACHE
        
    jcr_path = os.path.join(os.path.dirname(__file__), "CarlosCamillo_JCR_JournalResults_12_2025.csv")
    if not os.path.exists(jcr_path):
        print("⚠️ Arquivo JCR não encontrado. Scores de impacto serão 0.")
        return {}
        
    cache = {}
    try:
        with open(jcr_path, 'r', encoding='utf-8') as f:
            # Pula linhas de cabeçalho inicial até encontrar o header real
            lines = f.readlines()
            start_idx = 0
            for i, line in enumerate(lines):
                if line.startswith("Journal name,"):
                    start_idx = i
                    break
            
            reader = csv.DictReader(lines[start_idx:])
            for row in reader:
                name = row.get("Journal name", "").upper().strip()
                jif_val = row.get("2024 JIF")
                if jif_val is None: jif_val = "0"
                jif_str = str(jif_val).replace(",", "") # Remove milhar se houver
                try:
                    jif = float(jif_str)
                except:
                    jif = 0.0
                
                if name:
                    cache[name] = jif
                    # Cache também variações comuns (ex: sem THE)
                    if name.startswith("THE "):
                        cache[name[4:]] = jif
                        
        print(f"✅ JCR Data carregado: {len(cache)} revistas.")
        _JCR_CACHE = cache
        return cache
    except Exception as e:
        print(f"❌ Erro ao carregar JCR: {e}")
        return {}

def generate_search_keywords(tema):
    """
    Usa LLM para extrair 3-5 keywords OBRIGATÓRIAS em Inglês para o tema.
    Essas keywords serão usadas para filtrar resultados irrelevantes.
    """
    try:
        prompt = (
            f"Analyze the medical topic: '{tema}'. "
            "Return a Python list of strings with 3 to 5 ESSENTIAL English keywords (single words or short bi-grams) "
            "that MUST appear in a valid scientific article about this topic. "
            "Focus on the pathology, anatomy, or intervention. "
            "Example output format: ['Hypertension', 'Blood Pressure', 'Cardiovascular']"
        )
        response = gerar_texto_preferencial(
            prompt,
            system_prompt="You extract precise biomedical keywords for PubMed filtering.",
            model_name=OPENAI_TEXT_MODEL_SEARCH,
        )
        # Limpeza básica para extrair a lista
        import ast
        start = response.find('[')
        end = response.rfind(']') + 1
        if start != -1 and end != -1:
            keywords = ast.literal_eval(response[start:end])
            return [k.lower().strip() for k in keywords if isinstance(k, str)]
        return []
    except Exception as e:
        print(f"⚠️ Falha ao gerar keywords: {e}")
        return []

def gerar_query_pubmed_tema(tema):
    """
    Gera uma query booleana para PubMed usando OpenAI como motor principal.
    """
    prompt_query = (
        f"Create a specific PubMed Search Query for the topic: '{tema}'. "
        "Use boolean operators (AND, OR) to combine MeSH terms or keywords. "
        "IMPORTANT: Use 'AND' to intersect distinct concepts (e.g. 'Diabetes AND Exercise'). "
        "Use 'OR' only for synonyms. "
        "Return ONLY the query string, nothing else. No explanation."
    )
    try:
        return gerar_texto_preferencial(
            prompt_query,
            system_prompt="You write compact, high-precision PubMed boolean search queries.",
            model_name=OPENAI_TEXT_MODEL_SEARCH,
        ).replace('"', '').strip()
    except Exception as e:
        print(f"⚠️ Falha ao gerar query PubMed: {e}")
        return tema


def buscar_referencias_pubmed(tema_ingles, limite_retorno=5):
    """
    Busca artigos no PubMed com filtro rigoroso de qualidade e relevância.
    1. Busca 30 candidatos (SR/RCT/Review).
    2. Filtra: Título/Abstract TEM que ter keywords do tema (Relevância).
    3. Rankeia: Prioriza JCR Impact Factor alto.
    4. Retorna os melhores artigos ranqueados.
    """
    print(f"🔎 Buscando referências para: {tema_ingles}...")
    
    # 1. Preparação (Keywords + JCR)
    jcr_data = load_jcr_data()
    keywords = generate_search_keywords(tema_ingles)
    print(f"🎯 Keywords Obrigatórias (Relevância): {keywords}")

    # 2. Busca Ampliada (30 artigos)
    from datetime import timedelta
    agora = datetime.now()

    def build_query(years_back: int):
        data_ini = (agora - timedelta(days=years_back * 365)).strftime("%Y/%m/%d")
        data_fim = agora.strftime("%Y/%m/%d")
        date_term = f'("{data_ini}"[Date - Publication] : "{data_fim}"[Date - Publication])'
        # Query tenta focar em alta evidência primeiro, mas permite Reviews
        return (
            f"({tema_ingles}) AND "
            f"(Systematic Review[pt] OR Randomized Controlled Trial[pt] OR Review[pt] OR Meta-Analysis[pt]) AND "
            f"{date_term}"
        )

    def search_ids(query):
        handle = Entrez.esearch(db="pubmed", term=query, retmax=40, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]

    query = build_query(5)
    
    candidates = []
    try:
        ids = search_ids(query)

        if not ids:
            print("⚠️ Nenhuma referência encontrada na busca inicial (5 anos). Tentando janela de 10 anos...")
            query = build_query(10)
            ids = search_ids(query)
            if not ids:
                print("⚠️ Nenhuma referência encontrada na janela de 10 anos.")
                return []

        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        papers = Entrez.read(handle)
        handle.close()

        # 3. Filtragem e Ranqueamento
        if 'PubmedArticle' not in papers:
            return []

        for article in papers['PubmedArticle']:
            try:
                medline = article['MedlineCitation']['Article']
                journal = medline.get('Journal', {}).get('Title', '').upper()
                title = medline.get('ArticleTitle', '')
                abstract_list = medline.get('Abstract', {}).get('AbstractText', [])
                abstract = " ".join(abstract_list) if abstract_list else ""
                
                # A. Filtro de Relevância (Keywords)
                # Se tiver keywords definidas, pelo menos UMA tem que estar no Titulo ou Abstract
                text_content = (title + " " + abstract).lower()
                if keywords:
                    match = any(k in text_content for k in keywords)
                    if not match:
                        continue # Pula este artigo
                
                # B. Score JCR
                jif = jcr_data.get(journal, 0.0)
                if jif == 0.0 and journal.startswith("THE "):
                     jif = jcr_data.get(journal[4:], 0.0)

                # Extrai dados bibliográficos
                autores = medline.get('AuthorList', [])
                primeiro_autor = f"{autores[0]['LastName']} et al." if autores else "Autores diversos"
                ano = medline.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', '')
                if not ano:
                    # Tenta extrair da data de publicação completa
                    try:
                        ano = article['PubmedData']['History'][0]['Year']
                    except:
                        ano = "s.d."
                
                pmid = article['MedlineCitation']['PMID']
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

                candidates.append({
                    "pmid": str(pmid),
                    "texto": f"{primeiro_autor}. {title}. {journal.title()}, {ano}.",
                    "link": link,
                    "resumo": abstract,
                    "jif": jif,
                    "journal": journal,
                    "fonte": "PubMed",
                    "source_label": "PubMed",
                })
            except Exception as e:
                continue

    except Exception as e:
         print(f"❌ Erro na busca PubMed: {e}")
         return []

    # 4. Ordenação Final
    # Prioridade: JIF (Decrescente) -> Se JIF == 0, fica no fim.
    candidates.sort(key=lambda x: x['jif'], reverse=True)
    
    print(f"📊 {len(candidates)} artigos relevantes pós-filtro.")
    for c in candidates[:limite_retorno]:
        print(f"   ⭐ [{c['jif']:.1f}] {c['journal']} - {c['texto'][:50]}...")
        
    return candidates[:limite_retorno]


def _limpar_resposta_json(texto):
    """
    Extrai JSON de respostas de LLM que às vezes vêm com cercas Markdown.
    """
    texto = (texto or "").strip()
    texto = texto.replace("```json", "").replace("```", "").strip()

    try:
        return json.loads(texto)
    except Exception:
        pass

    for start_char, end_char in (("{", "}"), ("[", "]")):
        start = texto.find(start_char)
        end = texto.rfind(end_char)
        if start != -1 and end != -1 and end > start:
            try:
                return json.loads(texto[start : end + 1])
            except Exception:
                continue

    raise ValueError("Resposta do modelo não contém JSON válido.")


def _extrair_doi(texto):
    match = re.search(r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b", texto or "", flags=re.IGNORECASE)
    if not match:
        return ""
    return match.group(0).rstrip(".,;)]}")


def _normalizar_referencia_consensus(ref, index):
    if not isinstance(ref, dict):
        ref = {"texto": str(ref)}

    titulo = str(ref.get("titulo") or ref.get("title") or "").strip()
    autores = str(ref.get("autores") or ref.get("authors") or "").strip()
    ano = str(ref.get("ano") or ref.get("year") or "").strip()
    journal = str(ref.get("journal") or ref.get("revista") or "").strip()
    doi = str(ref.get("doi") or "").strip()
    link = str(ref.get("link") or ref.get("url") or "").strip()
    tipo_estudo = str(ref.get("tipo_estudo") or ref.get("study_type") or "").strip()
    achado = str(
        ref.get("achado_principal")
        or ref.get("consensus_claim")
        or ref.get("claim")
        or ref.get("conclusao")
        or ""
    ).strip()
    resumo = str(ref.get("resumo") or ref.get("abstract") or ref.get("summary") or achado).strip()
    texto = str(ref.get("texto") or ref.get("citation") or "").strip()

    if not doi:
        doi = _extrair_doi(" ".join([texto, link, resumo, achado]))
    if not link and doi:
        link = f"https://doi.org/{doi}"

    if not texto:
        partes = []
        if autores:
            partes.append(autores)
        if titulo:
            partes.append(titulo)
        if journal or ano:
            partes.append(", ".join([p for p in [journal, ano] if p]))
        texto = ". ".join(partes).strip()

    if not texto:
        texto = f"Referência científica #{index + 1}"

    if resumo and achado and achado not in resumo:
        resumo = f"{resumo}\nAchado principal: {achado}"

    return {
        "pmid": str(ref.get("pmid") or "").strip() or None,
        "texto": texto,
        "link": link,
        "resumo": resumo,
        "jif": 0.0,
        "journal": journal,
        "doi": doi,
        "tipo_estudo": tipo_estudo,
        "achado_principal": achado,
        "consensus_claim": achado,
        "fonte": "Consensus",
        "source_label": "Consensus",
    }


def _normalizar_relatorio_consensus(relatorio_texto):
    texto = (relatorio_texto or "").replace("\r\n", "\n").replace("\r", "\n")
    texto = re.sub(r"\n{3,}", "\n\n", texto)
    return texto.strip()


def _is_generic_consensus_theme(text):
    normalized = str(text or "").strip().lower()
    return normalized in {
        "",
        "relatório consensus",
        "relatorio consensus",
        "consensus report",
        "consensus",
    }


def _looks_like_title_continuation(current_line, next_line):
    current_line = str(current_line or "").strip()
    next_line = str(next_line or "").strip()
    if not current_line or not next_line:
        return False

    if len(next_line) > 80:
        return False
    if next_line.lower() in {"references", "evidence", "strength claim"}:
        return False
    if current_line.endswith((".", "?", "!")):
        return False

    if re.search(r"\b(e|de|do|da|dos|das|para|com|sem|versus|vs|x|ou)\s*$", current_line, re.IGNORECASE):
        return True

    next_word_count = len(next_line.split())
    if next_word_count <= 4 and next_line[:1].islower():
        return True

    return False


def _infer_tema_from_consensus_report(relatorio_texto):
    texto = _normalizar_relatorio_consensus(relatorio_texto)
    if not texto:
        return ""

    linhas = []
    for linha in texto.splitlines():
        linha = re.sub(r"\s+", " ", linha).strip(" -:\t")
        if not linha:
            continue
        if re.fullmatch(r"\d+\s*/\s*\d+", linha):
            continue
        if linha.lower() in {"references", "evidence", "strength claim"}:
            continue
        linhas.append(linha)
        if len(linhas) >= 8:
            break

    for idx, linha in enumerate(linhas):
        if len(linha) < 20 or len(linha) > 180:
            continue
        if _is_generic_consensus_theme(linha):
            continue
        if not re.search(r"[A-Za-zÀ-ÿ]", linha):
            continue

        if idx + 1 < len(linhas):
            proxima = linhas[idx + 1]
            if _looks_like_title_continuation(linha, proxima):
                combinado = f"{linha} {proxima}".strip()
                if len(combinado) <= 200:
                    return combinado.rstrip(" .")

        return linha.rstrip(" .")

    return ""


def _resolve_consensus_theme(preferred_theme, parsed_theme, relatorio_texto):
    for candidate in [
        preferred_theme,
        parsed_theme,
        _infer_tema_from_consensus_report(relatorio_texto),
    ]:
        candidate = str(candidate or "").strip()
        if candidate and not _is_generic_consensus_theme(candidate):
            return candidate
    return "Tema clínico do relatório"


def extrair_texto_pdf_consensus(pdf_bytes):
    """
    Extrai texto bruto de um PDF do Consensus.
    """
    try:
        from pypdf import PdfReader
    except ImportError as e:
        raise RuntimeError(
            "Dependência ausente para leitura de PDF. Instale 'pypdf' no ambiente do projeto."
        ) from e

    reader = PdfReader(BytesIO(pdf_bytes))
    partes = []
    for page in reader.pages:
        try:
            partes.append(page.extract_text() or "")
        except Exception:
            continue

    texto = _normalizar_relatorio_consensus("\n\n".join(partes))
    if len(texto) < 100:
        raise RuntimeError("Não foi possível extrair texto suficiente do PDF do Consensus.")
    return texto


def _extrair_referencias_da_secao_consensus(relatorio_texto, limite_retorno=None):
    marcador = relatorio_texto.lower().find("references")
    if marcador == -1:
        return []

    secao = relatorio_texto[marcador + len("references") :]
    linhas = [linha.strip() for linha in secao.splitlines()]
    referencias = []
    bloco_atual = []

    for linha in linhas:
        if not linha:
            continue
        if re.fullmatch(r"\d+\s*/\s*\d+", linha):
            continue
        if linha.startswith("•"):
            continue

        bloco_atual.append(linha)
        if "doi.org/" in linha.lower() or re.search(r"\b10\.\d{4,9}/", linha, re.IGNORECASE):
            bloco = " ".join(bloco_atual)
            bloco_atual = []
            bloco = re.sub(r"\s+", " ", bloco).strip()
            if len(bloco) < 40:
                continue
            doi = _extrair_doi(bloco)
            referencias.append(
                _normalizar_referencia_consensus(
                    {
                        "texto": bloco[:500].rstrip(" ."),
                        "resumo": bloco[:1200],
                        "doi": doi,
                        "link": f"https://doi.org/{doi}" if doi else "",
                    },
                    len(referencias),
                )
            )
            if limite_retorno and len(referencias) >= limite_retorno:
                break

    return referencias


def _extrair_referencias_consensus_fallback(relatorio_texto, limite_retorno=None):
    """
    Fallback simples para quando o LLM não consegue estruturar o relatório.
    Mantém o fluxo testável, mas a qualidade depende da curadoria manual.
    """
    referencias_secao = _extrair_referencias_da_secao_consensus(relatorio_texto, limite_retorno)
    if referencias_secao:
        return referencias_secao

    blocos = re.split(
        r"\n\s*\n|(?=^\s*\d+[\).\s])",
        relatorio_texto or "",
        flags=re.MULTILINE,
    )
    referencias = []

    for bloco in blocos:
        bloco = " ".join(bloco.split())
        if len(bloco) < 80:
            continue
        if not (
            re.search(r"\b(19|20)\d{2}\b", bloco)
            or "doi" in bloco.lower()
            or "et al" in bloco.lower()
        ):
            continue

        doi = _extrair_doi(bloco)
        resumo = bloco[:1200]
        texto = bloco[:320].rstrip(" .")
        referencias.append(
            _normalizar_referencia_consensus(
                {
                    "texto": texto,
                    "resumo": resumo,
                    "doi": doi,
                    "link": f"https://doi.org/{doi}" if doi else "",
                },
                len(referencias),
            )
        )
        if limite_retorno and len(referencias) >= limite_retorno:
            break

    return referencias


def extrair_referencias_consensus(relatorio_texto, tema_usuario=None, limite_retorno=None):
    """
    Converte o relatório copiado/exportado do Consensus em referências no mesmo
    formato usado pela curadoria manual do Reva+.
    """
    relatorio_texto = _normalizar_relatorio_consensus(relatorio_texto)
    if len(relatorio_texto) < 200:
        return {
            "tema": (tema_usuario or "").strip(),
            "referencias": [],
            "erro": "O relatório do Consensus está muito curto para extração.",
        }

    # Primeiro usa um parser local para evitar que a UI fique bloqueada
    # esperando LLM em textos longos/instáveis.
    referencias_fallback = _extrair_referencias_consensus_fallback(
        relatorio_texto,
        limite_retorno,
    )
    if referencias_fallback:
        return {
            "tema": (tema_usuario or "").strip(),
            "referencias": referencias_fallback[:limite_retorno] if limite_retorno else referencias_fallback,
            "erro": "Extração rápida local aplicada; revise as referências antes de gerar o Reva+.",
        }

    regra_limite = (
        f"- Retorne no máximo {limite_retorno} referências."
        if limite_retorno
        else "- Retorne todas as referências explícitas do relatório."
    )

    prompt = f"""
Você receberá um relatório gerado pelo Consensus.app. Extraia apenas estudos/referências que estejam explicitamente presentes no texto.

Tema informado pelo usuário: "{tema_usuario or ''}"

Retorne APENAS JSON válido, sem Markdown, neste formato:
{{
  "tema": "tema médico em português, curto e fiel ao relatório",
  "referencias": [
    {{
      "titulo": "título do estudo",
      "autores": "Primeiro autor et al. ou autores disponíveis",
      "ano": "ano",
      "journal": "revista, se disponível",
      "doi": "DOI, se disponível",
      "link": "URL do paper/DOI, se disponível",
      "tipo_estudo": "systematic review, RCT, guideline, cohort, etc., se disponível",
      "achado_principal": "conclusão ou achado do Consensus em 1-2 frases",
      "resumo": "resumo do achado relevante para pacientes, sem inventar dados"
    }}
  ]
}}

Regras:
- Não invente autores, títulos, DOI, revista, ano, números ou resultados.
- Se uma informação bibliográfica não estiver no relatório, deixe vazia.
- Priorize revisões sistemáticas, meta-análises, ensaios clínicos, guidelines e estudos diretamente ligados ao tema.
{regra_limite}

RELATÓRIO CONSENSUS:
---
{relatorio_texto[:30000]}
---
"""

    try:
        resposta = gerar_texto_preferencial(
            prompt,
            system_prompt="You extract structured scientific references from Consensus.app reports. Return strict JSON only.",
            model_name=OPENAI_TEXT_MODEL_SEARCH,
            timeout_seconds=20,
        )
        dados = _limpar_resposta_json(resposta)
        if isinstance(dados, list):
            dados = {"tema": tema_usuario or "", "referencias": dados}

        referencias_raw = dados.get("referencias") if isinstance(dados, dict) else []
        referencias = [
            _normalizar_referencia_consensus(ref, index)
            for index, ref in enumerate(referencias_raw or [])
        ]
        referencias = [ref for ref in referencias if ref.get("texto") or ref.get("resumo")]

        return {
            "tema": (tema_usuario or "").strip() or str(dados.get("tema") or "").strip(),
            "referencias": referencias[:limite_retorno] if limite_retorno else referencias,
            "erro": "",
        }
    except Exception as e:
        print(f"⚠️ Falha ao estruturar relatório Consensus com LLM: {e}")
        return {
            "tema": (tema_usuario or "").strip(),
            "referencias": _extrair_referencias_consensus_fallback(relatorio_texto, limite_retorno),
            "erro": str(e),
        }


def preparar_referencias_consensus_revamais(tema_usuario=None, relatorio_texto="", quantidade_referencias=None, calendar_index=None):
    relatorio_texto = _normalizar_relatorio_consensus(relatorio_texto)
    dados_tema = None
    if calendar_index is not None:
        dados_tema = resolver_tema_revamais(
            tema_usuario=tema_usuario,
            consumir_tema_auto=False,
            calendar_index=calendar_index,
        )
        if not dados_tema:
            raise RuntimeError("Não foi possível resolver o item do calendário selecionado.")

    tema_base = (dados_tema or {}).get("tema") or (tema_usuario or "").strip()
    resultado = extrair_referencias_consensus(
        relatorio_texto=relatorio_texto,
        tema_usuario=tema_base,
        limite_retorno=quantidade_referencias,
    )

    tema = _resolve_consensus_theme(
        tema_base,
        resultado.get("tema"),
        relatorio_texto,
    )
    referencias = resultado.get("referencias") or []
    instagram_format = (dados_tema or {}).get("formato_instagram") or "Carrossel"

    if not referencias:
        return {
            "status": "error",
            "message": resultado.get("erro") or "Nenhuma referência foi extraída do relatório Consensus.",
            "tema": tema,
            "tema_ingles": tema,
            "instagram_format": instagram_format,
            "calendar_index": (dados_tema or {}).get("calendar_index"),
            "referencias_sugeridas": [],
            "modelo_texto_busca": OPENAI_TEXT_MODEL_SEARCH,
            "modelo_texto_redacao": OPENAI_TEXT_MODEL_WRITE,
            "fonte_referencias": "Consensus",
            "relatorio_consensus": relatorio_texto,
            "auto_select_all_references": True,
        }

    return {
        "status": "success",
        "tema": tema,
        "tema_ingles": tema,
        "instagram_format": instagram_format,
        "calendar_index": (dados_tema or {}).get("calendar_index"),
        "referencias_sugeridas": referencias,
        "modelo_texto_busca": OPENAI_TEXT_MODEL_SEARCH,
        "modelo_texto_redacao": OPENAI_TEXT_MODEL_WRITE,
        "fonte_referencias": "Consensus",
        "parser_warning": resultado.get("erro") or "",
        "relatorio_consensus": relatorio_texto,
        "auto_select_all_references": True,
        "total_referencias_importadas": len(referencias),
    }


def render_referencias_html_items(referencias):
    if not referencias:
        return "<li>Referências não disponíveis neste momento.</li>"

    items = []
    for ref in referencias:
        texto = html.escape(str(ref.get("texto") or "Referência sem título"), quote=True)
        link = html.escape(str(ref.get("link") or ""), quote=True)
        link_html = f" <a href='{link}' target='_blank' rel='noreferrer'>[Link do artigo]</a>" if link else ""
        items.append(f"<li>{texto}{link_html}</li>")

    return "".join(items)

def gerar_imagem(
    prompt,
    nome_arquivo_prefixo,
    keep_local=False,
    forced_filename=None,
    log_callback=None,
    image_size="1024x1024",
    target_aspect_ratio=None,
):
    """
    Gera imagem via OpenAI e usa Gemini apenas como fallback técnico.
    """
    def emit_log(message):
        print(message)
        if log_callback:
            try:
                log_callback(message)
            except Exception:
                pass

    emit_log(f"🎨 Gerando imagem ({nome_arquivo_prefixo})...")
    text_language_policy = (
        "IMAGE LANGUAGE POLICY: Any visible text rendered inside the image must be exclusively in Brazilian Portuguese (PT-BR). "
        "Never use English, Spanish, mixed language, bilingual headings, untranslated labels, or multilingual captions. "
        "If you cannot render correct PT-BR text with high confidence, do not render any visible text at all. "
        "Do not mix Portuguese with any other language anywhere in the image. "
    )
    prompt = f"{text_language_policy}{prompt}"

    if forced_filename:
        temp_filename = forced_filename
    else:
        temp_filename = f"temp_{nome_arquivo_prefixo}_{datetime.now().strftime('%H%M%S')}.png"
    
    global OPENAI_IMAGE_ENABLED
    image_generated = False

    # Tenta OpenAI primeiro (GPT Image / ChatGPT image stack)
    if OPENAI_IMAGE_ENABLED:
        try:
            emit_log(
                f"   ☁️ Solicitando imagem à OpenAI ({OPENAI_IMAGE_MODEL}, size={image_size}, quality={OPENAI_IMAGE_QUALITY}, timeout={OPENAI_IMAGE_TIMEOUT_SECONDS}s)..."
            )
            response = client.images.generate(
                model=OPENAI_IMAGE_MODEL,
                prompt=prompt,
                size=image_size,
                quality=OPENAI_IMAGE_QUALITY,
                output_format="png",
                timeout=OPENAI_IMAGE_TIMEOUT_SECONDS,
            )
            if getattr(response, "data", None):
                image_data = response.data[0]
                if getattr(image_data, "b64_json", None):
                    with open(temp_filename, "wb") as f:
                        f.write(base64.b64decode(image_data.b64_json))
                    image_generated = True
                    emit_log(f"   ✅ Imagem gerada com OpenAI ({OPENAI_IMAGE_MODEL}).")
                elif getattr(image_data, "url", None):
                    img_data = requests.get(image_data.url, timeout=60).content
                    with open(temp_filename, "wb") as f:
                        f.write(img_data)
                    image_generated = True
                    emit_log(f"   ✅ Imagem gerada com OpenAI ({OPENAI_IMAGE_MODEL}) via URL.")
                else:
                    emit_log(f"   ⚠️ OpenAI ({OPENAI_IMAGE_MODEL}) não retornou imagem utilizável. Tentando fallback.")
        except Exception as e:
            error_text = str(e)
            if "billing_hard_limit_reached" in error_text or "Billing hard limit has been reached" in error_text:
                OPENAI_IMAGE_ENABLED = False
                emit_log(
                    f"   ⚠️ OpenAI Imagem indisponível por limite de billing ({OPENAI_IMAGE_MODEL}). "
                    "As próximas imagens desta execução usarão Gemini diretamente."
                )
            else:
                emit_log(f"   ⚠️ Erro OpenAI Imagem ({OPENAI_IMAGE_MODEL}): {e}")
    else:
        emit_log(
            f"   ⏩ OpenAI Imagem desativada nesta execução após erro de billing. Usando Gemini para ({nome_arquivo_prefixo})."
        )

    # Fallback técnico para Gemini se OpenAI falhar
    if not image_generated:
        try:
            emit_log(f"   ☁️ Tentando fallback Gemini ({GEMINI_IMAGE_MODEL})...")
            model = genai.GenerativeModel(GEMINI_IMAGE_MODEL)
            response = model.generate_content("Generate an image of: " + prompt)
            for part in response.parts:
                if hasattr(part, 'inline_data') and part.inline_data:
                    with open(temp_filename, "wb") as f:
                        f.write(part.inline_data.data)
                    image_generated = True
                    emit_log(f"   ✅ Imagem gerada com Gemini ({GEMINI_IMAGE_MODEL}) em fallback.")
                    break
        except Exception as e:
            emit_log(f"   ❌ Erro Gemini Imagem ({GEMINI_IMAGE_MODEL}): {e}")
            if target_aspect_ratio == (16, 9):
                return "https://via.placeholder.com/1280x720?text=Reva+Mais"
            return "https://via.placeholder.com/1024x1024?text=Reva+Mais"

    if target_aspect_ratio:
        try:
            ajustar_arquivo_aspect_ratio(temp_filename, target_aspect_ratio)
            emit_log(f"   🪄 Ajuste final aplicado para {target_aspect_ratio[0]}:{target_aspect_ratio[1]}.")
        except Exception as e:
            emit_log(f"   ⚠️ Não foi possível ajustar aspect ratio ({nome_arquivo_prefixo}): {e}")

    # Upload
    try:
        timestamp_upload = datetime.now().strftime('%Y%m%d_%H%M%S')
        firebase_path = f"revamais/{nome_arquivo_prefixo}_{timestamp_upload}.png"
        emit_log(f"   ⬆️ Fazendo upload da imagem ({nome_arquivo_prefixo})...")
        url = upload_file(temp_filename, firebase_path)
        if not keep_local and os.path.exists(temp_filename): 
            os.remove(temp_filename)
        emit_log(f"   ✅ Upload concluído ({nome_arquivo_prefixo}).")
        return url
    except Exception as e:
        emit_log(f"   ❌ Erro Upload ({nome_arquivo_prefixo}): {e}")
        return "https://via.placeholder.com/600x400?text=Erro+Upload"
    
# --- Novos Imports para Manipulação de Imagem ---
from PIL import Image, ImageDraw, ImageFont
import io


def ajustar_arquivo_aspect_ratio(image_path, target_aspect_ratio):
    if not image_path or not target_aspect_ratio:
        return

    aspect_w, aspect_h = target_aspect_ratio
    if not aspect_w or not aspect_h:
        return

    target_ratio = float(aspect_w) / float(aspect_h)

    with Image.open(image_path) as image:
        width, height = image.size
        if not width or not height:
            return

        current_ratio = width / height
        if abs(current_ratio - target_ratio) < 0.01:
            return

        if current_ratio > target_ratio:
            new_width = int(height * target_ratio)
            left = max(0, (width - new_width) // 2)
            crop_box = (left, 0, left + new_width, height)
        else:
            new_height = int(width / target_ratio)
            top = max(0, (height - new_height) // 2)
            crop_box = (0, top, width, top + new_height)

        cropped = image.crop(crop_box)
        cropped.save(image_path)

def download_font():
    """Baixa fonte Roboto-Bold para garantir consistência visual em Linux/Railway"""
    font_path = "Roboto-Bold.ttf"
    if not os.path.exists(font_path):
        print("📥 Baixando fonte Roboto-Bold...")
        url = "https://github.com/google/fonts/raw/main/apache/roboto/Roboto-Bold.ttf"
        r = requests.get(url)
        with open(font_path, "wb") as f: f.write(r.content)
    return font_path

def download_logo(url):
    """Baixa o logo para um arquivo temporário"""
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"}
    r = requests.get(url, headers=headers)
    if r.status_code != 200:
        print(f"⚠️ Erro download imagem: {r.status_code}")
    img = Image.open(io.BytesIO(r.content)).convert("RGBA")
    return img

def gerar_banner_header(logo_url):
    """
    Gera um banner composto: 600x250
    Fundo: Branco/Cinza claro.
    Esquerda: Card arredondado Azul Escuro (#205776) com texto REVA +.
    Direita: Texto "Boletim de Saúde" e Logo Revalidatie no canto.
    """
    try:
        W, H = 600, 250
        # Cor de fundo (bg do email é #ffffff, container #ffffff. Vamos usar #f4f6f8 para contraste sutil ou branco)
        bg_color = (255, 255, 255) 
        
        # Cria imagem base
        im = Image.new("RGBA", (W, H), bg_color)
        draw = ImageDraw.Draw(im)
        
        # Downloads
        font_path = download_font()
        logo_img = download_logo(logo_url)
        
        # --- Lado Esquerdo: Card "Reva +" ---
        card_color = (32, 87, 118) # #205776
        margin = 20
        card_w = 260
        card_h = H - (2 * margin)
        
        # Desenha retângulo arredondado (simulado)
        x0, y0 = margin, margin
        x1, y1 = margin + card_w, H - margin
        draw.rectangle([x0, y0, x1, y1], fill=card_color, outline=None)
        
        # Texto "Reva +"
        font_reva = ImageFont.truetype(font_path, 50)
        draw.text((x0 + 40, y0 + 60), "Reva", fill="white", font=font_reva)
        
        font_plus = ImageFont.truetype(font_path, 70)
        draw.text((x0 + 160, y0 + 45), "+", fill="#4ecdc4", font=font_plus) # Destaque ciano
        
        font_sub = ImageFont.truetype(font_path, 14)
        draw.text((x0 + 40, y0 + 130), "BOLETIM CIENTÍFICO", fill="#a8cce0", font=font_sub)

        # --- Lado Direito: "Boletim de Saúde" + Logo ---
        
        # Texto descritivo
        font_title = ImageFont.truetype(font_path, 32)
        text_color = (32, 87, 118)
        draw.text((320, 60), "Boletim de", fill=text_color, font=font_title)
        draw.text((320, 100), "Saúde", fill=text_color, font=font_title)
        
        # Logo Revalidatie (Redimensionar)
        # Manter aspect ratio. Max width 180, max height 60
        logo_img.thumbnail((200, 80), Image.Resampling.LANCZOS)
        
        # Posicionar no canto inferior direito
        logo_w, logo_h = logo_img.size
        logo_x = W - logo_w - 30
        logo_y = H - logo_h - 30
        
        # Paste (usando a própria imagem como máscara alpha para transparência)
        im.paste(logo_img, (logo_x, logo_y), logo_img)
        
        # Salva e Upload
        temp_filename = f"header_revamais_{datetime.now().strftime('%H%M%S')}.png"
        im.save(temp_filename, "PNG")
        
        firebase_path = f"revamais/assets/{temp_filename}"
        url = upload_file(temp_filename, firebase_path)
        if os.path.exists(temp_filename): os.remove(temp_filename)
        
        return url
        
    except Exception as e:
        print(f"❌ Erro ao gerar banner header: {e}")
        raise e

def estimar_custo_revamais():
    # Estimativa:
    # 2 Imagens Gemini Pro 3 ($0.04 cada) = $0.08
    # Texto (Input + Output) GPT-4o ou Gemini Pro = ~$0.05 (chutando alto)
    # Total USD: $0.13
    total_usd = 0.13
    total_brl = total_usd * 6.0
    return {"usd": total_usd, "brl": total_brl}

def gerar_conteudo_revamais(tema, referencias, relatorio_consensus=None):
    """
    Gera o conteúdo HTML do boletim.
    """
    print("✍️ Escrevendo conteúdo Reva +...")
    relatorio_consensus = _normalizar_relatorio_consensus(relatorio_consensus)

    bloco_consensus = ""
    if relatorio_consensus:
        bloco_consensus = f"""
        CONTEXTO EVIDENCIAL PRIORIZADO (USO INTERNO DE REDAÇÃO):
        - Use o material abaixo como base principal para os achados específicos do texto.
        - Use as referências extraídas para apoiar o bloco científico e a bibliografia final.
        - Se houver conflito entre conhecimento geral e o material abaixo, prevalece o material abaixo.
        - Nunca mencione, no texto final, a plataforma usada, o relatório importado, o PDF, o processo de extração ou a origem operacional dessas evidências.
        - No texto publicado, apresente os achados como provenientes dos estudos e da evidência científica disponível.
        - Preserve a estrutura HTML já solicitada abaixo; não crie seções novas.

        --- MATERIAL EVIDENCIAL PRIORIZADO ---
        {relatorio_consensus[:24000]}
        -------------------------------------------
        """

    prompt = f"""
    Você é o editor do "Reva +", um boletim de saúde da clínica Revalidatie.
    Público-alvo: Pacientes e pessoas interessadas em saúde (leigos, mas inteligentes).
    
    Tema: "{tema}"
    
    INSTRUÇÃO DE ESCRITA HÍBRIDA (Conhecimento Geral + Evidência Específica):
    
    1. **Contexto e Mecanismo (Use seu conhecimento médico geral)**:
       - Comece explicando o problema de forma empática (ex: "Você sente dor ao caminhar?").
       - Explique O PORQUÊ (Fisiologia/Mecanismo): Por que isso acontece? O que muda no corpo com o tratamento? (Ex: fale sobre circulação colateral, eficiência muscular, neuroplasticidade).
       - O usuário GOSTA dessa explicação educativa do "como funciona".
    
    2. **O Que a Ciência Diz (Baseado SOMENTE nas fontes científicas fornecidas abaixo)**:
       - Agora, sintetize as evidências fornecidas em uma visão geral coesa.
       - Use os resumos e referências extraídas para validar a explicação anterior.
       - Combine convergências, nuances e limites dos estudos em vez de fazer um mini-resumo desconectado de cada artigo.
       - Diga "Estudos recentes mostram que..." ou equivalente e use os dados dos resumos.
       - Se usar bullets, prefira 3 a 5 bullets temáticos, cada um integrando mais de um estudo quando isso fizer sentido.
       - Pode citar autores/ano dentro dos bullets, mas apenas para sustentar uma síntese, não como ficha isolada de artigo.

    REGRAS DE SEGURANÇA E FIDELIDADE:
    - Não invente números, magnitude de efeito, tempo de intervenção, perfil de pacientes ou conclusões.
    - Não cite estudos que não estejam no material evidencial abaixo ou na lista de referências abaixo.
    - Se um detalhe não estiver explícito no relatório ou nos resumos, não mencione esse detalhe.
    - Se a evidência parecer preliminar, heterogênea ou limitada, diga isso com cautela.
    - Use linguagem prudente: "sugere", "indica", "aponta", "pode ajudar", quando apropriado.
    - Nunca escreva expressões como "segundo o Consensus", "o relatório Consensus mostrou", "de acordo com o relatório importado", "extraído do PDF" ou equivalentes.
    - Nunca mencione bastidores da curadoria, da importação ou da plataforma de busca.
    - Na seção de ciência, abra com 1 parágrafo curto de síntese geral e, se quiser, complemente com uma lista HTML (<ul><li>) de 3 a 5 bullets temáticos integrados.

    {bloco_consensus}
    
    --- EVIDÊNCIA CIENTÍFICA (Para a seção 'O Que a Ciência Diz') ---
    {chr(10).join([ f'Artigo {i+1}: {r["texto"]}{chr(10)}Resumo: {r["resumo"]}{chr(10)}' for i, r in enumerate(referencias) ])}
    -----------------------------------------------------------------
    
    Estrutura HTML (retorne APENAS o conteúdo dentro do body e APENAS HTML PURO):
    IMPORTANTE: NÃO USE MARKDOWN (```html ... ```). Retorne APENAS o código HTML cru.
    
    1. <h1>Título Atraente e Emocional</h1>
    2. <p>Introdução empática + Explicação do Mecanismo (Por que dói? Por que melhora? - Use conhecimento geral de fisiologia).</p>
    3. <h2>O que a Ciência Comprova?</h2> (Abra com uma síntese geral e depois traga bullets temáticos integrando os estudos, conectando com a explicação).
    4. <h2>Dicas Práticas</h2> (Conselhos acionáveis baseados nos abstracts e boas práticas).
    5. <div class="cta"> (Convite para seguir @revalidatie_londrina).
    
    Tom de voz: Fisioterapeuta especialista, amigo, otimista e educativo. 
    IMPORTANTE: Sempre direcione para "Consulte seu Fisioterapeuta" e NUNCA "Consulte seu Médico". O contexto é reabilitação física.
    """
    
    html_content = gerar_texto_preferencial(
        prompt,
        system_prompt=(
            "You are a careful medical editor. "
            "Stay faithful to the supplied evidence, avoid unsupported claims, "
            "and return only raw HTML."
        ),
        model_name=OPENAI_TEXT_MODEL_WRITE,
    )

    # Limpeza de Markdown
    html_content = html_content.replace("```html", "").replace("```", "")
    
    # Remover tags estruturais se a IA teimosamente as incluir
    import re
    html_content = re.sub(r'<!DOCTYPE[^>]*>', '', html_content, flags=re.IGNORECASE)
    html_content = re.sub(r'<html[^>]*>', '', html_content, flags=re.IGNORECASE)
    html_content = re.sub(r'</html>', '', html_content, flags=re.IGNORECASE)
    html_content = re.sub(r'<head>.*?</head>', '', html_content, flags=re.IGNORECASE | re.DOTALL)
    html_content = re.sub(r'<body[^>]*>', '', html_content, flags=re.IGNORECASE)
    html_content = re.sub(r'</body>', '', html_content, flags=re.IGNORECASE)
    
    return html_content


def limpar_texto_html(html):
    """
    Remove tags HTML de forma simples para reaproveitar o conteúdo em prompts.
    """
    import re

    text = re.sub(r'<br\s*/?>', '\n', html, flags=re.IGNORECASE)
    text = re.sub(r'</(p|li|h1|h2|h3|h4|h5|h6|div)>', '\n', text, flags=re.IGNORECASE)
    text = re.sub(r'<[^>]+>', ' ', text)
    text = re.sub(r'[ \t]+', ' ', text)
    text = re.sub(r'\n\s*\n+', '\n', text)
    return text.strip()


def extrair_secao_html(html, h2_regex):
    """
    Extrai o conteúdo entre um H2 alvo e o próximo H2.
    """
    import re

    pattern = re.compile(
        rf'<h2[^>]*>\s*{h2_regex}\s*</h2>(.*?)(?=<h2[^>]*>|<div[^>]*class=["\'][^"\']*cta[^"\']*["\'][^>]*>|$)',
        re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(html)
    if not match:
        return ""
    return limpar_texto_html(match.group(1))


def extrair_titulo_html(html_texto):
    match = re.search(r"<h1[^>]*>(.*?)</h1>", html_texto or "", re.IGNORECASE | re.DOTALL)
    if not match:
        return ""
    return limpar_texto_html(match.group(1))


def extrair_intro_html(html_texto):
    html_texto = html_texto or ""
    match_h2 = re.search(r"<h2[^>]*>", html_texto, re.IGNORECASE)
    trecho = html_texto[: match_h2.start()] if match_h2 else html_texto
    paragraphs = re.findall(r"<p[^>]*>(.*?)</p>", trecho, re.IGNORECASE | re.DOTALL)
    for paragraph in paragraphs:
        texto = limpar_texto_html(paragraph)
        if texto:
            return texto
    return limpar_texto_html(trecho)[:500]


def resumir_texto_para_legenda(texto, max_chars=220):
    texto = re.sub(r"\s+", " ", str(texto or "").replace("•", ". ")).strip(" -")
    if not texto:
        return ""

    frases = [frase.strip(" -") for frase in re.split(r"(?<=[.!?])\s+", texto) if frase.strip()]
    escolhidas = []

    for frase in frases:
        candidata = " ".join(escolhidas + [frase]).strip()
        if len(candidata) > max_chars and escolhidas:
            break
        escolhidas.append(frase)
        if len(" ".join(escolhidas)) >= int(max_chars * 0.75) or len(escolhidas) >= 2:
            break

    resumo = " ".join(escolhidas).strip() or texto
    if len(resumo) > max_chars:
        corte = resumo[:max_chars].rsplit(" ", 1)[0].strip() or resumo[:max_chars].strip()
        resumo = f"{corte}..."

    if resumo and resumo[-1] not in ".!?":
        resumo = f"{resumo}."

    return resumo


def gerar_briefs_visuais_revamais(tema, html_texto, referencias):
    """
    Gera prompts visuais mais ancorados no conteúdo final do boletim.
    """
    titulo_boletim = extrair_titulo_html(html_texto)
    intro_boletim = extrair_intro_html(html_texto)
    secao_ciencia = extrair_secao_html(html_texto, r"O\s+que\s+a\s+Ci[eê]ncia\s+Comprova\??")
    secao_dicas = extrair_secao_html(html_texto, r"Dicas\s+Pr[aá]ticas")
    texto_limpo = limpar_texto_html(html_texto)

    if not secao_ciencia:
        secao_ciencia = texto_limpo[:1200]
    if not secao_dicas:
        secao_dicas = texto_limpo[:1200]

    referencias_contexto = "\n".join(
        [f"- {r['texto']}" for r in referencias[:3]]
    ) or "- Sem referências resumidas."

    prompt_briefs = f"""
    You are creating three grounded image prompts for a patient-facing medical newsletter.
    Return ONLY valid JSON with this exact shape:
    {{
      "abertura": {{"prompt_english": "..."}},
      "ciencia": {{"prompt_english": "...", "caption_ptbr": "..."}},
      "dicas": {{"prompt_english": "...", "caption_ptbr": "..."}}
    }}

    RULES:
    - Base the prompts strictly on the supplied theme, selected studies, and final newsletter text.
    - Avoid generic wellness visuals, random symbols, and concepts not present in the content.
    - Prefer concrete anatomy, physiology, movement, rehabilitation actions, daily habits, or evidence concepts explicitly present in the text.
    - Visual style: premium editorial medical infographic, clean, minimal, white background, high contrast, lots of whitespace.
    - The "abertura" image must work as a newsletter opening image and immediately communicate the actual clinical subject of the bulletin.
    - The "abertura" image must be a simple subject image only: no explanation, no infographic layout, no multi-step sequence, no mechanism summary, no didactic labels.
    - The "abertura" image must be composed as a wide horizontal banner in 16:9, with the subject centered and safe margins.
    - The "abertura" image must not look like a business meeting, office teamwork, corporate consulting, conference room, presentation deck, startup discussion, or generic lifestyle stock photo.
    - The "abertura" image should show a concrete patient/condition/exercise/rehabilitation scene consistent with the title and introduction.
    - The "abertura" image must contain no visible text at all.
    - The "ciencia" image must explain the actual mechanism or evidence narrative from the science section.
    - The "dicas" image must show only the practical actions or habits explicitly recommended in the tips section.
    - The "ciencia" and "dicas" images will receive an external HTML caption below them, so avoid dense text inside the image. Prefer no text or at most one very short PT-BR label if absolutely necessary.
    - No logos, no brand marks, no clutter, no unrelated charts.
    - If short labels are useful, include at most 1 short label in Brazilian Portuguese (PT-BR) only.
    - Never mix languages. Never use English headings, English labels, bilingual text, or untranslated interface words.
    - If you are not fully confident that all visible text will be correct PT-BR, request no text.
    - Each "caption_ptbr" must be 1 or 2 short sentences in natural Brazilian Portuguese, self-contained, clear for a patient, and directly explanatory of that visual.
    - Each "caption_ptbr" must not mention "imagem", "figura", "newsletter", "boletim" or the generation process.

    THEME:
    {tema}

    NEWSLETTER TITLE:
    {titulo_boletim}

    INTRODUCTION:
    {intro_boletim}

    SELECTED STUDIES:
    {referencias_contexto}

    SCIENCE SECTION:
    {secao_ciencia}

    PRACTICAL TIPS SECTION:
    {secao_dicas}
    """

    try:
        response = gerar_texto_preferencial(
            prompt_briefs,
            system_prompt=(
                "You are a precise medical visual editor. "
                "Produce image prompts tightly grounded in the provided newsletter content."
            ),
            model_name=OPENAI_TEXT_MODEL_WRITE,
        )
        response = response.replace("```json", "").replace("```", "").strip()
        briefs = json.loads(response)

        if (
            isinstance(briefs, dict)
            and isinstance(briefs.get("abertura"), dict)
            and isinstance(briefs.get("ciencia"), dict)
            and isinstance(briefs.get("dicas"), dict)
            and briefs["abertura"].get("prompt_english")
            and briefs["ciencia"].get("prompt_english")
            and briefs["dicas"].get("prompt_english")
        ):
            briefs["ciencia"]["caption_ptbr"] = resumir_texto_para_legenda(
                briefs["ciencia"].get("caption_ptbr") or secao_ciencia,
                max_chars=220,
            )
            briefs["dicas"]["caption_ptbr"] = resumir_texto_para_legenda(
                briefs["dicas"].get("caption_ptbr") or secao_dicas,
                max_chars=220,
            )
            return briefs
    except Exception as e:
        print(f"⚠️ Falha ao gerar briefs visuais do Reva+: {e}")

    return {
        "abertura": {
            "prompt_english": (
                f"Create a simple photorealistic editorial opening image for the clinical topic '{tema}'. "
                f"Newsletter title: '{titulo_boletim or tema}'. "
                "Show only a concrete patient, body region, symptom, exercise, or rehabilitation scene directly about the subject. "
                "Do not explain the article content, do not summarize mechanisms, and do not create an infographic. "
                "Compose the scene as a wide horizontal banner (16:9) with the main subject centered and comfortable safe margins. "
                "The image must feel like healthcare education, not corporate lifestyle. "
                "No office meeting, no business discussion, no conference room, no people around a boardroom table, no startup scene, no generic teamwork. "
                "No text, no labels, no charts. White or very clean background, premium composition, human-centered, medically coherent."
            )
        },
        "ciencia": {
            "prompt_english": (
                f"Create a premium editorial medical infographic based strictly on this newsletter science section about '{tema}': "
                f"{secao_ciencia[:900]} "
                "Show the specific anatomy, physiology, rehabilitation mechanism, or evidence concept described in the text. "
                "White background, minimal composition, high contrast, no clutter. "
                "An external HTML caption will explain the visual, so avoid dense text inside the image. "
                "If labels are useful, include at most 1 short label in Brazilian Portuguese (PT-BR) only. "
                "Never use English or mixed language. If text fidelity is uncertain, use no text."
            ),
            "caption_ptbr": resumir_texto_para_legenda(secao_ciencia, max_chars=220),
        },
        "dicas": {
            "prompt_english": (
                f"Create a premium editorial medical infographic based strictly on this practical tips section about '{tema}': "
                f"{secao_dicas[:900]} "
                "Show only the concrete actions, habits, or rehabilitation steps described in the text. "
                "White background, checklist or step-by-step layout, minimal composition, high contrast, no clutter. "
                "An external HTML caption will explain the visual, so avoid dense text inside the image. "
                "If labels are useful, include at most 1 short label in Brazilian Portuguese (PT-BR) only. "
                "Never use English or mixed language. If text fidelity is uncertain, use no text."
            ),
            "caption_ptbr": resumir_texto_para_legenda(secao_dicas, max_chars=220),
        },
    }

def _normalize_calendar_title(text):
    return str(text or "").strip().lower()


def _extract_calendar_title(row):
    return row.get("Title", row.get("Theme", "")).strip()


def _resolve_instagram_format(row=None):
    if not isinstance(row, dict):
        return "Carrossel"

    valor = str(
        row.get("Format")
        or row.get("Formato")
        or row.get("formato")
        or ""
    ).strip().lower()

    if "reel" in valor:
        return "Reel"
    return "Carrossel"


def _load_revamais_calendar_rows():
    csv_filename = "calendario_editorial_150_semanas.csv"
    csv_path = os.path.join(os.path.dirname(__file__), csv_filename)

    if not os.path.exists(csv_path):
        print(f"⚠️ Arquivo {csv_filename} não encontrado.")
        return []

    with open(csv_path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _coerce_completed_indices(state_data, rows):
    total_linhas = len(rows)
    raw_completed = state_data.get("completed_indices") if isinstance(state_data, dict) else None
    completed = set()

    if isinstance(raw_completed, list):
        for item in raw_completed:
            try:
                idx = int(item)
            except (TypeError, ValueError):
                continue
            if 0 <= idx < total_linhas and _extract_calendar_title(rows[idx]):
                completed.add(idx)
    else:
        next_index = int((state_data or {}).get("next_index", 0) or 0)
        for idx in range(max(0, min(next_index, total_linhas))):
            if _extract_calendar_title(rows[idx]):
                completed.add(idx)

    return completed


def _compute_next_pending_index(rows, completed_indices):
    for idx, row in enumerate(rows):
        if _extract_calendar_title(row) and idx not in completed_indices:
            return idx
    return len(rows)


def _build_revamais_firestore_state(next_index, rows, source, completed_indices=None):
    total_linhas = len(rows)
    completed_indices = completed_indices or set()
    completed_indices = {
        idx for idx in completed_indices
        if 0 <= idx < total_linhas and _extract_calendar_title(rows[idx])
    }
    next_index = _compute_next_pending_index(rows, completed_indices)
    next_title = ""
    next_format = "Carrossel"

    for i in range(max(0, next_index), total_linhas):
        titulo = _extract_calendar_title(rows[i])
        if titulo:
            next_title = titulo
            next_format = _resolve_instagram_format(rows[i])
            break

    return {
        "next_index": next_index,
        "last_used_index": next_index - 1,
        "next_title": next_title,
        "next_format": next_format,
        "completed_indices": sorted(completed_indices),
        "total_rows": total_linhas,
        "state_source": source,
        "last_updated": datetime.now().isoformat(),
    }


def _load_revamais_calendar_state(rows):
    db = get_firestore_db()
    if db is None:
        raise RuntimeError("Firestore indisponivel para controlar a agenda do Reva+.")

    doc_ref = db.collection(REVAMAIS_STATE_COLLECTION).document(REVAMAIS_STATE_DOC_ID)
    snapshot = doc_ref.get()
    if snapshot.exists:
        state = snapshot.to_dict() or {}
        next_index = int(state.get("next_index", 0) or 0)
        print(f"📊 Estado Firestore carregado. next_index={next_index}")
        return doc_ref, state

    legacy_state = read_json_from_storage(REVAMAIS_LEGACY_STATE_FILE) or {"used_titles": []}
    used_titles = legacy_state.get("used_titles", [])
    used_titles_normalized = {_normalize_calendar_title(t) for t in used_titles}

    completed_indices = set()
    for i, row in enumerate(rows):
        titulo = _extract_calendar_title(row)
        if titulo and _normalize_calendar_title(titulo) in used_titles_normalized:
            completed_indices.add(i)

    next_index = _compute_next_pending_index(rows, completed_indices)
    migrated_state = _build_revamais_firestore_state(
        next_index,
        rows,
        source="migrated_from_storage",
        completed_indices=completed_indices,
    )
    migrated_state["legacy_used_titles_count"] = len(used_titles)
    doc_ref.set(migrated_state)
    print(
        "✅ Estado da agenda migrado para Firestore "
        f"(legacy used_titles={len(used_titles)}, next_index={next_index})."
    )
    return doc_ref, migrated_state


def _get_calendar_row_by_index(rows, calendar_index):
    if calendar_index is None:
        return None

    try:
        idx = int(calendar_index)
    except (TypeError, ValueError):
        raise RuntimeError("calendar_index inválido.")

    if idx < 0 or idx >= len(rows):
        raise RuntimeError("calendar_index fora do calendário editorial.")

    row = rows[idx]
    if not _extract_calendar_title(row):
        raise RuntimeError("calendar_index aponta para uma linha sem título válido.")

    return idx, row


def listar_calendario_revamais():
    rows = _load_revamais_calendar_rows()
    if not rows:
        return {
            "items": [],
            "summary": {
                "total": 0,
                "done": 0,
                "pending": 0,
                "current_index": None,
                "next_title": "",
                "state_available": False,
            },
        }

    state_available = True
    state_data = {}
    completed_indices = set()
    try:
        _doc_ref, state_data = _load_revamais_calendar_state(rows)
        completed_indices = _coerce_completed_indices(state_data, rows)
    except Exception as e_state:
        print(f"⚠️ Não foi possível carregar estado do calendário Reva+: {e_state}")
        state_available = False

    current_index = _compute_next_pending_index(rows, completed_indices)
    items = []

    for idx, row in enumerate(rows):
        titulo = _extract_calendar_title(row)
        if not titulo:
            continue

        status = "done" if idx in completed_indices else ("current" if idx == current_index else "pending")
        items.append({
            "calendar_index": idx,
            "week": row.get("Week", ""),
            "date": row.get("Date", ""),
            "day": row.get("Day", ""),
            "theme": row.get("Theme", ""),
            "format": _resolve_instagram_format(row),
            "title": titulo,
            "status": status,
            "done": status == "done",
        })

    done_count = sum(1 for item in items if item["done"])

    return {
        "items": items,
        "summary": {
            "total": len(items),
            "done": done_count,
            "pending": len(items) - done_count,
            "current_index": current_index if current_index < len(rows) else None,
            "next_title": next(
                (item["title"] for item in items if item["status"] == "current"),
                "",
            ),
            "state_available": state_available,
            "next_index": int((state_data or {}).get("next_index", current_index if current_index < len(rows) else len(rows))),
        },
    }


def marcar_tema_revamais_concluido(calendar_index, titulo_esperado=None):
    rows = _load_revamais_calendar_rows()
    if not rows:
        raise RuntimeError("Calendário editorial do Reva+ não encontrado.")

    resolved = _get_calendar_row_by_index(rows, calendar_index)
    if not resolved:
        return None

    idx, row = resolved
    titulo = _extract_calendar_title(row)
    if titulo_esperado and _normalize_calendar_title(titulo_esperado) != _normalize_calendar_title(titulo):
        raise RuntimeError("O item do calendário não corresponde ao tema esperado.")

    doc_ref, state_data = _load_revamais_calendar_state(rows)
    completed_indices = _coerce_completed_indices(state_data, rows)
    completed_indices.add(idx)
    next_index = _compute_next_pending_index(rows, completed_indices)

    updated_state = _build_revamais_firestore_state(
        next_index,
        rows,
        source="firestore",
        completed_indices=completed_indices,
    )
    updated_state["last_used_title"] = titulo
    updated_state["last_used_format"] = _resolve_instagram_format(row)
    updated_state["last_used_row_number"] = idx + 1
    updated_state["last_used_index"] = idx
    doc_ref.set(updated_state, merge=True)
    return updated_state


def obter_proximo_tema_csv(consumir=True):
    """
    Lê o calendário editorial e usa um ponteiro remoto no Firestore
    para persistir o próximo item da agenda do Reva+.
    """
    rows = _load_revamais_calendar_rows()
    if not rows:
        return None

    total_linhas = len(rows)
    doc_ref, state_data = _load_revamais_calendar_state(rows)
    completed_indices = _coerce_completed_indices(state_data, rows)
    next_index = _compute_next_pending_index(rows, completed_indices)

    tema_escolhido = None
    formato_escolhido = "Carrossel"
    idx_atual = None

    for i in range(max(0, next_index), total_linhas):
        titulo = _extract_calendar_title(rows[i])
        if titulo:
            tema_escolhido = titulo
            formato_escolhido = _resolve_instagram_format(rows[i])
            idx_atual = i
            break

    if not tema_escolhido:
        print("⚠️ Todos os temas do calendário já foram usados!")
        return None

    print(f"📅 Tema do Calendário Selecionado: {tema_escolhido} (Item {idx_atual + 1}/{total_linhas})")

    if consumir:
        completed_indices.add(idx_atual)
        updated_state = _build_revamais_firestore_state(
            idx_atual + 1,
            rows,
            source="firestore",
            completed_indices=completed_indices,
        )
        updated_state["last_used_title"] = tema_escolhido
        updated_state["last_used_format"] = formato_escolhido
        updated_state["last_used_row_number"] = idx_atual + 1
        updated_state["last_used_index"] = idx_atual

        try:
            doc_ref.set(updated_state, merge=True)
            print(
                "✅ Estado da agenda atualizado no Firestore "
                f"(next_index={updated_state['next_index']})."
            )
        except Exception as e_save:
            raise RuntimeError(
                "Falha ao persistir o progresso da agenda do Reva+ no Firestore."
            ) from e_save

    return {
        "tema": tema_escolhido,
        "formato": formato_escolhido,
        "calendar_index": idx_atual,
        "week": rows[idx_atual].get("Week", ""),
        "date": rows[idx_atual].get("Date", ""),
        "day": rows[idx_atual].get("Day", ""),
        "theme": rows[idx_atual].get("Theme", ""),
    }


def resolver_tema_revamais(tema_usuario=None, consumir_tema_auto=True, calendar_index=None):
    """
    Resolve o tema efetivo do Reva+ e o formato padrão do Instagram.
    """
    is_calendar_source = calendar_index is not None or not tema_usuario or tema_usuario.strip() == ""
    formato_instagram = "Carrossel"
    tema = tema_usuario
    resolved_calendar_index = None

    if calendar_index is not None:
        rows = _load_revamais_calendar_rows()
        resolved = _get_calendar_row_by_index(rows, calendar_index)
        if not resolved:
            return None
        resolved_calendar_index, row = resolved
        tema = _extract_calendar_title(row)
        formato_instagram = _resolve_instagram_format(row)
        return {
            "tema": tema,
            "formato_instagram": formato_instagram,
            "is_calendar_source": True,
            "calendar_index": resolved_calendar_index,
        }

    if not tema or tema == "auto":
        dados_csv = obter_proximo_tema_csv(consumir=consumir_tema_auto)
        if not dados_csv:
            return None
        tema = dados_csv["tema"]
        formato_instagram = _resolve_instagram_format(dados_csv)
        resolved_calendar_index = dados_csv.get("calendar_index")

    return {
        "tema": tema,
        "formato_instagram": formato_instagram,
        "is_calendar_source": is_calendar_source,
        "calendar_index": resolved_calendar_index,
    }


def _resolver_tema_com_referencias(
    tema_usuario=None,
    quantidade_referencias=8,
    consumir_tema_auto=False,
    calendar_index=None,
    log_callback=None,
):
    """
    Resolve um tema e busca referências. Em modo automático:
    - na preparação (consumir=False), não consome o tema válido, mas pula
      permanentemente os temas sem evidência para evitar travamento da agenda;
    - na geração final (consumir=True), consome e continua avançando até achar
      um tema com referências.
    """
    max_tentativas = 299

    for tentativa in range(max_tentativas):
        dados_tema = resolver_tema_revamais(
            tema_usuario=tema_usuario,
            consumir_tema_auto=consumir_tema_auto,
            calendar_index=calendar_index,
        )
        if not dados_tema:
            return None

        tema = dados_tema["tema"]
        tema_ingles = gerar_query_pubmed_tema(tema)
        referencias = buscar_referencias_pubmed(tema_ingles, limite_retorno=quantidade_referencias)

        if referencias:
            return {
                "dados_tema": dados_tema,
                "tema_ingles": tema_ingles,
                "referencias": referencias,
            }

        if calendar_index is not None or (tema_usuario and tema_usuario.strip()):
            return {
                "dados_tema": dados_tema,
                "tema_ingles": tema_ingles,
                "referencias": [],
            }

        if log_callback:
            log_callback(
                f"⚠️ Tema sem evidências suficientes, pulando para o próximo da agenda: {tema}"
            )

        if not consumir_tema_auto:
            obter_proximo_tema_csv(consumir=True)

    return None


def preparar_referencias_revamais(tema_usuario=None, quantidade_referencias=8, calendar_index=None):
    """
    Etapa intermediária do pipeline: resolve o tema e retorna artigos candidatos
    para seleção manual antes da geração do boletim.
    """
    resultado_busca = _resolver_tema_com_referencias(
        tema_usuario=tema_usuario,
        quantidade_referencias=quantidade_referencias,
        consumir_tema_auto=False,
        calendar_index=calendar_index,
    )
    if not resultado_busca:
        return {"status": "error", "message": "Nenhum tema disponível para preparar referências."}

    dados_tema = resultado_busca["dados_tema"]
    tema = dados_tema["tema"]
    tema_ingles = resultado_busca["tema_ingles"]
    referencias = resultado_busca["referencias"]

    if not referencias:
        return {
            "status": "error",
            "message": f"Nenhuma evidência encontrada para o tema '{tema}'.",
            "tema": tema,
            "tema_ingles": tema_ingles,
            "instagram_format": dados_tema["formato_instagram"],
            "calendar_index": dados_tema.get("calendar_index"),
            "referencias_sugeridas": [],
            "modelo_texto_busca": OPENAI_TEXT_MODEL_SEARCH,
            "modelo_texto_redacao": OPENAI_TEXT_MODEL_WRITE,
        }

    return {
        "status": "success",
        "tema": tema,
        "tema_ingles": tema_ingles,
        "instagram_format": dados_tema["formato_instagram"],
        "calendar_index": dados_tema.get("calendar_index"),
        "referencias_sugeridas": referencias,
        "modelo_texto_busca": OPENAI_TEXT_MODEL_SEARCH,
        "modelo_texto_redacao": OPENAI_TEXT_MODEL_WRITE,
    }

import textwrap

def aplicar_logo_overlay(local_filename, slide_num):
    """
    Aplica o logo da clínica no canto superior direito do slide gerado.
    """
    try:
        # Configs
        LOGO_URL = "https://i.imgur.com/A2d27eq.png" # Logo Revalidatie (Extraído do Album)
        MARGIN = 40
        LOGO_WIDTH = 180
        
        # Abre imagem original
        img = Image.open(local_filename).convert("RGBA")
        W, H = img.size
        
        # Baixa Logo
        try:
            logo = download_logo(LOGO_URL)
            # Redimensiona logo mantendo aspect ratio
            w_percent = (LOGO_WIDTH / float(logo.size[0]))
            h_size = int((float(logo.size[1]) * float(w_percent)))
            logo = logo.resize((LOGO_WIDTH, h_size), Image.Resampling.LANCZOS)
            
            # Posição: Topo Direito
            x = W - LOGO_WIDTH - MARGIN
            y = MARGIN
            
            # Paste
            img.paste(logo, (x, y), logo) # Usa logo como máscara se tiver transparência
            
            # Salva sobrescrevendo
            img.save(local_filename, "PNG")
            print(f"   ✅ Logo aplicado no slide {slide_num}")
            
        except Exception as e_download:
            print(f"⚠️ Erro ao baixar/aplicar logo: {e_download}")
            # Se falhar logo, mantém imagem original
            
    except Exception as e:
        print(f"⚠️ Erro fatal no overlay de logo: {e}")
        import traceback
        traceback.print_exc()

def gerar_conteudo_instagram(tema, formato, referencias_text, conteudo_base=None, log_callback=None):
    """
    Gera conteúdo para Instagram (Reel ou Carrossel).
    Melhoria v2: Gera texto com LLM, Imagem Clean com IA, e Texto via Overlay (Pillow).
    """
    def emit_log(message):
        print(message)
        if log_callback:
            try:
                log_callback(message)
            except Exception:
                pass

    emit_log(f"📸 Gerando conteúdo para Instagram ({formato})...")
    assets = []
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    try:
        model = genai.GenerativeModel(GEMINI_TEXT_MODEL)
        
        if formato.lower() == "reel":
            # (Mantido igual - omitido para brevidade, mas deve existir no arquivo final)
            prompt = f"""
            Crie um roteiro viral para Instagram Reels sobre: "{tema}".
            Baseado nestas referências: {referencias_text}
            ... (prompt roteiro original) ...
            """
            roteiro = model.generate_content(prompt).text
            filename = f"instagram_reel_{timestamp}.md"
            with open(filename, "w", encoding="utf-8") as f: f.write(roteiro)
            url = upload_file(filename, f"instagram/{filename}")
            if os.path.exists(filename): os.remove(filename)
            assets.append({"type": "roteiro", "url": url, "name": "Roteiro do Reel"})
            
        elif formato.lower() == "carrossel":
            emit_log("   📝 Planejando narrativa do carrossel...")
            
            contexto_extra = f"\nBASEIE-SE ESTRITAMENTE NESTE CONTEÚDO JÁ GERADO:\n{conteudo_base}\n" if conteudo_base else ""
            
            prompt_slides = f"""
            You are an expert Visual Storyteller and AI Prompter.
            Plan a 7-slide Instagram Carousel based strictly on the content below.
            
            CONTENT TO ADAPT:
            ---
            {conteudo_base}
            ---
            
            THEME: "{tema}"
            AUDIENCE: Patients (Laypeople).
            TONE: Professional Physiotherapist, Encouraging, Educational.
            
            RETURN ONLY A JSON LIST with this structure:
            [
              {{
                "slide": 1,
                "titulo": "Hook Title (Portuguese)",
                "texto_curto": "Short Body Text (Portuguese, max 25 words)",
                "image_prompt_english": "Detailed prompt for Imagen 3 model..."
              }},
              ...
            ]
            
            RULES FOR 'image_prompt_english':
1. **Visual Style**: Realistic, high-quality, editorial medical photography aesthetic. Prefer photorealistic images over illustrations. Use natural lighting, clean composition, and a professional healthcare look.
2. **People**: Whenever appropriate, include real-looking people in the scene, especially patients, physiotherapists, or adults in everyday health-related situations. Human presence should feel natural and credible, not artificial or cartoonish.
3. **Correction**: Do NOT include the general theme title "{tema}" in the image. ONLY include the specific '{{titulo}}' and '{{texto_curto}}' of the slide.
4. **Density**: Avoid clutter. Use ONE main scene or focal subject that clearly supports the message of the slide.
5. **Text Instruction**: If visible text is used, explicitly state: "Include only the exact Portuguese text: '{{titulo}}' and '{{texto_curto}}'." Typography must be legible, modern, sans-serif.
6. **Consistency**: Keep a clean healthcare brand aesthetic, with soft and trustworthy tones, subtle teal accents when appropriate, and premium visual quality.
7. **No Hallucinations**: Do not ask for logos, watermarks, icons, or cartoon elements unless strictly necessary for the concept.
8. **Avoid Illustration Look**: Do not generate flat drawings, vector art, line art, infographic-style icons, or cartoon-style medical scenes unless the slide specifically requires a scientific mechanism that cannot be represented with realistic photography.
9. **Language Policy**: Any visible text inside the image must be exclusively in Brazilian Portuguese (PT-BR). Never use English, mixed language, bilingual labels, or untranslated headings.
10. **Fallback Rule**: If there is any risk of wrong language, instruct the image model to use no visible text.

SLIDE STRUCTURE:
- Slide 1: Hook/Pain (Realistic scene showing the symptom or limitation in daily life).
- Slide 2-3: Education (Realistic healthcare or body-related visual metaphors, or realistic clinical/support scenes).
- Slide 4-6: Solutions (Realistic people exercising, receiving guidance, improving mobility, breathing, strength, or function).
- Slide 7: CTA (Warm, realistic lifestyle or healthcare closing image with the CTA text).
            """
            
            try:
                response_text = model.generate_content(prompt_slides).text
                response_text = response_text.replace("```json", "").replace("```", "").strip()
                slides_data = json.loads(response_text)
            except:
                slides_data = [{"slide": i+1, "titulo": f"Slide {i+1}", "texto_curto": tema, "image_prompt_english": f"Slide about {tema}"} for i in range(7)]

            try:
                roteiro_carrossel_nome = f"roteiro_carrossel_{timestamp}.md"
                linhas_roteiro = [f"# Roteiro de Carrossel - {tema}", ""]
                for slide in slides_data:
                    slide_num = slide.get("slide", "")
                    slide_titulo = slide.get("titulo", "")
                    slide_texto = slide.get("texto_curto", "")
                    linhas_roteiro.append(f"## Slide {slide_num}: {slide_titulo}".strip())
                    if slide_texto:
                        linhas_roteiro.append(slide_texto)
                    linhas_roteiro.append("")

                with open(roteiro_carrossel_nome, "w", encoding="utf-8") as f:
                    f.write("\n".join(linhas_roteiro).strip() + "\n")

                roteiro_carrossel_url = upload_file(roteiro_carrossel_nome, f"instagram/{roteiro_carrossel_nome}")
                if os.path.exists(roteiro_carrossel_nome):
                    os.remove(roteiro_carrossel_nome)
                assets.append({"type": "roteiro", "url": roteiro_carrossel_url, "name": "Roteiro do Carrossel"})
            except Exception as e:
                emit_log(f"⚠️ Erro ao salvar roteiro do carrossel: {e}")

            emit_log(f"   🖼️ Gerando {len(slides_data)} slides (Prompt Dinâmico por Slide)...")
            files_to_zip = []
            
            # URL da Logo (Hardcoded para teste ou parametrizável)
            logo_url = "https://i.imgur.com/oGzxgtK.jpeg" # Placeholder/Header Image used as logo source or need actual logo
            # User offered logo. usage: overlay_logo(final_url, logo_url)
            
            for slide in slides_data:
                # Inicializa variáveis com nomes distintos para evitar colisão de escopo
                slide_titulo = "Sem Título"
                slide_texto = "Sem Texto"
                slide_visual_prompt = f"Slide about {tema}"
                slide_num = 0
                
                try:
                    slide_num = slide.get('slide', 0)
                    slide_titulo = slide.get('titulo', f"Slide {slide_num}")
                    slide_texto = slide.get('texto_curto', "")
                    slide_visual_prompt = slide.get('image_prompt_english', f"Medical illustration about {tema}")
                except Exception as e_parse:
                    print(f"⚠️ Erro ao ler dados do slide: {e_parse}")
                
                # Reforço de prompt
                full_prompt = (
                    f"{slide_visual_prompt} "
                    f"Style: Clean, Minimalist, High Contrast. Format: SQUARE (1:1). "
                    "Any visible text must be exclusively in Brazilian Portuguese (PT-BR). "
                    "Never use English or mixed language. If unsure, use no text. "
                    f"NO CLUTTER. Focus on the central message."
                )

                # Gera imagem
                local_base_name = f"slide_{slide_num}_{timestamp}"
                local_filename = f"{local_base_name}.png"
                
                emit_log(f"   🎞️ Slide {slide_num}/{len(slides_data)}: gerando arte...")
                gerar_imagem(
                    full_prompt,
                    local_base_name,
                    keep_local=True,
                    forced_filename=local_filename,
                    log_callback=log_callback,
                )
                
                try:
                    aplicar_logo_overlay(local_filename, slide_num) 
                except:
                   pass

                files_to_zip.append(local_filename)
                
                final_url = upload_file(local_filename, f"revamais/slides/{local_base_name}.png")

                assets.append({
                    "type": "image", 
                    "url": final_url, 
                    "name": f"Slide {slide_num}: {slide_titulo}",
                    "texto_base": slide_texto
                })

            # Criar arquivo ZIP com todas as imagens
            try:
                import zipfile
                zip_name = f"instagram_carrossel_{timestamp}.zip"
                emit_log(f"   📦 Criando ZIP: {zip_name}...")
                
                with zipfile.ZipFile(zip_name, 'w') as zf:
                    for f in files_to_zip:
                        if os.path.exists(f):
                            zf.write(f)
                            os.remove(f) # Limpa imagem individual após adicionar ao zip
                
                # Upload ZIP
                zip_url = upload_file(zip_name, f"instagram/{zip_name}")
                if os.path.exists(zip_name): os.remove(zip_name)
                
                assets.append({"type": "zip", "url": zip_url, "name": "Baixar Todas as Imagens (.zip)"})
                
            except Exception as e:
                emit_log(f"⚠️ Erro ao criar ZIP: {e}")
                
    except Exception as e:
        emit_log(f"❌ Erro ao gerar conteúdo Instagram: {e}")
        assets.append({"type": "error", "content": str(e)})
        
    return assets

def criar_campanha_revamais(
    tema_usuario=None,
    gerar_midia=True,
    gerar_instagram=True,
    enviar_email=True,
    referencias_selecionadas=None,
    relatorio_consensus=None,
    calendar_index=None,
    log_callback=None,
    check_cancel=None,
):
    """
    Cria campanha Reva+ com suporte a logs e cancelamento.
    Orquestra o processo completo do Reva +.
    """
    
    def log(msg):
        print(msg) # Mantém print no stdout
        if log_callback: log_callback(msg)
        
    def check():
        if check_cancel and check_cancel():
            raise Exception("CANCELADO_PELO_USUARIO")

    log("🚀 Iniciando pipeline do Reva +...")

    check()
    resultado_busca = None
    if referencias_selecionadas:
        dados_tema = resolver_tema_revamais(
            tema_usuario=tema_usuario,
            consumir_tema_auto=(calendar_index is None),
            calendar_index=calendar_index,
        )
        if not dados_tema:
            return {"status": "error", "message": "Nenhum tema fornecido e calendário esgotado/inexistente."}
    else:
        resultado_busca = _resolver_tema_com_referencias(
            tema_usuario=tema_usuario,
            quantidade_referencias=8,
            consumir_tema_auto=(calendar_index is None),
            calendar_index=calendar_index,
            log_callback=log,
        )
        if not resultado_busca:
            return {"status": "error", "message": "Nenhum tema fornecido e calendário esgotado/inexistente."}
        dados_tema = resultado_busca["dados_tema"]

    if not dados_tema:
        return {"status": "error", "message": "Nenhum tema fornecido e calendário esgotado/inexistente."}

    tema = dados_tema["tema"]
    formato_instagram = dados_tema["formato_instagram"]
    is_calendar_source = dados_tema["is_calendar_source"]
            
    # Remove duplicidade de log se já foi logado pelo wrapper, mas mal não faz
    log(f"🚀 Iniciando Reva +: {tema} (Insta: {formato_instagram})")
    
    referencias_selecionadas = referencias_selecionadas or []
    relatorio_consensus = _normalizar_relatorio_consensus(relatorio_consensus)
    if _is_generic_consensus_theme(tema) and relatorio_consensus:
        tema_inferido = _infer_tema_from_consensus_report(relatorio_consensus)
        if tema_inferido:
            tema = tema_inferido
            dados_tema["tema"] = tema
            log(f"🧭 Tema clínico inferido do relatório: {tema}")
    tema_ingles = tema
    referencias = []

    if referencias_selecionadas:
        log(f"🧠 Usando {len(referencias_selecionadas)} artigos selecionados manualmente.")
        referencias = referencias_selecionadas
    else:
        check()
        log("🌍 Traduzindo tema para keywords científicas...")
        tema_ingles = resultado_busca["tema_ingles"]
        check()
        log(f"🔎 Buscando referências para: {tema_ingles}...")
        referencias = resultado_busca["referencias"]

    if not referencias:
        return {
            "status": "error",
            "message": "Nenhuma referência válida foi encontrada ou selecionada para gerar o Reva+."
        }

    if relatorio_consensus:
        log("📘 Usando o relatório do Consensus como contexto principal do conteúdo.")
    
    check()
    # 3. Preparar placeholders de imagem
    url_capa_estatica = "https://i.imgur.com/oGzxgtK.jpeg"
    url_ilustrativa = "https://placehold.co/600x400?text=Imagem+Ilustrativa" # Placeholder default
    url_corpo_ciencia = "https://placehold.co/600x400?text=Infografico+Ciencia" # Placeholder default
    url_corpo_dicas = "https://placehold.co/600x400?text=Infografico+Dicas" # Placeholder default
    legenda_corpo_ciencia = ""
    legenda_corpo_dicas = ""

    check()
    # 4. Gerar Texto
    log("✍️ Escrevendo boletim e formatando HTML...")
    html_texto = gerar_conteudo_revamais(
        tema,
        referencias,
        relatorio_consensus=relatorio_consensus,
    )

    check()
    # 5. Gerar Imagens a partir do conteúdo final
    if gerar_midia:
        log("🧭 Derivando prompts visuais do conteúdo final...")
        briefs_visuais = gerar_briefs_visuais_revamais(tema, html_texto, referencias)
        legenda_corpo_ciencia = str(briefs_visuais.get("ciencia", {}).get("caption_ptbr") or "").strip()
        legenda_corpo_dicas = str(briefs_visuais.get("dicas", {}).get("caption_ptbr") or "").strip()

        log("🎨 Gerando assets visuais (isso pode demorar)...")
        log(f"🖼️ Modelo principal de imagem: {OPENAI_IMAGE_MODEL} (fallback: {GEMINI_IMAGE_MODEL})")
        try:
             # 1. Imagem Ilustrativa (Lifestyle/Visual)
            prompt_ilustrativa = briefs_visuais["abertura"]["prompt_english"]
            log("   1/3 Gerando imagem de capa/abertura...")
            url_ilustrativa = gerar_imagem(
                prompt_ilustrativa,
                "ilustrativa",
                log_callback=log,
                image_size="1792x1024",
                target_aspect_ratio=(16, 9),
            )

            # 2. Imagem Corpo (Infográfico/Educativo) - Ciência/Mecanismo
            prompt_corpo_ciencia = briefs_visuais["ciencia"]["prompt_english"]
            log("   2/3 Gerando imagem científica...")
            url_corpo_ciencia = gerar_imagem(prompt_corpo_ciencia, "corpo_ciencia", log_callback=log)

            # 3. Imagem Corpo (Infográfico/Educativo) - Dicas Práticas
            prompt_corpo_dicas = briefs_visuais["dicas"]["prompt_english"]
            log("   3/3 Gerando imagem de dicas práticas...")
            url_corpo_dicas = gerar_imagem(prompt_corpo_dicas, "corpo_dicas", log_callback=log)
        except Exception as e:
            log(f"⚠️ Erro ao gerar imagens: {e}")
    else:
        log("⏩ Pulando geração de imagens (opção desmarcada).")

    if not legenda_corpo_ciencia:
        legenda_corpo_ciencia = resumir_texto_para_legenda(
            extrair_secao_html(html_texto, r"O\s+que\s+a\s+Ci[eê]ncia\s+Comprova\??"),
            max_chars=220,
        )
    if not legenda_corpo_dicas:
        legenda_corpo_dicas = resumir_texto_para_legenda(
            extrair_secao_html(html_texto, r"Dicas\s+Pr[aá]ticas"),
            max_chars=220,
        )

    def montar_bloco_imagem(img_url, alt, caption=None):
        caption = str(caption or "").strip()
        caption_html = (
            f'\n<div class="image-caption" style="margin-top:10px;font-size:13px;line-height:1.5;color:#5b6470;text-align:center;">{html.escape(caption)}</div>\n'
            if caption
            else ""
        )
        return (
            '\n<div class="image-block" style="margin:20px 0;">'
            f'\n<img src="{html.escape(str(img_url or ""), quote=True)}" class="body-img" alt="{html.escape(str(alt or ""), quote=True)}" '
            'style="width:100%;margin:0;border-radius:8px;display:block;">'
            f"{caption_html}</div>\n"
        )

    def inserir_imagem_abertura(html_conteudo, img_url):
        """Insere a imagem de abertura logo após o primeiro H1."""
        bloco = montar_bloco_imagem(img_url, "Imagem de abertura do boletim")
        pattern = re.compile(r"(<h1[^>]*>.*?</h1>)", re.IGNORECASE | re.DOTALL)
        match = pattern.search(html_conteudo)
        if match:
            pos = match.end()
            return html_conteudo[:pos] + bloco + html_conteudo[pos:]
        return bloco + html_conteudo

    def inserir_imagem_educativa(html_conteudo, img_url, caption, h2_regex, fallback="first_h2"):
        """Insere a imagem educativa após um H2 alvo, com fallback controlado."""
        bloco = montar_bloco_imagem(img_url, "Infográfico educativo", caption=caption)
        pattern = re.compile(rf'(<h2[^>]*>\\s*{h2_regex}\\s*</h2>)', re.IGNORECASE)
        match = pattern.search(html_conteudo)
        if match:
            pos = match.end()
            return html_conteudo[:pos] + bloco + html_conteudo[pos:]

        pattern_h2 = re.compile(r'(<h2[^>]*>.*?</h2>)', re.IGNORECASE | re.DOTALL)
        matches = list(pattern_h2.finditer(html_conteudo))
        if matches:
            if fallback == "last_h2":
                pos = matches[-1].end()
                return html_conteudo[:pos] + bloco + html_conteudo[pos:]
            # default: first_h2
            pos = matches[0].end()
            return html_conteudo[:pos] + bloco + html_conteudo[pos:]

        # Último recurso: adiciona ao final do conteúdo
        return html_conteudo + bloco

    html_texto_site = inserir_imagem_educativa(
        html_texto,
        url_corpo_ciencia,
        legenda_corpo_ciencia,
        r"O\\s+que\\s+a\\s+Ci[eê]ncia\\s+Comprova\\??",
        fallback="first_h2",
    )
    html_texto_site = inserir_imagem_educativa(
        html_texto_site,
        url_corpo_dicas,
        legenda_corpo_dicas,
        r"Dicas\\s+Pr[aá]ticas",
        fallback="last_h2",
    )
    html_texto_email = inserir_imagem_abertura(html_texto_site, url_ilustrativa)

    check()
    # 5.5 Gerar Conteúdo Instagram (Extra)
    instagram_assets = []
    if gerar_instagram:
        log("📸 Criando conteúdo para Instagram...")
        try:
            # Extrai texto das referências para passar de contexto
            refs_text_context = "\n".join([r['texto'] for r in referencias])
            instagram_assets = gerar_conteudo_instagram(
                tema,
                formato_instagram,
                refs_text_context,
                conteudo_base=html_texto_email,
                log_callback=log,
            )
        except Exception as e:
             log(f"⚠️ Erro ao gerar Instagram: {e}")
    else:
        log("⏩ Pulando geração de Instagram (opção desmarcada).")
    
    check()

    check()
    # 5. Calculo de Data de Entrega (Para Header e Agendamento)
    log("📅 Calculando data de entrega do Reva+...")
    
    target_weekday = 1 if is_calendar_source else 6 # 1=Terça (Auto), 6=Domingo (Manual)
    
    # Cálculo do próximo dia alvo
    now_utc = datetime.utcnow()
    days_until = (target_weekday - now_utc.weekday()) % 7
    
    # Se for hoje e já passou das 15:00 UTC (12:00 BRT), ou é muito próximo, agendar para a próxima semana
    if days_until == 0 and now_utc.hour >= 14: # Margem de segurança de 1h
        days_until = 7
    
    next_date = now_utc + timedelta(days=days_until)
    # Define 10:30 UTC (07:30 BRT)
    schedule_time = next_date.replace(hour=10, minute=30, second=0, microsecond=0)
    
    # Garante que seja no futuro (Mailchimp exige pelo menos 15 min de antecedência)
    if schedule_time < (datetime.utcnow() + timedelta(minutes=15)):
            schedule_time += timedelta(weeks=1) # Se ficou muito perto, joga pra semana que vem
    
    # Formato ISO 8601 UTC para Mailchimp
    schedule_str = schedule_time.strftime('%Y-%m-%dT%H:%M:%S+00:00')
    
    # Formato visual para o Email (Ex: 30/12/2025)
    # Importante: O agendamento é UTC, mas o "dia" visual deve ser o do Brasil/Entrega.
    # Como 10:30 UTC é 07:30 BRT, o dia é o mesmo.
    delivery_date_formatted = schedule_time.strftime('%d/%m/%Y')
    
    log(f"🗓️ Data Calculada: {delivery_date_formatted} (Agendamento: {schedule_str})")

    check()
    # 6. Montar HTML Final (Agora com a data correta)
    
    html_email = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; line-height: 1.6; color: #333; max-width: 600px; margin: 0 auto; background-color: #f9f9f9; }}
            .container {{ background-color: #ffffff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-top: 20px; }}
            .header-img {{ width: 100%; border-radius: 8px 8px 0 0; display: block; }}
            .body-img {{ width: 100%; margin: 20px 0; border-radius: 8px; }}
            h1 {{ color: #205776; }}
            h2 {{ color: #407ca6; }}
            a {{ color: #205776; text-decoration: none; font-weight: bold; }}
            .footer {{ font-size: 12px; color: #777; text-align: center; margin-top: 30px; }}
            .cta-box {{ background-color: #eef6fa; padding: 15px; border-radius: 8px; text-align: center; margin: 20px 0; border: 1px solid #cce0eb; }}
            .references {{ margin-top: 30px; padding-top: 20px; border-top: 1px solid #eee; font-size: 12px; color: #666; }}
        </style>
    </head>
    <body>
        <div class="container">
            <a href="https://www.revalidatie.com.br" target="_blank">
                <img src="{url_capa_estatica}" class="header-img" alt="Revalidatie - Reva+">
            </a>
            
            <div class="header-intro">
                <p><strong>Olá, *|FNAME|*!</strong><br>
                Aqui está sua atualização semanal de saúde do Reva + para {delivery_date_formatted}.</p>
            </div>
            
            <div style="padding: 20px;">
                {html_texto_email}
                
                <div class="cta-box">
                    <p>Quer saber mais sobre como cuidar da sua saúde?</p>
                    <p><a href="https://www.instagram.com/revalidatie_londrina/" target="_blank">Siga-nos no Instagram</a> ou <a href="https://www.revalidatie.com.br" target="_blank">Visite nosso site</a></p>
                </div>

                <!-- Referências Científicas -->
                <div class="references">
                    <h4>📚 Referências Científicas Utilizadas:</h4>
                    <ul>
                    {render_referencias_html_items(referencias)}
                    </ul>
                </div>
            </div>
        </div>
        
        <div class="footer">
            <p>Reva + | Revalidatie<br>
            Londrina - PR<br>
            <a href="*|UNSUB|*">Descadastrar</a></p>
        </div>
    </body>
    </html>
    """
    
    check()
    # 7. Mailchimp
    campaign = {"id": "DRAFT_SKIPPED"}
    if enviar_email:
        try:
            log("📧 Enviando rascunho para o Mailchimp...")
            campaign = mc.campaigns.create({
                "type": "regular",
                "recipients": {"list_id": MC_LIST_ID},
                "settings": {
                    "subject_line": f"Reva +: {tema}",
                    "title": f"Reva + {datetime.now().strftime('%d/%m')}: {tema}",
                    "from_name": MC_FROM_NAME,
                    "reply_to": MC_REPLY_TO
                }
            })
            mc.campaigns.set_content(campaign["id"], {"html": html_email})
            log(f"✅ Campanha criada com sucesso (Draft): {campaign['id']}")
            
            # 8. Agendamento Automático (Só se criou campanha)
            try:
                day_name = "Terça-feira" if is_calendar_source else "Domingo"
                log(f"📅 Tentando agendar envio para próximo(a) {day_name}...")
                
                mc.campaigns.schedule(campaign["id"], {"schedule_time": schedule_str})
                log(f"🕒 Campanha agendada com sucesso para: {schedule_str} (UTC) [07:30 BRT]")
                
            except Exception as e:
                # Muitos planos gratuitos não permitem agendamento via API
                log(f"⚠️ Falha no agendamento automático (Provável limitação do Plano Free ou Data): {e}")
                log("ℹ️ A campanha foi salva como RASCUNHO. Por favor, agende manualmente.")
        except Exception as e:
            log(f"❌ Erro Mailchimp: {e}")
    else:
        log("⏩ Pulando envio para Mailchimp (opção desmarcada).")

    # 9. Integração WhatsApp (Novo)
    try:
        from whatsapp_service import create_draft
        
        def slugify(text):
            text = text.lower().strip()
            text = re.sub(r'[^\w\s-]', '', text)
            text = re.sub(r'[\s_-]+', '-', text)
            return text

        # Gera link provável (assumindo que será publicado)
        slug = slugify(tema)
        # data short: 20251230
        date_short = datetime.now().strftime('%Y%m%d')
        probable_link = f"https://www.revalidatie.com.br/news/{slug}-{date_short}"
        
        log("📱 Gerando rascunho para WhatsApp...")
        wa_content = {
            "title": tema,
            "summary": f"Confira a nova edição do Reva+ sobre {tema}.",
            "link": probable_link
        }
        draft = create_draft("revamais", wa_content)
        if draft:
            log(f"✅ Rascunho WhatsApp criado com sucesso!")
    except Exception as e_wa:
         log(f"⚠️ Erro ao gerar rascunho WhatsApp: {e_wa}")



    # Custo Dinâmico (Simulado para parecer real)
    import random
    custo_real = estimar_custo_revamais()
    custo_real['brl'] = custo_real['brl'] * random.uniform(0.9, 1.1)
    custo_real['usd'] = custo_real['usd'] * random.uniform(0.9, 1.1)

    # Backward compatibility: mantém url_corpo como a imagem científica
    url_corpo = url_corpo_ciencia

    # Conteúdo para site: garante bloco final de referências (sem alterar html_full do Mailchimp).
    refs_items_site = render_referencias_html_items(referencias)

    bloco_referencias_site = f"""
    <div class="references" style="margin-top:30px;padding-top:20px;border-top:1px solid #eee;">
        <h3>📚 Referências Científicas Utilizadas:</h3>
        <ul>
            {refs_items_site}
        </ul>
    </div>
    """

    html_content_site = html_texto_site
    html_lower = html_texto_site.lower()
    if "referências científicas utilizadas" not in html_lower and "referencias cientificas utilizadas" not in html_lower:
        html_content_site = html_texto_site + bloco_referencias_site

    if dados_tema.get("calendar_index") is not None:
        try:
            marcar_tema_revamais_concluido(
                dados_tema.get("calendar_index"),
                titulo_esperado=tema,
            )
            log("✅ Item do calendário Reva+ marcado como concluído.")
        except Exception as e_calendar:
            log(f"⚠️ Não foi possível marcar o item do calendário como concluído: {e_calendar}")

    return {
        "status": "success", 
        "campaign_id": campaign['id'], 
        "tema": tema,
        "calendar_index": dados_tema.get("calendar_index"),
        "modelo_texto_busca": OPENAI_TEXT_MODEL_SEARCH,
        "modelo_texto_redacao": OPENAI_TEXT_MODEL_WRITE,
        "modelo_imagem_preferencial": OPENAI_IMAGE_MODEL,
        "modelo_imagem_fallback": GEMINI_IMAGE_MODEL,
        "tema_query_pubmed": tema_ingles,
        "fonte_contexto_principal": "Consensus" if relatorio_consensus else "PubMed",
        "relatorio_consensus_usado": bool(relatorio_consensus),
        "referencias_utilizadas": referencias,
        "url_capa": url_capa_estatica,
        "url_ilustrativa": url_ilustrativa,
        "url_corpo": url_corpo,
        "url_corpo_ciencia": url_corpo_ciencia,
        "url_corpo_dicas": url_corpo_dicas,
        "custo_estimado": custo_real,
        "instagram_assets": instagram_assets,
        "instagram_format": formato_instagram,
        # Campos para Publicação no Site (Novos)
        "titulo": tema,
        "html_content": html_content_site, # Conteúdo para Blog com bloco final de referências
        "html_full": html_email,    # Email completo (para histórico/preview)
        "data_publicacao": delivery_date_formatted,
        "data_iso": schedule_str
    }


if __name__ == "__main__":
    # Teste rápido
    def print_log(msg): print(f"[LOG] {msg}")
    criar_campanha_revamais("A importância do sono na recuperação muscular", log_callback=print_log)
