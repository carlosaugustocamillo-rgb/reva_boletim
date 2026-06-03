"""
boletim_service.py

Serviço principal que:
- Busca artigos no PubMed (várias consultas)
- Traduz resumos para português usando OpenAI
- Gera boletim principal (HTML) e boletim detalhado (texto)
- Gera roteiro e áudio do podcast (TTS OpenAI + pydub)
- Cria e agenda campanha no Mailchimp
- Salva tudo em um diretório base (PODCAST_BASE_PATH)

CONFIGURAÇÃO POR VARIÁVEIS DE AMBIENTE (obrigatórias):
- OPENAI_API_KEY
- ENTREZ_EMAIL
- MC_API_KEY
- MC_SERVER
- MC_LIST_ID
- MC_FROM_NAME
- MC_REPLY_TO

OPCIONAIS:
- PODCAST_BASE_PATH  (padrão: "/data/podcast")
- INTRO_FILENAME     (padrão: "intro_guto.mp3")
"""

import os
import json
import copy
import uuid
import textwrap
import re
from datetime import datetime, timedelta

import pytz
import requests
from Bio import Entrez
from openai import OpenAI
from pydub import AudioSegment
from mailchimp_marketing import Client
from mailchimp_marketing.api_client import ApiClientError
from dotenv import load_dotenv

# Carrega variáveis de ambiente do arquivo .env
# Carrega variáveis de ambiente do arquivo .env
load_dotenv()

import google.generativeai as genai

# Configuração Gemini
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)


# ======================================================================
# CONFIGURAÇÕES GERAIS
# ======================================================================

# Diretório base: sempre uma pasta "data" dentro do próprio projeto
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(PROJECT_ROOT, "data")
DATA_DIR = BASE_DIR # Alias para compatibilidade
os.makedirs(BASE_DIR, exist_ok=True)

# Subpasta para áudios temporários
AUDIO_DIR = os.path.join(BASE_DIR, "audios")
os.makedirs(AUDIO_DIR, exist_ok=True)

# Arquivo de intro (deixe intro_guto.mp3 dentro da pasta data/)
INTRO_FILENAME = "intro_guto.mp3"
INTRO_PATH = os.path.join(BASE_DIR, INTRO_FILENAME)


# --- OpenAI ---
OPENAI_API_KEY = os.environ["OPENAI_API_KEY"]
client = OpenAI(api_key=OPENAI_API_KEY)

# --- PubMed / Entrez ---
Entrez.email = os.environ["ENTREZ_EMAIL"]

# --- Mailchimp ---
MC_API_KEY = os.environ["MC_API_KEY"]
MC_SERVER = os.environ["MC_SERVER"]
MC_LIST_ID = os.environ["MC_LIST_ID"]
MC_FROM_NAME = os.environ["MC_FROM_NAME"]
MC_REPLY_TO = os.environ["MC_REPLY_TO"]

mc = Client()
mc.set_config({"api_key": MC_API_KEY, "server": MC_SERVER})

from elevenlabs.client import ElevenLabs
from elevenlabs import VoiceSettings

# ... (imports)

# --- ElevenLabs ---
def _bool_env(name, default=False):
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def _int_env(name, default=None):
    raw = os.environ.get(name)
    if raw in (None, ""):
        return default
    try:
        return int(raw)
    except Exception:
        return default


ELEVENLABS_API_KEY = os.environ.get("ELEVENLABS_API_KEY")
ELEVEN_VOICE_ID_HOST = os.environ.get("ELEVEN_VOICE_ID_HOST", "p5oveq8dCbyBIAaD6gzR")
ELEVEN_VOICE_ID_COHOST = os.environ.get("ELEVEN_VOICE_ID_COHOST", "tnSpp4vdxKPjI9w0GnoV")
ELEVEN_AUDIO_DIALOGUE_ENABLED = _bool_env("ELEVEN_AUDIO_DIALOGUE_ENABLED", True)
ELEVEN_AUDIO_DIALOGUE_MODEL = os.environ.get("ELEVEN_AUDIO_DIALOGUE_MODEL", "eleven_v3")
ELEVEN_AUDIO_FALLBACK_MODEL = os.environ.get("ELEVEN_AUDIO_FALLBACK_MODEL", "eleven_multilingual_v2")
ELEVEN_AUDIO_LANGUAGE_CODE = os.environ.get("ELEVEN_AUDIO_LANGUAGE_CODE", "pt")
ELEVEN_AUDIO_DIALOGUE_MAX_CHARS = _int_env("ELEVEN_AUDIO_DIALOGUE_MAX_CHARS", 1800)
ELEVEN_AUDIO_DIALOGUE_SEED = _int_env("ELEVEN_AUDIO_DIALOGUE_SEED")
elevenlabs_client = None
if ELEVENLABS_API_KEY:
    elevenlabs_client = ElevenLabs(api_key=ELEVENLABS_API_KEY)
    print(f"✅ ElevenLabs configurado com vozes:")
    print(f"   HOST: {ELEVEN_VOICE_ID_HOST}")
    print(f"   COHOST: {ELEVEN_VOICE_ID_COHOST}")


# ======================================================================
# FUNÇÕES AUXILIARES (copiadas/adaptadas do seu Colab)
# ======================================================================

import time

def buscar_ids(query):
    """
    Busca TODOS os IDs no PubMed para uma query dos últimos 7 dias.
    """
    tz = pytz.timezone("America/Sao_Paulo")
    agora = datetime.now(tz)
    
    # Últimos 7 dias (Sábado a Sexta), assumindo que hoje é Sexta-feira
    data_final = agora.date()
    data_inicial = data_final - timedelta(days=6)

    mindate = data_inicial.strftime("%Y/%m/%d")
    maxdate = data_final.strftime("%Y/%m/%d")
    print(f"🔎 Buscando artigos de {mindate} até {maxdate} (Sábado a Sexta - 7 dias)")
    print(f"Query: {query}")

    # 1) Primeiro esearch: só pra saber o COUNT
    time.sleep(2) # Rate limit
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=0,
            mindate=mindate,
            maxdate=maxdate,
            datetype="pdat"
        )
    except Exception as e:
        print(f"⚠️ Erro no esearch (tentando de novo em 5s): {e}")
        time.sleep(5)
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=0,
            mindate=mindate,
            maxdate=maxdate,
            datetype="pdat"
        )
    record = Entrez.read(handle)
    total_encontrado = int(record.get("Count", 0))

    if total_encontrado == 0:
        print("Nenhum artigo encontrado para essa query nesse período.")
        return []

    print(f"Encontrados {total_encontrado} artigos. Buscando TODOS os IDs...")

    # 2) Segundo esearch: traz todos os IDs (sem limite artificial)
    time.sleep(2) # Rate limit
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=total_encontrado,
        mindate=mindate,
        maxdate=maxdate,
        datetype="pdat"
    )
    record = Entrez.read(handle)
    return record["IdList"]


def buscar_info_estruturada(ids):
    if not ids:
        return []
    time.sleep(2) # Rate limit
    try:
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="xml", retmode="xml")
    except Exception as e:
        print(f"⚠️ Erro no efetch (tentando de novo em 5s): {e}")
        time.sleep(5)
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="xml", retmode="xml")
    records = Entrez.read(handle)
    artigos = []
    for article in records['PubmedArticle']:
        artigo = {}
        art_data = article['MedlineCitation']['Article']
        artigo['titulo'] = art_data['ArticleTitle']
        artigo['autores'] = [
            a.get('LastName', '') + ' ' + a.get('Initials', '')
            for a in art_data.get('AuthorList', [])
        ]
        
        # Extrai tipos de publicação
        artigo['tipos'] = [str(pt) for pt in art_data.get('PublicationTypeList', [])]
        
        journal_info = art_data['Journal']['JournalIssue']
        artigo['journal'] = art_data['Journal']['Title']
        artigo['ano'] = journal_info['PubDate'].get('Year', 's/ano')
        artigo['volume'] = journal_info.get('Volume', 's/vol')
        artigo['issue'] = journal_info.get('Issue', 's/issue')
        artigo['paginas'] = art_data.get('Pagination', {}).get('MedlinePgn', 's/páginas')
        artigo['pmid'] = article['MedlineCitation']['PMID']
        artigo['doi'] = None
        for id_ in article['PubmedData']['ArticleIdList']:
            if id_.attributes.get('IdType') == 'doi':
                artigo['doi'] = str(id_)
        # AbstractText pode vir como lista de strings ou blocos com Label
        abstract_data = art_data.get('Abstract', {}).get('AbstractText', [])
        resumo_parts = []
        
        if isinstance(abstract_data, list):
            for item in abstract_data:
                if hasattr(item, 'attributes') and 'Label' in item.attributes:
                    label = item.attributes['Label']
                    text = str(item)
                    resumo_parts.append(f"{label}: {text}")
                else:
                    resumo_parts.append(str(item))
            abstract = "\n\n".join(resumo_parts)
        else:
            abstract = str(abstract_data)
            
        artigo['resumo_original'] = abstract or ""
        artigos.append(artigo)
    return artigos


def traduzir_resumo(texto):
    """
    Tradução literal do resumo do inglês para o português do Brasil,
    sem resumir, sem interpretar e sem acrescentar informações.
    """
    if not texto.strip():
        return ""

    prompt = f"""
    Traduza o texto abaixo do inglês para o português do Brasil.

    REGRAS:
    - NÃO resuma, NÃO reestruture, NÃO interprete.
    - Mantenha TODAS as informações presentes no texto original.
    - Preserve o sentido de cada frase, mas não mude a ordem das ideias.
    - Preserve todos os números, porcentagens e termos técnicos exatamente como estão.
    - Não adicione comentários, títulos, seções (como "Objetivos", "Resultados") ou qualquer texto extra.
    - Apenas devolva o texto traduzido, em parágrafos, na mesma ordem do original.

    Texto original (inglês):
    {texto}

    Tradução literal para português do Brasil:
    """

    try:
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        resposta = model.generate_content(prompt)
        return resposta.text.strip()
    except Exception as e:
        print(f"⚠️ Erro Gemini na tradução: {e}. Tentando fallback OpenAI...")
        resposta = client.chat.completions.create(
            model="gpt-4o", # Ajustado para modelo válido
            messages=[
                {
                    "role": "system",
                    "content": "Você é um tradutor científico que faz traduções literais, sem resumir ou interpretar."
                },
                {"role": "user", "content": prompt}
            ],
            temperature=0.0,
        )
        return resposta.choices[0].message.content.strip()


def resumo_para_podcast(titulo, resumo_pt, primeiro_autor, idx=0, is_last=False):
    """
    Gera um roteiro de podcast em formato de CONVERSA entre dois apresentadores,
    com base no RESUMO TRADUZIDO. Retorna uma lista de dicionários com speaker e text.
    """
    import random
    
    # Frases de transição variadas para evitar repetição
    frases_transicao = [
        "Bom, vamos em frente... nesse próximo estudo...",
        "Seguindo nossa pauta de hoje...",
        "Mudando um pouco de assunto, mas ainda dentro da nossa área...",
        "O próximo artigo traz um ponto interessante sobre...",
        "Agora, olha só que curioso esse próximo estudo...",
        "Avançando aqui, temos um trabalho sobre..."
    ]
    
    # Define se é o primeiro estudo ou não para ajustar a transição
    if idx == 0:
        contexto_inicial = "Este é o PRIMEIRO estudo do episódio. O HOST DEVE começar com uma frase de ânimo tipo 'Vamos lá então, pessoal! O nosso 1o estudo de hoje fala sobre...' ou similar."
    else:
        # Escolhe uma frase aleatória (baseada no índice para garantir variação se for re-executado, ou random mesmo)
        frase_escolhida = frases_transicao[idx % len(frases_transicao)] # Usa módulo para ciclar se acabarem as frases
        contexto_inicial = f"Este é o estudo número {idx + 1}. O HOST DEVE começar com uma transição natural. Sugestão: '{frase_escolhida}' (ou similar, mas varie o vocabulário)."
    
    # Define se deve ter despedida no final
    if is_last:
        contexto_final = """Este é o ÚLTIMO estudo. Após discutir o estudo, ENCERRE o episódio SCRIPTADO EXATAMENTE ASSIM (pode adaptar levemente, mas mantenha a essência):
        HOST: "E com isso a gente encerra mais um episódio do RevaCast Weekly. Obrigado pela audiência, e a gente se vê no próximo episódio! Até lá!"
        COHOST: "Até mais pessoal, até a próxima!" """
    else:
        contexto_final = "NÃO finalize o podcast. Deixe a conversa aberta para o próximo estudo."
    
    prompt = f"""
Atue como um roteirista sênior de podcast sobre CIÊNCIA.
Crie um diálogo NATURAL e CONVERSACIONAL entre o HOST (Ivo) e a COHOST (Manu).

TOM DE VOZ: "Professional Casual".
Imagine dois colegas médicos/cientistas conversando no corredor ou num café. Eles são amigos, mas estão discutindo ciência séria.

REGRAS DE ESTILO (CRÍTICO):
1. **Naturalidade**: Use frases fluidas, mas evite gírias excessivas.
2. **Marcadores de Conversa**: Use "Então...", "Pois é...", mas VARIE. Não comece todo estudo igual.
3. **Sem Formalidade Rígida**: Fale como colegas.
4. **Reações Inteligentes**: A Manu deve fazer perguntas pertinentes.
5. **Transição e Encerramento (OBRIGATÓRIO)**:
   - {contexto_inicial}
   - {contexto_final}
6. **Entonação via Pontuação**: Use reticências (...) e exclamações (!).

ESTRUTURA TÉCNICA (Rigor Obrigatório):
- Apresente o estudo: "{titulo}" ({primeiro_autor}).
- Metodologia: Explique o desenho de forma clara.
- Resultados: Diga os números exatos (P-valor, IC), pois o público é médico/técnico.
- Conclusão: Implicação prática.

FORMATO:
Retorne APENAS um JSON array: [{{"speaker": "HOST", "text": "..."}}, {{"speaker": "COHOST", "text": "..."}}]
Use a notação para entonação se a IA de voz suportar, mas foque no TEXTO ser expressivo.

Contexto do estudo:
Título: {titulo}
Primeiro autor: {primeiro_autor}

Resumo traduzido:
{resumo_pt}

FORMATO DE RETORNO (JSON array):
Retorne APENAS um array JSON válido.
"""

    resposta = client.chat.completions.create(
        model="gpt-4o", # ou gpt-5.1
        messages=[
            {
                "role": "system",
                "content": "Você é um roteirista de podcast. Retorne SEMPRE e SOMENTE um JSON array válido com o diálogo."
            },
            {"role": "user", "content": prompt}
        ],
        temperature=0.7,
    )
    
    import json
    try:
        conteudo = resposta.choices[0].message.content.strip()
        # Remove markdown code blocks se existirem
        if conteudo.startswith('```'):
            conteudo = conteudo.split('```')[1]
            if conteudo.startswith('json'):
                conteudo = conteudo[4:]
        conteudo = conteudo.strip()
        
        # Tenta parsear como JSON
        dialogo = json.loads(conteudo)
        
        # Formata abreviações em cada fala
        if isinstance(dialogo, list):
            for fala in dialogo:
                if 'text' in fala:
                    fala['text'] = formatar_abreviacoes(fala['text'])
        
        # Se é um objeto com uma chave, extrai a lista
        if isinstance(dialogo, dict):
            dialogo = dialogo.get('dialogue', dialogo.get('dialog', dialogo.get('conversation', [])))
        
        return dialogo if isinstance(dialogo, list) else []
    except Exception as e:
        print(f"Erro ao parsear diálogo JSON: {e}")
        # Fallback: retorna texto simples como HOST
        return [{"speaker": "HOST", "text": conteudo}]


def normalizar_dialogo(dialogo):
    """
    Normaliza o diálogo para o formato:
    [{"speaker": "HOST|COHOST", "text": "..."}]
    """
    if isinstance(dialogo, dict):
        dialogo = dialogo.get('dialogue', dialogo.get('dialog', dialogo.get('conversation', [])))

    if not isinstance(dialogo, list):
        return []

    dialogo_normalizado = []
    for fala in dialogo:
        speaker = "HOST"
        text = ""

        if isinstance(fala, dict):
            if 'text' in fala:
                speaker = str(fala.get('speaker', 'HOST')).strip().upper()
                text = str(fala.get('text', '')).strip()
            elif 'HOST' in fala:
                speaker = "HOST"
                text = str(fala.get('HOST', '')).strip()
            elif 'COHOST' in fala:
                speaker = "COHOST"
                text = str(fala.get('COHOST', '')).strip()
        elif isinstance(fala, str):
            text = fala.strip()

        if not text:
            continue

        if speaker not in ("HOST", "COHOST"):
            speaker = "HOST"

        dialogo_normalizado.append({
            "speaker": speaker,
            "text": formatar_abreviacoes(text)
        })

    return dialogo_normalizado


def parse_dialogo_json(conteudo):
    """
    Faz parse seguro de texto JSON retornado por LLM para diálogo de podcast.
    """
    if not conteudo:
        return []

    conteudo = conteudo.strip()
    if conteudo.startswith("```"):
        partes = conteudo.split("```")
        if len(partes) >= 2:
            conteudo = partes[1]
            if conteudo.startswith("json"):
                conteudo = conteudo[4:]
        conteudo = conteudo.strip()

    try:
        parsed = json.loads(conteudo)
    except Exception:
        return []

    return normalizar_dialogo(parsed)


def parse_roteiro_completo_json(conteudo):
    """
    Parseia retorno JSON de revisão de ROTEIRO COMPLETO.
    Esperado:
    [
      {"study_index": 1, "title": "...", "dialogue": [...]},
      ...
    ]
    """
    if not conteudo:
        return []

    conteudo = conteudo.strip()
    if conteudo.startswith("```"):
        partes = conteudo.split("```")
        if len(partes) >= 2:
            conteudo = partes[1]
            if conteudo.startswith("json"):
                conteudo = conteudo[4:]
        conteudo = conteudo.strip()

    try:
        parsed = json.loads(conteudo)
    except Exception:
        return []

    if isinstance(parsed, dict):
        parsed = parsed.get("studies", parsed.get("roteiro", []))

    if not isinstance(parsed, list):
        return []

    saida = []
    for item in parsed:
        if not isinstance(item, dict):
            continue
        dialogo = item.get("dialogue", item.get("dialogo", item.get("script", [])))
        saida.append(normalizar_dialogo(dialogo))

    return saida


def revisar_roteiro_completo_para_podcast(roteiros_audio, titulos_estudos):
    """
    Etapa 8.5 revisada:
    Lê TODO o roteiro (todos os estudos em sequência) e ajusta fluidez global,
    principalmente as transições entre estudos.
    """
    roteiro_base = [normalizar_dialogo(d) for d in roteiros_audio if isinstance(d, list)]
    if not roteiro_base:
        return roteiro_base

    payload = []
    for idx, dialogo in enumerate(roteiro_base):
        titulo_atual = titulos_estudos[idx] if idx < len(titulos_estudos) else f"Estudo {idx+1}"
        proximo_titulo = (
            titulos_estudos[idx + 1]
            if idx + 1 < len(titulos_estudos)
            else "Encerramento do episódio"
        )
        payload.append({
            "study_index": idx + 1,
            "title": titulo_atual,
            "next_title": proximo_titulo,
            "dialogue": dialogo
        })

    prompt = f"""
Você é um roteirista profissional especializado em podcasts científicos conversacionais.

Vou te fornecer o roteiro gerado na ETAPA 8 (todos os estudos em sequência). Sua tarefa é produzir a ETAPA 8.5:
refinar o texto para fluidez, transições, ritmo e naturalidade, sem alterar o conteúdo científico.

REGRAS CLARAS:
1) NÃO resuma e NÃO remova fatos, números, resultados ou conclusões.
2) NÃO invente novas evidências e NÃO adicione dados.
3) NÃO use frases mecânicas de relatório, por exemplo:
   - "E no próximo estudo a gente discute..."
   - "O título é..."
   - "O objetivo foi..."
   - "Eles fizeram uma busca de..."
4) Reescreva em tom de conversa real de podcast:
   - reações contextualizadas,
   - perguntas naturais,
   - pausas conceituais curtas ("isso importa porque...", "olha que ponto..."),
   - transições temáticas orgânicas e variadas.
5) Transições entre estudos devem ser coerentes com `next_title`, sem citar tema fora do estudo atual e do próximo.
6) Tom profissional, educativo e conversacional (sem academicismo duro e sem informalidade excessiva).
7) No estudo final, manter encerramento caloroso de podcast.
8) Não use cabeçalhos dentro das falas; a conversa deve soar contínua.
9) Manter ORDEM e QUANTIDADE de estudos.

FORMATO DE SAÍDA (OBRIGATÓRIO):
Retorne SOMENTE JSON válido:
[
  {{
    "study_index": 1,
    "title": "...",
    "dialogue": [{{"speaker":"HOST","text":"..."}},{{"speaker":"COHOST","text":"..."}}]
  }}
]

Roteiro ETAPA 8 para revisar:
{json.dumps(payload, ensure_ascii=False, indent=2)}
"""

    # 1) Tenta Gemini
    if GEMINI_API_KEY:
        try:
            model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
            resposta = model.generate_content(prompt)
            revisado = parse_roteiro_completo_json(getattr(resposta, "text", ""))
            if len(revisado) == len(roteiro_base):
                return revisado
            print("⚠️ Revisão completa Gemini inválida/incompleta. Tentando OpenAI...")
        except Exception as e:
            print(f"⚠️ Erro Gemini na revisão completa: {e}. Tentando fallback OpenAI...")

    # 2) Fallback OpenAI
    try:
        resposta = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {
                    "role": "system",
                    "content": (
                        "Você revisa roteiros completos de podcast científico sem alucinar. "
                        "Retorne somente JSON válido."
                    ),
                },
                {"role": "user", "content": prompt},
            ],
            temperature=0.3,
        )
        conteudo = resposta.choices[0].message.content.strip()
        revisado = parse_roteiro_completo_json(conteudo)
        if len(revisado) == len(roteiro_base):
            return revisado
        print("⚠️ Revisão completa OpenAI inválida/incompleta. Mantendo roteiro original.")
    except Exception as e:
        print(f"⚠️ Erro OpenAI na revisão completa: {e}. Mantendo roteiro original.")

    return roteiro_base


def _strip_code_fence(texto):
    if not texto:
        return ""
    saida = texto.strip()
    if saida.startswith("```"):
        partes = saida.split("```")
        if len(partes) >= 2:
            saida = partes[1]
            if saida.startswith("text"):
                saida = saida[4:]
            elif saida.startswith("markdown"):
                saida = saida[8:]
    return saida.strip()


def _roteiro_para_texto_continuo(roteiros_audio):
    linhas = []
    for dialogo in roteiros_audio or []:
        for fala in normalizar_dialogo(dialogo):
            speaker = fala.get("speaker", "HOST")
            texto = (fala.get("text") or "").strip()
            if texto:
                linhas.append(f"{speaker}: {texto}")
    return "\n".join(linhas)


def _brief_spotify_fallback(titulos_estudos, data_ref):
    linhas = [
        f"Neste episódio do RevaCast Weekly ({data_ref}), Ivo e Manu discutem as evidências mais recentes em oncologia com foco em intervenções não farmacológicas, reabilitação e qualidade de vida.",
        "",
        "Os temas abordados incluem:"
    ]
    for titulo in titulos_estudos[:8]:
        linhas.append(f"- {titulo}")
    linhas.extend([
        "",
        "Um episódio com leitura crítica da evidência e foco em aplicabilidade clínica no cuidado centrado no paciente."
    ])
    return "\n".join(linhas)


def gerar_brief_spotify(roteiros_audio, titulos_estudos, data_ref):
    """
    Gera um brief/resumo do episódio para descrição no Spotify.
    Usa Gemini primeiro e OpenAI como fallback.
    """
    roteiro_texto = _roteiro_para_texto_continuo(roteiros_audio)
    if not roteiro_texto.strip():
        return ""

    prompt = f"""
Você é editor de texto para podcasts científicos em português.

Tarefa:
Produza um resumo/brief para descrição de episódio no Spotify com base EXCLUSIVA no roteiro abaixo.

Regras:
1) Não invente dados, números, desfechos ou conclusões.
2) Não remova os principais achados científicos discutidos.
3) Mantenha tom profissional, educativo e conversacional.
4) Estrutura desejada:
   - 1 parágrafo de abertura (2 a 4 frases).
   - Bloco "Entre os destaques da semana:" com bullets curtos (1 por estudo).
   - 1 fechamento curto convidando a ouvir o episódio.
5) Não usar tabelas, não usar JSON e não usar cabeçalhos técnicos.
6) Texto final em português do Brasil, pronto para colar no Spotify.

Data de referência: {data_ref}
Títulos dos estudos:
{json.dumps(titulos_estudos, ensure_ascii=False, indent=2)}

Roteiro do episódio:
{roteiro_texto}
"""

    # 1) Gemini
    if GEMINI_API_KEY:
        try:
            model = genai.GenerativeModel("gemini-2.5-flash-preview-09-2025")
            resposta = model.generate_content(prompt)
            texto = _strip_code_fence(getattr(resposta, "text", ""))
            if len(texto) >= 180:
                return texto
            print("⚠️ Brief Spotify Gemini curto/inválido. Tentando OpenAI...")
        except Exception as e:
            print(f"⚠️ Erro Gemini ao gerar brief Spotify: {e}. Tentando OpenAI...")

    # 2) OpenAI
    try:
        resposta = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {
                    "role": "system",
                    "content": "Você escreve descrições de episódios científicos para Spotify sem alucinar.",
                },
                {"role": "user", "content": prompt},
            ],
            temperature=0.35,
        )
        texto = _strip_code_fence(resposta.choices[0].message.content.strip())
        if len(texto) >= 180:
            return texto
        print("⚠️ Brief Spotify OpenAI curto/inválido. Usando fallback local.")
    except Exception as e:
        print(f"⚠️ Erro OpenAI ao gerar brief Spotify: {e}. Usando fallback local.")

    return _brief_spotify_fallback(titulos_estudos=titulos_estudos, data_ref=data_ref)


def _eh_fala_transicao(texto):
    if not texto:
        return False
    t = texto.lower()
    pistas = [
        "próximo estudo", "proximo estudo", "no próximo", "no proximo",
        "a seguir", "vamos para", "vamos pro", "vamos ao próximo", "vamos ao proximo",
        "próximo tema", "proximo tema", "puxando esse fio", "na sequência",
        "na sequencia", "gancho natural", "abre caminho", "próximo bloco",
        "proximo bloco", "vamos ver o próximo artigo", "vamos ver o proximo artigo"
    ]
    return any(p in t for p in pistas)


def _eh_transicao_mecanica(texto):
    if not texto:
        return False
    t = texto.lower()
    mecanicas = [
        "e no próximo estudo a gente discute",
        "e no proximo estudo a gente discute",
        "o título é",
        "o titulo é",
        "o objetivo foi",
        "eles fizeram uma busca de",
    ]
    return any(m in t for m in mecanicas)


def _resumir_titulo_para_transicao(titulo):
    t = " ".join((titulo or "").split()).strip().strip(".")
    if not t:
        return "o próximo estudo"
    if ":" in t:
        t = t.split(":", 1)[0].strip()
    if len(t) > 110:
        t = t[:107].rstrip(" ,;:-") + "..."
    return t


def _tema_pt_do_titulo(titulo):
    t = (titulo or "").lower()

    if "peripheral neuropathy" in t or "chemotherapy-induced" in t:
        return "neuropatia periférica induzida por quimioterapia"
    if "advanced lung cancer" in t:
        return "exercício em pessoas com câncer de pulmão avançado"
    if "social support" in t or "older adult survivors of cancer" in t:
        return "apoio social ao exercício em sobreviventes de câncer mais velhos"
    if "pilates" in t and "dance" in t:
        return "pilates e dança na funcionalidade de membros superiores após cirurgia de câncer de mama"
    if "badunjin" in t or "baduanjin" in t:
        return "efeitos do Badunjin sentado no manejo da carga de sintomas"

    titulo_curto = _resumir_titulo_para_transicao(titulo)
    return titulo_curto[0].lower() + titulo_curto[1:] if titulo_curto else "o próximo tema"


def _construir_transicao_organica(proximo_titulo, idx):
    tema = _tema_pt_do_titulo(proximo_titulo).strip().rstrip(".")
    templates = [
        "Isso conversa diretamente com o próximo tema: {tema}.",
        "Puxando esse fio, a gente entra agora em {tema}.",
        "Na sequência, o foco passa para {tema}.",
        "E tem um gancho natural aqui para o que vem agora: {tema}.",
        "Esse ponto abre caminho para o próximo bloco, sobre {tema}.",
    ]
    tpl = templates[idx % len(templates)]
    return tpl.format(tema=tema)


def suavizar_frases_mecanicas(roteiros_audio):
    """
    Suaviza frases de relatório que possam escapar da revisão LLM.
    """
    if not roteiros_audio:
        return roteiros_audio

    regras = [
        (r"^\s*O objetivo foi\b", "A pergunta central desse estudo foi"),
        (r"^\s*Eles fizeram uma busca de\b", "Eles fizeram uma varredura de"),
        (r"^\s*O título é\b", "Esse estudo se chama"),
    ]

    saida = copy.deepcopy(roteiros_audio)
    for i, dialogo in enumerate(saida):
        dialogo_novo = []
        for fala in normalizar_dialogo(dialogo):
            txt = fala.get("text", "")
            for padrao, repl in regras:
                txt = re.sub(padrao, repl, txt, flags=re.IGNORECASE)
            dialogo_novo.append({"speaker": fala.get("speaker", "HOST"), "text": txt})
        saida[i] = dialogo_novo
    return saida


def forcar_transicoes_ancoradas_no_proximo_titulo(roteiros_audio, titulos_estudos):
    """
    Ajuste de segurança pós-LLM:
    garante que cada transição entre estudos aponte para o título real do próximo estudo.
    """
    if not roteiros_audio:
        return roteiros_audio

    resultado = copy.deepcopy(roteiros_audio)
    for idx in range(len(resultado) - 1):
        dialogo = normalizar_dialogo(resultado[idx])
        proximo_titulo = titulos_estudos[idx + 1] if idx + 1 < len(titulos_estudos) else "o próximo estudo"
        texto_transicao = _construir_transicao_organica(proximo_titulo, idx)

        # Remove transições redundantes ou mecânicas no fim do estudo (HOST/COHOST)
        indices_transicao = []
        for j in range(len(dialogo) - 1, max(-1, len(dialogo) - 8), -1):
            if j >= 0 and (_eh_fala_transicao(dialogo[j].get("text", "")) or _eh_transicao_mecanica(dialogo[j].get("text", ""))):
                indices_transicao.append(j)

        # Mantém só uma transição final e orgânica
        for j in sorted(indices_transicao, reverse=True):
            if j >= 0:
                dialogo.pop(j)

        # Adiciona uma única transição final, sempre no HOST (curta e temática)
        # Só adiciona se não for idêntica à última fala HOST já existente.
        ultima_host = ""
        for fala in reversed(dialogo):
            if fala.get("speaker") == "HOST":
                ultima_host = (fala.get("text") or "").strip()
                break
        if ultima_host != texto_transicao.strip():
            dialogo.append({"speaker": "HOST", "text": texto_transicao})

        # Segurança extra: remove duplicidade consecutiva idêntica no final.
        if len(dialogo) >= 2:
            a = dialogo[-1]
            b = dialogo[-2]
            if (
                a.get("speaker") == b.get("speaker") == "HOST"
                and (a.get("text") or "").strip() == (b.get("text") or "").strip()
            ):
                dialogo.pop()

        resultado[idx] = normalizar_dialogo(dialogo)

    return resultado


def salvar_roteiro_txt(caminho, roteiros_audio, titulo):
    """Salva roteiro em formato texto legível."""
    with open(caminho, "w", encoding="utf-8") as f:
        f.write(f"{titulo}\n")
        f.write("=" * 60 + "\n\n")
        for estudo_idx, dialogo in enumerate(roteiros_audio):
            f.write(f"\n{'='*60}\nESTUDO {estudo_idx + 1}\n{'='*60}\n\n")
            if isinstance(dialogo, list):
                for fala in dialogo:
                    f.write(f"{fala.get('speaker')}: {fala.get('text')}\n\n")
            else:
                f.write(f"HOST: {dialogo}\n\n")


def formatar_abreviacoes(texto):
    """
    Formata abreviações comuns adicionando pontos entre as letras
    para o ElevenLabs pronunciar corretamente.
    """
    abreviacoes = {
        r'\bDPOC\b': 'D.P.O.C.',
        r'\bDPI\b': 'D.P.I.',
        r'\bVEF1\b': 'V.E.F.1',
        r'\bVEF\b': 'V.E.F.',
        r'\bCVF\b': 'C.V.F.',
        r'\bFEV1\b': 'F.E.V.1',
        r'\bIMC\b': 'I.M.C.',
        r'\bOMS\b': 'O.M.S.',
        r'\bUSA\b': 'U.S.A.',
        r'\bEUA\b': 'E.U.A.',
        r'\bUK\b': 'U.K.',
        r'\bCOPD\b': 'C.O.P.D.',
    }
    
    import re
    for abrev, formatada in abreviacoes.items():
        texto = re.sub(abrev, formatada, texto, flags=re.IGNORECASE)
    
    return texto


TTS_TRANSICOES_CURTAS = [
    "Isso abre caminho para o proximo estudo da pauta.",
    "Com esse gancho, a gente vai para o proximo estudo.",
    "Na sequencia, a gente passa para o proximo estudo.",
]


def _encurtar_transicao_para_audio(texto):
    base = (texto or "").strip()
    if not base:
        return ""
    idx = sum(ord(c) for c in base) % len(TTS_TRANSICOES_CURTAS)
    return TTS_TRANSICOES_CURTAS[idx]


def preparar_texto_para_audio(texto, speaker="HOST"):
    """
    Prepara o texto somente para TTS.
    Nao altera email, HTML nem o roteiro salvo.
    """
    txt = " ".join(str(texto or "").split()).strip()
    if not txt:
        return ""

    # Transicoes longas, especialmente quando carregam titulo em ingles,
    # soam artificiais e quebram o fluxo no TTS.
    if (_eh_fala_transicao(txt) or _eh_transicao_mecanica(txt)) and (":" in txt or len(txt) > 90):
        return _encurtar_transicao_para_audio(txt)

    substituicoes = [
        (r"\b6MWD\b", "teste de caminhada de seis minutos"),
        (r"\b6MWT\b", "teste de caminhada de seis minutos"),
        (r"\bRCT\b", "ensaio clinico randomizado"),
        (r"\brandomised controlled trial\b", "ensaio clinico randomizado"),
        (r"\brandomized controlled trial\b", "ensaio clinico randomizado"),
        (r"\bUPR\b", "U.P.R."),
        (r"\bIRE1\b", "I.R.E.1"),
        (r"\bPERK\b", "P.E.R.K."),
        (r"\bATF6\b", "A.T.F.6"),
        (r"\bDPOC\b", "D.P.O.C."),
        (r"\bCIPN\b", "C.I.P.N."),
        (r"\bVO2peak\b", "V.O.2 peak"),
        (r"\(neo-\)adjuvant\b", "neo ou adjuvante"),
        (r"\(neo-\)adjuvante\b", "neo ou adjuvante"),
        (r"\(Review\)", "(revisao)"),
    ]
    for padrao, repl in substituicoes:
        txt = re.sub(padrao, repl, txt, flags=re.IGNORECASE)

    if speaker == "COHOST":
        txt = re.sub(r"^Com certeza,\s*Ivo!\s*", "Com certeza, Ivo. ", txt)

    return txt.strip()


def normalizar_dialogo_para_audio(dialogo):
    dialogo_norm = []
    for fala in normalizar_dialogo(dialogo):
        speaker = fala.get("speaker", "HOST")
        text = preparar_texto_para_audio(fala.get("text", ""), speaker=speaker)
        if text:
            dialogo_norm.append({"speaker": speaker, "text": text})
    return dialogo_norm


def _agrupar_dialogo_para_v3(dialogo, limite_chars=1800):
    grupos = []
    atual = []
    chars_atuais = 0

    for fala in dialogo:
        texto = fala.get("text", "")
        tamanho = len(texto)
        if atual and chars_atuais + tamanho > limite_chars:
            grupos.append(atual)
            atual = []
            chars_atuais = 0
        atual.append(fala)
        chars_atuais += tamanho

    if atual:
        grupos.append(atual)

    return grupos


def gerar_dialogo_com_eleven_v3(dialogo, caminho_saida, seed=None):
    """
    Usa o endpoint Text to Dialogue do Eleven v3.
    """
    if not ELEVENLABS_API_KEY:
        raise ValueError("Sem chave ElevenLabs configurada.")

    dialogo_audio = normalizar_dialogo_para_audio(dialogo)
    if not dialogo_audio:
        raise ValueError("Dialogo vazio para audio.")

    grupos = _agrupar_dialogo_para_v3(
        dialogo_audio,
        limite_chars=max(500, ELEVEN_AUDIO_DIALOGUE_MAX_CHARS),
    )
    headers = {
        "xi-api-key": ELEVENLABS_API_KEY,
        "Content-Type": "application/json",
        "Accept": "audio/mpeg",
    }
    endpoint = "https://api.elevenlabs.io/v1/text-to-dialogue?output_format=mp3_44100_128"
    caminhos_partes = []

    for idx, grupo in enumerate(grupos):
        payload = {
            "model_id": ELEVEN_AUDIO_DIALOGUE_MODEL,
            "language_code": ELEVEN_AUDIO_LANGUAGE_CODE,
            "apply_text_normalization": "auto",
            "inputs": [
                {
                    "text": fala["text"],
                    "voice_id": ELEVEN_VOICE_ID_HOST if fala["speaker"] == "HOST" else ELEVEN_VOICE_ID_COHOST,
                }
                for fala in grupo
            ],
        }
        if seed is not None:
            payload["seed"] = (seed + idx) % 4294967295

        response = requests.post(endpoint, headers=headers, json=payload, timeout=180)
        if response.status_code >= 400:
            detalhe = response.text[:500] if response.text else "sem detalhe"
            raise RuntimeError(f"Eleven v3 dialogue falhou ({response.status_code}): {detalhe}")

        if len(grupos) == 1:
            caminho_parte = caminho_saida
        else:
            raiz, ext = os.path.splitext(caminho_saida)
            caminho_parte = f"{raiz}_parte{idx+1}{ext}"

        with open(caminho_parte, "wb") as f:
            f.write(response.content)
        caminhos_partes.append(caminho_parte)

    if len(caminhos_partes) == 1:
        return caminho_saida, dialogo_audio

    combinado = AudioSegment.empty()
    for idx, caminho in enumerate(caminhos_partes):
        combinado += AudioSegment.from_file(caminho, format="mp3")
        if idx < len(caminhos_partes) - 1:
            combinado += AudioSegment.silent(duration=350)

    combinado.export(caminho_saida, format="mp3")
    return caminho_saida, dialogo_audio


def dividir_texto(texto, limite=4096):
    # simples quebra por tamanho aproximado de caracteres
    return textwrap.wrap(texto, width=limite, break_long_words=False, break_on_hyphens=False)


def artigo_tem_exercicio_no_resumo(artigo):
    resumo = artigo.get('resumo_original', '').lower()
    titulo = artigo.get('titulo', '').lower()
    texto_completo = f"{titulo} {resumo}"
    
    # Termos OBRIGATÓRIOS (pelo menos um destes deve estar presente)
    termos_relevantes = [
        'exercise', 'exercício', 'exercising', 
        'rehabilitation', 'reabilitação',
        'physical therapy', 'fisioterapia',
        'training', 'treinamento', # Mantendo 'training' mas com caution (ex: "resistance training")
        'physical activity', 'atividade física',
        'fitness', 'aptidão',
        'gait', 'caminhada', 'walking'
    ]
    
    # Termos de EXCLUSÃO (se tiver isso, provavelmente não é rehab clínica ou é inconclusivo)
    # Ex: "training set" (IA), "training cohort" (estudos de coorte sem intervenção)
    termos_exclusao = [
        'training set', 'validation set', 'training cohort', 
        'medical training', 'resident training', 'surgeon training',
        'exercise of authority', 'exercise of power'
    ]
    
    if any(excl in texto_completo for excl in termos_exclusao):
        return False
        
    return any(term in texto_completo for term in termos_relevantes)


def formatar_artigo_para_html(artigo, resumo_traduzido):
    """Formata metadados + resumo traduzido em HTML simples (sem bullet inteligente)."""

    info_journal = artigo.get('journal', 'N/A')
    info_ano = artigo.get('ano', 'N/A')
    info_volume = f"{artigo.get('volume', '')}" if artigo.get('volume') and 's/vol' not in str(artigo.get('volume')) else ""
    info_issue = f"({artigo.get('issue', '')})" if artigo.get('issue') and 's/issue' not in str(artigo.get('issue')) else ""
    info_paginas = f":{artigo.get('paginas', '')}" if artigo.get('paginas') and 's/páginas' not in str(artigo.get('paginas')) else ""

    publicacao_completa = f"{info_journal}, {info_ano};{info_volume}{info_issue}{info_paginas}"

    autores = ', '.join(artigo['autores']) if artigo['autores'] else "Autores não informados"

    link_doi_url = f"https://doi.org/{artigo['doi']}" if artigo['doi'] else "DOI não disponível"
    link_pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{artigo['pmid']}/"

    links_formatados = f"<strong>DOI:</strong> {link_doi_url}"
    if link_doi_url != "DOI não disponível":
        links_formatados = (
            f'<strong>DOI:</strong> '
            f'<a href="{link_doi_url}" target="_blank" '
            f'style="color:#065e77;text-decoration:underline;">{link_doi_url}</a>'
        )

    links_formatados += (
        f' | <strong>PubMed:</strong> '
        f'<a href="{link_pubmed_url}" target="_blank" '
        f'style="color:#065e77;text-decoration:underline;">{link_pubmed_url}</a>'
    )

    html_meta = f"""
<ul style="margin-left: 0; padding-left: 20px; list-style-type: disc; font-family: Helvetica, Arial, sans-serif; color: #333;">
    <li style="margin-bottom: 5px;"><strong>Título:</strong> {artigo['titulo']}</li>
    <li style="margin-bottom: 5px;"><strong>Autores:</strong> {autores}</li>
    <li style="margin-bottom: 5px;"><strong>Publicação:</strong> {publicacao_completa}</li>
    <li style="margin-bottom: 5px;">{links_formatados}</li>
</ul>
"""

    # Quebra de linhas básica para o resumo traduzido
    texto = resumo_traduzido.replace("\r\n", "\n").replace("\r", "\n")
    paragrafos = [p.strip() for p in texto.split("\n") if p.strip()]
    html_paragrafos = "".join(
        f'<p style="font-family: Helvetica, Arial, sans-serif; color: #333; '
        f'line-height:1.5; margin: 0 0 10px 0;">{p}</p>'
        for p in paragrafos
    )
    if not html_paragrafos:
        html_paragrafos = (
            '<p style="font-family: Helvetica, Arial, sans-serif; color: #333;">'
            'Resumo não disponível.</p>'
        )

    html_final = f"""
<div style="margin-bottom: 25px; padding-bottom: 25px; border-bottom: 1px solid #cccccc;">
    {html_meta}
    {html_paragrafos}
</div>
"""
    return html_final


# ======================================================================
# CONSULTAS (iguais às do seu notebook)
# ======================================================================

CONSULTAS_PRINCIPAIS = [
    'DPOC: ("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab]) AND ("Pulmonary Disease, Chronic Obstructive"[Mesh] OR COPD[tiab] OR "chronic obstructive pulmonary disease"[tiab])',
    'Doenças Intersticiais: ("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab]) AND ("Lung Diseases, Interstitial"[Mesh] OR "interstitial lung disease"[tiab] OR "interstitial lung diseases"[tiab] OR ILD[tiab])',
    'Asma: ("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab]) AND ("Asthma"[Mesh] OR asthma[tiab])',
    'Fibrose Cística: ("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab]) AND ("Cystic Fibrosis"[Mesh] OR "cystic fibrosis"[tiab])',
    'Câncer: ("exercise"[tiab] OR "Exercise"[Mesh]) AND ("Lung Neoplasms"[Mesh] OR "Breast Neoplasms"[Mesh] OR "Colorectal Neoplasms"[Mesh] OR "lung cancer"[tiab] OR "breast cancer"[tiab] OR "colorectal cancer"[tiab] OR "colon cancer"[tiab]) AND (randomized controlled trial[pt] OR controlled clinical trial[pt] OR meta-analysis[pt] OR systematic[sb] OR review[pt])'
]

CONSULTAS_DETALHADAS = {
    "DPOC": (
        '('
        '("Pulmonary Disease, Chronic Obstructive"[Mesh] OR COPD[tiab] OR "chronic obstructive pulmonary disease"[tiab])'
        ' AND '
        '("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab])'
        ' AND '
        '('
        ' randomized controlled trial[pt]'
        ' OR controlled clinical trial[pt]'
        ' OR meta-analysis[pt]'
        ' OR systematic[sb]'
        ' OR review[pt]'
        ')'
        ' AND humans[Mesh]'
        ')'
    ),

    "Asma": (
        '('
        '("Asthma"[Mesh] OR asthma[tiab])'
        ' AND '
        '("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab])'
        ' AND '
        '('
        ' randomized controlled trial[pt]'
        ' OR controlled clinical trial[pt]'
        ' OR meta-analysis[pt]'
        ' OR systematic[sb]'
        ' OR review[pt]'
        ')'
        ' AND humans[Mesh]'
        ')'
    ),

    "Fibrose Cística": (
        '('
        '("Cystic Fibrosis"[Mesh] OR "cystic fibrosis"[tiab])'
        ' AND '
        '("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab])'
        ' AND '
        '('
        ' randomized controlled trial[pt]'
        ' OR controlled clinical trial[pt]'
        ' OR meta-analysis[pt]'
        ' OR systematic[sb]'
        ' OR review[pt]'
        ')'
        ' AND humans[Mesh]'
        ')'
    ),

    "Doenças Intersticiais": (
        '('
        '("Lung Diseases, Interstitial"[Mesh] OR "interstitial lung disease"[tiab] OR "interstitial lung diseases"[tiab] OR ILD[tiab])'
        ' AND '
        '("pulmonary rehabilitation"[tiab] OR "Exercise"[tiab])'
        ' AND '
        '('
        ' randomized controlled trial[pt]'
        ' OR controlled clinical trial[pt]'
        ' OR meta-analysis[pt]'
        ' OR systematic[sb]'
        ' OR review[pt]'
        ')'
        ' AND humans[Mesh]'
        ')'
    ),

    "Câncer": (
        '('
        '("Neoplasms"[Mesh] OR cancer[tiab] OR cancers[tiab])'
        ' AND '
        '("exercise"[tiab])'
        ' AND '
        '('
        ' randomized controlled trial[pt]'
        ' OR controlled clinical trial[pt]'
        ' OR meta-analysis[pt]'
        ' OR systematic[sb]'
        ' OR review[pt]'
        ')'
        ' AND humans[Mesh]'
        ')'
    )
}


# ======================================================================
# MAILCHIMP – TEMPLATE HTML BASE
# ======================================================================

TEMPLATE_HTML_BASE = """
<!DOCTYPE html>
<html lang="pt-br">
<head><meta charset="UTF-8"><title>Boletim Científico Semanal | Revalidatie</title></head>
<body style="background:#fcfdff;margin:0;padding:0;">
  <div style="text-align:center;font-size:14px;margin-top:20px;"><a href="*|ARCHIVE|*" style="color:#065e77;text-decoration:underline;">Ver este e-mail no seu navegador</a></div>
  <div style="text-align:center;margin:30px 0 10px 0;"><a href="https://www.revalidatie.com.br" target="_blank"><img src="https://i.imgur.com/6FIUeHX.png" alt="Logo Revalidatie" style="max-width:270px;width:100%;height:auto;"></a></div>
  <div style="width:100%;max-width:900px;margin:auto;text-align:center;"><img src="https://i.imgur.com/cJdqW3l.png" alt="Boletim Científico Semanal" style="width:100%;max-width:900px;height:auto;border-radius:16px;"></div>
  <div style="width:100%;max-width:900px;margin:22px auto 10px auto;text-align:center;">
    <div style="display:inline-block;text-align:center;background:#ffffff;border-radius:18px;padding:16px 20px;box-shadow:0 4px 10px rgba(0,0,0,0.06);">
      <a href="https://open.spotify.com/show/5NcU5h7u11n5WJDqPS2ZYb?si=Ng42jHrZQK-tQdygp1jf6Q" target="_blank" style="text-decoration:none;">
        <img src="https://i.imgur.com/1b57Ych.png" alt="RevaCast Weekly no Spotify" style="max-width:140px;width:45%;min-width:110px;height:auto;border-radius:26px;display:block;margin:0 auto 10px auto;">
        <div style="font-family:Helvetica,Arial,sans-serif;color:#205776;font-size:1.02em;font-weight:bold;margin-bottom:6px;">
          Prefere ouvir esse boletim?
        </div>
        <div style="font-family:Helvetica,Arial,sans-serif;color:#556;font-size:0.96em;margin-bottom:12px;">
          Clique abaixo para acessar o episódio no Spotify.
        </div>
        <span style="display:inline-block;background:#205776;color:#fff;padding:9px 26px;border-radius:999px;font-size:0.98em;font-weight:bold;">
          Ouvir no Spotify
        </span>
      </a>
    </div>
  </div>
  <table align="center" border="0" cellpadding="0" cellspacing="0" width="92%" style="max-width:760px; margin:auto; background:#fff;">
    <tr><td style="padding: 36px 20px 0 20px; text-align: center;"><h1 style="margin:0 0 10px 0;font-size:2.4em;color:#205776;font-family:Helvetica,Arial,sans-serif;font-weight:bold;">Olá, *|FNAME|*!</h1></td></tr>
    <tr>
      <td style="padding:0 20px 36px 20px; text-align:left;">
        <span style="color:#407ca6;font-size:1.07em;font-family:Helvetica,Arial,sans-serif;">Segue abaixo os principais destaques da semana na literatura científica.</span>
        <div style="height:18px;"></div>
        <div style="font-size:1.1em; color:#111;font-family:Helvetica,Arial,sans-serif;">{conteudo_aqui}</div>
        <div style="height:34px;"></div>
      </td>
    </tr>
  </table>
  <div style="background:#222C36;color:#fff;padding:36px 0 24px 0;text-align:center;font-size:1.08em;font-family:Helvetica,Arial,sans-serif;">
    <div style="margin-bottom:12px;"><img src="https://i.imgur.com/6FIUeHX.png" alt="Revalidatie" style="max-width:180px;width:100%;height:auto;"></div>
    <div style="margin-bottom:15px;">Copyright (C ) 2025 Revalidatie. Todos os direitos reservados.
      <br>Você está recebendo este email porque se inscreveu para receber atualizações científicas da Revalidatie.</div>
    <div class="disclaimer" style="color:#ddd;font-size:0.96em;">Quer alterar como recebe estes emails?
      <a href="*|UPDATE_PROFILE|*" style="color:#5beaff;">Atualizar preferências</a> ou <a href="*|UNSUB|*" style="color:#5beaff;">descadastrar</a></div>
  </div>
</body>
</html>
"""


# ======================================================================
# FUNÇÃO PRINCIPAL: rodar_boletim()
# ======================================================================

def rodar_boletim(opcoes=None):
    """
    Executa o pipeline conforme as opções selecionadas.
    opcoes: dict com chaves booleanas:
      - 'resumos': Busca artigos e gera textos (Principal e Detalhado)
      - 'roteiro': Gera o roteiro do podcast (texto)
      - 'revisao_roteiro': Revisa fluidez do roteiro (etapa 8.5)
      - 'brief_spotify': Gera resumo final para descrição do episódio no Spotify
      - 'audio': Gera o áudio (ElevenLabs)
      - 'mailchimp': Cria e agenda campanha
      - 'firebase': Upload e RSS
    """
    if opcoes is None:
        opcoes = {
            'resumos': True,
            'roteiro': True,
            'revisao_roteiro': True,
            'brief_spotify': True,
            'audio': True,
            'mailchimp': True,
            'firebase': True
        }

    yield f"🚀 Iniciando pipeline com opções: {opcoes}"

    hoje = datetime.today().strftime('%Y-%m-%d')

    # Paths principais
    boletim_path = os.path.join(BASE_DIR, f"boletim_pubmed_{hoje}.txt")
    revisao_path = os.path.join(BASE_DIR, f"boletim_para_revisao_{hoje}.txt")
    boletim_detalhado_path = os.path.join(BASE_DIR, f"boletim_detalhado_{hoje}.txt")
    roteiro_path = os.path.join(BASE_DIR, f"roteiro_podcast_{hoje}.txt")
    base_episodio_name = f"episodio_boletim_{hoje}"
    episodio_filename = f"{base_episodio_name}.mp3"
    episodio_path = os.path.join(AUDIO_DIR, episodio_filename) # Usando AUDIO_DIR para manter organizado
    brief_spotify_path = os.path.join(BASE_DIR, f"brief_spotify_{hoje}.txt")
    if not os.path.exists(AUDIO_DIR): os.makedirs(AUDIO_DIR, exist_ok=True)
    
    # Versionamento: se já existe, cria _1, _2...
    counter = 1
    while os.path.exists(episodio_path):
        episodio_filename = f"{base_episodio_name}_{counter}.mp3"
        episodio_path = os.path.join(AUDIO_DIR, episodio_filename)
        counter += 1

    # Variáveis de estado para passar entre etapas
    roteiros_audio = []
    roteiros_audio_original = []
    titulos_podcast = []
    brief_spotify_text = ""
    total_chars_elevenlabs = 0
    
    # ------------------------------------------------------------------
    # 1) BOLETIM PRINCIPAL & DETALHADO (RESUMOS)
    # ------------------------------------------------------------------
    if opcoes.get('resumos'):
        yield "🔎 1/5: Buscando artigos no PubMed e gerando Resumos..."
        
        # --- BOLETIM PRINCIPAL ---
        boletim_final = (
            '<p style="font-family: Helvetica, Arial, sans-serif; color: #333; line-height: 1.5;">'
            "Os dados a seguir mostram os estudos publicados no PubMed na última semana "
            "(sábado a sexta), com os resumos traduzidos literalmente do inglês para o português."
            "</p>"
        )
        boletim_revisao = (
            "Os artigos a seguir foram encontrados no PubMed, mas não continham resumo disponível "
            "ou houve falha técnica na tradução. Eles requerem revisão manual.\n\n"
        )
        artigos_vistos = set()

        for consulta in CONSULTAS_PRINCIPAIS:
            titulo_secao, query = consulta.split(": ", 1)
            yield f"   - Processando seção: {titulo_secao}..."
            
            # 1. Busca
            yield f"   🔎 Buscando artigos no PubMed para '{titulo_secao}'..."
            ids = buscar_ids(query)
            
            if not ids:
                print(f"   Nenhum artigo encontrado para '{titulo_secao}'.")
                continue
                
            # Pega detalhes
            yield f"   📄 Baixando detalhes de {len(ids)} artigos para '{titulo_secao}'..."
            artigos = buscar_info_estruturada(ids)
            artigos_unicos = [a for a in artigos if a['pmid'] not in artigos_vistos]

            if not artigos_unicos:
                continue

            cabecalho_secao_html = f"""
<h2 style="font-family: Helvetica, Arial, sans-serif; font-size: 22px; font-weight: bold; color: #205776; border-bottom: 2px solid #205776; padding-bottom: 5px; margin-top: 30px; margin-bottom: 15px;">
    {titulo_secao}
</h2>
"""
            boletim_final += cabecalho_secao_html
            boletim_revisao += f"### {titulo_secao}\n\n"

            for art in artigos_unicos:
                artigos_vistos.add(art['pmid'])
                autores_str = ', '.join(art['autores']) if art['autores'] else "Autores não informados"
                link_pubmed = f"https://pubmed.ncbi.nlm.nih.gov/{art['pmid']}/"
                link_doi = f"https://doi.org/{art['doi']}" if art['doi'] else "DOI não disponível"

                info_basica_artigo = f"""* {art['titulo']}
* {autores_str}
* {art['journal']}, {art['ano']}; {art['volume']}({art['issue']} ): {art['paginas']}
* DOI: {link_doi}
* PubMed: {link_pubmed}
"""
                resumo_original = art.get('resumo_original', '').strip()

                if not resumo_original:
                    boletim_revisao += (
                        info_basica_artigo
                        + "\nMotivo da falha: Resumo ausente no PubMed.\n\n---\n\n"
                    )
                    continue

                try:
                    resumo_traduzido = traduzir_resumo(resumo_original)
                    html_formatado = formatar_artigo_para_html(art, resumo_traduzido)
                    boletim_final += html_formatado
                except Exception as e:
                    print(f"❌ Erro na tradução para o PMID {art['pmid']}: {e}")
                    boletim_revisao += (
                        info_basica_artigo
                        + f"\nMotivo da falha: Erro na chamada da API de tradução - {e}\n\n---\n\n"
                    )

        boletim_final += (
            '<p style="font-family: Helvetica, Arial, sans-serif; color: #333; line-height: 1.5; margin-top: 30px;">'
            "Espero que estas traduções sejam úteis. Siga nosso podcast para mais!<br>"
            "Abraços,<br>Guto"
            "</p>"
        )

        with open(boletim_path, "w", encoding="utf-8") as f:
            f.write(boletim_final)
        print(f"✅ Boletim principal salvo como: {boletim_path}")

        with open(revisao_path, "w", encoding="utf-8") as f:
            f.write(boletim_revisao)
        print(f"✅ Boletim para revisão salvo como: {revisao_path}")

        # --- BOLETIM DETALHADO ---
        yield "📝 Gerando Boletim Detalhado..."
        boletim_detalhado = (
            "Boletim detalhado com os estudos mais relevantes sobre exercício em diferentes "
            "condições crônicas, publicados na última semana (sábado a sexta). "
            "Os resumos abaixo são traduções literais do PubMed.\n\n"
        )
        
        todos_artigos_relevantes = []
        artigos_vistos_podcast = set()

        for tema, query in CONSULTAS_DETALHADAS.items():
            boletim_detalhado += f"## {tema}\n\n"
            ids = buscar_ids(query)
            artigos = buscar_info_estruturada(ids)

            artigos_relevantes = [a for a in artigos if artigo_tem_exercicio_no_resumo(a)]
            if not artigos_relevantes:
                boletim_detalhado += "Nenhum estudo relevante encontrado nesta semana.\n\n"
                continue

            for idx, art in enumerate(artigos_relevantes):
                autores_str = ', '.join(art['autores']) if art['autores'] else "Autores não informados"
                link_pubmed = f"https://pubmed.ncbi.nlm.nih.gov/{art['pmid']}/"
                link_doi = f"https://doi.org/{art['doi']}" if art['doi'] else "DOI não disponível"
                resumo_original = art.get('resumo_original', '').strip()

                if not resumo_original:
                    continue # Simplificado para brevidade

                try:
                    resumo_traduzido = traduzir_resumo(resumo_original)
                except Exception as e:
                    continue

                boletim_detalhado += f"""* {art['titulo']}
* {autores_str}
* {art['journal']}, {art['ano']}; {art['volume']}({art['issue']}): {art['paginas']}
* DOI: {link_doi}
* PubMed: {link_pubmed}

{resumo_traduzido}

---

"""
                # [FIX]: Garante que o mesmo artigo não seja adicionado duas vezes (ex: se caiu em DPOC e ASMA)
                if art['pmid'] not in artigos_vistos_podcast:
                    artigos_vistos_podcast.add(art['pmid'])
                    primeiro_autor = art['autores'][0] if art['autores'] else "Autor não identificado"
                    todos_artigos_relevantes.append({
                        'titulo': art['titulo'],
                        'resumo_traduzido': resumo_traduzido,
                        'primeiro_autor': primeiro_autor,
                        'tipos': art.get('tipos', []), # Importante passar metadados
                        'resumo_original': art.get('resumo_original', '')
                    })
        
        boletim_detalhado += "\nCompartilhe com colegas. RevaCast Pesquisa Detalhada!"

        with open(boletim_detalhado_path, "w", encoding="utf-8") as f:
            f.write(boletim_detalhado)
        print(f"✅ Boletim detalhado salvo como: {boletim_detalhado_path}")
        
        # --- ROTEIRO (Parte do passo de texto, mas opcional) ---
        if opcoes.get('roteiro'):
            yield "📝 Gerando Roteiros de Podcast..."
            
            # Filtra para o Podcast: Apenas RCT, Systematic Review, Meta-Analysis, Guidelines
            # E que tenham resumo traduzido disponível
            tipos_podcast = [
                'Randomized Controlled Trial', 
                'Systematic Review', 
                'Meta-Analysis', 
                'Practice Guideline',
                'Guideline',
                'Clinical Trial',
                'Review',
                'Observational Study',
                'Comparative Study'
            ]
            
            artigos_podcast = []
            for art in todos_artigos_relevantes:
                if not art.get('resumo_traduzido'): continue
                
                # Verifica se tem algum dos tipos aceitos NOS METADADOS
                tipos_artigo = art.get('tipos', [])
                eh_alta_evidencia = any(t in tipos_artigo for t in tipos_podcast)
                
                # FALLBACK: Se não achou nos metadados (comum em artigos muito recentes),
                # procura palavras-chave no TÍTULO ou RESUMO ORIGINAL
                if not eh_alta_evidencia:
                    texto_completo = (art.get('titulo', '') + ' ' + art.get('resumo_original', '')).lower()
                    termos_chave = [
                        'randomized', 'randomised', 'controlled trial', 
                        'systematic review', 'meta-analysis', 'guideline',
                        'consensus', 'position statement'
                    ]
                    eh_alta_evidencia = any(termo in texto_completo for termo in termos_chave)
                
                if eh_alta_evidencia:
                    artigos_podcast.append(art)
            
            # Ordenação Inteligente por Nível de Evidência
            def get_evidence_score(art):
                tipos = [t.lower() for t in art.get('tipos', [])]
                texto = (art.get('titulo', '') + ' ' + art.get('resumo_original', '')).lower()
                
                # Peso 3: Guidelines e Consensos
                if any(t in tipos for t in ['practice guideline', 'guideline', 'consensus development conference']) or \
                   any(w in texto for w in ['guideline', 'consensus statement', 'position paper']):
                    return 3
                
                # Peso 2: Meta-análises e Revisões Sistemáticas
                if any(t in tipos for t in ['meta-analysis', 'systematic review']) or \
                   any(w in texto for w in ['meta-analysis', 'systematic review']):
                    return 2
                
                # Peso 1: RCTs
                if any(t in tipos for t in ['randomized controlled trial', 'controlled clinical trial']) or \
                   any(w in texto for w in ['randomized', 'randomised', 'controlled trial']):
                    return 1
                
                return 0

            # Ordena: Maior score primeiro. Desempate pela ordem original (que é data descrescente no PubMed)
            artigos_podcast.sort(key=get_evidence_score, reverse=True)

            # LIMITA A 6 ESTUDOS (aprox 20-25 min de áudio)
            LIMIT_PODCAST = 6
            artigos_cortados = []
            
            if len(artigos_podcast) > LIMIT_PODCAST:
                artigos_cortados = artigos_podcast[LIMIT_PODCAST:]
                artigos_podcast = artigos_podcast[:LIMIT_PODCAST]

            # RELATÓRIO DE CURADORIA
            print("\n" + "="*50)
            print("🎙️ RELATÓRIO DE CURADORIA DO PODCAST")
            print("="*50)
            
            print(f"\n✅ SELECIONADOS ({len(artigos_podcast)}):")
            for art in artigos_podcast:
                tipos = ", ".join(art.get('tipos', [])[:2]) # Mostra só os 2 primeiros tipos
                print(f"   - {art.get('titulo', '')[:80]}... [{tipos}]")
                
            if artigos_cortados:
                print(f"\n❌ CORTADOS PELO LIMITE ({len(artigos_cortados)}):")
                for art in artigos_cortados:
                    tipos = ", ".join(art.get('tipos', [])[:2])
                    print(f"   - {art.get('titulo', '')[:80]}... [{tipos}]")
            
            ignored_count = len(todos_artigos_relevantes) - len(artigos_podcast) - len(artigos_cortados)
            print(f"\n⚠️ IGNORADOS (Baixa evidência/Outros): {ignored_count} estudos.")
            print("="*50 + "\n")

            # FALLBACK DE SEGURANÇA: Se não sobrou nada (muito restrito), pega os top 3 gerais
            if not artigos_podcast and todos_artigos_relevantes:
                print("⚠️ Nenhum estudo de alta evidência encontrado. Usando fallback (Top 3 gerais).")
                artigos_podcast = todos_artigos_relevantes[:3]

            roteiros_audio = []
            titulos_podcast = []
            if artigos_podcast:
                yield f"🎙️ Gerando roteiro para {len(artigos_podcast)} estudos selecionados..."
                
                for idx, art in enumerate(artigos_podcast):
                    is_last = (idx == len(artigos_podcast) - 1)
                    yield f"   - Roteirizando estudo {idx+1}/{len(artigos_podcast)}: {art.get('titulo', 'Sem título')[:30]}..."
                    
                    autores_list = art.get('autores', [])
                    primeiro_autor = autores_list[0] if autores_list else "Autor desconhecido"
                    
                    dialogo = resumo_para_podcast(
                        titulo=art.get('titulo', 'Sem título'),
                        resumo_pt=art.get('resumo_traduzido', ''),
                        primeiro_autor=primeiro_autor,
                        idx=idx,
                        is_last=is_last
                    )
                    roteiros_audio.append(normalizar_dialogo(dialogo))
                    titulos_podcast.append(art.get('titulo', f'Estudo {idx+1}'))

            # Guarda versão da etapa 8 (antes da revisão de fluidez)
            roteiros_audio_original = copy.deepcopy(roteiros_audio)

            # ------------------------------------------------------------------
            # 8.5) REVISÃO DE FLUIDEZ DO ROTEIRO (GEMINI/OPENAI)
            # ------------------------------------------------------------------
            if roteiros_audio and opcoes.get('revisao_roteiro', True):
                yield "🧠 2.5/5: Revisando roteiro COMPLETO para soar mais podcast..."
                roteiros_audio = revisar_roteiro_completo_para_podcast(
                    roteiros_audio=roteiros_audio,
                    titulos_estudos=titulos_podcast
                )
                yield "🧹 Suavizando frases mecânicas residuais..."
                roteiros_audio = suavizar_frases_mecanicas(roteiros_audio)
                yield "🧭 Ajustando transições orgânicas entre estudos..."
                roteiros_audio = forcar_transicoes_ancoradas_no_proximo_titulo(
                    roteiros_audio=roteiros_audio,
                    titulos_estudos=titulos_podcast
                )
            elif roteiros_audio:
                yield "⏭️ Etapa 8.5 desativada (usando roteiro bruto da etapa 8)."

            # Salva TXT original (etapa 8) e revisado (etapa 8.5)
            roteiro_original_txt = os.path.join(BASE_DIR, f"roteiro_podcast_{hoje}_etapa8_original.txt")
            roteiro_revisado_txt = os.path.join(BASE_DIR, f"roteiro_podcast_{hoje}_etapa8_5_revisado.txt")

            salvar_roteiro_txt(
                caminho=roteiro_original_txt,
                roteiros_audio=roteiros_audio_original,
                titulo="ROTEIRO COMPLETO DO PODCAST - RevaCast Weekly (ETAPA 8 - ORIGINAL)"
            )
            salvar_roteiro_txt(
                caminho=roteiro_revisado_txt,
                roteiros_audio=roteiros_audio,
                titulo="ROTEIRO COMPLETO DO PODCAST - RevaCast Weekly (ETAPA 8.5 - REVISADO)"
            )

            # Mantém compatibilidade: arquivo "roteiro_path" sempre recebe a versão revisada (8.5)
            salvar_roteiro_txt(
                caminho=roteiro_path,
                roteiros_audio=roteiros_audio,
                titulo="ROTEIRO COMPLETO DO PODCAST - RevaCast Weekly"
            )
            print(f"✅ Roteiro etapa 8 salvo em: {roteiro_original_txt}")
            print(f"✅ Roteiro etapa 8.5 salvo em: {roteiro_revisado_txt}")
            print(f"✅ Roteiro ativo para áudio salvo em: {roteiro_path}")
            
            # Salva também em JSON estruturado (etapa 8, etapa 8.5 e alias compatível)
            roteiro_json_dir = os.path.join(BASE_DIR, "roteiros")
            roteiro_json_path = os.path.join(roteiro_json_dir, f"roteiro_estruturado_{hoje}.json")
            roteiro_json_etapa8 = os.path.join(roteiro_json_dir, f"roteiro_estruturado_{hoje}_etapa8.json")
            roteiro_json_etapa85 = os.path.join(roteiro_json_dir, f"roteiro_estruturado_{hoje}_etapa8_5.json")
            try:
                os.makedirs(roteiro_json_dir, exist_ok=True)
                with open(roteiro_json_etapa8, "w", encoding="utf-8") as f:
                    json.dump(roteiros_audio_original, f, indent=2, ensure_ascii=False)
                with open(roteiro_json_etapa85, "w", encoding="utf-8") as f:
                    json.dump(roteiros_audio, f, indent=2, ensure_ascii=False)
                with open(roteiro_json_path, "w", encoding="utf-8") as f:
                    json.dump(roteiros_audio, f, indent=2, ensure_ascii=False)

                print(f"✅ Roteiro JSON etapa 8 salvo em: {roteiro_json_etapa8}")
                print(f"✅ Roteiro JSON etapa 8.5 salvo em: {roteiro_json_etapa85}")
                print(f"✅ Roteiro JSON ativo salvo em: {roteiro_json_path}")
                
                # Upload para Firestore se habilitado
                if opcoes.get('firebase'):
                    from firebase_service import save_firestore_document
                    doc_data = {
                        "date": hoje,
                        "script": roteiros_audio,
                        "script_original": roteiros_audio_original,
                        "created_at": datetime.now().isoformat()
                    }
                    save_firestore_document("roteiros", f"roteiro_{hoje}", doc_data)
                    
            except Exception as e:
                print(f"❌ Erro ao salvar JSON do roteiro: {e}")

            if roteiros_audio and opcoes.get('brief_spotify', True):
                yield "🧾 2.8/5: Gerando brief do episódio para Spotify..."
                brief_spotify_text = gerar_brief_spotify(
                    roteiros_audio=roteiros_audio,
                    titulos_estudos=titulos_podcast,
                    data_ref=hoje
                ).strip()

                if brief_spotify_text:
                    with open(brief_spotify_path, "w", encoding="utf-8") as f:
                        f.write(f"Resumo Spotify - RevaCast Weekly ({hoje})\n")
                        f.write("=" * 50 + "\n\n")
                        f.write(brief_spotify_text + "\n")
                    yield f"✅ Brief Spotify salvo: {brief_spotify_path}"
                else:
                    yield "⚠️ Brief Spotify não foi gerado (texto vazio)."
        else:
            yield "⏭️ Pulando geração de Roteiro."

    else:
        yield "⏭️ Pulando busca e geração de Resumos (usando arquivos existentes se houver)."
        # Tenta carregar roteiro existente se necessário para o áudio
        if opcoes.get('audio') and os.path.exists(roteiro_path):
             # Lógica simplificada: ler o roteiro do arquivo seria complexo de parsear de volta para JSON.
             # Por enquanto, assumimos que se pulou resumos, não tem roteiro em memória.
             # Se o usuário quiser gerar áudio sem gerar roteiro, precisaria carregar do disco.
             # Para simplificar: Se pulou resumos/roteiro, não gera áudio NOVO, apenas usa existente.
             pass

    # ------------------------------------------------------------------
    # 3) GERAÇÃO DO ÁUDIO
    # ------------------------------------------------------------------
    # Safety check: Se Firebase estiver ativado mas sem credenciais, não gera áudio
    if opcoes.get('audio') and opcoes.get('firebase'):
        from firebase_service import is_firebase_ready
        if not is_firebase_ready():
            yield "⚠️ AVISO: Firebase não configurado (falta credenciais). Desativando áudio para economizar créditos."
            opcoes['audio'] = False

    if opcoes.get('audio'):
        yield "🎙️ 3/5: Verificando áudio..."
        total_chars_elevenlabs = 0
        
        if os.path.exists(episodio_path):
             yield f"⚠️ Áudio já encontrado! (Isso não deve acontecer com versionamento)"
             
        elif roteiros_audio:
            yield "🎙️ Gerando Áudio (ElevenLabs)..."
            audio_paths = []
            
            # --- GERAÇÃO DA ABERTURA FALADA (IVO E MANU) ---
            import random
            yield "   - Gerando apresentação dos hosts..."
            
            POOL_ABERTURAS = [
                # Opção 1 (Sugerida)
                [
                    {"speaker": "HOST", "text": "Olááá! Sejam muito bem-vindos a mais um episódio do RevaCast Weekly. Eu sou o Ivo..."},
                    {"speaker": "COHOST", "text": "E eu sou a Manu."},
                    {"speaker": "HOST", "text": "Juntos a gente dá uma passeada na literatura científica da última semana. Bora lá, Manu?"},
                    {"speaker": "COHOST", "text": "Bora!"}
                ],
                # Opção 2 (Mais direta)
                [
                    {"speaker": "HOST", "text": "Fala pessoal! Começando mais uma edição do RevaCast Weekly, o seu resumo de ciência. Aqui é o Ivo."},
                    {"speaker": "COHOST", "text": "E aqui é a Manu. Tudo pronto para as atualizações desta semana."},
                    {"speaker": "HOST", "text": "Exato. Separamos os artigos mais importantes para discutir. Vamos nessa?"},
                    {"speaker": "COHOST", "text": "Com certeza, vamos lá!"}
                ],
                # Opção 3 (Mais energética)
                [
                    {"speaker": "HOST", "text": "Sejam bem-vindos ao RevaCast Weekly! Eu sou o Ivo e já estou com os estudos na mão."},
                    {"speaker": "COHOST", "text": "Oi gente, eu sou a Manu! Vamos descomplicar as evidências da semana?"},
                    {"speaker": "HOST", "text": "É isso aí. Sem enrolação, vamos ver o que saiu de novo."},
                    {"speaker": "COHOST", "text": "Partiu!"}
                ]
            ]
            
            abertura_escolhida = random.choice(POOL_ABERTURAS)
            audios_abertura_temp = []
            path_abertura_final = os.path.join(AUDIO_DIR, f"intro_falada_{hoje}.mp3")
            
            try:
                usou_dialogue_v3 = False

                if ELEVEN_AUDIO_DIALOGUE_ENABLED:
                    try:
                        gerar_dialogo_com_eleven_v3(
                            abertura_escolhida,
                            path_abertura_final,
                            seed=ELEVEN_AUDIO_DIALOGUE_SEED,
                        )
                        usou_dialogue_v3 = True
                    except Exception as e_v3_intro:
                        print(f"⚠️ Eleven v3 na abertura falhou, usando fallback: {e_v3_intro}")

                if not usou_dialogue_v3:
                    for idx_intro, fala in enumerate(normalizar_dialogo_para_audio(abertura_escolhida)):
                        voice_id = ELEVEN_VOICE_ID_HOST if fala['speaker'] == 'HOST' else ELEVEN_VOICE_ID_COHOST
                        audio_gen = elevenlabs_client.text_to_speech.convert(
                            voice_id=voice_id,
                            text=fala['text'],
                            model_id=ELEVEN_AUDIO_FALLBACK_MODEL,
                            voice_settings=VoiceSettings(
                                stability=0.75,
                                similarity_boost=0.75,
                                style=0.0,
                                use_speaker_boost=True,
                            )
                        )
                        path_intro = os.path.join(AUDIO_DIR, f"temp_intro_{idx_intro}.mp3")
                        with open(path_intro, "wb") as f:
                            for chunk in audio_gen:
                                f.write(chunk)
                        audios_abertura_temp.append(path_intro)

                    abertura_combinada = AudioSegment.empty()
                    for p in audios_abertura_temp:
                        abertura_combinada += AudioSegment.from_file(p) + AudioSegment.silent(duration=300)
                    abertura_combinada.export(path_abertura_final, format="mp3")

            except Exception as e:
                print(f"Erro ao gerar abertura: {e}")

            # --- FIM DA ABERTURA ---
            
            # Gera ID único para essa execução para não misturar arquivos temp
            run_id = str(uuid.uuid4())[:8]
            
            # Se for string (erro de parse anterior), tenta corrigir ou pula
            if isinstance(roteiros_audio, str):
                try:
                    roteiros_audio = json.loads(roteiros_audio)
                except:
                    yield "⚠️ Erro ao ler roteiro JSON (formato inválido)."
                    return
            
            # Adiciona o caminho da abertura NA LISTA DE CAMINHOS se existir
            # Adiciona o caminho da abertura NA LISTA DE CAMINHOS se existir
            if path_abertura_final and os.path.exists(path_abertura_final):
                audio_paths.append(path_abertura_final)

            for estudo_idx, dialogo in enumerate(roteiros_audio):
                if not isinstance(dialogo, list): continue
                
                yield f"   - Sintetizando estudo {estudo_idx+1}..."
                
                # Caminho único para este estudo nesta execução
                caminho_estudo = os.path.join(DATA_DIR, "audios", f"temp_{run_id}_estudo{estudo_idx+1}.mp3")
                
                try:
                    if not elevenlabs_client: raise ValueError("Sem chave ElevenLabs")

                    dialogo_audio = normalizar_dialogo_para_audio(dialogo)
                    total_chars_elevenlabs += sum(len(fala.get("text", "")) for fala in dialogo_audio)

                    usou_dialogue_v3 = False
                    if ELEVEN_AUDIO_DIALOGUE_ENABLED:
                        try:
                            caminho_gerado, _ = gerar_dialogo_com_eleven_v3(
                                dialogo_audio,
                                caminho_estudo,
                                seed=(
                                    ELEVEN_AUDIO_DIALOGUE_SEED + estudo_idx
                                    if ELEVEN_AUDIO_DIALOGUE_SEED is not None
                                    else None
                                ),
                            )
                            audio_paths.append(caminho_gerado)
                            usou_dialogue_v3 = True

                            if opcoes.get('firebase'):
                                try:
                                    from firebase_service import upload_file
                                    dest_blob = f"audios_raw/{hoje}/estudo{estudo_idx+1}_dialogue_v3.mp3"
                                    upload_file(caminho_gerado, dest_blob)
                                except Exception as e_upload:
                                    print(f"⚠️ Erro upload audio v3: {e_upload}")
                        except Exception as e_v3:
                            print(f"⚠️ Eleven v3 falhou no estudo {estudo_idx+1}, usando fallback: {e_v3}")

                    if usou_dialogue_v3:
                        continue

                    estudo_audios = []
                    for fala_idx, fala in enumerate(dialogo_audio):
                        # Normaliza formato do JSON (pode vir {"HOST": "texto"} ou {"speaker": "HOST", "text": "texto"})
                        speaker_clean = "HOST"
                        text = ""
                        
                        if 'text' in fala:
                            # Formato estruturado
                            text = fala.get('text', '')
                            speaker_clean = fala.get('speaker', 'HOST').strip().upper()
                        else:
                            # Formato chave-valor
                            if 'HOST' in fala:
                                speaker_clean = 'HOST'
                                text = fala['HOST']
                            elif 'COHOST' in fala:
                                speaker_clean = 'COHOST'
                                text = fala['COHOST']
                            
                        if not text: 
                            print(f"⚠️ Fala vazia ou formato desconhecido no índice {fala_idx}: {fala}")
                            continue
                        
                        # DEBUG: Ver quem está falando
                        # print(f"   [FALA {fala_idx+1}] Speaker: '{speaker_clean}'")
                        
                        # DEBUG: Ver quem está falando
                        print(f"   [FALA {fala_idx+1}] Speaker JSON: '{fala.get('speaker')}' -> Clean: '{speaker_clean}'")
                        
                        voice_id = ELEVEN_VOICE_ID_HOST if speaker_clean == 'HOST' else ELEVEN_VOICE_ID_COHOST
                        
                        # Configurações ajustadas para fala mais rápida/dinâmica
                        audio_generator = elevenlabs_client.text_to_speech.convert(
                            voice_id=voice_id,
                            text=text,
                            model_id=ELEVEN_AUDIO_FALLBACK_MODEL,
                            output_format="mp3_44100_128",
                            voice_settings=VoiceSettings(
                                stability=0.50,       # "More Emotion" setting
                                similarity_boost=0.75, 
                                style=0.20,           # "With Style" setting approved by user
                                use_speaker_boost=True
                            )
                        )
                        
                        segmento_path = os.path.join(DATA_DIR, "audios", f"temp_{run_id}_e{estudo_idx}_f{fala_idx}.mp3")
                        # Assuming 'save' is a function that writes the audio generator to a file
                        # If 'save' is not defined, this will cause a NameError.
                        # For now, I'll assume it's defined elsewhere or replace with the original file writing logic.
                        # Given the instruction, I'll use 'save' as provided.
                        # If 'save' is not available, the original code's way of writing to file is:
                        with open(segmento_path, "wb") as f:
                            for chunk in audio_generator: f.write(chunk)
                        estudo_audios.append(segmento_path)
                        
                        # Upload individual para Firebase Storage (se habilitado)
                        if opcoes.get('firebase'):
                            try:
                                from firebase_service import upload_file
                                dest_blob = f"audios_raw/{hoje}/temp_estudo{estudo_idx+1}_fala{fala_idx+1}.mp3"
                                upload_file(segmento_path, dest_blob)
                            except Exception as e_upload:
                                print(f"⚠️ Erro upload audio temp: {e_upload}")
                    
                    if estudo_audios:
                        estudo_combinado = AudioSegment.empty()
                        for idx, audio_path in enumerate(estudo_audios):
                            estudo_combinado += AudioSegment.from_file(audio_path, format="mp3")
                            if idx < len(estudo_audios) - 1: estudo_combinado += AudioSegment.silent(duration=300)
                        
                        # Transição
                        # Assuming POOL_TRANSICTIONS is defined elsewhere
                        # transicao = random.choice(POOL_TRANSICTIONS)
                        # ... (logica transicao simplificada ou removida se nao tiver audio pronto)
                        
                        estudo_combinado.export(caminho_estudo, format="mp3")
                        audio_paths.append(caminho_estudo)
                        
                        # Mantendo arquivos temporários para permitir edição/remixagem posterior
                        # for temp_path in estudo_audios: os.remove(temp_path)
                except Exception as e:
                    print(f"Erro audio: {e}")

            # Intro e mixagem final
            if audio_paths:
                yield "   - Montando episódio final..."
                print(f"📂 Arquivos para mixagem ({len(audio_paths)}):")
                for p in audio_paths:
                    sz = os.path.getsize(p) if os.path.exists(p) else 0
                    print(f"   -> {p} ({sz/1024:.1f} KB)")
                    
                try:
                    # Carrega intro musical fixa (intro_guto.mp3 ou similar) se existir, senão silêncio
                    if os.path.exists(INTRO_PATH):
                        intro_audio = AudioSegment.from_file(INTRO_PATH, format="mp3")
                    else:
                        intro_audio = AudioSegment.silent(duration=1000)

                    episodio = intro_audio + AudioSegment.silent(duration=1000)
                    
                    for caminho in audio_paths:
                        try:
                            segmento = AudioSegment.from_file(caminho, format="mp3")
                            episodio += segmento + AudioSegment.silent(duration=2500)
                        except Exception as e_seg:
                            print(f"⚠️ Erro ao adicionar segmento {caminho}: {e_seg}")
                    
                    episodio.export(episodio_path, format="mp3")
                    print(f"🎧 Episódio salvo: {episodio_path}")
                    yield f"💰 Consumo ElevenLabs: {total_chars_elevenlabs} chars."
                    
                    if opcoes.get('firebase'):
                        # Upload do episódio final
                        from firebase_service import upload_file, update_podcast_feed
                        duracao = len(episodio) / 1000.0
                        tamanho = os.path.getsize(episodio_path)
                        
                        audio_url = upload_file(episodio_path, f"episodios/{os.path.basename(episodio_path)}")
                        if audio_url:
                            # Update Feed
                            rss_url = update_podcast_feed(
                                audio_url, 
                                f"RevaCast Weekly - {hoje}", 
                                f"Episódio gerado automaticamente em {hoje}. Destaques da semana.",
                                datetime.now(pytz.utc),
                                duracao,
                                tamanho
                            )
                            yield {"audio_url": audio_url, "rss_url": rss_url}
                            
                except Exception as e:
                    print(f"Erro montagem final: {e}")
                    raise e
        else:
            yield "⚠️ Sem roteiro novo para gerar áudio."
    else:
        yield "⏭️ Pulando geração de Áudio."

    # ------------------------------------------------------------------
    # 4) MAILCHIMP
    # ------------------------------------------------------------------
    campaign_id = None
    mailchimp_status = "skipped"
    mailchimp_error = None
    mailchimp_schedule_time = None

    if opcoes.get('mailchimp'):
        yield "📧 4/5: Mailchimp..."
        if os.path.exists(boletim_path):
            # ... (Lógica Mailchimp mantida) ...
            try:
                assunto = "Boletim Científico Semanal | RevaCast"
                with open(boletim_path, "r", encoding="utf-8") as f: conteudo_boletim = f.read()
                # O conteúdo já está em HTML formatado, não precisamos de replacements destrutivos
                html_final_completo = TEMPLATE_HTML_BASE.format(conteudo_aqui=conteudo_boletim)

                campaign = mc.campaigns.create({
                    "type": "regular", "recipients": {"list_id": MC_LIST_ID},
                    "settings": {"subject_line": assunto, "title": f"Boletim {hoje}", "from_name": MC_FROM_NAME, "reply_to": MC_REPLY_TO}
                })
                campaign_id = campaign["id"]
                mc.campaigns.set_content(campaign_id, {"html": html_final_completo})
                
                brasilia_tz = pytz.timezone("America/Sao_Paulo")
                agora_brasilia = datetime.now(brasilia_tz)
                dias_ate_sabado = (5 - agora_brasilia.weekday() + 7) % 7
                if dias_ate_sabado == 0: dias_ate_sabado = 7
                data_envio = (agora_brasilia + timedelta(days=dias_ate_sabado)).replace(hour=7, minute=30, second=0, microsecond=0)
                
                mc.campaigns.schedule(campaign_id, {"schedule_time": data_envio.isoformat()})
                mailchimp_status = "scheduled"
                mailchimp_schedule_time = data_envio.isoformat()
                yield f"✅ Campanha agendada para {data_envio}"
            except Exception as e:
                mailchimp_status = "error"
                mailchimp_error = str(e)
                yield f"❌ Erro Mailchimp: {e}"
        else:
            yield "⚠️ Arquivo do boletim não encontrado para envio."
    else:
        yield "⏭️ Pulando Mailchimp."

    # ------------------------------------------------------------------
    # 5) UPLOAD E RSS
    # ------------------------------------------------------------------
    rss_url = None
    audio_url = None

    if opcoes.get('firebase'):
        yield "☁️ 5/5: Firebase Upload..."
        
        # Tenta subir o episódio completo
        if episodio_path and os.path.exists(episodio_path):
            try:
                from firebase_service import upload_file, update_podcast_feed
                filename = os.path.basename(episodio_path)
                audio_url = upload_file(episodio_path, f"episodios/{filename}")
                if audio_url:
                    audio_segment = AudioSegment.from_file(episodio_path)
                    duracao_seg = len(audio_segment) / 1000
                    tamanho_bytes = os.path.getsize(episodio_path)
                    
                    rss_url = update_podcast_feed(
                        episodio_audio_url=audio_url,
                        episodio_titulo=f"Boletim {hoje}",
                        episodio_descricao="Resumo semanal das evidências científicas.",
                        data_pub=datetime.now(pytz.timezone("America/Sao_Paulo")),
                        duracao_segundos=duracao_seg,
                        tamanho_bytes=tamanho_bytes
                    )
                    yield f"✅ Upload concluído! RSS atualizado: {rss_url}"
            except Exception as e:
                yield f"❌ Erro no upload: {e}"
        else:
            yield "⚠️ Arquivo final do episódio não encontrado (erro na mixagem?). Tentando resgatar partes..."
            # Lógica de Resgate: Procura por arquivos parciais dos estudos
            try:
                from firebase_service import upload_file
                import glob
                # Procura por estudo*_completo.mp3 na pasta de áudio
                partes = glob.glob(os.path.join(AUDIO_DIR, "estudo*_completo.mp3"))
                if partes:
                    yield f"🚑 Resgatando {len(partes)} arquivos parciais..."
                    for parte in partes:
                        fname = os.path.basename(parte)
                        url_parte = upload_file(parte, f"resgate/{hoje}/{fname}")
                        yield f"   - Salvo: {fname}"
                else:
                    yield "❌ Nenhum arquivo parcial encontrado para resgate."
            except Exception as e_rescue:
                yield f"❌ Erro no resgate: {e_rescue}"




    # ------------------------------------------------------------------
    # 5.5) INTEGRAÇÃO WHATSAPP (NOVO)
    # ------------------------------------------------------------------
    if opcoes.get('firebase'):
         try:
            from whatsapp_service import create_draft
            yield "📱 Gerando rascunho WhatsApp (Weekly)..."
            
            # Link prioritário: RSS > Audio > Site
            link_final = rss_url if rss_url else (audio_url if audio_url else "https://www.revalidatie.com.br")
            
            wa_content = {
                "title": f"RevaCast Weekly - {hoje}",
                "summary": f"Novo episódio no ar! Confira os destaques da semana em {link_final}",
                "link": link_final
            }
            
            draft = create_draft("revacast_weekly", wa_content)
            if draft:
                 yield f"✅ Rascunho WhatsApp criado com sucesso!"
         except Exception as e_wa:
             yield f"⚠️ Erro ao gerar rascunho WhatsApp: {e_wa}"

    # ------------------------------------------------------------------
    # 6) RETORNO
    # ------------------------------------------------------------------
    # ======================================================================
    # RELATÓRIO DE CUSTOS ESTIMADO
    # ======================================================================
    # Preços aproximados (Dez 2024)
    # Gemini 1.5 Flash: ~$0.075 / 1M tokens (input), $0.30 / 1M tokens (output) -> Quase zero para esse uso
    # OpenAI GPT-4o: $2.50 / 1M (input), $10.00 / 1M (output)
    # ElevenLabs: ~$0.30 / 1000 chars (Creator tier aprox) ou variavel
    # Vamos usar valores conservadores para estimativa
    
    # Estimativa de tokens (muito aproximada, pois não pegamos o usage exato de cada call no código atual sem refatorar tudo)
    # Assumindo média de 3000 tokens input / 1000 tokens output para todo o fluxo de texto (OpenAI + Gemini)
    custo_texto_estimado = 0.15 # $0.15 fixo como "teto" para texto (GPT-4o + Gemini)
    
    # ElevenLabs é o mais caro e mensurável por caracteres
    custo_audio_estimado = (total_chars_elevenlabs / 1000) * 0.30 # $0.30 por 1k chars (exemplo, ajuste conforme seu plano)
    
    custo_total = custo_texto_estimado + custo_audio_estimado
    
    moeda_cambio = 6.0 # 1 USD = 6 BRL (margem de segurança)
    custo_brl = custo_total * moeda_cambio

    resultado = {
        "data_referencia": hoje,
        "boletim_path": boletim_path,
        "episodio_path": episodio_path,
        "brief_spotify_path": brief_spotify_path if os.path.exists(brief_spotify_path) else None,
        "brief_spotify_text": brief_spotify_text if brief_spotify_text else None,
        "brief_spotify_download_url": f"/baixar-brief/{hoje}" if os.path.exists(brief_spotify_path) else None,
        "audio_url": audio_url,
        "rss_url": rss_url,
        "mailchimp": {"status": mailchimp_status, "error": mailchimp_error},
        "custos": {
            "elevenlabs_chars": total_chars_elevenlabs,
            "estimativa_usd": round(custo_total, 2),
            "estimativa_brl": round(custo_brl, 2),
            "detalhe": f"Texto: ~$0.15 | Áudio: ~${custo_audio_estimado:.2f} ({total_chars_elevenlabs} chars)"
        }
    }

    print("✅ Pipeline finalizado.")
    print(f"💰 Custo Estimado: R$ {custo_brl:.2f}")
    yield resultado


# Permite testar localmente: python boletim_service.py
if __name__ == "__main__":
    for msg in rodar_boletim():
        print(msg)

# ------------------------------------------------------------------------------
# FERRAMENTA DE RECUPERAÇÃO / REMIXAGEM
# ------------------------------------------------------------------------------
from firebase_admin import storage

def remixar_audio_from_firebase():
    """Baixa segmentos do dia atual do Firebase e remonta o episódio."""
    print("🔄 Iniciando Remixagem via Firebase...")
    try:
        bucket = storage.bucket()
        hoje = datetime.now().strftime("%Y-%m-%d")
        
        # Procura na pasta do dia
        prefix = f"audios_raw/{hoje}/"
        blobs = list(bucket.list_blobs(prefix=prefix))
        
        valid_blobs = [b for b in blobs if b.name.endswith(".mp3")]
        
        if not valid_blobs:
             return {"success": False, "message": f"Nenhum áudio encontrado em {prefix}"}
             
        print(f"Baixando {len(valid_blobs)} segmentos de {prefix}...")
        
        local_dir = os.path.join(DATA_DIR, "remix_temp")
        os.makedirs(local_dir, exist_ok=True)
        
        downloaded_files = []
        for blob in valid_blobs:
            local_path = os.path.join(local_dir, os.path.basename(blob.name))
            blob.download_to_filename(local_path)
            downloaded_files.append(local_path)
            
        # Ordenação inteligente
        def sort_key(fname):
            name = os.path.basename(fname)
            try:
                import re
                # Busca 'estudoX' e 'falaY'
                # Exemplo: temp_UUID_e1_f3.mp3 ou temp_estudo1_fala3.mp3
                estudo_match = re.search(r'(?:estudo|e)(\d+)', name)
                fala_match = re.search(r'(?:fala|f)(\d+)', name)
                
                e_idx = int(estudo_match.group(1)) if estudo_match else 999
                f_idx = int(fala_match.group(1)) if fala_match else 0
                
                if "abertura" in name or "intro" in name: return (-1, 0)
                return (e_idx, f_idx)
            except:
                return (999, 999)

        downloaded_files.sort(key=sort_key)
        
        print(f"Montando episódio com {len(downloaded_files)} arquivos ordenados...")
        episodio = AudioSegment.silent(duration=500)
        
        # Opcional: Intro musical fixa se existir
        # path_intro = os.path.join(DATA_DIR, "intro.mp3")
        # if os.path.exists(path_intro): episodio += AudioSegment.from_file(path_intro)
        
        for fpath in downloaded_files:
            seg = AudioSegment.from_file(fpath)
            episodio += seg + AudioSegment.silent(duration=350)
            
        final_path = os.path.join(DATA_DIR, "audios", f"episodio_remixado_{hoje}.mp3")
        os.makedirs(os.path.dirname(final_path), exist_ok=True)
        episodio.export(final_path, format="mp3")
        
        # Upload
        from firebase_service import upload_file, update_podcast_feed
        duracao = len(episodio) / 1000.0
        tamanho = os.path.getsize(final_path)
        
        audio_url = upload_file(final_path, f"episodios/{os.path.basename(final_path)}")
        rss_url = ""
        if audio_url:
            rss_url = update_podcast_feed(
                audio_url, 
                f"RevaCast Weekly - {hoje} (Remix)", 
                f"Episódio recuperado/remixado em {hoje}.",
                datetime.now(pytz.utc),
                duracao,
                tamanho
            )
        
        # Cleanup
        try:
            for f in downloaded_files: os.remove(f)
            os.rmdir(local_dir)
        except: pass
        
        return {
            "success": True, 
            "audio_url": audio_url, 
            "rss_url": rss_url, 
            "message": f"Remixado {len(downloaded_files)} arquivos. Duração: {duracao:.1f}s"
        }

    except Exception as e:
        print(f"Erro remix: {e}")
        return {"success": False, "error": str(e)}
