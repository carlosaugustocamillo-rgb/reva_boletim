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
import textwrap
from datetime import datetime, timedelta

import pytz
from Bio import Entrez
from openai import OpenAI
from pydub import AudioSegment
from mailchimp_marketing import Client
from mailchimp_marketing.api_client import ApiClientError
from dotenv import load_dotenv

# Carrega variáveis de ambiente do arquivo .env
load_dotenv()


# ======================================================================
# CONFIGURAÇÕES GERAIS
# ======================================================================

# Diretório base: sempre uma pasta "data" dentro do próprio projeto
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.join(PROJECT_ROOT, "data")
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

# ... (imports)

# --- ElevenLabs ---
ELEVENLABS_API_KEY = os.environ.get("ELEVENLABS_API_KEY")
ELEVEN_VOICE_ID_HOST = os.environ.get("ELEVEN_VOICE_ID_HOST", "WWL28Z00upcD5SGFqY2n")
ELEVEN_VOICE_ID_COHOST = os.environ.get("ELEVEN_VOICE_ID_COHOST", "x3mAOLD9WzlmrFCwA1S3")
elevenlabs_client = None
if ELEVENLABS_API_KEY:
    elevenlabs_client = ElevenLabs(api_key=ELEVENLABS_API_KEY)
    print(f"✅ ElevenLabs configurado com vozes:")
    print(f"   HOST: {ELEVEN_VOICE_ID_HOST}")
    print(f"   COHOST: {ELEVEN_VOICE_ID_COHOST}")


# ======================================================================
# FUNÇÕES AUXILIARES (copiadas/adaptadas do seu Colab)
# ======================================================================

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
        # AbstractText pode vir como lista de strings ou blocos – aqui pego tudo concatenado
        abstract = art_data.get('Abstract', {}).get('AbstractText', [])
        if isinstance(abstract, list):
            abstract = " ".join(str(p) for p in abstract)
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

    resposta = client.chat.completions.create(
        model="gpt-5.1",
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
    
    # Define se é o primeiro estudo ou não para ajustar a transição
    if idx == 0:
        contexto_inicial = "Este é o PRIMEIRO estudo do episódio. Comece direto apresentando o estudo, sem saudações."
    else:
        contexto_inicial = f"Este é o estudo número {idx + 1}. NÃO faça saudações. Comece DIRETO com uma transição natural tipo 'Agora vamos falar de outro estudo...' ou 'O próximo artigo trata de...'."
    
    # Define se deve ter despedida no final
    if is_last:
        contexto_final = "Este é o ÚLTIMO estudo do episódio. Após discutir o estudo, FINALIZE o podcast com uma despedida calorosa. Agradeça os ouvintes e diga 'até a próxima!'"
    else:
        contexto_final = "NÃO finalize o podcast. Deixe a conversa aberta para o próximo estudo."
    
    prompt = f"""
Você é um roteirista sênior do RevaCast Weekly. Seu objetivo é transformar um resumo científico em uma CONVERSA EXTREMAMENTE NATURAL e ENGAJANTE, no estilo "NotebookLM Audio Overview".

O tom deve ser de dois colegas/amigos apaixonados por ciência conversando no corredor ou num café. Nada de "palestra" ou "leitura de texto".

PERSONAGENS:
- HOST: Especialista, entusiasta, traz a novidade com energia ("Cara, olha que incrível isso aqui").
- COHOST: Curioso, inteligente, reage com emoção ("Sério?", "Não acredito!", "Uau"), faz perguntas pertinentes e tenta conectar os pontos.

REGRAS DE ESTILO (CRUCIAIS):
1. USE LINGUAGEM FALADA REAL: Use marcadores de discurso ("Então", "Sabe...", "Olha só", "Pois é").
2. REAÇÕES GENUÍNAS: Se o resultado for bom, o COHOST deve ficar impressionado. Se for polêmico, deve ficar surpreso.
3. FLUIDEZ: Uma fala deve "enganchar" na outra. Evite perguntas e respostas robóticas (tipo ping-pong).
4. SEM FORMALIDADES: Evite "O estudo concluiu que". Prefira "O que eles descobriram foi...", "O mais louco é que...".
5. EXPLICAR O "PORQUÊ": Não jogue apenas dados. Explique o impacto prático disso.

REGRAS DE SEGURANÇA (CRÍTICO - RISCO DE ALUCINAÇÃO):
- CRIATIVIDADE APENAS NO TOM E NA CONVERSA.
- RIGOR ABSOLUTO NOS DADOS: Números, p-values, tamanhos de amostra e resultados devem ser COPIADOS EXATAMENTE do resumo.
- PROIBIDO INVENTAR: Se o resumo não diz "melhorou 20%", NÃO invente "melhorou muito" ou "20%". Diga apenas "houve melhora significativa".
- Se não souber um dado, o COHOST deve perguntar e o HOST deve dizer "O resumo não especifica esse detalhe".

REGRAS TÉCNICAS:
- {contexto_inicial}
- {contexto_final}
- SEM SAUDAÇÕES ("Olá", "Bem-vindos").
- SEM APRESENTAÇÕES ("Eu sou o Guto").
- Duração: 6 a 10 trocas de fala (para dar tempo de aprofundar um pouco).
- Use as abreviações formatadas: D.P.O.C., V.E.F.1, I.M.C.

Contexto do estudo:
Título: {titulo}
Primeiro autor: {primeiro_autor}

Resumo traduzido:
{resumo_pt}

FORMATO DE RETORNO (JSON array):
[
  {{"speaker": "HOST", "text": "Cara, você não vai acreditar nesse estudo sobre DPOC que saiu..."}},
  {{"speaker": "COHOST", "text": "Sério? O que tem de novo? É sobre reabilitação?"}},
  {{"speaker": "HOST", "text": "Exatamente! Eles pegaram um grupo enorme e..."}}
]
"""

    resposta = client.chat.completions.create(
        model="gpt-5.1",
        messages=[
            {
                "role": "system",
                "content": "Você é um roteirista de podcast. Retorne SEMPRE e SOMENTE um JSON array válido com o diálogo. NUNCA use saudações ou apresentações."
            },
            {"role": "user", "content": prompt}
        ],
        temperature=0.85,
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


def dividir_texto(texto, limite=4096):
    # simples quebra por tamanho aproximado de caracteres
    return textwrap.wrap(texto, width=limite, break_long_words=False, break_on_hyphens=False)


def artigo_tem_exercicio_no_resumo(artigo):
    resumo = artigo.get('resumo_original', '').lower()
    termos = ['exercise', 'exercício', 'exercícios', 'exercising', 'training', 'atividade física']
    return any(term in resumo for term in termos)


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
    'Câncer: ("exercise"[tiab]) AND ("Neoplasms"[Mesh] OR cancer[tiab] OR cancers[tiab])'
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

def rodar_boletim():
    """
    Executa TODO o pipeline:
    - Boletim principal (HTML para e-mail)
    - Boletim de revisão
    - Boletim detalhado
    - Roteiro + áudio do podcast
    - Cria e agenda campanha no Mailchimp para sábado 7:30 (Brasília)

    Retorna um gerador que emite mensagens de status e, por fim, o dicionário de resultado.
    """
    yield "🚀 Iniciando pipeline do boletim científico..."

    hoje = datetime.today().strftime('%Y-%m-%d')

    # Paths principais
    boletim_path = os.path.join(BASE_DIR, f"boletim_pubmed_{hoje}.txt")
    revisao_path = os.path.join(BASE_DIR, f"boletim_para_revisao_{hoje}.txt")
    boletim_detalhado_path = os.path.join(BASE_DIR, f"boletim_detalhado_{hoje}.txt")
    episodio_path = os.path.join(BASE_DIR, f"episodio_boletim_{hoje}.mp3")

    # ------------------------------------------------------------------
    # 1) BOLETIM PRINCIPAL
    # ------------------------------------------------------------------
    yield "🔎 1/5: Buscando artigos no PubMed e gerando Boletim Principal..."
    
    boletim_final = (
        "Os dados a seguir mostram os estudos publicados no PubMed na última semana "
        "(sábado a sexta), com os resumos traduzidos literalmente do inglês para o português.\n\n"
    )
    boletim_revisao = (
        "Os artigos a seguir foram encontrados no PubMed, mas não continham resumo disponível "
        "ou houve falha técnica na tradução. Eles requerem revisão manual.\n\n"
    )
    artigos_vistos = set()

    for consulta in CONSULTAS_PRINCIPAIS:
        titulo_secao, query = consulta.split(": ", 1)
        yield f"   - Processando seção: {titulo_secao}..."
        ids = buscar_ids(query)
        artigos = buscar_info_estruturada(ids)
        artigos_unicos = [a for a in artigos if a['pmid'] not in artigos_vistos]

        if not artigos_unicos:
            continue

        cabecalho_secao_html = f"""
<h2 style="font-family: Helvetica, Arial, sans-serif; color: #205776; border-bottom: 2px solid #205776; padding-bottom: 5px; margin-top: 30px;">
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

    boletim_final += "\nEspero que estas traduções sejam úteis. Siga nosso podcast para mais!\nAbraços,\nGuto"

    # Salvar boletins principal e revisão
    with open(boletim_path, "w", encoding="utf-8") as f:
        f.write(boletim_final)
    print(f"✅ Boletim principal salvo como: {boletim_path}")

    with open(revisao_path, "w", encoding="utf-8") as f:
        f.write(boletim_revisao)
    print(f"✅ Boletim para revisão salvo como: {revisao_path}")

    # ------------------------------------------------------------------
    # 2) BOLETIM DETALHADO + ROTEIROS
    # ------------------------------------------------------------------
    yield "📝 2/5: Gerando Boletim Detalhado e Roteiros de Podcast..."
    
    boletim_detalhado = (
        "Boletim detalhado com os estudos mais relevantes sobre exercício em diferentes "
        "condições crônicas, publicados na última semana (sábado a sexta). "
        "Os resumos abaixo são traduções literais do PubMed.\n\n"
    )
    roteiros_audio = []
    
    # Primeiro, coleta todos os artigos relevantes
    todos_artigos_relevantes = []

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
                boletim_detalhado += f"""* {art['titulo']}
* {autores_str}
* {art['journal']}, {art['ano']}; {art['volume']}({art['issue']}): {art['paginas']}
* DOI: {link_doi}
* PubMed: {link_pubmed}

Resumo ausente no PubMed. Recomenda-se leitura direta do artigo.

---

"""
                continue

            try:
                resumo_traduzido = traduzir_resumo(resumo_original)
            except Exception as e:
                boletim_detalhado += f"""* {art['titulo']}
* {autores_str}
* {art['journal']}, {art['ano']}; {art['volume']}({art['issue']}): {art['paginas']}
* DOI: {link_doi}
* PubMed: {link_pubmed}

Falha na tradução automática do resumo ({e}). Recomenda-se revisão manual.

---

"""
                continue

            boletim_detalhado += f"""* {art['titulo']}
* {autores_str}
* {art['journal']}, {art['ano']}; {art['volume']}({art['issue']}): {art['paginas']}
* DOI: {link_doi}
* PubMed: {link_pubmed}

{resumo_traduzido}

---

"""
            primeiro_autor = art['autores'][0] if art['autores'] else "Autor não identificado"
            todos_artigos_relevantes.append({
                'titulo': art['titulo'],
                'resumo_traduzido': resumo_traduzido,
                'primeiro_autor': primeiro_autor
            })
    
    # Agora gera os roteiros, sabendo qual é o último
    total_artigos = len(todos_artigos_relevantes)
    for idx, artigo_info in enumerate(todos_artigos_relevantes):
        is_last = (idx == total_artigos - 1)
        yield f"   - Gerando roteiro para estudo {idx+1}/{total_artigos}..."
        roteiro = resumo_para_podcast(
            artigo_info['titulo'], 
            artigo_info['resumo_traduzido'], 
            artigo_info['primeiro_autor'], 
            idx=idx,
            is_last=is_last
        )
        roteiros_audio.append(roteiro)
    
    # Salva o roteiro completo como arquivo de texto
    roteiro_path = os.path.join(BASE_DIR, f"roteiro_podcast_{hoje}.txt")
    with open(roteiro_path, "w", encoding="utf-8") as f:
        f.write("ROTEIRO COMPLETO DO PODCAST - RevaCast Weekly\n")
        f.write("=" * 60 + "\n\n")
        
        for estudo_idx, dialogo in enumerate(roteiros_audio):
            f.write(f"\n{'='*60}\n")
            f.write(f"ESTUDO {estudo_idx + 1}\n")
            f.write(f"{'='*60}\n\n")
            
            if isinstance(dialogo, list):
                for fala_idx, fala in enumerate(dialogo):
                    speaker = fala.get('speaker', 'HOST')
                    text = fala.get('text', '')
                    f.write(f"{speaker}: {text}\n\n")
            else:
                f.write(f"HOST: {dialogo}\n\n")
    
    print(f"✅ Roteiro completo salvo como: {roteiro_path}")

    boletim_detalhado += "\nCompartilhe com colegas. RevaCast Pesquisa Detalhada!"

    with open(boletim_detalhado_path, "w", encoding="utf-8") as f:
        f.write(boletim_detalhado)
    print(f"✅ Boletim detalhado salvo como: {boletim_detalhado_path}")

    # ------------------------------------------------------------------
    # 3) GERAÇÃO DO ÁUDIO EM FORMATO DE CONVERSA (se houver roteiros)
    # ------------------------------------------------------------------
    yield "🎙️ 3/5: Verificando áudio..."
    
    total_chars_elevenlabs = 0
    
    if os.path.exists(episodio_path):
        yield f"⚠️ Áudio já encontrado para hoje! Pulando geração do ElevenLabs para economizar seus créditos."
        print(f"Áudio existente mantido: {episodio_path}")
    elif roteiros_audio:
        yield "🎙️ Gerando Áudio (ElevenLabs) - Isso pode demorar..."
        audio_paths = []
        
        # Para cada estudo, gerar o diálogo completo como um único áudio conversacional
        for estudo_idx, dialogo in enumerate(roteiros_audio):
            if not isinstance(dialogo, list):
                print(f"Aviso: roteiro do estudo {estudo_idx+1} não é uma lista. Pulando.")
                continue
            
            try:
                if not elevenlabs_client:
                    raise ValueError("A chave da API ELEVENLABS_API_KEY não foi configurada.")
                
                # Preparar o roteiro conversacional completo
                roteiro_texto = ""
                for fala in dialogo:
                    speaker = fala.get('speaker', 'HOST')
                    text = fala.get('text', '')
                    if text:
                        # Marca quem está falando usando tags especiais
                        speaker_name = "apresentador 1" if speaker == 'HOST' else "apresentador 2"
                        roteiro_texto += f"{speaker_name}: {text}\n\n"
                
                if not roteiro_texto:
                    continue
                
                # Gera um único áudio conversacional para todo o estudo
                # Alternando automaticamente entre as vozes
                yield f"   - Sintetizando vozes para estudo {estudo_idx+1}..."
                print(f"🎙️ Gerando áudio conversacional para estudo {estudo_idx+1}...")
                
                # Como a API do ElevenLabs não tem suporte nativo para conversação com múltiplas vozes
                # em uma única chamada, vamos gerar com pausas menores entre as falas
                estudo_audios = []
                for fala_idx, fala in enumerate(dialogo):
                    speaker = fala.get('speaker', 'HOST')
                    text = fala.get('text', '')
                    
                    if not text:
                        continue
                    
                    # Contabiliza caracteres para estimativa de custo
                    total_chars_elevenlabs += len(text)
                    
                    # Escolhe a voz baseado no speaker
                    # Garante que espaços extras não atrapalhem a verificação
                    speaker_clean = speaker.strip().upper()
                    voice_id = ELEVEN_VOICE_ID_HOST if speaker_clean == 'HOST' else ELEVEN_VOICE_ID_COHOST
                    
                    print(f"   → {speaker} ({len(text)} chars): usando voice_id = {voice_id}")
                    
                    from elevenlabs import VoiceSettings
                    
                    # Trocando de volta para o multilingual_v2 para garantir a fidelidade da voz (HOST correto)
                    # Ajustando configurações para mais entonação e velocidade natural
                    audio_generator = elevenlabs_client.text_to_speech.convert(
                        voice_id=voice_id,
                        text=text,
                        model_id="eleven_multilingual_v2",
                        voice_settings=VoiceSettings(
                            stability=0.50,  # Reduzido: permite mais variação de tom e velocidade (menos robótico)
                            similarity_boost=0.80,  # Mantém a identidade da voz
                            style=0.55,  # Aumentado: adiciona a "entonação" e expressividade pedida
                            use_speaker_boost=True
                        )
                    )
                    
                    caminho_temp = os.path.join(AUDIO_DIR, f"temp_estudo{estudo_idx+1}_fala{fala_idx+1}.mp3")
                    
                    with open(caminho_temp, "wb") as f:
                        for chunk in audio_generator:
                            f.write(chunk)
                    
                    estudo_audios.append(caminho_temp)
                
                # Combina todos os áudios do estudo com pausas curtas (300ms)
                if estudo_audios:
                    estudo_combinado = AudioSegment.empty()
                    for idx, audio_path in enumerate(estudo_audios):
                        estudo_combinado += AudioSegment.from_file(audio_path, format="mp3")
                        # Pausa curta entre falas (exceto na última)
                        if idx < len(estudo_audios) - 1:
                            estudo_combinado += AudioSegment.silent(duration=300)
                    
                    caminho_estudo = os.path.join(AUDIO_DIR, f"estudo{estudo_idx+1}_completo.mp3")
                    estudo_combinado.export(caminho_estudo, format="mp3")
                    audio_paths.append(caminho_estudo)
                    
                    # Remove arquivos temporários
                    for temp_path in estudo_audios:
                        try:
                            os.remove(temp_path)
                        except:
                            pass
                    
                    print(f"✅ Áudio conversacional do estudo {estudo_idx+1} gerado: {caminho_estudo}")
                    
            except Exception as e:
                print(f"Erro ao gerar áudio do estudo {estudo_idx+1}: {e}")

        # Carregar intro
        yield "   - Montando episódio final com intro..."
        try:
            intro_audio = AudioSegment.from_file(INTRO_PATH, format="mp3")
        except Exception as e:
            raise FileNotFoundError(
                f"Não foi possível encontrar sua introdução em {INTRO_PATH}: {e}"
            )

        episodio = intro_audio + AudioSegment.silent(duration=1000)
        for caminho in audio_paths:
            try:
                episodio += AudioSegment.from_file(caminho, format="mp3") + AudioSegment.silent(duration=2500)
            except Exception as e:
                print(f"Erro ao juntar {caminho}: {e}")

        episodio.export(episodio_path, format="mp3")
        print(f"\n🎧 Episódio final salvo em: {episodio_path}")
        yield f"💰 Consumo ElevenLabs: {total_chars_elevenlabs} caracteres usados neste episódio."
    else:
        print("Nenhum roteiro para gerar áudio.")
        episodio_path = None

    # ------------------------------------------------------------------
    # 4) MAILCHIMP – CRIAR E AGENDAR CAMPANHA
    # ------------------------------------------------------------------
    yield "📧 4/5: Criando e agendando campanha no Mailchimp..."
    
    assunto = "Boletim Científico Semanal | RevaCast"

    with open(boletim_path, "r", encoding="utf-8") as f:
        conteudo_boletim = f.read()

    conteudo_boletim = conteudo_boletim.replace('"', '&quot;')
    conteudo_boletim = conteudo_boletim.replace("'", '&apos;')
    conteudo_boletim = conteudo_boletim.replace("\n\n", "</p><p>").replace("\n", "<br>")
    conteudo_boletim = f"<p>{conteudo_boletim}</p>"

    html_final_completo = TEMPLATE_HTML_BASE.format(conteudo_aqui=conteudo_boletim)

    campaign_id = None
    mailchimp_status = "not_started"
    mailchimp_error = None
    mailchimp_schedule_time = None

    try:
        print("MAILCHIMP PASSO 1: Criando a 'casca' da campanha...")
        campaign = mc.campaigns.create({
            "type": "regular",
            "recipients": {"list_id": MC_LIST_ID},
            "settings": {
                "subject_line": assunto,
                "title": f"Boletim {hoje} (Agendado)",
                "from_name": MC_FROM_NAME,
                "reply_to": MC_REPLY_TO
            }
        })
        campaign_id = campaign["id"]
        print(f"✅ 'Casca' da campanha criada com sucesso: {campaign_id}")

        print("\nMAILCHIMP PASSO 2: Inserindo o conteúdo HTML na campanha...")
        mc.campaigns.set_content(campaign_id, {"html": html_final_completo})
        print("✅ Conteúdo HTML enviado para a campanha.")

        print("\nMAILCHIMP PASSO 3: Calculando data e agendando o envio para 7:30 AM (Horário de Brasília)...")

        brasilia_tz = pytz.timezone("America/Sao_Paulo")
        agora_brasilia = datetime.now(brasilia_tz)

        # próximo sábado
        dias_ate_sabado = (5 - agora_brasilia.weekday() + 7) % 7
        if dias_ate_sabado == 0:
            dias_ate_sabado = 7

        data_envio_brasilia = (agora_brasilia + timedelta(days=dias_ate_sabado)).replace(
            hour=7, minute=30, second=0, microsecond=0
        )

        # Mailchimp aceita ISO com timezone; o tz da datetime já é America/Sao_Paulo
        data_envio_iso = data_envio_brasilia.isoformat()

        print(f"Horário calculado em Brasília: {data_envio_brasilia.strftime('%Y-%m-%d %H:%M:%S %Z%z')}")
        print(f"🚀 Tentando agendar a campanha para: {data_envio_iso}")

        mc.campaigns.schedule(campaign_id, {"schedule_time": data_envio_iso})

        print("\n🏆🎉 SUCESSO! Campanha agendada e pronta para o envio às 7:30 (Horário de Brasília)!")
        mailchimp_status = "scheduled"
        mailchimp_schedule_time = data_envio_iso

    except ApiClientError as error:
        error_detail = error.text
        if campaign_id:
            try:
                status_check = mc.campaigns.get(campaign_id)
                error_detail += f" | Status atual da campanha: {status_check.get('status')}"
            except ApiClientError:
                error_detail += " | Não foi possível obter o status da campanha."

        print(f"\n❌ ERRO MAILCHIMP: {error_detail}")
        mailchimp_status = "error"
        mailchimp_error = error_detail

    # ------------------------------------------------------------------
    # 5) UPLOAD E RSS (FIREBASE)
    # ------------------------------------------------------------------
    yield "☁️ 5/5: Fazendo upload para Firebase e atualizando RSS..."
    
    rss_url = None
    audio_url = None
    
    if episodio_path and os.path.exists(episodio_path):
        try:
            from firebase_service import upload_file, update_podcast_feed
            
            # 1. Upload do Áudio
            filename = os.path.basename(episodio_path)
            audio_url = upload_file(episodio_path, f"episodios/{filename}")
            
            if audio_url:
                # 2. Metadados para o RSS
                audio_segment = AudioSegment.from_file(episodio_path)
                duration_sec = len(audio_segment) / 1000.0
                file_size = os.path.getsize(episodio_path)
                
                tz = pytz.timezone("America/Sao_Paulo")
                pub_date = datetime.now(tz)
                
                # Título e Descrição do episódio
                titulo_ep = f"RevaCast Weekly - {hoje}"
                descricao_ep = f"Resumo semanal dos artigos científicos. Confira o boletim completo em nosso site."
                
                # 3. Atualiza RSS
                rss_url = update_podcast_feed(
                    audio_url, 
                    titulo_ep, 
                    descricao_ep, 
                    pub_date, 
                    duration_sec, 
                    file_size
                )
                print(f"📡 Feed RSS atualizado: {rss_url}")
        except Exception as e:
            print(f"❌ Erro na etapa de Upload/RSS: {e}")

    # ------------------------------------------------------------------
    # 6) RETORNO PARA O AGENTE / API
    # ------------------------------------------------------------------
    resultado = {
        "data_referencia": hoje,
        "boletim_path": boletim_path,
        "revisao_path": revisao_path,
        "boletim_detalhado_path": boletim_detalhado_path,
        "episodio_path": episodio_path,
        "audio_url": audio_url,
        "rss_url": rss_url,
        "mailchimp": {
            "campaign_id": campaign_id,
            "status": mailchimp_status,
            "schedule_time": mailchimp_schedule_time,
            "error": mailchimp_error,
        },
    }

    print("✅ Pipeline do boletim finalizado.")
    yield resultado


# Permite testar localmente: python boletim_service.py
if __name__ == "__main__":
    res = rodar_boletim()
    print("\nRESULTADO FINAL:")
    print(res)
