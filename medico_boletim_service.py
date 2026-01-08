"""
medico_boletim_service.py

Service autônomo para gerar o Boletim Médico Semanal (foco em tratamento/intervenção, excluindo exercício).
Funcionalidades:
- Busca no PubMed com filtros específicos (Asma, DPOC, etc - Foco Farmacológico/Cirúrgico).
- Filtra estritamente por Journals Q1 (baseado em CSV local).
- Traduz resumos usando LLM.
- Gera HTML e cria campanha no Mailchimp.

NÃO MODIFIQUE ARQUIVOS EXISTENTES DO PROJETO.
"""

import os
import csv
import time
import textwrap
from datetime import datetime, timedelta
import pytz

try:
    from Bio import Entrez
except ImportError:
    print("⚠️ BioPython não instalado. Execute: pip install biopython")

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

import google.generativeai as genai
from openai import OpenAI

try:
    from mailchimp_marketing import Client
except ImportError:
    pass

# ======================================================================
# CONFIGURAÇÕES -- Carrega do .env
# ======================================================================

OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL", "email@example.com")
MC_API_KEY = os.environ.get("MC_API_KEY")
MC_SERVER = os.environ.get("MC_SERVER")
# Tenta pegar um ID específico para o Médico, senão usa o padrão
MC_LIST_ID = os.environ.get("MC_MEDICO_LIST_ID") or os.environ.get("MC_LIST_ID")
MC_FROM_NAME = os.environ.get("MC_FROM_NAME", "Revalidatie Médico")
MC_REPLY_TO = os.environ.get("MC_REPLY_TO", "contato@revalidatie.com.br")

# Configura APIs
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)

Entrez.email = ENTREZ_EMAIL

# CSV de Journals Q1
CSV_JCR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "CarlosCamillo_JCR_JournalResults_12_2025.csv")

# ======================================================================
# 1. CARREGAMENTO E FILTRAGEM DE JOURNALS (Q1 e Q2)
# ======================================================================

def load_journal_rankings(csv_path):
    """
    Carrega o CSV do JCR e retorna um dict: { "JOURNAL NAME": "Q1" or "Q2" ... }
    Ignora Q3/Q4.
    """
    journal_map = {}
    if not os.path.exists(csv_path):
        print(f"❌ Erro: Arquivo CSV não encontrado em {csv_path}")
        return journal_map

    print(f"📂 Carregando Rankings de: {csv_path}")
    try:
        with open(csv_path, mode='r', encoding='utf-8') as f:
            # Skip preamble lines until we find the header
            lines = f.readlines()
            start_line = 0
            for i, line in enumerate(lines):
                if line.strip().startswith("Journal name") or line.strip().startswith('"Journal name"'):
                    start_line = i
                    break
            
            # Re-join from header onwards to feed DictReader
            content = lines[start_line:]
            
            reader = csv.DictReader(content)
            count_q1 = 0
            count_q2 = 0
            
            for row in reader:
                journal_name = (row.get("Journal name") or "").strip().upper()
                quartile = (row.get("JIF Quartile") or "").strip().upper()
                category = (row.get("Category") or "").strip().upper()
                
                # FILTRO DE CATEGORIA: Exclui Reabilitação conforme pedido do usuário (07/01/2026)
                if "REHABILITATION" in category:
                    continue

                # PRIORIZA Q1:
                # Se o journal já existe no mapa:
                # - Se já é Q1, mantém Q1 (ignora se o atual for Q2).
                # - Se é Q2 e o atual é Q1, atualiza para Q1.
                # - Se não existe, adiciona.
                
                if journal_name:
                    current_val = journal_map.get(journal_name)
                    
                    if quartile == "Q1":
                        if current_val != "Q1":
                            journal_map[journal_name] = "Q1"
                            count_q1 += 1
                            if current_val == "Q2": count_q2 -= 1 # Ajusta contagem se fez upgrade
                    
                    elif quartile == "Q2":
                        if journal_name not in journal_map:
                            journal_map[journal_name] = "Q2"
                            count_q2 += 1
                        
            print(f"✅ Carregados {count_q1} journals Q1 e {count_q2} journals Q2.")
    except Exception as e:
        print(f"❌ Erro ao ler CSV: {e}")
    
    return journal_map

# Carregamento Global
JOURNAL_RANKING_MAP = {}

# ======================================================================
# 2. DEFINIÇÃO DAS CONSULTAS (QUERIES) MÉDICAS
# ======================================================================
# Estratégia: (Doença) AND (Tratamento/Cirurgia) AND (Alta Evidência) NOT (Exercício/Reabilitação)
# Subheadings úteis do MeSH: /drug therapy, /surgery, /therapy
# Exclusão: "Exercise"[Mesh], "Rehabilitation"[Mesh], "Physical Therapy Modalities"[Mesh]

FILTRO_EVIDENCIA = '(randomized controlled trial[pt] OR systematic review[pt] OR meta-analysis[pt] OR practice guideline[pt])'
FILTRO_EXCLUSAO = 'NOT ("Exercise"[Mesh] OR "Exercise Therapy"[Mesh] OR "Rehabilitation"[Mesh] OR "Physical Therapy Modalities"[Mesh] OR "exercise"[tiab] OR "training"[tiab] OR "rehabilitation"[tiab])'

# Montando queries específicas
def build_query(core_term):
    return f'({core_term}) AND {FILTRO_EVIDENCIA} {FILTRO_EXCLUSAO} AND humans[Mesh]'

CONSULTAS_MEDICAS = {
    "Asma": build_query('("Asthma/drug therapy"[Mesh] OR "Asthma/surgery"[Mesh] OR "Asthma/therapy"[Mesh])'),
    
    "DPOC": build_query('("Pulmonary Disease, Chronic Obstructive/drug therapy"[Mesh] OR "Pulmonary Disease, Chronic Obstructive/surgery"[Mesh] OR "Pulmonary Disease, Chronic Obstructive/therapy"[Mesh])'),
    
    "Doenças Intersticiais": build_query('("Lung Diseases, Interstitial/drug therapy"[Mesh] OR "Lung Diseases, Interstitial/surgery"[Mesh] OR "Lung Diseases, Interstitial/therapy"[Mesh] OR "Idiopathic Pulmonary Fibrosis/drug therapy"[Mesh])'),
    
    "Hipertensão Pulmonar": build_query('("Hypertension, Pulmonary/drug therapy"[Mesh] OR "Hypertension, Pulmonary/surgery"[Mesh] OR "Hypertension, Pulmonary/therapy"[Mesh])'),
    
    "Bronquiectasias": build_query('("Bronchiectasis/drug therapy"[Mesh] OR "Bronchiectasis/surgery"[Mesh] OR "Bronchiectasis/therapy"[Mesh])'),
    
    "Broncoscopia/Intervencionista": build_query('("Bronchoscopy"[Mesh] OR "Procedures and Techniques"[Mesh]) AND ("Lung Diseases"[Mesh])'),

    "Infecções Respiratórias": build_query('("Respiratory Tract Infections/drug therapy"[Mesh] OR "Pneumonia/drug therapy"[Mesh])')
}

# ======================================================================
# 3. FUNÇÕES DE BUSCA E TRADUÇÃO (Reimplementadas para isolamento)
# ======================================================================

def buscar_ids(query, dias_atras=7):
    """Busca IDs dos últimos X dias."""
    tz = pytz.timezone("America/Sao_Paulo")
    agora = datetime.now(tz)
    data_final = agora.date()
    data_inicial = data_final - timedelta(days=dias_atras)
    
    mindate = data_inicial.strftime("%Y/%m/%d")
    maxdate = data_final.strftime("%Y/%m/%d")
    
    print(f"   🔎 Buscando (PubMed): {mindate} a {maxdate}")
    
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=50, # Limite razoável por tema para não estourar
            mindate=mindate,
            maxdate=maxdate,
            datetype="pdat"
        )
        record = Entrez.read(handle)
        return record["IdList"]
    except Exception as e:
        print(f"   ⚠️ Erro na busca: {e}")
        return []

def buscar_detalhes(ids):
    if not ids: return []
    try:
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="xml", retmode="xml")
        records = Entrez.read(handle)
        artigos = []
        
        for article in records['PubmedArticle']:
            art_data = article['MedlineCitation']['Article']
            pmid = article['MedlineCitation']['PMID']
            
            # Dados básicos
            titulo = art_data.get('ArticleTitle', '')
            journal = art_data.get('Journal', {}).get('Title', '')
            
            # Autores
            autores_list = art_data.get('AuthorList', [])
            autores = [f"{a.get('LastName','')}" for a in autores_list if 'LastName' in a]
            
            # DOI
            doi = ''
            for eid in article['PubmedData']['ArticleIdList']:
                if eid.attributes.get('IdType') == 'doi':
                    doi = str(eid)
                    break
            
            # Abstract
            abstract_text = ""
            abs_data = art_data.get('Abstract', {}).get('AbstractText', [])
            if isinstance(abs_data, list):
                # Concatena partes (ex: BACKGROUND: ..., METHODS: ...)
                texto_parts = []
                for item in abs_data:
                    if hasattr(item, 'attributes') and 'Label' in item.attributes:
                        texto_parts.append(f"{item.attributes['Label']}: {item}")
                    else:
                        texto_parts.append(str(item))
                abstract_text = "\n\n".join(texto_parts)
            else:
                abstract_text = str(abs_data)

            # Metadata Data
            journal_info = art_data['Journal']['JournalIssue']
            ano = journal_info['PubDate'].get('Year', '2025')
            
            artigos.append({
                'pmid': pmid,
                'titulo': titulo,
                'journal': journal,
                'autores': autores,
                'doi': doi,
                'abstract': abstract_text,
                'ano': ano
            })
            
        return artigos
    except Exception as e:
        print(f"   ⚠️ Erro ao baixar detalhes: {e}")
        return []

def traduzir_resumo_medico(texto):
    """
    Tradução técnica para médicos. 
    Usa Gemini Flash ou OpenAI como fallback.
    """
    if not texto.strip(): return "Resumo não disponível."
    
    prompt = f"""
    Você é um assistente especialista em pneumologia.
    Traduza o seguinte resumo científico para o Português (Brasil).
    Público alvo: Médicos especialistas.
    Mantenha termos técnicos precisos.
    Não resuma, faça uma tradução fiel e fluida.
    
    Texto Original:
    {texto}
    """
    
    try:
        # Tenta Gemini primeiro (mais rápido/barato)
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        resp = model.generate_content(prompt)
        return resp.text.strip()
    except Exception as e:
        # Fallback OpenAI
        # print(f"   (Gemini falhou, usando GPT): {e}")
        try:
            client = OpenAI(api_key=OPENAI_API_KEY)
            resp = client.chat.completions.create(
                model="gpt-4o",
                messages=[{"role": "user", "content": prompt}],
                temperature=0.1
            )
            return resp.choices[0].message.content.strip()
        except Exception as e2:
            return f"Erro na tradução: {e2}"

# ======================================================================
# 4. GERAÇÃO DE HTML (MAILCHIMP)
# ======================================================================

TEMPLATE_HTML_MEDICO = """
<!DOCTYPE html>
<html lang="pt-br">
<head><meta charset="UTF-8"><title>Boletim Médico Semanal | Revalidatie</title></head>
<body style="background:#fcfdff;margin:0;padding:0;">
  <div style="text-align:center;font-size:14px;margin-top:20px;"><a href="*|ARCHIVE|*" style="color:#065e77;text-decoration:underline;">Ver este e-mail no seu navegador</a></div>
  
  <div style="text-align:center;margin:30px 0 10px 0;"><a href="https://www.revalidatie.com.br" target="_blank"><img src="https://i.imgur.com/6FIUeHX.png" alt="Logo Revalidatie" style="max-width:270px;width:100%;height:auto;"></a></div>
  
  <!-- CABEÇALHO NOVO (Usuario pediu imgur/5e12tRk) -->
  <div style="width:100%;max-width:900px;margin:auto;text-align:center;"><img src="https://i.imgur.com/5e12tRk.png" alt="Reva Medico Weekly" style="width:100%;max-width:900px;height:auto;border-radius:16px;"></div>

  <table align="center" border="0" cellpadding="0" cellspacing="0" width="92%" style="max-width:760px; margin:auto; background:#fff;">
    <tr><td style="padding: 36px 20px 0 20px; text-align: center;"><h1 style="margin:0 0 10px 0;font-size:2.4em;color:#205776;font-family:Helvetica,Arial,sans-serif;font-weight:bold;">Olá, *|FNAME|*!</h1></td></tr>
    <tr>
      <td style="padding:0 20px 36px 20px; text-align:left;">
        <span style="color:#407ca6;font-size:1.07em;font-family:Helvetica,Arial,sans-serif;">Confira a seleção semanal de artigos publicados nos periódicos de maior impacto (Q1 e Q2) em Pneumologia e Cirurgia Torácica.</span>
        <div style="height:18px;"></div>
        <div style="font-size:1.1em; color:#111;font-family:Helvetica,Arial,sans-serif;">{conteudo}</div>
        <div style="height:34px;"></div>
      </td>
    </tr>
  </table>
  
  <div style="background:#222C36;color:#fff;padding:36px 0 24px 0;text-align:center;font-size:1.08em;font-family:Helvetica,Arial,sans-serif;">
    <div style="margin-bottom:12px;"><img src="https://i.imgur.com/6FIUeHX.png" alt="Revalidatie" style="max-width:180px;width:100%;height:auto;"></div>
    <div style="margin-bottom:15px;">Copyright (C ) 2025 Revalidatie. Todos os direitos reservados.
      <br>Você está recebendo este email porque se inscreveu na lista Médica da Revalidatie.</div>
    <div class="disclaimer" style="color:#ddd;font-size:0.96em;">Quer alterar como recebe estes emails?
      <a href="*|UPDATE_PROFILE|*" style="color:#5beaff;">Atualizar preferências</a> ou <a href="*|UNSUB|*" style="color:#5beaff;">descadastrar</a></div>
  </div>
</body>
</html>
"""

def formatar_artigo_html(artigo, resumo_traduzido, quartil):
    # Formata metadados
    info_journal = artigo.get('journal', 'N/A')
    info_ano = artigo.get('ano', 'N/A')
    
    # Badge do Quartil
    q_badge = f"<span style='background:#f0f0f0; color:#444; padding:2px 6px; border-radius:4px; font-size:0.8em; margin-left:5px; border:1px solid #ccc;'>{quartil}</span>"

    autores = ', '.join(artigo['autores']) if artigo['autores'] else "Autores não informados"
    
    link_doi_url = f"https://doi.org/{artigo['doi']}" if artigo['doi'] else "#"
    link_pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{artigo['pmid']}/"

    links_finais = f'<strong>DOI:</strong> <a href="{link_doi_url}" target="_blank" style="color:#065e77;text-decoration:underline;">{link_doi_url}</a>'
    links_finais += f' | <strong>PubMed:</strong> <a href="{link_pubmed_url}" target="_blank" style="color:#065e77;text-decoration:underline;">{link_pubmed_url}</a>'

    # Lista de metadados estilo RevaCast
    html_meta = f"""
    <ul style="margin-left: 0; padding-left: 20px; list-style-type: disc; font-family: Helvetica, Arial, sans-serif; color: #333;">
        <li style="margin-bottom: 5px;"><strong>Título:</strong> <a href="{link_pubmed_url}" style="text-decoration:none; color:#205776; font-weight:bold;">{artigo['titulo']}</a></li>
        <li style="margin-bottom: 5px;"><strong>Journal:</strong> {info_journal} {q_badge} ({info_ano})</li>
        <li style="margin-bottom: 5px;"><strong>Autores:</strong> {autores}</li>
        <li style="margin-bottom: 5px;">{links_finais}</li>
    </ul>
    """

    # Formata Resumo (parágrafos)
    texto = resumo_traduzido.replace("\r\n", "\n").replace("\r", "\n")
    paragrafos = [p.strip() for p in texto.split("\n") if p.strip()]
    if not paragrafos: 
        html_paragrafos = '<p>Resumo não disponível.</p>'
    else:
        html_paragrafos = "".join(
            f'<p style="font-family: Helvetica, Arial, sans-serif; color: #333; line-height:1.5; margin: 0 0 10px 0;">{p}</p>'
            for p in paragrafos
        )

    return f"""
    <div style="margin-bottom: 30px; padding-bottom: 25px; border-bottom: 1px solid #cccccc;">
        {html_meta}
        <div style="background:#f8f9fa; padding:15px; border-radius:6px; border-left:4px solid #205776;">
             {html_paragrafos}
        </div>
    </div>
    """

# ======================================================================
# 5. LÓGICA PRINCIPAL
# ======================================================================

def rodar_boletim_medico(dry_run=False, output_file="boletim_medico_preview.html", log_callback=print):
    """
    Executa o pipeline completo com estratégia de fallback.
    1. Carrega Mapa de Quartis (Q1, Q2).
    2. Coleta TODOS os candidatos por tema (Meta-data only).
    3. Separa em listas Q1 e Q2.
    4. Se Total Q1 > 0 -> Processa apenas Q1.
    5. Se Total Q1 == 0 -> Processa Q1 + Q2.
    6. Traduz e Gera HTML.
    """
    log_callback("\n🚀 INICIANDO REVA MED WEEKLY PIPELINE")
    
    # 1. Load Rankings
    global JOURNAL_RANKING_MAP
    JOURNAL_RANKING_MAP = load_journal_rankings(CSV_JCR_PATH)
    
    # Armazena candidatos: { 'Tema': {'q1': [art...], 'q2': [art...]} }
    candidatos_por_tema = {}
    total_q1_encontrados = 0
    
    log_callback("\n🕵️  COLETANDO CANDIDATOS (Passo 1/2)...")
    for tema, query in CONSULTAS_MEDICAS.items():
        log_callback(f"   > Tema: {tema}")
        candidatos_por_tema[tema] = {'q1': [], 'q2': []}
        
        ids = buscar_ids(query)
        if not ids: continue
        
        artigos = buscar_detalhes(ids)
        for art in artigos:
            journal_upper = art['journal'].upper().strip()
            
            # Identifica Quartil
            quartil = JOURNAL_RANKING_MAP.get(journal_upper)
            
            if not quartil:
                # Debug leve
                # log_callback(f"     [Ignorado - Sem Quartil] {journal_upper}")
                continue
            
            # FILTRO EXERCICIO (Segurança extra)
            texto_chk = (art['titulo'] + " " + art['abstract']).lower()
            termos_banidos = ["physical therapy", "physiotherapy", "rehabilitation program", "exercise training", "aerobic training"]
            if any(termo in texto_chk for termo in termos_banidos):
                # log_callback(f"     [Ignorado - Termo Banido] {art['titulo'][:30]}...")
                continue
            
            if quartil == 'Q1':
                candidatos_por_tema[tema]['q1'].append(art)
                total_q1_encontrados += 1
            elif quartil == 'Q2':
                candidatos_por_tema[tema]['q2'].append(art)
    
    # DECISÃO DE ESCOPO
    usar_q2 = False
    if total_q1_encontrados == 0:
        log_callback("\n⚠️  NENHUM ARTIGO Q1 ENCONTRADO! Ativando o modo de contingência (Incluindo Q2)...")
        usar_q2 = True
    else:
        log_callback(f"\n✅ Encontrados {total_q1_encontrados} artigos Q1. Processando apenas Q1.")

    # PROCESSAMENTO E TRADUÇÃO
    log_callback("\n🧠 TRADUZINDO E GERANDO BOLETIM (Passo 2/2)...")
    html_corpo = ""
    total_final = 0
    
    for tema, grupos in candidatos_por_tema.items():
        # Lista final para este tema
        lista_arts = grupos['q1']
        if usar_q2:
            lista_arts.extend(grupos['q2'])
            
        if not lista_arts: continue
        
        html_corpo += f"<h2 style='color: #d9534f; margin-top: 40px; border-bottom: 2px solid #d9534f;'>{tema}</h2>"
        
        for art in lista_arts:
            q_label = JOURNAL_RANKING_MAP.get(art['journal'].upper().strip(), 'Q?')
            log_callback(f"     - [{q_label}] {art['titulo'][:60]}...")
            
            resumo_pt = traduzir_resumo_medico(art['abstract'])
            html_corpo += formatar_artigo_html(art, resumo_pt, q_label)
            total_final += 1

    # Finaliza HTML
    if total_final == 0:
        log_callback("\n❌ Nenhum artigo selecionado em nenhum tema (Mesmo com Q2 se ativado).")
        return
        
    html_final = TEMPLATE_HTML_MEDICO.format(conteudo=html_corpo)
    
    # Salva arquivos locais
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_final)
    log_callback(f"\n💾 Arquivo HTML gerado: {output_file}")
    

    # Mailchimp
    if not dry_run and MC_API_KEY:
        log_callback("\n📧 Criando e AGENDANDO campanha no Mailchimp...")
        try:
            mc = Client()
            mc.set_config({"api_key": MC_API_KEY, "server": MC_SERVER})
            
            hoje_str = datetime.now().strftime("%d/%m/%Y")
            subject = f"Revalida Medical - {hoje_str}"
            if usar_q2:
                 subject += " (Edição Ampliada)"
            else:
                 subject += " (Destaques Q1)"

            campaign = mc.campaigns.create({
                "type": "regular",
                "recipients": {"list_id": MC_LIST_ID},
                "settings": {
                    "subject_line": subject,
                    "title": f"Medical Weekly {hoje_str}",
                    "from_name": MC_FROM_NAME,
                    "reply_to": MC_REPLY_TO
                }
            })
            cid = campaign['id']
            mc.campaigns.set_content(cid, {"html": html_final})
            
            # AGENDAMENTO
            schedule_time = get_next_friday_730()
            try:
                mc.campaigns.schedule(cid, {"schedule_time": schedule_time})
                log_callback(f"✅ Campanha criada e AGENDADA com sucesso para {schedule_time}! ID: {cid}")
            except Exception as e_sched:
                log_callback(f"⚠️ Campanha criada (ID: {cid}), mas falha no agendamento (ficou como Rascunho): {e_sched}")

        except Exception as e:
            log_callback(f"❌ Erro Mailchimp: {e}")
    elif dry_run:
        log_callback("\n🛑 Dry run: Mailchimp ignorado.")

def get_next_friday_730():
    """Retorna a próxima Sexta-feira às 07:30 AM (UTC)."""
    # Mailchimp schedules in UTC usually, let's just stick to ISO format. 
    # But usually API expects UTC. Let's assume server time is Brazil (-3).
    # To simplicity, we will calculate next Friday 07:30 Local and convert, or simply assume local if MC timezone settings apply.
    # Mailchimp API "schedule_time" requires UTC datetime in ISO 8601 format (e.g. 2026-06-28T14:47:00+00:00)
    
    tz = pytz.timezone("America/Sao_Paulo")
    now = datetime.now(tz)
    
    # 0 = Seg, 4 = Sex
    days_ahead = (4 - now.weekday() + 7) % 7
    # Se hoje é sexta e já passou das 7:30, joga pra próxima semana
    if days_ahead == 0 and (now.hour > 7 or (now.hour == 7 and now.minute >= 30)):
        days_ahead = 7
        
    next_friday = now + timedelta(days=days_ahead)
    target_time = next_friday.replace(hour=7, minute=30, second=0, microsecond=0)
    
    # Convert to UTC for Mailchimp API
    target_utc = target_time.astimezone(pytz.utc)
    return target_utc.isoformat()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    rodar_boletim_medico(dry_run=args.dry_run)
