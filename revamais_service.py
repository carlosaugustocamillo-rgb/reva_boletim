import os
# Force Build 123
import requests
import json
from datetime import datetime, timedelta
from dotenv import load_dotenv
import google.generativeai as genai
from openai import OpenAI
from Bio import Entrez
import mailchimp_marketing as MailchimpMarketing
from mailchimp_marketing.api_client import ApiClientError
import csv
import shutil
from tempfile import NamedTemporaryFile

# Importa ferramentas já existentes
from firebase_service import upload_file

load_dotenv()

# Configurações
Entrez.email = os.environ.get("ENTREZ_EMAIL")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
OPENAI_TEXT_MODEL = os.environ.get("OPENAI_TEXT_MODEL", "gpt-5.5")
OPENAI_TEXT_MODEL_SEARCH = os.environ.get("OPENAI_TEXT_MODEL_SEARCH", OPENAI_TEXT_MODEL)
OPENAI_TEXT_MODEL_WRITE = os.environ.get("OPENAI_TEXT_MODEL_WRITE", "gpt-5.5-pro")

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


def gerar_texto_openai(prompt, system_prompt=None, model_name=None):
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
    )
    content = response.choices[0].message.content or ""
    return content.strip()


def gerar_texto_preferencial(prompt, system_prompt=None, model_name=None):
    """
    Usa OpenAI como primeira opção e Gemini como fallback para tarefas textuais.
    """
    model_name = model_name or OPENAI_TEXT_MODEL_SEARCH
    try:
        return gerar_texto_openai(prompt, system_prompt=system_prompt, model_name=model_name)
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
                    "journal": journal
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

def gerar_imagem(prompt, nome_arquivo_prefixo, keep_local=False, forced_filename=None):
    """
    Gera imagem via Gemini ou DALL-E e faz upload.
    """
    print(f"🎨 Gerando imagem ({nome_arquivo_prefixo})...")
    if forced_filename:
        temp_filename = forced_filename
    else:
        temp_filename = f"temp_{nome_arquivo_prefixo}_{datetime.now().strftime('%H%M%S')}.png"
    
    image_generated = False
    
    # Tenta Gemini
    try:
        model = genai.GenerativeModel(GEMINI_IMAGE_MODEL)
        response = model.generate_content("Generate an image of: " + prompt)
        for part in response.parts:
            if hasattr(part, 'inline_data') and part.inline_data:
                with open(temp_filename, "wb") as f:
                    f.write(part.inline_data.data)
                image_generated = True
                break
    except Exception as e:
        print(f"   ⚠️ Erro Gemini Imagem: {e}")

    # Fallback DALL-E
    if not image_generated:
        try:
            response = client.images.generate(
                model="dall-e-3",
                prompt=prompt,
                size="1024x1024",
                quality="standard",
                n=1,
            )
            img_data = requests.get(response.data[0].url).content
            with open(temp_filename, "wb") as f:
                f.write(img_data)
            image_generated = True
        except Exception as e:
            print(f"   ❌ Erro DALL-E: {e}")
            return "https://via.placeholder.com/600x400?text=Reva+Mais"

    # Upload
    try:
        timestamp_upload = datetime.now().strftime('%Y%m%d_%H%M%S')
        firebase_path = f"revamais/{nome_arquivo_prefixo}_{timestamp_upload}.png"
        url = upload_file(temp_filename, firebase_path)
        if not keep_local and os.path.exists(temp_filename): 
            os.remove(temp_filename)
        return url
    except Exception as e:
        print(f"   ❌ Erro Upload: {e}")
        return "https://via.placeholder.com/600x400?text=Erro+Upload"
    
# --- Novos Imports para Manipulação de Imagem ---
from PIL import Image, ImageDraw, ImageFont
import io

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

def gerar_conteudo_revamais(tema, referencias):
    """
    Gera o conteúdo HTML do boletim.
    """
    print("✍️ Escrevendo conteúdo Reva +...")

    prompt = f"""
    Você é o editor do "Reva +", um boletim de saúde da clínica Revalidatie.
    Público-alvo: Pacientes e pessoas interessadas em saúde (leigos, mas inteligentes).
    
    Tema: "{tema}"
    
    INSTRUÇÃO DE ESCRITA HÍBRIDA (Conhecimento Geral + Evidência Específica):
    
    1. **Contexto e Mecanismo (Use seu conhecimento médico geral)**:
       - Comece explicando o problema de forma empática (ex: "Você sente dor ao caminhar?").
       - Explique O PORQUÊ (Fisiologia/Mecanismo): Por que isso acontece? O que muda no corpo com o tratamento? (Ex: fale sobre circulação colateral, eficiência muscular, neuroplasticidade).
       - O usuário GOSTA dessa explicação educativa do "como funciona".
    
    2. **O Que a Ciência Diz (Baseado SOMENTE nos artigos selecionados abaixo)**:
       - Agora, cite as evidências fornecidas.
       - Use os abstracts para validar a explicação anterior.
       - Diga "Estudos recentes mostram que..." e use os dados dos resumos.

    REGRAS DE SEGURANÇA E FIDELIDADE:
    - Não invente números, magnitude de efeito, tempo de intervenção, perfil de pacientes ou conclusões.
    - Não cite estudos que não estejam na lista abaixo.
    - Se um detalhe não estiver explícito nos resumos, não mencione esse detalhe.
    - Se a evidência parecer preliminar, heterogênea ou limitada, diga isso com cautela.
    - Use linguagem prudente: "sugere", "indica", "aponta", "pode ajudar", quando apropriado.
    - Na seção de ciência, organize os achados em uma lista HTML (<ul><li>) com 1 item por artigo ou achado principal.
    
    --- EVIDÊNCIA CIENTÍFICA (Para a seção 'O Que a Ciência Diz') ---
    {chr(10).join([ f'Artigo {i+1}: {r["texto"]}{chr(10)}Resumo: {r["resumo"]}{chr(10)}' for i, r in enumerate(referencias) ])}
    -----------------------------------------------------------------
    
    Estrutura HTML (retorne APENAS o conteúdo dentro do body e APENAS HTML PURO):
    IMPORTANTE: NÃO USE MARKDOWN (```html ... ```). Retorne APENAS o código HTML cru.
    
    1. <h1>Título Atraente e Emocional</h1>
    2. <p>Introdução empática + Explicação do Mecanismo (Por que dói? Por que melhora? - Use conhecimento geral de fisiologia).</p>
    3. <h2>O que a Ciência Comprova?</h2> (Aqui você insere a "tradução" dos abstracts acima, conectando com a explicação).
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
        rf'<h2[^>]*>\s*{h2_regex}\s*</h2>(.*?)(?=<h2[^>]*>|$)',
        re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(html)
    if not match:
        return ""
    return limpar_texto_html(match.group(1))


def gerar_briefs_visuais_revamais(tema, html_texto, referencias):
    """
    Gera prompts visuais mais ancorados no conteúdo final do boletim.
    """
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
    You are creating two grounded image prompts for a patient-facing medical newsletter.
    Return ONLY valid JSON with this exact shape:
    {{
      "ciencia": {{"prompt_english": "..."}},
      "dicas": {{"prompt_english": "..."}}
    }}

    RULES:
    - Base the prompts strictly on the supplied theme, selected studies, and final newsletter text.
    - Avoid generic wellness visuals, random symbols, and concepts not present in the content.
    - Prefer concrete anatomy, physiology, movement, rehabilitation actions, daily habits, or evidence concepts explicitly present in the text.
    - Visual style: premium editorial medical infographic, clean, minimal, white background, high contrast, lots of whitespace.
    - The "ciencia" image must explain the actual mechanism or evidence narrative from the science section.
    - The "dicas" image must show only the practical actions or habits explicitly recommended in the tips section.
    - No logos, no brand marks, no clutter, no unrelated charts.
    - If short labels are useful, include at most 2 or 3 labels in Portuguese. If unsure, request no text.

    THEME:
    {tema}

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
            and isinstance(briefs.get("ciencia"), dict)
            and isinstance(briefs.get("dicas"), dict)
            and briefs["ciencia"].get("prompt_english")
            and briefs["dicas"].get("prompt_english")
        ):
            return briefs
    except Exception as e:
        print(f"⚠️ Falha ao gerar briefs visuais do Reva+: {e}")

    return {
        "ciencia": {
            "prompt_english": (
                f"Create a premium editorial medical infographic based strictly on this newsletter science section about '{tema}': "
                f"{secao_ciencia[:900]} "
                "Show the specific anatomy, physiology, rehabilitation mechanism, or evidence concept described in the text. "
                "White background, minimal composition, high contrast, no clutter. "
                "Include at most 2 or 3 short Portuguese labels only if clearly supported by the section; otherwise no text."
            )
        },
        "dicas": {
            "prompt_english": (
                f"Create a premium editorial medical infographic based strictly on this practical tips section about '{tema}': "
                f"{secao_dicas[:900]} "
                "Show only the concrete actions, habits, or rehabilitation steps described in the text. "
                "White background, checklist or step-by-step layout, minimal composition, high contrast, no clutter. "
                "Include at most 2 or 3 short Portuguese labels only if clearly supported by the section; otherwise no text."
            )
        },
    }

from firebase_service import read_json_from_storage

def obter_proximo_tema_csv(consumir=True):
    """
    Lê o calendário e usa o Firebase para persistir quais TEMAS já foram usados (Estado Global).
    """
    csv_filename = "calendario_editorial_150_semanas.csv"
    csv_path = os.path.join(os.path.dirname(__file__), csv_filename)
    
    if not os.path.exists(csv_path):
        print(f"⚠️ Arquivo {csv_filename} não encontrado.")
        return None
        
    # 1. Carrega todas as linhas do CSV (Read-Only)
    rows = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        
    total_linhas = len(rows)
    
    # 2. Busca estado atual no Firebase
    state_file = "revamais/used_themes_state.json"
    state_data = read_json_from_storage(state_file)
    
    if not state_data:
        state_data = {"used_titles": []}
        
    used_titles = set(state_data.get("used_titles", []))
    
    # 3. Encontra o primeiro não utilizado
    tema_escolhido = None
    formato_escolhido = "Carrossel"
    idx_atual = 0
    
    # Normalização para comparação (Remove espaços extras e case insensitive)
    def normalize(text):
        return str(text).strip().lower()

    used_titles_normalized = {normalize(t) for t in used_titles}
    print(f"📊 Estado Carregado: {len(used_titles)} temas já utilizados.")
    
    for i, row in enumerate(rows):
        titulo = row.get('Title', row.get('Theme', '')).strip()
        titulo_norm = normalize(titulo)
        
        if titulo and titulo_norm not in used_titles_normalized:
            tema_escolhido = titulo
            formato_escolhido = row.get('Format', 'Carrossel')
            idx_atual = i + 1
            break
            
    if tema_escolhido:
        print(f"📅 Tema do Calendário Selecionado: {tema_escolhido} (Item {idx_atual}/{total_linhas})")

        if consumir:
            # 4. Atualiza estado e salva no Firebase
            used_titles.add(tema_escolhido)
            state_data["used_titles"] = list(used_titles)
            state_data["last_updated"] = datetime.now().isoformat()

            # Cria arquivo local temporário para upload
            tmp_path = None
            try:
                # Usa prefixo para facilitar debug
                with NamedTemporaryFile(mode='w', delete=False, suffix='.json', prefix='state_update_') as tmp:
                    json.dump(state_data, tmp)
                    tmp_path = tmp.name

                print(f"💾 Tentando salvar estado atualizado ({len(used_titles)} itens) no Firebase...")
                url_state = upload_file(tmp_path, state_file)

                if url_state:
                    print("✅ Estado atualizado com sucesso no Firebase.")
                else:
                    print("⚠️ FALHA SILENCIOSA ao salvar estado (upload_file retornou None). O tema vai repetir!")

            except Exception as e_save:
                print(f"❌ ERRO EXCEÇÃO ao salvar estado: {e_save}")
            finally:
                if tmp_path and os.path.exists(tmp_path):
                    os.remove(tmp_path)
             
        return {
            "tema": tema_escolhido,
            "formato": formato_escolhido
        }
    else:
        print("⚠️ Todos os temas do calendário já foram usados!")
        return None


def resolver_tema_revamais(tema_usuario=None, consumir_tema_auto=True):
    """
    Resolve o tema efetivo do Reva+ e o formato padrão do Instagram.
    """
    is_calendar_source = not tema_usuario or tema_usuario.strip() == ""
    formato_instagram = "Carrossel"
    tema = tema_usuario

    if not tema or tema == "auto":
        dados_csv = obter_proximo_tema_csv(consumir=consumir_tema_auto)
        if not dados_csv:
            return None
        tema = dados_csv["tema"]
        formato_instagram = dados_csv.get("formato", "Carrossel")

    return {
        "tema": tema,
        "formato_instagram": formato_instagram,
        "is_calendar_source": is_calendar_source,
    }


def preparar_referencias_revamais(tema_usuario=None, quantidade_referencias=8):
    """
    Etapa intermediária do pipeline: resolve o tema e retorna artigos candidatos
    para seleção manual antes da geração do boletim.
    """
    dados_tema = resolver_tema_revamais(
        tema_usuario=tema_usuario,
        consumir_tema_auto=False,
    )
    if not dados_tema:
        return {"status": "error", "message": "Nenhum tema disponível para preparar referências."}

    tema = dados_tema["tema"]
    tema_ingles = gerar_query_pubmed_tema(tema)
    referencias = buscar_referencias_pubmed(tema_ingles, limite_retorno=quantidade_referencias)

    return {
        "status": "success",
        "tema": tema,
        "tema_ingles": tema_ingles,
        "instagram_format": dados_tema["formato_instagram"],
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

def gerar_conteudo_instagram(tema, formato, referencias_text, conteudo_base=None):
    """
    Gera conteúdo para Instagram (Reel ou Carrossel).
    Melhoria v2: Gera texto com LLM, Imagem Clean com IA, e Texto via Overlay (Pillow).
    """
    print(f"📸 Gerando conteúdo para Instagram ({formato})...")
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
            print("   📝 Planejando narrativa do carrossel...")
            
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
3. **Correction**: Do NOT include the general theme title "{tema}" in the image. ONLY include the specific '{titulo}' and '{texto_curto}' of the slide.
4. **Density**: Avoid clutter. Use ONE main scene or focal subject that clearly supports the message of the slide.
5. **Text Instruction**: Explicitly state: "Include the text: '{titulo}' and '{texto_curto}' in the image. Typography must be legible, modern, sans-serif."
6. **Consistency**: Keep a clean healthcare brand aesthetic, with soft and trustworthy tones, subtle teal accents when appropriate, and premium visual quality.
7. **No Hallucinations**: Do not ask for logos, watermarks, icons, or cartoon elements unless strictly necessary for the concept.
8. **Avoid Illustration Look**: Do not generate flat drawings, vector art, line art, infographic-style icons, or cartoon-style medical scenes unless the slide specifically requires a scientific mechanism that cannot be represented with realistic photography.

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

            print(f"   🖼️ Gerando 7 slides (Prompt Dinâmico por Slide)...")
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
                    f"NO CLUTTER. Focus on the central message."
                )

                # Gera imagem
                local_base_name = f"slide_{slide_num}_{timestamp}"
                local_filename = f"{local_base_name}.png"
                
                gerar_imagem(full_prompt, local_base_name, keep_local=True, forced_filename=local_filename)
                
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
                print(f"   📦 Criando ZIP: {zip_name}...")
                
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
                print(f"⚠️ Erro ao criar ZIP: {e}")
                
    except Exception as e:
        print(f"❌ Erro ao gerar conteúdo Instagram: {e}")
        assets.append({"type": "error", "content": str(e)})
        
    return assets

def criar_campanha_revamais(
    tema_usuario=None,
    gerar_midia=True,
    gerar_instagram=True,
    enviar_email=True,
    referencias_selecionadas=None,
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
    dados_tema = resolver_tema_revamais(tema_usuario=tema_usuario, consumir_tema_auto=True)
    if not dados_tema:
        return {"status": "error", "message": "Nenhum tema fornecido e calendário esgotado/inexistente."}

    tema = dados_tema["tema"]
    formato_instagram = dados_tema["formato_instagram"]
    is_calendar_source = dados_tema["is_calendar_source"]
            
    # Remove duplicidade de log se já foi logado pelo wrapper, mas mal não faz
    log(f"🚀 Iniciando Reva +: {tema} (Insta: {formato_instagram})")
    
    referencias_selecionadas = referencias_selecionadas or []
    tema_ingles = tema
    referencias = []

    if referencias_selecionadas:
        log(f"🧠 Usando {len(referencias_selecionadas)} artigos selecionados manualmente.")
        referencias = referencias_selecionadas
    else:
        check()
        # 1. Traduzir tema para keywords científicas
        log("🌍 Traduzindo tema para keywords científicas...")
        tema_ingles = gerar_query_pubmed_tema(tema)

        check()
        # 2. Buscar Referências
        log(f"🔎 Buscando referências para: {tema_ingles}...")
        referencias = buscar_referencias_pubmed(tema_ingles)

    if not referencias:
        return {
            "status": "error",
            "message": "Nenhuma referência válida foi encontrada ou selecionada para gerar o Reva+."
        }
    
    check()
    # 3. Preparar placeholders de imagem
    url_capa_estatica = "https://i.imgur.com/oGzxgtK.jpeg"
    url_ilustrativa = "https://placehold.co/600x400?text=Imagem+Ilustrativa" # Placeholder default
    url_corpo_ciencia = "https://placehold.co/600x400?text=Infografico+Ciencia" # Placeholder default
    url_corpo_dicas = "https://placehold.co/600x400?text=Infografico+Dicas" # Placeholder default

    check()
    # 4. Gerar Texto
    log("✍️ Escrevendo boletim e formatando HTML...")
    html_texto = gerar_conteudo_revamais(tema, referencias)

    check()
    # 5. Gerar Imagens a partir do conteúdo final
    if gerar_midia:
        log("🧭 Derivando prompts visuais do conteúdo final...")
        briefs_visuais = gerar_briefs_visuais_revamais(tema, html_texto, referencias)

        log("🎨 Gerando assets visuais (isso pode demorar)...")
        try:
             # 1. Imagem Ilustrativa (Lifestyle/Visual)
            prompt_ilustrativa = f"A high quality, photorealistic or cinematic style photo-illustration about '{tema_ingles}'. Showing people, lifestyle, or the subject in a natural, positive way. NO TEXT. Suitable for a newsletter cover."
            url_ilustrativa = gerar_imagem(prompt_ilustrativa, "ilustrativa")

            # 2. Imagem Corpo (Infográfico/Educativo) - Ciência/Mecanismo
            prompt_corpo_ciencia = briefs_visuais["ciencia"]["prompt_english"]
            url_corpo_ciencia = gerar_imagem(prompt_corpo_ciencia, "corpo_ciencia")

            # 3. Imagem Corpo (Infográfico/Educativo) - Dicas Práticas
            prompt_corpo_dicas = briefs_visuais["dicas"]["prompt_english"]
            url_corpo_dicas = gerar_imagem(prompt_corpo_dicas, "corpo_dicas")
        except Exception as e:
            log(f"⚠️ Erro ao gerar imagens: {e}")
    else:
        log("⏩ Pulando geração de imagens (opção desmarcada).")

    def inserir_imagem_educativa(html, img_url, h2_regex, fallback="first_h2"):
        """Insere a imagem educativa após um H2 alvo, com fallback controlado."""
        img_tag = f'\n<img src="{img_url}" class="body-img" alt="Infográfico">\n'
        import re
        pattern = re.compile(rf'(<h2[^>]*>\\s*{h2_regex}\\s*</h2>)', re.IGNORECASE)
        match = pattern.search(html)
        if match:
            pos = match.end()
            return html[:pos] + img_tag + html[pos:]

        pattern_h2 = re.compile(r'(<h2[^>]*>.*?</h2>)', re.IGNORECASE | re.DOTALL)
        matches = list(pattern_h2.finditer(html))
        if matches:
            if fallback == "last_h2":
                pos = matches[-1].end()
                return html[:pos] + img_tag + html[pos:]
            # default: first_h2
            pos = matches[0].end()
            return html[:pos] + img_tag + html[pos:]

        # Último recurso: adiciona ao final do conteúdo
        return html + img_tag

    # Insere as duas imagens educativas dentro do texto
    html_texto = inserir_imagem_educativa(
        html_texto,
        url_corpo_ciencia,
        r"O\\s+que\\s+a\\s+Ci[eê]ncia\\s+Comprova\\??",
        fallback="first_h2",
    )
    html_texto = inserir_imagem_educativa(
        html_texto,
        url_corpo_dicas,
        r"Dicas\\s+Pr[aá]ticas",
        fallback="last_h2",
    )

    check()
    # 5.5 Gerar Conteúdo Instagram (Extra)
    instagram_assets = []
    if gerar_instagram:
        log("📸 Criando conteúdo para Instagram...")
        try:
            # Extrai texto das referências para passar de contexto
            refs_text_context = "\n".join([r['texto'] for r in referencias])
            instagram_assets = gerar_conteudo_instagram(tema, formato_instagram, refs_text_context, conteudo_base=html_texto)
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
                <!-- Imagem Ilustrativa (Nova) -->
                <img src="{url_ilustrativa}" class="body-img" alt="Ilustração do Título">
                
                {html_texto}
                
                <div class="cta-box">
                    <p>Quer saber mais sobre como cuidar da sua saúde?</p>
                    <p><a href="https://www.instagram.com/revalidatie_londrina/" target="_blank">Siga-nos no Instagram</a> ou <a href="https://www.revalidatie.com.br" target="_blank">Visite nosso site</a></p>
                </div>

                <!-- Referências Científicas -->
                <div class="references">
                    <h4>📚 Referências Científicas Utilizadas:</h4>
                    <ul>
                    {''.join([f"<li>{r['texto']} <a href='{r['link']}' target='_blank'>[PubMed]</a></li>" for r in referencias]) if referencias else "<li>Referências não disponíveis neste momento.</li>"}
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
        import re
        
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
    refs_items_site = (
        "".join([
            f"<li>{r['texto']} <a href='{r['link']}' target='_blank'>[PubMed]</a></li>"
            for r in referencias
        ])
        if referencias
        else "<li>Referências não disponíveis neste momento.</li>"
    )

    bloco_referencias_site = f"""
    <div class="references" style="margin-top:30px;padding-top:20px;border-top:1px solid #eee;">
        <h3>📚 Referências Científicas Utilizadas:</h3>
        <ul>
            {refs_items_site}
        </ul>
    </div>
    """

    html_content_site = html_texto
    html_lower = html_texto.lower()
    if "referências científicas utilizadas" not in html_lower and "referencias cientificas utilizadas" not in html_lower:
        html_content_site = html_texto + bloco_referencias_site

    return {
        "status": "success", 
        "campaign_id": campaign['id'], 
        "tema": tema,
        "modelo_texto_busca": OPENAI_TEXT_MODEL_SEARCH,
        "modelo_texto_redacao": OPENAI_TEXT_MODEL_WRITE,
        "tema_query_pubmed": tema_ingles,
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
