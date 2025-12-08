import os
import requests
import json
from datetime import datetime
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

# Configuração Gemini
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)

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

def buscar_referencias_pubmed(tema_ingles):
    """
    Busca 3 artigos recentes e relevantes no PubMed para servir de referência.
    Filtra por Revisões Sistemáticas e RCTs em revistas de alto impacto (Core Clinical Journals).
    """
    print(f"🔎 Buscando referências de alta qualidade para: {tema_ingles}...")
    
    # Filtro por tipos de estudo nobres e revistas 'Core' (jsubsetaim é um bom proxy para alto impacto)
    query = (
        f"({tema_ingles}) AND "
        f"(Systematic Review[pt] OR Randomized Controlled Trial[pt]) AND "
        f"jsubsetaim[text] AND "
        f"(last 5 years[dp])"
    )
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=3, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        
        if not ids:
            print("⚠️ Nenhuma referência 'Core' encontrada. Tentando busca mais ampla...")
            # Fallback sem o filtro de 'jsubsetaim' mas mantendo SR/RCT
            query_fallback = (
                f"({tema_ingles}) AND "
                f"(Systematic Review[pt] OR Randomized Controlled Trial[pt]) AND "
                f"(last 10 years[dp])" # Aumentei para 10 anos
            )
            # Remove jsubsetaim (revistas core) para ter mais chances
            handle = Entrez.esearch(db="pubmed", term=query_fallback, retmax=3, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            ids = record["IdList"]

        if not ids:
            return []

        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        referencias = []
        for article in records['PubmedArticle']:
            titulo = article['MedlineCitation']['Article']['ArticleTitle']
            autores_list = article['MedlineCitation']['Article'].get('AuthorList', [])
            primeiro_autor = f"{autores_list[0].get('LastName', '')} et al." if autores_list else "Autor desconhecido"
            journal = article['MedlineCitation']['Article']['Journal']['Title']
            ano = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate'].get('Year', 's/d')
            pmid = article['MedlineCitation']['PMID']
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            
            referencias.append({
                "texto": f"{primeiro_autor}. {titulo}. {journal}, {ano}.",
                "link": link
            })
            
        return referencias

    except Exception as e:
        print(f"⚠️ Erro ao buscar referências: {e}")
        return []

def gerar_imagem(prompt, nome_arquivo_prefixo):
    """
    Gera imagem via Gemini ou DALL-E e faz upload.
    """
    print(f"🎨 Gerando imagem ({nome_arquivo_prefixo})...")
    temp_filename = f"temp_{nome_arquivo_prefixo}_{datetime.now().strftime('%H%M%S')}.png"
    image_generated = False
    
    # Tenta Gemini
    try:
        model = genai.GenerativeModel('gemini-3-pro-image-preview')
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
        firebase_path = f"revamais/{nome_arquivo_prefixo}_{datetime.now().strftime('%Y-%m-%d')}.png"
        url = upload_file(temp_filename, firebase_path)
        if os.path.exists(temp_filename): os.remove(temp_filename)
        return url
    except Exception as e:
        print(f"   ❌ Erro Upload: {e}")
        return "https://via.placeholder.com/600x400?text=Erro+Upload"

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
    
    refs_html = ""
    if referencias:
        refs_html = "<h3>📚 Referências Científicas (Nível de Evidência Elevado):</h3><ul>"
        for ref in referencias:
            refs_html += f"<li>{ref['texto']} <a href='{ref['link']}' target='_blank'>[PubMed]</a></li>"
        refs_html += "</ul>"

    prompt = f"""
    Você é o editor do "Reva +", um boletim de saúde da clínica Revalidatie.
    Público-alvo: Pacientes e pessoas interessadas em saúde (leigos, mas inteligentes).
    
    Tema: "{tema}"
    
    Escreva um artigo curto, informativo e acolhedor.
    
    Estrutura HTML (retorne APENAS o conteúdo dentro do body, sem tags html/body):
    1. <h1>Título Atraente</h1>
    2. <p>Introdução conectando com o dia a dia.</p>
    3. <h2>O que a ciência diz?</h2> (Explicação simples, sem jargões difíceis).
    4. <h2>Dicas Práticas</h2> (Lista bullet points com conselhos aplicáveis).
    5. <div class="cta"> (Um texto convidando para seguir no Instagram @revalidatie ou visitar o site).
    
    Use tom otimista e baseado em evidências.
    """
    
    try:
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        response = model.generate_content(prompt)
        html_content = response.text
    except:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "user", "content": prompt}]
        )
        html_content = response.choices[0].message.content

    # Limpeza básica de markdown se houver
    html_content = html_content.replace("```html", "").replace("```", "")
    
    return html_content + "<hr>" + refs_html

    return html_content + "<hr>" + refs_html

def obter_proximo_tema_csv():
    """
    Lê o arquivo 'calendario_editorial_150_semanas.csv', pega o próximo tema
    não utilizado, marca como usado e retorna o Título.
    """
    csv_filename = "calendario_editorial_150_semanas.csv"
    csv_path = os.path.join(os.path.dirname(__file__), csv_filename)
    
    if not os.path.exists(csv_path):
        print(f"⚠️ Arquivo {csv_filename} não encontrado.")
        return None
        
    tema_escolhido = None
    rows = []
    headers = []
    
    # Lê todo o arquivo
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        # Garante que coluna 'Used' exista nos headers
        if 'Used' not in headers:
            headers.append('Used')
        
        for row in reader:
            rows.append(row)
            
    # Procura o primeiro não usado
    idx_escolhido = -1
    for i, row in enumerate(rows):
        is_used = (row.get('Used') or '').strip()
        if not is_used:
            tema_escolhido = row.get('Title', row.get('Theme', 'Tema Genérico'))
            formato_escolhido = row.get('Format', 'Carrossel')
            # Marca como usado
            rows[i]['Used'] = datetime.now().isoformat()
            idx_escolhido = i
            break
    
    if tema_escolhido:
        # Salva o arquivo atualizado
        with open(csv_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)
        print(f"📅 Tema do Calendário Selecionado: {tema_escolhido}")
        return {
            "tema": tema_escolhido,
            "formato": formato_escolhido
        }
    else:
        print("⚠️ Todos os temas do calendário já foram usados!")
        return None

def gerar_conteudo_instagram(tema, formato, referencias_text):
    """
    Gera conteúdo para Instagram baseado no formato (Reel ou Carrossel).
    Retorna lista de assets (URLs ou Texto).
    """
    print(f"📸 Gerando conteúdo para Instagram ({formato})...")
    assets = []
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    try:
        if formato.lower() == "reel":
            # Gera Roteiro
            prompt = f"""
            Crie um roteiro viral para Instagram Reels sobre: "{tema}".
            Baseado nestas referências: {referencias_text}
            
            Estrutura:
            1. Gancho (0-3s): Algo impactante visual e auditivo.
            2. Desenvolvimento (30-45s): Explicação rápida e dinâmica.
            3. Call to Action (CTA): Convite para ler a legenda ou seguir.
            4. Legenda sugerida para o post (com hashtags).
            
            Formato: Markdown.
            """
            model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
            roteiro = model.generate_content(prompt).text
            
            # Salva MD e sobe pro Firebase
            filename = f"instagram_reel_{timestamp}.md"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(roteiro)
            
            url = upload_file(filename, f"instagram/{filename}")
            if os.path.exists(filename): os.remove(filename)
            
            assets.append({"type": "roteiro", "url": url, "name": "Roteiro do Reel"})
            
        elif formato.lower() == "carrossel":
            # Gera 7 Imagens (Capa + Conteúdo Profundo + CTA)
            slides_prompts = [
                f"Slide 1 (Cover): Title '{tema}'. Clean, bold typography, medical background. Text in Portuguese.",
                f"Slide 2: Introduction/Context about '{tema}'. Short text in Portuguese.",
                f"Slide 3: Scientific Explanation (What science says). Educational diagram. Text in Portuguese.",
                f"Slide 4: Deep Dive/Mechanism (Why it happens). Detailed infographic style. Text in Portuguese.",
                f"Slide 5: Practical Tip #1 about '{tema}'. Iconography style. Text in Portuguese.",
                f"Slide 6: Practical Tip #2 or Checklist. Clean list style. Text in Portuguese.",
                f"Slide 7 (CTA): Conclusion & Text 'Gostou? Siga @revalidatie'. Text in Portuguese."
            ]
            
            print(f"   🖼️ Gerando 7 slides para carrossel...")
            for i, p_prompt in enumerate(slides_prompts):
                full_prompt = f"Create an Instagram Carousel Slide (1080x1080). {p_prompt}. High quality, professional health clinic branding (Teal/Lavender colors)."
                
                # Tenta gerar com Gemini 3 (reusa função existente mas adaptada ou chama direto)
                # Vamos chamar gerar_imagem mas precisamos garantir que use o Gemini 3 para texto
                # A função gerar_imagem já usa 'gemini-3-pro-image-preview'.
                
                slider_name = f"slide_{i+1}_{timestamp}"
                url = gerar_imagem(full_prompt, slider_name)
                assets.append({"type": "image", "url": url, "name": f"Slide {i+1}"})
                
    except Exception as e:
        print(f"❌ Erro ao gerar conteúdo Instagram: {e}")
        assets.append({"type": "error", "content": str(e)})
        
    return assets

def criar_campanha_revamais(tema=None, log_callback=None, check_cancel=None):
    """
    Cria campanha Reva+ com suporte a logs e cancelamento.
    log_callback: função(str) para enviar logs.
    check_cancel: função() -> bool que retorna True se deve cancelar.
    """
    def log(msg):
        print(msg) # Mantém print no stdout
        if log_callback: log_callback(msg)
        
    def check():
        if check_cancel and check_cancel():
            raise Exception("CANCELADO_PELO_USUARIO")

    formato_instagram = "Carrossel" # Default

    check()
    # Se não veio tema, tenta pegar do CSV
    if not tema or tema == "auto":
        log("🤖 Modo Automático: Buscando tema no calendário editorial...")
        dados_csv = obter_proximo_tema_csv()
        if not dados_csv:
            return {"status": "error", "message": "Nenhum tema fornecido e calendário esgotado/inexistente."}
        tema = dados_csv['tema']
        formato_instagram = dados_csv.get('formato', 'Carrossel')
            
    log(f"🚀 Iniciando Reva +: {tema} (Insta: {formato_instagram})")
    
    check()
    # 1. Traduzir tema para busca
    try:
        log("🌍 Traduzindo tema para busca científica...")
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        tema_ingles = model.generate_content(f"Translate to English for medical search: {tema}").text.strip()
    except:
        tema_ingles = tema

    check()
    # 2. Buscar Referências
    log(f"🔎 Buscando referências para: {tema_ingles}...")
    referencias = buscar_referencias_pubmed(tema_ingles)
    
    check()
    # 3. Gerar Imagens
    log("🎨 Gerando assets visuais (isso pode demorar)...")
    
    # Capa: SEM TEXTO, apenas visual
    prompt_capa = f"A welcoming, modern health illustration about '{tema_ingles}'. Soft lighting, professional photography style or high quality 3D render. Minimalist. NO TEXT. NO WORDS."
    url_capa = gerar_imagem(prompt_capa, "capa")
    check()
    
    # Imagem Corpo: Se tiver texto, DEVE SER EM PORTUGUÊS
    prompt_corpo = f"An educational infographic or diagram about '{tema_ingles}'. Clean lines, easy to understand, white background. IMPORTANT: Any text or labels MUST BE IN PORTUGUESE (PT-BR). If you cannot generate correct Portuguese text, do not include any text."
    url_corpo = gerar_imagem(prompt_corpo, "corpo")
    
    check()
    # 4. Gerar Texto
    log("✍️ Escrevendo boletim e formatando HTML...")
    html_texto = gerar_conteudo_revamais(tema, referencias)

    check()
    # 4.5 Gerar Conteúdo Instagram (Extra)
    log("📸 Criando conteúdo para Instagram...")
    # Extrai texto das referências para passar de contexto
    refs_text_context = "\n".join([r['texto'] for r in referencias])
    instagram_assets = gerar_conteudo_instagram(tema, formato_instagram, refs_text_context)
    
    check()
    # 5. Montar HTML Final
    # Logo da Revalidatie (usado no boletim_service.py)
    # Se o usuário definir uma no .env, usa. Senão, usa a padrão.
    DEFAULT_LOGO = "https://i.imgur.com/in5MW0g.png"
    logo_url = os.environ.get("REVA_LOGO_URL", DEFAULT_LOGO)
    
    today_formatted = datetime.now().strftime('%d/%m/%Y')
    
    html_email = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; line-height: 1.6; color: #333; max-width: 600px; margin: 0 auto; background-color: #f9f9f9; }}
            .container {{ background-color: #ffffff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-top: 20px; }}
            .logo-header {{ text-align: center; margin-bottom: 20px; }}
            .logo-header img {{ max-width: 200px; height: auto; }}
            .header-intro {{ text-align: center; margin-bottom: 20px; font-size: 14px; color: #555; }}
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
            <div class="logo-header">
                <a href="https://www.revalidatie.com.br" target="_blank">
                    <img src="{logo_url}" alt="Revalidatie">
                </a>
            </div>
            
            <div class="header-intro">
                <p><strong>Olá, *|FNAME|*!</strong><br>
                Aqui está sua atualização semanal de saúde do Reva + para {today_formatted}.</p>
            </div>
            
            <img src="{url_capa}" class="header-img" alt="Capa">
            
            <div style="padding: 20px;">
                {html_texto}
                
                <img src="{url_corpo}" class="body-img" alt="Ilustração">
                
                <div class="cta-box">
                    <p>Quer saber mais sobre como cuidar da sua saúde?</p>
                    <p><a href="https://www.instagram.com/revalidatie_londrina/" target="_blank">Siga-nos no Instagram @revalidatie_londrina</a> ou <a href="https://www.revalidatie.com.br" target="_blank">Visite nosso site</a></p>
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
    # 6. Mailchimp
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
        log(f"✅ Campanha criada com sucesso: {campaign['id']}")
        
        return {
            "status": "success", 
            "campaign_id": campaign['id'], 
            "url_capa": url_capa,
            "url_corpo": url_corpo,
            "custo_estimado": estimar_custo_revamais(),
            "instagram_assets": instagram_assets,
            "instagram_format": formato_instagram
        }
    except Exception as e:
        log(f"❌ Erro Mailchimp: {e}")
        return {"status": "error", "message": str(e)}

if __name__ == "__main__":
    # Teste rápido
    def print_log(msg): print(f"[LOG] {msg}")
    criar_campanha_revamais("A importância do sono na recuperação muscular", log_callback=print_log)
