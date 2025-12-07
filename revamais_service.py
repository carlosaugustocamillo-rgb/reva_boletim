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
# Tenta pegar ID específico para Reva+, senão usa o padrão (Weekly)
MC_LIST_ID = os.environ.get("MC_LIST_ID_REVAMAIS") or os.environ.get("MC_LIST_ID")
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

def criar_campanha_revamais(tema):
    print(f"🚀 Iniciando Reva +: {tema}")
    
    # 1. Traduzir tema para busca
    try:
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        tema_ingles = model.generate_content(f"Translate to English for medical search: {tema}").text.strip()
    except:
        tema_ingles = tema

    # 2. Buscar Referências
    referencias = buscar_referencias_pubmed(tema_ingles)
    
    # 3. Gerar Imagens
    # Capa: SEM TEXTO, apenas visual
    prompt_capa = f"A welcoming, modern health illustration about '{tema_ingles}'. Soft lighting, professional photography style or high quality 3D render. Minimalist. NO TEXT. NO WORDS."
    url_capa = gerar_imagem(prompt_capa, "capa")
    
    # Imagem Corpo: Se tiver texto, DEVE SER EM PORTUGUÊS
    prompt_corpo = f"An educational infographic or diagram about '{tema_ingles}'. Clean lines, easy to understand, white background. IMPORTANT: Any text or labels MUST BE IN PORTUGUESE (PT-BR). If you cannot generate correct Portuguese text, do not include any text."
    url_corpo = gerar_imagem(prompt_corpo, "corpo")
    
    # 4. Gerar Texto
    html_texto = gerar_conteudo_revamais(tema, referencias)
    
    # 5. Montar HTML Final
    # Logo da Revalidatie (usado no boletim_service.py)
    # Logo da Revalidatie
    # Se o usuário definir uma no .env, usa. Senão, usa a padrão.
    DEFAULT_LOGO = "https://i.imgur.com/6FIUeHX.png"
    logo_url = os.environ.get("REVA_LOGO_URL", DEFAULT_LOGO)
    
    html_email = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; line-height: 1.6; color: #333; max-width: 600px; margin: 0 auto; background-color: #f9f9f9; }}
            .container {{ background-color: #ffffff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-top: 20px; }}
            .logo-header {{ text-align: center; margin-bottom: 20px; }}
            .logo-header img {{ max-width: 200px; height: auto; }}
            .header-img {{ width: 100%; border-radius: 8px 8px 0 0; display: block; }}
            .body-img {{ width: 100%; margin: 20px 0; border-radius: 8px; }}
            h1 {{ color: #205776; }}
            h2 {{ color: #407ca6; }}
            a {{ color: #205776; text-decoration: none; font-weight: bold; }}
            .footer {{ font-size: 12px; color: #777; text-align: center; margin-top: 30px; }}
            .cta-box {{ background-color: #eef6fa; padding: 15px; border-radius: 8px; text-align: center; margin: 20px 0; border: 1px solid #cce0eb; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="logo-header">
                <a href="https://www.revalidatie.com.br" target="_blank">
                    <img src="{logo_url}" alt="Revalidatie">
                </a>
            </div>
            
            <img src="{url_capa}" class="header-img" alt="Capa">
            
            <div style="padding: 20px;">
                {html_texto}
                
                <img src="{url_corpo}" class="body-img" alt="Ilustração">
                
                <div class="cta-box">
                    <p>Quer saber mais sobre como cuidar da sua saúde?</p>
                    <p><a href="https://www.instagram.com/revalidatie/" target="_blank">Siga-nos no Instagram</a> ou <a href="https://www.revalidatie.com.br" target="_blank">Visite nosso site</a></p>
                </div>
            </div>
        </div>
        
        <div class="footer">
            <p>Reva + | Revalidatie<br>
            <a href="*|UNSUB|*">Descadastrar</a></p>
        </div>
    </body>
    </html>
    """
    
    # 6. Mailchimp
    try:
        print("📧 Criando rascunho no Mailchimp...")
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
        print(f"✅ Campanha criada: {campaign['id']}")
        return {
            "status": "success", 
            "campaign_id": campaign['id'], 
            "url_capa": url_capa,
            "url_corpo": url_corpo,
            "custo_estimado": estimar_custo_revamais()
        }
        }
    except Exception as e:
        print(f"❌ Erro Mailchimp: {e}")
        return {"status": "error", "message": str(e)}

if __name__ == "__main__":
    # Teste rápido
    criar_campanha_revamais("A importância do sono na recuperação muscular")
