import os
import requests
import json
from datetime import datetime, timedelta
import pytz
from dotenv import load_dotenv

load_dotenv()

from openai import OpenAI
from Bio import Entrez
import mailchimp_marketing as MailchimpMarketing
from mailchimp_marketing.api_client import ApiClientError

# Importa ferramentas já existentes
from firebase_service import upload_file

import google.generativeai as genai

# Configurações
Entrez.email = os.environ.get("ENTREZ_EMAIL")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

# Configuração Gemini
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)
else:
    print("⚠️ GEMINI_API_KEY não encontrada. O script pode falhar se tentar usar o Gemini.")

# Configuração Mailchimp
mc = MailchimpMarketing.Client()
mc.set_config({
    "api_key": os.environ.get("MC_API_KEY"),
    "server": os.environ.get("MC_SERVER")
})
MC_LIST_ID = os.environ.get("MC_LIST_ID") # Se for uma lista nova, altere no .env ou passe como argumento
MC_FROM_NAME = os.environ.get("MC_FROM_NAME", "Revalidatie")
MC_REPLY_TO = os.environ.get("MC_REPLY_TO", "contato@revalidatie.com.br")

def buscar_evidencia_cientifica(tema_pesquisa):
    """
    Realiza uma 'Deep Research' no PubMed buscando Guidelines e Revisões Sistemáticas
    sobre o tema para embasar a resposta.
    """
    print(f"🔎 Investigando no PubMed: {tema_pesquisa}...")
    
    # Query focada em segurança e diretrizes
    query = f"({tema_pesquisa}) AND (Practice Guideline[pt] OR Systematic Review[pt] OR Review[pt]) AND (safety[tiab] OR recommendations[tiab] OR contraindications[tiab])"
    
    try:
        # 1. Busca IDs (focando em relevância, não data)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        
        if not ids:
            return "Nenhuma diretriz específica encontrada no PubMed. O texto será baseado em conhecimento clínico geral."

        # 2. Busca Detalhes (Resumos)
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        contexto_cientifico = ""
        for article in records['PubmedArticle']:
            titulo = article['MedlineCitation']['Article']['ArticleTitle']
            try:
                abstract_list = article['MedlineCitation']['Article']['Abstract']['AbstractText']
                resumo = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
            except KeyError:
                resumo = "Resumo não disponível."
            
            contexto_cientifico += f"- Estudo: {titulo}\n- Resumo: {resumo}\n\n"
            
        return contexto_cientifico

    except Exception as e:
        print(f"Erro na busca do PubMed: {e}")
        return "Erro ao buscar evidências. O texto será baseado no conhecimento do modelo."

def gerar_imagem_capa(tema):
    """
    Gera uma imagem via Gemini (gemini-3-pro-image-preview) ou DALL-E 3 (fallback) e faz upload para o Firebase.
    """
    print("🎨 Gerando imagem de capa...")
    
    prompt_imagem = f"""
    A clean, modern, abstract medical illustration representing: {tema}.
    Style: Minimalist, soft colors (blue, white, light orange), professional vector art style.
    IMPORTANT: NO TEXT, NO LETTERS, NO NUMBERS, NO WORDS inside the image.
    Just pure visual symbolism and shapes.
    Suitable for a health newsletter header.
    """
    
    temp_filename = f"temp_cover_{datetime.now().strftime('%H%M%S')}.png"
    image_generated = False
    
    # Tenta usar Gemini primeiro
    try:
        print("   👉 Tentando Gemini (gemini-3-pro-image-preview)...")
        model = genai.GenerativeModel('gemini-3-pro-image-preview')
        response = model.generate_content("Generate an image of: " + prompt_imagem)
        
        for part in response.parts:
            if hasattr(part, 'inline_data') and part.inline_data:
                img_data = part.inline_data.data
                with open(temp_filename, "wb") as f:
                    f.write(img_data)
                image_generated = True
                print("   ✅ Imagem gerada com Gemini!")
                break
                
        if not image_generated:
            print("   ⚠️ Gemini não retornou imagem. Tentando DALL-E 3...")
            
    except Exception as e:
        print(f"   ❌ Erro Gemini Imagem: {e}")
        print("   👉 Tentando DALL-E 3 (fallback)...")

    # Fallback para DALL-E 3 se Gemini falhar
    if not image_generated:
        try:
            print("   🎨 Usando DALL-E 3...")
            response = client.images.generate(
                model="dall-e-3",
                prompt=prompt_imagem,
                size="1024x1024",
                quality="standard",
                n=1,
            )
            image_url = response.data[0].url
            img_data = requests.get(image_url).content
            with open(temp_filename, "wb") as f:
                f.write(img_data)
            image_generated = True
            
        except Exception as e:
            print(f"❌ Erro ao gerar imagem (DALL-E 3): {e}")
            return "https://via.placeholder.com/600x300?text=RevaCast+Weekly" # Fallback final

    try:
        # Upload para Firebase
        firebase_path = f"campanhas/cover_{datetime.now().strftime('%Y-%m-%d')}.png"
        public_url = upload_file(temp_filename, firebase_path)
        
        # Limpa temp
        if os.path.exists(temp_filename):
            os.remove(temp_filename)
        
        return public_url
        
    except Exception as e:
        print(f"❌ Erro no upload: {e}")
        return "https://via.placeholder.com/600x300?text=RevaCast+Weekly"

def gerar_email_educativo(tema, contexto_cientifico):
    """
    Usa o GPT-4 para escrever o e-mail baseado na evidência encontrada.
    """
    print("✍️ Escrevendo e-mail educativo...")
    
    prompt = f"""
    Você é um especialista em fisiologia do exercício e comunicação em saúde da 'Revalidatie'.
    
    TEMA DO E-MAIL: "{tema}"
    
    EVIDÊNCIA CIENTÍFICA ENCONTRADA (Use isso como base de verdade):
    {contexto_cientifico}
    
    OBJETIVO:
    Escrever um e-mail para pacientes e leigos explicando esse tema de forma simples, segura e prática.
    
    ESTRUTURA DO E-MAIL (HTML):
    1. Título (H1): Uma pergunta ou afirmação chamativa.
    2. Introdução: Conecte com a dor/dúvida do leitor.
    3. A Ciência Simplificada: Explique o que acontece no corpo (baseado nos abstracts acima).
    4. O Veredito Prático: Pode ou não pode? Quais os cuidados? (Ex: verificar cetonas).
    5. Call to Action (CTA): Convide para conhecer a Revalidatie ou ouvir o podcast.
    
    TOM DE VOZ:
    Empático, autoridade, claro, sem "mediquês" desnecessário. Use negrito para destacar o importante.
    
    FORMATO:
    Retorne APENAS o código HTML do corpo do e-mail (sem tags <html> ou <body>, apenas o conteúdo interno: <h1>, <p>, <ul>, etc).
    """
    
    try:
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        response = model.generate_content(prompt)
        return response.text
    except Exception as e:
        print(f"❌ Erro Gemini: {e}. Tentando fallback para GPT-4...")
        # Fallback para GPT-4 se Gemini falhar
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "Você é um redator de saúde sênior."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7
        )
        return response.choices[0].message.content

def criar_campanha_tematica(tema_usuario):
    """
    Orquestra todo o processo.
    """
    print(f"🚀 Iniciando campanha temática: {tema_usuario}")
    
    # 1. Deep Research
    # Traduzimos o tema para inglês para buscar no PubMed, mas mantemos o contexto em PT
    try:
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        tema_ingles_gpt = model.generate_content(f"Translate this medical topic to English keywords for PubMed search: {tema_usuario}").text.strip()
    except Exception:
        tema_ingles_gpt = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "user", "content": f"Translate this medical topic to English keywords for PubMed search: {tema_usuario}"}]
        ).choices[0].message.content
    
    evidencia = buscar_evidencia_cientifica(tema_ingles_gpt)
    
    # 2. Gerar Imagem
    url_imagem = gerar_imagem_capa(tema_ingles_gpt)
    
    # 3. Gerar Texto
    corpo_email_html = gerar_email_educativo(tema_usuario, evidencia)
    
    # 4. Montar HTML Final
    html_completo = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{ font-family: Helvetica, Arial, sans-serif; line-height: 1.6; color: #333; max-width: 600px; margin: 0 auto; }}
            .header-img {{ width: 100%; max-width: 600px; border-radius: 8px; margin-bottom: 20px; }}
            .content {{ padding: 20px; }}
            .cta-button {{ display: inline-block; background-color: #205776; color: white; padding: 12px 24px; text-decoration: none; border-radius: 5px; font-weight: bold; margin-top: 20px; }}
            .footer {{ margin-top: 40px; font-size: 12px; color: #888; text-align: center; border-top: 1px solid #eee; padding-top: 20px; }}
        </style>
    </head>
    <body>
        <img src="{url_imagem}" alt="Capa do E-mail" class="header-img">
        
        <div class="content">
            {corpo_email_html}
            
            <div style="text-align: center;">
                <a href="https://www.revalidatie.com.br" class="cta-button">Agendar Avaliação na Revalidatie</a>
            </div>
        </div>
        
        <div class="footer">
            <p>Revalidatie - Reabilitação Pulmonar e Exercício<br>
            Baseado em evidências científicas.</p>
            <p><a href="*|UNSUB|*">Descadastrar</a></p>
        </div>
    </body>
    </html>
    """
    
    # 5. Mailchimp
    try:
        print("📧 Criando campanha no Mailchimp...")
        campaign = mc.campaigns.create({
            "type": "regular",
            "recipients": {"list_id": MC_LIST_ID},
            "settings": {
                "subject_line": f"RevaCast: {tema_usuario}",
                "title": f"Campanha Temática: {tema_usuario[:20]}...",
                "from_name": MC_FROM_NAME,
                "reply_to": MC_REPLY_TO
            }
        })
        
        mc.campaigns.set_content(campaign["id"], {"html": html_completo})
        
        # Agendar para amanhã às 8:00 (exemplo)
        # Para teste imediato, vamos apenas criar como rascunho
        print(f"✅ Campanha criada com sucesso! ID: {campaign['id']}")
        print("⚠️ A campanha foi salva como RASCUNHO. Acesse o Mailchimp para revisar e enviar.")
        
        return {"status": "success", "campaign_id": campaign['id'], "image": url_imagem}
        
    except ApiClientError as error:
        print(f"❌ Erro Mailchimp: {error.text}")
        return {"status": "error", "detail": error.text}

if __name__ == "__main__":
    # Teste direto
    tema = "Benefícios do HIIT (Treino Intervalado) no tratamento da Hipertensão Arterial"
    criar_campanha_tematica(tema)
