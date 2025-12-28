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
    # Usando intervalo de datas explícito para maior robustez
    try:
        from datetime import timedelta
        agora = datetime.now()
        data_ini = (agora - timedelta(days=5*365)).strftime("%Y/%m/%d")
        data_fim = agora.strftime("%Y/%m/%d")
        date_term = f'("{data_ini}"[Date - Publication] : "{data_fim}"[Date - Publication])'
    except:
        date_term = "last 5 years[dp]"

    query_core = (
        f"({tema_ingles}) AND "
        f"(Systematic Review[pt] OR Randomized Controlled Trial[pt]) AND "
        f"jsubsetaim[text] AND "
        f"{date_term}"
    )
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query_core, retmax=3, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        
        if not ids:
            print("⚠️ Nenhuma referência 'Core' encontrada. Tentando busca mais ampla...")
            # Fallback 1: Sem jsubsetaim, mas ainda SR/RCT e 5 anos
            query_fallback_1 = (
                f"({tema_ingles}) AND "
                f"(Systematic Review[pt] OR Randomized Controlled Trial[pt]) AND "
                f"{date_term}"
            )
            handle = Entrez.esearch(db="pubmed", term=query_fallback_1, retmax=3, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            ids = record["IdList"]

        if not ids:
            print("⚠️ Ainda sem referências. Tentando busca por Review geral...")
            # Fallback 2: Review geral
            query_fallback_2 = (
                f"({tema_ingles}) AND (Review[pt]) AND {date_term}"
            )
            handle = Entrez.esearch(db="pubmed", term=query_fallback_2, retmax=3, sort="relevance")
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
            
            # Extração do Abstract
            abstract_text = ""
            if 'Abstract' in article['MedlineCitation']['Article']:
                abstract_data = article['MedlineCitation']['Article']['Abstract'].get('AbstractText', [])
                if isinstance(abstract_data, list):
                    abstract_text = " ".join([str(item) for item in abstract_data])
                else:
                    abstract_text = str(abstract_data)
            
            referencias.append({
                "texto": f"{primeiro_autor}. {titulo}. {journal}, {ano}.",
                "link": link,
                "resumo": abstract_text
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
    
    INSTRUÇÃO CRÍTICA:
    Abaixo estão 3 resumos científicos (abstracts) recentes sobre o tema.
    Você DEVE escrever o artigo baseando-se EXCLUSIVAMENTE nas evidências apresentadas nestes resumos.
    Não alucine informações. Se os resumos não cobrirem algo, generalize com cautela, mas priorize os dados abaixo.
    
    --- EVIDÊNCIA CIENTÍFICA ---
    {chr(10).join([ f'Artigo {i+1}: {r["texto"]}{chr(10)}Resumo: {r["resumo"]}{chr(10)}' for i, r in enumerate(referencias) ])}
    ----------------------------
    
    Escreva um artigo curto, informativo e acolhedor.
    
    Estrutura HTML (retorne APENAS o conteúdo dentro do body, sem tags html/body):
    1. <h1>Título Atraente</h1>
    2. <p>Introdução conectando com o dia a dia.</p>
    3. <h2>O que a ciência diz?</h2> (Explicação baseada nos abstracts acima. Cite "estudos recentes mostram..." sem ser técnico demais).
    4. <h2>Dicas Práticas</h2> (Conclusões práticas extraídas dos abstracts).
    5. <div class="cta"> (Um texto convidando para seguir no Instagram @revalidatie_londrina ou visitar o site).
    
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
    
    return html_content

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

import textwrap

def gerar_slide_instagram_composto(background_url, titulo, texto, slide_num):
    """
    Baixa o background gerado por IA, e desenha o texto por cima com Pillow.
    Garante legibilidade e branding.
    """
    try:
        # Configurações
        W, H = 1080, 1080
        font_path = download_font() # Garante fonte
        
        # Baixa imagem de fundo
        bg_img = download_logo(background_url).convert("RGBA").resize((W, H))
        
        # Cria overlay para escurecer e melhorar leitura
        overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)
        
        # Gradiente ou bloco sólido: Vamos usar um bloco semi-transparente na parte inferior/meio
        # Style: Clean Modern. 
        # Title at top (with dark bg behind it maybe?)
        # Body text centered or bottom.
        
        # Desenha um retângulo branco semi-transparente cobrindo a area do texto
        # Margens
        margin_x = 80
        margin_y = 150
        
        # --- TITLE ---
        # Caixa do Tiulo - Azul Escuro
        title_bg_color = (32, 87, 118, 240) # #205776 com alpha
        draw.rectangle([0, 80, W, 280], fill=title_bg_color)
        
        font_title = ImageFont.truetype(font_path, 70)
        
        # Wrap title
        title_lines = textwrap.wrap(titulo, width=25) # Ajustar width conforme font size
        current_h = 110
        for line in title_lines:
            # Centraliza texto horizontalmente
            bbox = draw.textbbox((0, 0), line, font=font_title)
            text_w = bbox[2] - bbox[0]
            draw.text(((W - text_w) / 2, current_h), line, font=font_title, fill="white")
            current_h += 80

        # --- BODY TEXT ---
        # Se tiver texto
        if texto:
            font_body = ImageFont.truetype(font_path, 45)
            # Caixa de texto corpo - Branco suave
            body_bg_color = (255, 255, 255, 220)
            
            # Area pro texto: Y=350 ate Y=900
            draw.rectangle([40, 350, W-40, 900], fill=body_bg_color, outline=None, width=0)
            
            body_lines = textwrap.wrap(texto, width=40)
            current_h = 400
            for line in body_lines:
                bbox = draw.textbbox((0, 0), line, font=font_body)
                text_w = bbox[2] - bbox[0]
                draw.text(((W - text_w) / 2, current_h), line, font=font_body, fill=(50, 50, 50))
                current_h += 55

        # --- FOOTER / BRANDING ---
        # Pequena barra embaixo
        draw.rectangle([0, H-60, W, H], fill=(32, 87, 118, 255))
        font_footer = ImageFont.truetype(font_path, 24)
        footer_text = "@revalidatie_londrina | Boletim Reva +"
        draw.text((margin_x, H-45), footer_text, fill="white", font=font_footer)
        
        # Paginação
        draw.text((W - margin_x - 50, H-45), f"{slide_num}/7", fill="white", font=font_footer)

        # Composite
        final_img = Image.alpha_composite(bg_img, overlay)
        
        # Salva e Upload
        temp_filename = f"slide_comp_{slide_num}_{datetime.now().strftime('%H%M%S')}.png"
        final_img.save(temp_filename, "PNG")
        
        firebase_path = f"revamais/slides/{temp_filename}"
        url = upload_file(temp_filename, firebase_path)
        try:
            if os.path.exists(temp_filename): os.remove(temp_filename)
        except:
            pass
        
        return url
        
    except Exception as e:
        print(f"⚠️ Erro overlay slide: {e}")
        import traceback
        traceback.print_exc()
        return background_url # Fallback para imagem original sem texto overlay

def gerar_conteudo_instagram(tema, formato, referencias_text):
    """
    Gera conteúdo para Instagram (Reel ou Carrossel).
    Melhoria v2: Gera texto com LLM, Imagem Clean com IA, e Texto via Overlay (Pillow).
    """
    print(f"📸 Gerando conteúdo para Instagram ({formato})...")
    assets = []
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    try:
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        
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
            prompt_slides = f"""
            Crie o planejamento de um Carrossel Educativo (7 slides) para Instagram.
            Tema: "{tema}"
            Tom: Educativo, sério porém acessível, foco em fisioterapia/saúde.
            
            Retorne APENAS um JSON:
            [
              {{"slide": 1, "titulo": "...", "texto_curto": "..."}},
              ...
            ]
            
            Regras para os Textos:
            - Títulos: Máximo 5 palavras. Impactantes.
            - Texto Curto (Corpo): Máximo 30 palavras. Resumido, direto ao ponto.
            - Slide 1: Título é a dor/problema, Texto é uma provocação.
            - Slide 2-3: Explicação científica simplificada.
            - Slide 4-6: Dicas, Exercícios ou Soluções.
            - Slide 7: Título "Gostou?", Texto "Siga @revalidatie_londrina para mais dicas.".
            """
            
            try:
                response_text = model.generate_content(prompt_slides).text
                response_text = response_text.replace("```json", "").replace("```", "").strip()
                slides_data = json.loads(response_text)
            except:
                # Fallback structure
                slides_data = [{"slide": i+1, "titulo": f"Slide {i+1}", "texto_curto": tema} for i in range(7)]

            print(f"   🖼️ Gerando 7 slides COMPOSTOS (IA + Texto Overlay)...")
            for slide in slides_data:
                slide_num = slide.get('slide')
                titulo = slide.get('titulo')
                texto = slide.get('texto_curto')
                
                # Prompt APENAS VISUAL (Sem texto)
                visual_prompt = (
                    f"A professional medical illustration background for a slide about: '{titulo}'. "
                    f"Context: {texto}. "
                    f"Style: Clean, minimalist, white/teal background, abstract shapes or anatomical diagram. "
                    f"IMPORTANT: NO TEXT. NO LABELS. EMPTY SPACE IN THE CENTER for text overlay."
                )

                # Gera background limpo
                temp_name = f"bg_slide_{slide_num}_{timestamp}"
                bg_url = gerar_imagem(visual_prompt, temp_name)
                
                # Composição (Texto overlay)
                final_url = gerar_slide_instagram_composto(bg_url, titulo, texto, slide_num)
                
                assets.append({
                    "type": "image", 
                    "url": final_url, 
                    "name": f"Slide {slide_num}: {titulo}",
                    "texto_base": texto
                })
                
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
    # Identifica fonte do tema para agendamento
    is_calendar_source = not tema or tema.strip() == ""
    
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
        log("🌍 Traduzindo tema para keywords científicas...")
        model = genai.GenerativeModel('gemini-2.5-flash-preview-09-2025')
        tema_ingles = model.generate_content(f"Translate this topic to English medical keywords for PubMed search. Return ONLY the keywords, no sentences: {tema}").text.strip()
    except:
        tema_ingles = tema

    check()
    # 2. Buscar Referências
    log(f"🔎 Buscando referências para: {tema_ingles}...")
    referencias = buscar_referencias_pubmed(tema_ingles)
    
    check()
    # 3. Gerar Imagens
    log("🎨 Gerando assets visuais (isso pode demorar)...")
    
    # Gerar Banner Composto (Header)
    try:
        log("🖼️ Criando banner de cabeçalho personalizado...")
        # This section is removed as per instruction.
    except Exception as e:
        log(f"⚠️ Falha ao criar banner composto: {e}. Usando logo padrão.")
    
    # Capa: Estática (solicitada pelo usuário)
    url_capa_estatica = "https://i.imgur.com/oGzxgtK.jpeg"
    
    # 1. Imagem Ilustrativa (Lifestyle/Visual)
    # Ex: Pessoa na esteira, feliz, etc. Sem texto.
    prompt_ilustrativa = f"A high quality, photorealistic or cinematic style photo-illustration about '{tema_ingles}'. Showing people, lifestyle, or the subject in a natural, positive way. NO TEXT. Suitable for a newsletter cover."
    url_ilustrativa = gerar_imagem(prompt_ilustrativa, "ilustrativa")

    # 2. Imagem Corpo (Infográfico/Educativo)
    # Se tiver texto, DEVE SER EM PORTUGUÊS
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
    
    today_formatted = datetime.now().strftime('%d/%m/%Y')
    
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
                Aqui está sua atualização semanal de saúde do Reva + para {today_formatted}.</p>
            </div>
            
            <div style="padding: 20px;">
                <!-- Imagem Ilustrativa (Nova) -->
                <img src="{url_ilustrativa}" class="body-img" alt="Ilustração do Título">
                
                {html_texto}

                <!-- Imagem Educativa (Existente) -->
                <img src="{url_corpo}" class="body-img" alt="Infográfico">
                
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
        log(f"✅ Campanha criada com sucesso (Draft): {campaign['id']}")
        
        # 7. Agendamento Automático (Terça-feira 12:00 BRT)
        try:
            target_weekday = 1 if is_calendar_source else 6 # 1=Terça, 6=Domingo
            day_name = "Terça-feira" if is_calendar_source else "Domingo"
            
            log(f"📅 Tentando agendar envio para próximo(a) {day_name} (12:00 BRT)...")
            
            # Cálculo do próximo dia alvo
            now_utc = datetime.utcnow()
            days_until = (target_weekday - now_utc.weekday()) % 7
            
            # Se for hoje e já passou das 15:00 UTC (12:00 BRT), agendar para a próxima semana
            if days_until == 0 and now_utc.hour >= 15:
                days_until = 7
            elif days_until == 0 and now_utc.hour < 15:
                 # Se for hoje e ainda não passou, agenda pra hoje mesmo (0 dias)
                 pass
            
            next_date = now_utc + timedelta(days=days_until)
            # Define 15:00 UTC (12:00 BRT)
            schedule_time = next_date.replace(hour=15, minute=0, second=0, microsecond=0)
            
            # Formato ISO 8601 UTC
            schedule_str = schedule_time.strftime('%Y-%m-%dT%H:%M:%S+00:00')
            
            mc.campaigns.schedule(campaign["id"], {"schedule_time": schedule_str})
            log(f"🕒 Campanha agendada com sucesso para: {schedule_str} (UTC)")
            
        except Exception as e:
            # Muitos planos gratuitos não permitem agendamento via API
            log(f"⚠️ Não foi possível agendar automatiamente: {e}")
            log("ℹ️ A campanha foi salva como RASCUNHO. Por favor, agende ou envie manualmente pelo painel do Mailchimp.")

        return {
            "status": "success", 
            "campaign_id": campaign['id'], 
            "url_capa": url_capa_estatica,
            "url_ilustrativa": url_ilustrativa,
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
