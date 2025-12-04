import os
import firebase_admin
from firebase_admin import credentials, storage
from datetime import datetime
import xml.etree.ElementTree as ET
from email.utils import formatdate
from dotenv import load_dotenv

# Carrega variáveis de ambiente
load_dotenv()

# Configuração Firebase
if not firebase_admin._apps:
    cred_path = os.environ.get("FIREBASE_CREDENTIALS_JSON", "firebase_credentials.json")
    if os.path.exists(cred_path):
        cred = credentials.Certificate(cred_path)
        firebase_admin.initialize_app(cred, {
            'storageBucket': os.environ.get("FIREBASE_BUCKET_NAME")
        })
    else:
        print("❌ Erro: Credenciais do Firebase não encontradas.")
        exit(1)

BUCKET_NAME = os.environ.get("FIREBASE_BUCKET_NAME")
bucket = storage.bucket()

def force_update_feed():
    print("🔄 Iniciando atualização forçada do Feed RSS...")
    
    # 1. Lista todos os episódios na pasta 'episodios/'
    blobs = bucket.list_blobs(prefix="episodios/")
    episodios = []
    
    for blob in blobs:
        if blob.name.endswith(".mp3"):
            # Tenta extrair data do nome (episodio_boletim_YYYY-MM-DD.mp3)
            try:
                data_str = blob.name.split("_")[-1].replace(".mp3", "")
                data_obj = datetime.strptime(data_str, "%Y-%m-%d")
            except:
                # Se falhar, usa a data de criação do arquivo
                data_obj = blob.time_created
            
            # Torna o blob público para pegar o link
            try:
                blob.make_public()
            except:
                pass
                
            episodios.append({
                'titulo': f"Boletim Científico - {data_obj.strftime('%d/%m/%Y')}",
                'data': data_obj,
                'url': blob.public_url,
                'tamanho': blob.size,
                'descricao': f"Destaques da semana na ciência da reabilitação. (Data: {data_obj.strftime('%d/%m/%Y')})"
            })
    
    # Ordena por data (mais recente primeiro)
    episodios.sort(key=lambda x: x['data'], reverse=True)
    
    print(f"✅ Encontrados {len(episodios)} episódios.")
    
    # 2. Gera o XML do Feed
    rss = ET.Element("rss", version="2.0", xmlns_itunes="http://www.itunes.com/dtds/podcast-1.0.dtd")
    channel = ET.SubElement(rss, "channel")
    
    ET.SubElement(channel, "title").text = "RevaCast Weekly"
    ET.SubElement(channel, "description").text = "O seu boletim semanal de ciência da reabilitação."
    ET.SubElement(channel, "language").text = "pt-br"
    ET.SubElement(channel, "link").text = "https://www.revalidatie.com.br"
    
    # Adiciona episódios
    for ep in episodios:
        item = ET.SubElement(channel, "item")
        ET.SubElement(item, "title").text = ep['titulo']
        ET.SubElement(item, "description").text = ep['descricao']
        ET.SubElement(item, "pubDate").text = formatdate(ep['data'].timestamp())
        
        enclosure = ET.SubElement(item, "enclosure")
        enclosure.set("url", ep['url'])
        enclosure.set("type", "audio/mpeg")
        enclosure.set("length", str(ep['tamanho']))
        
        ET.SubElement(item, "guid").text = ep['url']
    
    # 3. Salva e Upload
    feed_path = "feed_temp.xml"
    tree = ET.ElementTree(rss)
    tree.write(feed_path, encoding="UTF-8", xml_declaration=True)
    
    blob_feed = bucket.blob("feed.xml")
    blob_feed.upload_from_filename(feed_path)
    blob_feed.make_public()
    
    print(f"🚀 Feed atualizado com sucesso!")
    print(f"🔗 URL do Feed: {blob_feed.public_url}")
    
    # Limpa temp
    if os.path.exists(feed_path):
        os.remove(feed_path)

if __name__ == "__main__":
    force_update_feed()
