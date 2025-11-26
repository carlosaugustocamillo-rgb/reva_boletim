import os
import firebase_admin
from firebase_admin import credentials, storage
from podgen import Podcast, Episode, Media, Person, Category
from datetime import datetime, timedelta
import pytz

# Configuração
# Certifique-se de ter o arquivo firebase_credentials.json na raiz
CREDENTIALS_PATH = "firebase_credentials.json"
BUCKET_NAME = os.environ.get("FIREBASE_BUCKET_NAME", "seu-bucket.appspot.com")

# Inicializa Firebase
if not firebase_admin._apps:
    if os.path.exists(CREDENTIALS_PATH):
        cred = credentials.Certificate(CREDENTIALS_PATH)
        firebase_admin.initialize_app(cred, {
            'storageBucket': BUCKET_NAME
        })
    else:
        print(f"⚠️ Aviso: Arquivo {CREDENTIALS_PATH} não encontrado. Uploads falharão.")

def upload_file(local_path, destination_blob_name):
    """Faz upload de um arquivo para o Firebase Storage e retorna a URL pública."""
    if not firebase_admin._apps:
        print("❌ Firebase não inicializado.")
        return None

    bucket = storage.bucket()
    blob = bucket.blob(destination_blob_name)
    
    print(f"⬆️ Iniciando upload de {local_path} para {destination_blob_name}...")
    blob.upload_from_filename(local_path)
    blob.make_public()
    
    print(f"✅ Upload concluído: {blob.public_url}")
    return blob.public_url

def update_podcast_feed(episodio_audio_url, episodio_titulo, episodio_descricao, data_pub, duracao_segundos, tamanho_bytes):
    """
    Gera ou atualiza o feed RSS do podcast e faz upload para o Firebase.
    """
    rss_filename = "podcast_feed.xml"
    local_rss_path = os.path.join("data", rss_filename)
    
    # Configurações do Podcast (Idealmente viriam de variáveis de ambiente)
    p = Podcast(
        name="RevaCast Weekly",
        description="Resumo semanal das evidências científicas em reabilitação pulmonar e exercício físico.",
        website="https://www.revalidatie.com.br",
        explicit=False,
        image="https://i.imgur.com/1b57Ych.png", # Use uma URL válida da sua capa
        authors=[Person("Revalidatie", "contato@revalidatie.com.br")],
        language="pt-br",
        category=Category("Health & Fitness", "Medicine"),
        owner=Person("Revalidatie", "contato@revalidatie.com.br"),
    )
    
    # Tenta carregar feed existente (lógica simplificada: recria o feed com o novo episódio + placeholder para antigos se necessário)
    # Para um sistema real, idealmente baixaríamos o XML atual, faríamos parse e adicionaríamos o novo.
    # Como estamos começando, vamos assumir que o script mantém um histórico local ou recria.
    # A lib podgen não faz "append" em XML existente facilmente sem parsear. 
    # VAMOS SIMPLIFICAR: Vamos criar um feed com o episódio atual. 
    # **IMPORTANTE**: Para produção, você precisa de um banco de dados dos episódios passados 
    # ou ler o XML antigo. Vou implementar a leitura do XML antigo se existir localmente.
    
    if os.path.exists(local_rss_path):
        try:
            p = Podcast.from_file(local_rss_path)
        except Exception as e:
            print(f"⚠️ Erro ao ler feed existente, criando novo: {e}")

    # Adiciona o novo episódio
    ep = Episode()
    ep.title = episodio_titulo
    ep.summary = episodio_descricao
    ep.publication_date = data_pub # datetime com timezone
    ep.media = Media(episodio_audio_url, size=tamanho_bytes, type="audio/mpeg", duration=timedelta(seconds=duracao_segundos))
    
    p.add_episode(ep)
    
    # Salva localmente
    p.rss_file(local_rss_path)
    
    # Upload do RSS atualizado
    rss_url = upload_file(local_rss_path, rss_filename)
    return rss_url
