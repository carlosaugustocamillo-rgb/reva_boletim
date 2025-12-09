import firebase_admin
from firebase_admin import credentials, storage
import os

CREDENTIALS_PATH = "firebase_credentials.json"
BUCKET_NAME = "revalidatie-website.firebasestorage.app"
DOWNLOAD_DIR = "downloads_resgate"

def download_rescue():
    if not firebase_admin._apps:
        if os.path.exists(CREDENTIALS_PATH):
            cred = credentials.Certificate(CREDENTIALS_PATH)
            firebase_admin.initialize_app(cred, {'storageBucket': BUCKET_NAME})
        else:
            print("❌ Credenciais não encontradas.")
            return

    if not os.path.exists(DOWNLOAD_DIR):
        os.makedirs(DOWNLOAD_DIR)

    bucket = storage.bucket()
    print("🔎 Listando arquivos...")
    blobs = list(bucket.list_blobs(prefix="resgate/"))
    
    count = 0
    for b in blobs:
        if b.name.endswith(".mp3"):
            # Cria subpastas locais se necessário (ex: downloads_resgate/2025-11-29/)
            local_path = os.path.join(DOWNLOAD_DIR, b.name.replace("resgate/", ""))
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            
            print(f"⬇️ Baixando {b.name}...")
            b.download_to_filename(local_path)
            count += 1

    print(f"\n✅ {count} arquivos baixados para a pasta '{DOWNLOAD_DIR}'")

if __name__ == "__main__":
    download_rescue()
