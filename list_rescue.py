import firebase_admin
from firebase_admin import credentials, storage
import os

CREDENTIALS_PATH = "firebase_credentials.json"
BUCKET_NAME = "revalidatie-website.firebasestorage.app"

def list_rescue_urls():
    if not firebase_admin._apps:
        if os.path.exists(CREDENTIALS_PATH):
            cred = credentials.Certificate(CREDENTIALS_PATH)
            firebase_admin.initialize_app(cred, {'storageBucket': BUCKET_NAME})
        else:
            print("❌ Credenciais não encontradas.")
            return

    bucket = storage.bucket()
    blobs = list(bucket.list_blobs(prefix="resgate/"))
    
    print("\n🔗 Links dos arquivos de resgate:")
    for b in blobs:
        if b.name.endswith(".mp3"):
            # Garante que é público (caso não seja)
            # b.make_public() 
            # Ou usa o link de media do Firebase se make_public falhar por permissão
            print(f"- {b.name}: {b.public_url}")

if __name__ == "__main__":
    list_rescue_urls()
