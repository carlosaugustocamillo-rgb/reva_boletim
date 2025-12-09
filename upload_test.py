import os
import firebase_admin
from firebase_admin import credentials, storage

CREDENTIALS_PATH = "firebase_credentials.json"
BUCKET_NAME = "revalidatie-website.firebasestorage.app"

def upload_test():
    if not firebase_admin._apps:
        cred = credentials.Certificate(CREDENTIALS_PATH)
        firebase_admin.initialize_app(cred, {'storageBucket': BUCKET_NAME})

    bucket = storage.bucket()
    blob = bucket.blob("testes/teste_host_isolado.mp3")
    blob.upload_from_filename("teste_host_isolado.mp3")
    blob.make_public()
    print(f"URL: {blob.public_url}")

if __name__ == "__main__":
    upload_test()
