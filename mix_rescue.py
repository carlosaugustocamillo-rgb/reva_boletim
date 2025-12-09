import os
import firebase_admin
from firebase_admin import credentials, storage
from pydub import AudioSegment
import io

# Configurações
CREDENTIALS_PATH = "firebase_credentials.json"
# Nome do bucket padrão (fallback do código original)
BUCKET_NAME = "revalidatie-website.firebasestorage.app" 
INTRO_PATH = "data/intro_guto.mp3"

def mix_rescue():
    print("🚀 Iniciando resgate de episódio...")
    
    # Inicializa Firebase
    if not firebase_admin._apps:
        if os.path.exists(CREDENTIALS_PATH):
            cred = credentials.Certificate(CREDENTIALS_PATH)
            firebase_admin.initialize_app(cred, {'storageBucket': BUCKET_NAME})
        else:
            print(f"❌ Erro: {CREDENTIALS_PATH} não encontrado.")
            return

    bucket = storage.bucket()

    # Lista arquivos na pasta resgate/
    print("🔎 Procurando arquivos em resgate/...")
    blobs = list(bucket.list_blobs(prefix="resgate/"))
    
    # Agrupa por data (resgate/YYYY-MM-DD/...)
    folders = set()
    for b in blobs:
        parts = b.name.split('/')
        if len(parts) >= 3:
            folders.add(parts[1])
    
    if not folders:
        print("❌ Nenhuma pasta encontrada em resgate/.")
        return

    # Pega a data mais recente
    latest_date = sorted(list(folders))[-1]
    print(f"📅 Data encontrada: {latest_date}")

    # Filtra blobs dessa data e que sejam mp3
    parts_blobs = [b for b in blobs if f"resgate/{latest_date}/" in b.name and b.name.endswith(".mp3")]
    
    # Ordena por número do estudo
    def get_study_num(blob):
        try:
            # Ex: resgate/2025-11-29/estudo1_completo.mp3
            fname = blob.name.split('/')[-1]
            num = fname.replace('estudo', '').replace('_completo.mp3', '')
            return int(num)
        except:
            return 999

    parts_blobs.sort(key=get_study_num)

    if not parts_blobs:
        print("❌ Nenhum arquivo de áudio encontrado para essa data.")
        return

    print(f"📥 Baixando {len(parts_blobs)} partes...")
    
    audio_segments = []
    
    # Adiciona Intro
    if os.path.exists(INTRO_PATH):
        print("🎵 Adicionando Intro...")
        try:
            intro = AudioSegment.from_file(INTRO_PATH)
            audio_segments.append(intro)
            audio_segments.append(AudioSegment.silent(duration=1000))
        except Exception as e:
            print(f"⚠️ Erro ao carregar intro: {e}")
    else:
        print(f"⚠️ Intro não encontrada em {INTRO_PATH}. O episódio ficará sem intro.")

    # Baixa e adiciona partes
    for b in parts_blobs:
        print(f"   - Processando {b.name}...")
        blob_content = b.download_as_bytes()
        seg = AudioSegment.from_file(io.BytesIO(blob_content), format="mp3")
        audio_segments.append(seg)
        audio_segments.append(AudioSegment.silent(duration=2000)) # 2s silêncio entre estudos

    # Mixagem
    print("🎛️ Mixando episódio final...")
    final_audio = AudioSegment.empty()
    for seg in audio_segments:
        final_audio += seg

    # Exporta localmente
    output_filename = f"episodio_boletim_{latest_date}_RESGATE.mp3"
    final_audio.export(output_filename, format="mp3")
    print(f"💾 Salvo localmente: {output_filename}")

    # Upload
    print("☁️ Fazendo upload para episodios/...")
    blob_dest = bucket.blob(f"episodios/{output_filename}")
    blob_dest.upload_from_filename(output_filename)
    blob_dest.make_public()
    
    print("-" * 30)
    print(f"✅ SUCESSO! Episódio resgatado disponível em:")
    print(blob_dest.public_url)
    print("-" * 30)

if __name__ == "__main__":
    mix_rescue()
