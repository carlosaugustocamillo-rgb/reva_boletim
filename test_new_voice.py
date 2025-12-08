
import os
import time
from dotenv import load_dotenv
from elevenlabs.client import ElevenLabs
from elevenlabs import VoiceSettings
from pydub import AudioSegment

# Carregar ambiente
load_dotenv()

# Configurações
API_KEY = os.environ.get("ELEVENLABS_API_KEY")
HOST_ID = "k0eIjUFv1GH1DnJFJK46"
COHOST_ID = "x3mAOLD9WzlmrFCwA1S3"

print(f"🔹 Teste rápido de áudio com:")
print(f"   HOST: {HOST_ID}")
print(f"   COHOST: {COHOST_ID}")

client = ElevenLabs(api_key=API_KEY)

# Lista de falas para o teste
dialogo = [
    {"speaker": "HOST", "text": "Olá! Estamos testando esta nova voz para o RevaCast. Espero que a qualidade esteja boa."},
    {"speaker": "COHOST", "text": "Oi! Estou te ouvindo bem. O tom parece profissional e claro para o podcast."},
    {"speaker": "HOST", "text": "Excelente notícia. Vamos encerrar este teste agora. Obrigado!"}
]

audios_gerados = []
output_dir = "data/teste_rapido"
os.makedirs(output_dir, exist_ok=True)

print("\n🎙️ Gerando áudios...")

for i, fala in enumerate(dialogo):
    speaker = fala["speaker"]
    text = fala["text"]
    voice_id = HOST_ID if speaker == "HOST" else COHOST_ID
    
    print(f"   Gerando fala {i+1} ({speaker})...")
    
    try:
        audio_generator = client.text_to_speech.convert(
            voice_id=voice_id,
            text=text,
            model_id="eleven_multilingual_v2",
             voice_settings=VoiceSettings(
                stability=0.75,
                similarity_boost=0.85,
                style=0.0,
                use_speaker_boost=True
            )
        )
        
        filename = f"{output_dir}/fala_{i+1}_{speaker}.mp3"
        with open(filename, "wb") as f:
            for chunk in audio_generator:
                f.write(chunk)
        
        audios_gerados.append(filename)
        print(f"   ✅ Salvo: {filename}")
        
    except Exception as e:
        print(f"   ❌ Erro ao gerar: {e}")

# Tentar juntar com introdução se existir
print("\n🎧 Montando episódio de teste...")
final_audio = AudioSegment.empty()

# Tenta carregar intro
intro_path = "data/intro_guto.mp3"
if os.path.exists(intro_path):
    print("   + Adicionando introdução")
    final_audio += AudioSegment.from_file(intro_path, format="mp3")
    final_audio += AudioSegment.silent(duration=1000)
else:
    print("   (Introdução não encontrada, pulando)")

# Adiciona falas
for audio_path in audios_gerados:
    if os.path.exists(audio_path):
        final_audio += AudioSegment.from_file(audio_path, format="mp3")
        final_audio += AudioSegment.silent(duration=500)

# Salva final
final_path = "data/episodio_teste_novo_host.mp3"
final_audio.export(final_path, format="mp3")

print(f"\n✅ TESTE CONCLUÍDO!")
print(f"Arquivo final gerado em: {final_path}")
print("Verifique este arquivo para validar a nova voz.")
