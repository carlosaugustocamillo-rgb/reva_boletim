"""
Script de teste para verificar a velocidade do áudio com pydub
"""
import os
import time
from dotenv import load_dotenv
from elevenlabs.client import ElevenLabs
from elevenlabs import VoiceSettings
from pydub import AudioSegment

load_dotenv()

ELEVENLABS_API_KEY = os.environ.get("ELEVENLABS_API_KEY")
ELEVEN_VOICE_ID_HOST = os.environ.get("ELEVEN_VOICE_ID_HOST")
ELEVEN_VOICE_ID_COHOST = os.environ.get("ELEVEN_VOICE_ID_COHOST")

client = ElevenLabs(api_key=ELEVENLABS_API_KEY)

# Configurações de teste
SPEED_FACTOR = 1.30  # Aumentar velocidade em 30%

def gerar_e_acelerar(voice_id, nome, texto, filename):
    print(f"\n🎙️ Gerando áudio para {nome}...")
    start = time.time()
    
    # 1. Gerar áudio normal
    audio_generator = client.text_to_speech.convert(
        voice_id=voice_id,
        text=texto,
        model_id="eleven_multilingual_v2",
        voice_settings=VoiceSettings(
            stability=0.40,       # Menor estabilidade = mais dinâmico/rápido
            similarity_boost=0.80, 
            style=0.60,           # Mais estilo = mais expressivo
            use_speaker_boost=True
        )
    )
    
    temp_file = f"temp_{filename}"
    with open(temp_file, "wb") as f:
        for chunk in audio_generator:
            f.write(chunk)
            
    # 2. Acelerar com pydub
    print(f"   ⚡ Acelerando áudio em {int((SPEED_FACTOR-1)*100)}%...")
    audio = AudioSegment.from_file(temp_file)
    
    # Acelera mudando o frame rate (sem alterar o pitch de forma complexa, mas pydub simples altera pitch)
    # Para manter o pitch, precisamos de um algoritmo mais complexo ou aceitar pequena mudança.
    # O método speedup do pydub tenta manter o pitch cortando pedacinhos.
    audio_rapido = audio.speedup(playback_speed=SPEED_FACTOR)
    
    audio_rapido.export(filename, format="mp3")
    
    # Limpa temp
    os.remove(temp_file)
    print(f"   ✅ Salvo: {filename} ({time.time()-start:.2f}s)")

# Textos de teste
texto_ivo = "Olá, eu sou o Ivo. Estamos testando uma nova velocidade de fala para o RevaCast, para que o conteúdo fique mais dinâmico."
texto_manu = "Oi, aqui é a Manu. Eu também vou falar um pouco mais rápido, assim conseguimos passar mais informação em menos tempo."

# Executa
gerar_e_acelerar(ELEVEN_VOICE_ID_HOST, "Ivo (Host)", texto_ivo, "teste_ivo_rapido.mp3")
gerar_e_acelerar(ELEVEN_VOICE_ID_COHOST, "Manu (Cohost)", texto_manu, "teste_manu_rapido.mp3")

print("\n🎧 Teste concluído! Verifique os arquivos 'teste_ivo_rapido.mp3' e 'teste_manu_rapido.mp3'.")
