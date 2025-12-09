import os
from elevenlabs.client import ElevenLabs
from elevenlabs import VoiceSettings
from dotenv import load_dotenv

# Carrega variáveis de ambiente
load_dotenv()

API_KEY = os.environ.get("ELEVENLABS_API_KEY")
VOICE_ID = "p5oveq8dCbyBIAaD6gzR" # Forçando ID novo para teste

def test_host():
    if not API_KEY:
        print("❌ Erro: ELEVENLABS_API_KEY não encontrada no .env")
        return

    client = ElevenLabs(api_key=API_KEY)
    
    text = "Olá, este é um teste rápido para verificar a qualidade da voz do Host no modelo Turbo 2.5."
    
    print(f"🎙️ Gerando áudio para Host ({VOICE_ID})...")
    print("⚙️ Config: Model=eleven_turbo_v2_5, Stability=0.40, Similarity=0.80, Style=0.60")

    try:
        audio_generator = client.text_to_speech.convert(
            voice_id=VOICE_ID,
            text=text,
            model_id="eleven_turbo_v2_5",
            voice_settings=VoiceSettings(
                stability=0.40,
                similarity_boost=0.80,
                style=0.60,
                use_speaker_boost=True
            )
        )
        
        output_file = "teste_host_isolado.mp3"
        with open(output_file, "wb") as f:
            for chunk in audio_generator:
                f.write(chunk)
                
        print(f"✅ Áudio salvo em: {output_file}")
        
    except Exception as e:
        print(f"❌ Erro na geração: {e}")

if __name__ == "__main__":
    test_host()
