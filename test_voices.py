"""
Script de teste para verificar se os IDs de voz do ElevenLabs são válidos
"""
import os
from dotenv import load_dotenv
from elevenlabs.client import ElevenLabs

load_dotenv()

ELEVENLABS_API_KEY = os.environ.get("ELEVENLABS_API_KEY")
ELEVEN_VOICE_ID_HOST = os.environ.get("ELEVEN_VOICE_ID_HOST")
ELEVEN_VOICE_ID_COHOST = os.environ.get("ELEVEN_VOICE_ID_COHOST")

print("=" * 60)
print("TESTE DE VOZES DO ELEVENLABS")
print("=" * 60)
print(f"\nAPI Key: {ELEVENLABS_API_KEY[:20]}...")
print(f"HOST Voice ID: {ELEVEN_VOICE_ID_HOST}")
print(f"COHOST Voice ID: {ELEVEN_VOICE_ID_COHOST}")
print("\n" + "=" * 60)

client = ElevenLabs(api_key=ELEVENLABS_API_KEY)

# Lista todas as vozes disponíveis na sua conta
print("\n📋 VOZES DISPONÍVEIS NA SUA CONTA:")
print("=" * 60)

try:
    voices = client.voices.get_all()
    
    if hasattr(voices, 'voices'):
        for voice in voices.voices:
            print(f"\n✓ Nome: {voice.name}")
            print(f"  ID: {voice.voice_id}")
            print(f"  Categoria: {voice.category if hasattr(voice, 'category') else 'N/A'}")
            
            # Verifica se é uma das vozes configuradas
            if voice.voice_id == ELEVEN_VOICE_ID_HOST:
                print(f"  ⭐ Esta é a voz HOST configurada!")
            if voice.voice_id == ELEVEN_VOICE_ID_COHOST:
                print(f"  ⭐ Esta é a voz COHOST configurada!")
    else:
        print("Nenhuma voz encontrada")
        
except Exception as e:
    print(f"❌ Erro ao buscar vozes: {e}")

print("\n" + "=" * 60)
print("\n🔍 VERIFICANDO VOICE IDs CONFIGURADOS:")
print("=" * 60)

# Testa se os IDs estão válidos
def testar_voz(voice_id, nome):
    try:
        print(f"\n{nome} ({voice_id}):")
        # Tenta gerar um áudio pequeno de teste
        audio = client.text_to_speech.convert(
            voice_id=voice_id,
            text="Teste",
            model_id="eleven_multilingual_v2"
        )
        # Consome o gerador
        for _ in audio:
            pass
        print(f"  ✅ ID VÁLIDO - Voz encontrada e funcionando!")
        return True
    except Exception as e:
        print(f"  ❌ ID INVÁLIDO - Erro: {e}")
        return False

host_ok = testar_voz(ELEVEN_VOICE_ID_HOST, "HOST")
cohost_ok = testar_voz(ELEVEN_VOICE_ID_COHOST, "COHOST")

print("\n" + "=" * 60)
print("RESULTADO FINAL:")
print("=" * 60)
if host_ok and cohost_ok:
    print("✅ Ambas as vozes estão configuradas corretamente!")
else:
    print("❌ Há problemas com os IDs de voz configurados.")
    print("   Verifique os IDs no seu .env e compare com a lista acima.")
