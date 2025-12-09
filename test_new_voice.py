
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
HOST_ID = "p5oveq8dCbyBIAaD6gzR"
COHOST_ID = "tnSpp4vdxKPjI9w0GnoV"

print(f"🔹 Teste rápido de áudio com:")
print(f"   HOST: {HOST_ID}")
print(f"   COHOST: {COHOST_ID}")

client = ElevenLabs(api_key=API_KEY)

# Lista de falas para o teste (5 falas cada)
dialogo = [
    {"speaker": "HOST", "text": "Vamos lá então, pessoal! O nosso primeiro estudo de hoje fala sobre os efeitos do exercício na DPOC."},
    {"speaker": "COHOST", "text": "Isso mesmo. E os resultados parecem bem promissores, não é? O que eles descobriram exatamente?"},
    {"speaker": "HOST", "text": "Pois é. Eles viram que a reabilitação pulmonar aumentou significativamente a capacidade funcional dos pacientes."},
    {"speaker": "COHOST", "text": "Interessante! E quanto tempo durou essa intervenção? Foi algo curto ou de longo prazo?"},
    {"speaker": "HOST", "text": "Foram 12 semanas de treinamento resistido e aeróbico, três vezes por semana."},
    {"speaker": "COHOST", "text": "Bastante tempo. E imagino que a qualidade de vida também tenha melhorado com esse protocolo."},
    {"speaker": "HOST", "text": "Exatamente. Os questionários de qualidade de vida mostraram pontuações muito superiores ao grupo controle."},
    {"speaker": "COHOST", "text": "Muito bom saber disso. Fica a dica clínica importante para nossos ouvintes fisioterapeutas."},
    {"speaker": "HOST", "text": "Com certeza. Bom, vamos em frente... nesse próximo estudo a gente vai falar sobre asma grave."},
    {"speaker": "COHOST", "text": "Opa, esse também é um tema quente. Vamos ver o que tem de novo."}
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
