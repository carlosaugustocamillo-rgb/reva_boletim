"""
Script de diagnóstico para testar a geração de áudio isoladamente e capturar erros.
"""
import os
import sys
from dotenv import load_dotenv

# Carrega variáveis
load_dotenv()

print("🔍 Iniciando diagnóstico de áudio...")

try:
    # Tenta importar o módulo principal
    print("1. Importando boletim_service...")
    from boletim_service import rodar_boletim, resumo_para_podcast
    print("   ✅ Importação bem sucedida!")
except Exception as e:
    print(f"   ❌ Erro na importação: {e}")
    sys.exit(1)

try:
    # Tenta simular uma geração de áudio curta
    print("\n2. Simulando geração de áudio...")
    
    # Mock de um roteiro simples
    roteiro_teste = [
        {"speaker": "HOST", "text": "Teste de som, um dois três."},
        {"speaker": "COHOST", "text": "Confirmando teste de som."}
    ]
    
    # Precisamos "enganar" a função rodar_boletim para ela achar que já tem roteiro
    # Ou podemos chamar a parte do áudio diretamente se refatorarmos, mas vamos tentar rodar o fluxo
    # com opções restritas.
    
    # Como rodar_boletim é um gerador, precisamos iterar
    print("   Executando pipeline (apenas áudio)...")
    
    # Hack: Vamos injetar o roteiro diretamente na variável local da função? Não dá.
    # Vamos criar um arquivo de roteiro fake para ele ler.
    from datetime import datetime
    hoje = datetime.today().strftime('%Y-%m-%d')
    base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    roteiro_path = os.path.join(base_dir, f"roteiro_podcast_{hoje}.txt")
    
    with open(roteiro_path, "w", encoding="utf-8") as f:
        f.write("ROTEIRO TESTE\n\nHOST: Teste de som.\n\nCOHOST: Teste confirmado.")
        
    # Agora rodamos com opção 'audio' apenas. 
    # OBS: O código atual do boletim_service lê o roteiro da memória (roteiros_audio) se gerado na hora,
    # ou PULA se não tiver. Ele não lê do arquivo txt para gerar o áudio (eu notei isso na leitura anterior).
    # Ele diz: "Se pulou resumos, não tem roteiro em memória... pass".
    
    # Então precisamos testar a função de geração de áudio DIRETAMENTE.
    # Mas ela está dentro de rodar_boletim.
    
    # Vamos testar apenas a importação e configuração do ElevenLabs por enquanto,
    # pois se o erro for de sintaxe ou import, já teria aparecido no passo 1.
    
    from elevenlabs.client import ElevenLabs
    from elevenlabs import VoiceSettings
    
    api_key = os.environ.get("ELEVENLABS_API_KEY")
    if not api_key:
        print("   ❌ ELEVENLABS_API_KEY não encontrada!")
    else:
        client = ElevenLabs(api_key=api_key)
        print("   ✅ Cliente ElevenLabs inicializado.")
        
        print("   Tentando gerar áudio com VoiceSettings...")
        audio = client.text_to_speech.convert(
            voice_id="WWL28Z00upcD5SGFqY2n", # Host
            text="Teste de diagnóstico.",
            model_id="eleven_multilingual_v2",
            voice_settings=VoiceSettings(stability=0.4, similarity_boost=0.8, style=0.6, use_speaker_boost=True)
        )
        # Consome
        for _ in audio: pass
        print("   ✅ Geração de áudio com VoiceSettings funcionou!")

except Exception as e:
    print(f"\n❌ ERRO DETECTADO DURANTE A EXECUÇÃO:\n{e}")
    import traceback
    traceback.print_exc()

print("\n🏁 Diagnóstico concluído.")
