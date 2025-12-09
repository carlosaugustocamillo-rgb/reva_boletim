import requests
import time

URL = "https://revaboletim-production.up.railway.app/teste-audio-curto"

print(f"🚀 Chamando endpoint de teste no Railway: {URL}")
try:
    start = time.time()
    response = requests.get(URL, timeout=120) # Timeout longo pois gera áudio
    end = time.time()
    
    if response.status_code == 200:
        data = response.json()
        print("\n✅ SUCESSO!")
        print(f"⏱️ Tempo: {end - start:.2f}s")
        print(f"🔗 URL do Áudio: {data.get('url')}")
        print(f"🗣️ Host Voice: {data.get('host_voice')}")
        print(f"🗣️ Cohost Voice: {data.get('cohost_voice')}")
    else:
        print(f"\n❌ Erro {response.status_code}:")
        print(response.text)

except Exception as e:
    print(f"\n❌ Falha na requisição: {e}")
