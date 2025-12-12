from revamais_service import criar_campanha_revamais
import sys

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)

# Tema de teste
tema = "Benefícios da caminhada para idosos"

print(f"🧪 Iniciando teste do Reva + com o tema: '{tema}'")

try:
    # Pass a simple flush callback to ensure logs appear immediately
    def log_callback(msg):
        print(f"[LOG] {msg}")
        sys.stdout.flush()

    resultado = criar_campanha_revamais(tema, log_callback=log_callback)
    print("\n✅ Teste finalizado com sucesso!")
    print(resultado)
except Exception as e:
    print(f"\n❌ Erro no teste: {e}")
    import traceback
    traceback.print_exc()
