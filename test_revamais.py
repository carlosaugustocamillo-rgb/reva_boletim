from revamais_service import criar_campanha_revamais
import sys

# Tema de teste
tema = "Benefícios da caminhada para idosos"

print(f"🧪 Iniciando teste do Reva + com o tema: '{tema}'")

try:
    resultado = criar_campanha_revamais(tema)
    print("\n✅ Teste finalizado com sucesso!")
    print(resultado)
except Exception as e:
    print(f"\n❌ Erro no teste: {e}")
    import traceback
    traceback.print_exc()
