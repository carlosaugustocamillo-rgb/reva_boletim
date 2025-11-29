"""
Script para fazer upload dos testes de áudio para o Firebase e retornar URLs públicas
"""
import os
from firebase_service import upload_file

# Arquivos gerados pelo teste anterior
files = ["teste_ivo_rapido.mp3", "teste_manu_rapido.mp3"]

print("📤 Iniciando upload dos testes para o Firebase...")

for filename in files:
    if os.path.exists(filename):
        destination_blob = f"testes_velocidade/{filename}"
        try:
            url = upload_file(filename, destination_blob)
            print(f"\n✅ {filename} disponível em:\n👉 {url}")
        except Exception as e:
            print(f"❌ Erro ao subir {filename}: {e}")
    else:
        print(f"⚠️ Arquivo {filename} não encontrado localmente.")
