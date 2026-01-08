#!/usr/bin/env python3
"""
run_medico_boletim.py

Script de linha de comando para rodar manualmente o Boletim Médico Semanal (Q1 Only).
Uso:
    python run_medico_boletim.py [--dry-run]

--dry-run: Gera o HTML localmente mas NÃO envia para o Mailchimp.
"""

import argparse
from medico_boletim_service import rodar_boletim_medico

def main():
    parser = argparse.ArgumentParser(description="Gera o Boletim Médico Semanal (Q1 Only).")
    parser.add_argument("--dry-run", action="store_true", help="Executa sem enviar para o Mailchimp (apenas gera HTML).")
    parser.add_argument("--output", default="boletim_medico_preview.html", help="Nome do arquivo HTML de saída.")
    
    args = parser.parse_args()
    
    print(f"🔬 Iniciando execução...")
    if args.dry_run:
        print("⚠️  MODO DRY-RUN ATIVADO: Nenhuma campanha será criada no Mailchimp.")
        
    try:
        rodar_boletim_medico(dry_run=args.dry_run, output_file=args.output)
    except KeyboardInterrupt:
        print("\n🛑 Interrompido pelo usuário.")
    except Exception as e:
        print(f"\n❌ Erro fatal: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
