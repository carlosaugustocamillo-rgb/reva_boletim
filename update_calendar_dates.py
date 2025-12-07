import csv
from datetime import datetime, timedelta

def get_next_weekday(start_date, weekday):
    """
    Get the next date for the given weekday (0=Monday, 1=Tuesday, ...).
    If start_date is already that weekday, returns start_date (or next week depending on logic, here we want next occurrence).
    """
    days_ahead = weekday - start_date.weekday()
    if days_ahead <= 0: # Target day already happened this week
        days_ahead += 7
    return start_date + timedelta(days=days_ahead)

def update_csv_dates():
    csv_file = 'calendario_editorial_150_semanas.csv'
    
    # Próxima Terça-feira a partir de hoje (07/12/2025 - Domingo)
    # Terça = 1
    hoje = datetime.now()
    # Para teste/garantia, vamos fixar a data base
    # hoje = datetime(2025, 12, 7) 
    
    proxima_terca = get_next_weekday(hoje, 1) # 1 = Tuesday
    proxima_sexta = get_next_weekday(hoje, 4) # 4 = Friday
    
    # Se a próxima sexta for ANTES da próxima terça (ex: hoje é quarta), ajusta
    # Nesse caso (Domingo), Terça vem antes (Terca dia 9, Sexta dia 12). OK.
    
    current_date = proxima_terca
    
    rows = []
    # Use utf-8-sig to handle potential BOM from Excel
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)
            
    # Lógica de atualização
    # O arquivo original alterna Terça e Sexta. Vamos manter isso.
    # Semana 1: Terça, Sexta
    # Semana 2: Terça, Sexta...
    
    # Data inicial para o loop
    # Vamos assumir que a primeira linha É a próxima terça.
    # E vamos recalcular todas as datas subsequentes baseadas nisso.
    
    date_cursor = proxima_terca
    
    print(f"📅 Atualizando datas iniciando em: {date_cursor.strftime('%Y-%m-%d')} (Terça)")
    
    for i, row in enumerate(rows):
        # Define se é Terça ou Sexta baseado no indice ou na coluna 'Day' original
        # O arquivo original tem coluna 'Day'. Vamos confiar nela ou alternar?
        # Vamos olhar a coluna 'Day' para ser mais seguro.
        
        dia_semana_str = row.get('Day', '').lower()
        
        # Se for a primeira linha, já definimos como proxima_terca
        if i == 0:
            new_date = proxima_terca
        else:
            # Lógica simples: se a linha anterior era Terça, essa deve ser Sexta (mesma semana)
            # Se a linha anterior era Sexta, essa deve ser Terça (próxima semana)
            
            # Mas vamos recalcular tudo sequencialmente para garantir
            # Padrão observado: Terça -> Sexta -> Terça (prox semana)
            
            # Se o indice é PAR (0, 2, 4...), é Terça (início da semana)
            # Se o indice é IMPAR (1, 3, 5...), é Sexta
            
            if i % 2 == 0: # Terça
                # Soma 4 dias em relação a Sexta anterior? Não.
                # Terça (0) -> Sexta (1) -> Terça (2)
                # Data = Data da Sexta anterior + 4 dias (Sexta->Sábado->Domingo->Segunda->Terça)
                # Ou Data inicial + (semanas * 7)
                
                weeks_passed = i // 2
                new_date = proxima_terca + timedelta(weeks=weeks_passed)
                row['Day'] = 'Terça'
            else: # Sexta
                # Sexta é sempre 3 dias depois da Terça DA MESMA semana
                # Data = Data da Terça dessa semana + 3 dias
                weeks_passed = i // 2
                terca_da_semana = proxima_terca + timedelta(weeks=weeks_passed)
                new_date = terca_da_semana + timedelta(days=3)
                row['Day'] = 'Sexta'

        row['Date'] = new_date.strftime('%Y-%m-%d')
        
        # Atualiza o numero da semana também, só pra garantir
        row['Week'] = (i // 2) + 1

    # Salva
    with open(csv_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        
    print(f"✅ {len(rows)} datas atualizadas com sucesso!")

if __name__ == "__main__":
    update_csv_dates()
