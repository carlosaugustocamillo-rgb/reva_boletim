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


def update_planned_rows(rows, first_tuesday):
    planned_rows = [
        row for row in rows
        if str(row.get('Status', 'Planejado')).strip().casefold() == 'planejado'
    ]

    for planned_index, row in enumerate(planned_rows):
        weeks_passed = planned_index // 2
        days_after_tuesday = 0 if planned_index % 2 == 0 else 3
        new_date = first_tuesday + timedelta(
            weeks=weeks_passed,
            days=days_after_tuesday,
        )
        row['Date'] = new_date.strftime('%Y-%m-%d')
        row['Day'] = 'Terça' if planned_index % 2 == 0 else 'Sexta'
        row['Week'] = new_date.isocalendar().week

    return len(planned_rows)

def update_csv_dates():
    csv_file = 'calendario_editorial_150_semanas.csv'
    hoje = datetime.now()
    proxima_terca = get_next_weekday(hoje, 1)

    rows = []
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    updated_count = update_planned_rows(rows, proxima_terca)
    print(
        f"📅 Atualizando {updated_count} pautas planejadas a partir de "
        f"{proxima_terca.strftime('%Y-%m-%d')} (Terça)"
    )

    # Salva
    with open(csv_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        
    print(f"✅ {updated_count} datas planejadas atualizadas; histórico preservado.")

if __name__ == "__main__":
    update_csv_dates()
