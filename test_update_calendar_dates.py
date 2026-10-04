import unittest
from datetime import datetime

from update_calendar_dates import update_planned_rows


class UpdateCalendarDatesTest(unittest.TestCase):
    def test_updates_only_planned_rows(self):
        rows = [
            {
                "ID": "mc-published",
                "Week": "1",
                "Date": "2026-01-01",
                "Day": "Quinta",
                "Status": "Concluído",
            },
            {"ID": "planned-1", "Status": "Planejado"},
            {"ID": "planned-2", "Status": "Planejado"},
            {"ID": "planned-3", "Status": "Planejado"},
        ]

        updated = update_planned_rows(rows, datetime(2026, 10, 6))

        self.assertEqual(updated, 3)
        self.assertEqual(rows[0]["Date"], "2026-01-01")
        self.assertEqual(rows[0]["Day"], "Quinta")
        self.assertEqual(rows[1]["Date"], "2026-10-06")
        self.assertEqual(rows[1]["Day"], "Terça")
        self.assertEqual(rows[2]["Date"], "2026-10-09")
        self.assertEqual(rows[2]["Day"], "Sexta")
        self.assertEqual(rows[3]["Date"], "2026-10-13")
        self.assertEqual(rows[3]["Day"], "Terça")


if __name__ == "__main__":
    unittest.main()
