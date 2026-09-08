import unittest
import ast
from datetime import date
from pathlib import Path


def load_calendar_helpers():
    source = Path("revamais_service.py").read_text(encoding="utf-8")
    module = ast.parse(source)
    needed = {
        "_extract_calendar_title",
        "_parse_calendar_date",
        "_compute_next_pending_index",
    }
    funcs = [
        node for node in module.body
        if isinstance(node, ast.FunctionDef) and node.name in needed
    ]
    namespace = {
        "date": date,
        "_today_for_calendar": lambda: date(2026, 9, 8),
    }
    exec(compile(ast.Module(body=funcs, type_ignores=[]), "calendar_helpers", "exec"), namespace)
    return namespace["_compute_next_pending_index"]


class RevaMaisCalendarSelectionTest(unittest.TestCase):
    def setUp(self):
        self.compute_next_pending_index = load_calendar_helpers()

    def test_prefers_pending_item_from_today_over_old_backlog(self):
        rows = [
            {"Title": "Old pending", "Date": "2026-04-21"},
            {"Title": "Already done today", "Date": "2026-09-08"},
            {"Title": "Today pending", "Date": "2026-09-08"},
            {"Title": "Future pending", "Date": "2026-09-11"},
        ]

        current = self.compute_next_pending_index(
            rows,
            completed_indices={1},
            today=date(2026, 9, 8),
        )

        self.assertEqual(current, 2)

    def test_falls_back_to_old_backlog_when_no_current_or_future_pending_exists(self):
        rows = [
            {"Title": "Old done", "Date": "2026-04-21"},
            {"Title": "Old pending", "Date": "2026-04-24"},
        ]

        current = self.compute_next_pending_index(
            rows,
            completed_indices={0},
            today=date(2026, 9, 8),
        )

        self.assertEqual(current, 1)

    def test_skips_rows_without_valid_titles_or_dates_for_date_priority(self):
        rows = [
            {"Title": "", "Date": "2026-09-08"},
            {"Title": "Bad date backlog", "Date": "invalid"},
            {"Title": "Future pending", "Date": "2026-09-11"},
        ]

        current = self.compute_next_pending_index(
            rows,
            completed_indices=set(),
            today=date(2026, 9, 8),
        )

        self.assertEqual(current, 2)


if __name__ == "__main__":
    unittest.main()
