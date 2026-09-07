import unittest

from connected_papers import build_manual_report, parse_bibtex


SAMPLE = """
@article{anchor,
title = {Effects of downhill walking},
year = {2020},
url = {https://www.semanticscholar.org/paper/anchor},
abstract = {Main abstract with {nested} braces.},
author = {C. Camillo and C. Osadnik},
journal = {European Respiratory Journal},
doi = {10.1234/anchor},
pmid = {32444407},
}

@article{related,
title = {Earlier eccentric training},
year = {2017},
url = {https://www.semanticscholar.org/paper/related},
abstract = {A useful abstract with a comma, and multiple lines.},
author = {A. Researcher and B. Scientist},
journal = {Chest},
doi = {10.1234/related},
pmid = {12345678},
}

@article{doi_only,
title = {DOI only},
year = {2015},
abstract = {Short abstract.},
author = {C. Author},
doi = {10.1234/doi-only},
journal = {null},
}
"""


class ConnectedPapersTests(unittest.TestCase):
    def test_parse_anchor_and_candidates(self):
        result = parse_bibtex(SAMPLE)
        self.assertEqual(result["total_entradas"], 3)
        self.assertEqual(result["ancora"]["pmid"], "32444407")
        self.assertEqual(len(result["candidatos"]), 2)
        self.assertEqual(result["candidatos"][0]["autores"], ["A. Researcher", "B. Scientist"])
        self.assertIn("nested", result["ancora"]["resumo_original"])
        self.assertEqual(result["candidatos"][1]["journal"], "")

    def test_malformed_and_empty_files_fail_closed(self):
        with self.assertRaises(ValueError):
            parse_bibtex("")
        with self.assertRaises(ValueError):
            parse_bibtex("@article{broken, title = {missing")

    def test_manual_report_matches_podcast_contract(self):
        parsed = parse_bibtex(SAMPLE)
        anchor = parsed["ancora"]
        reference = parsed["candidatos"][0]
        report = build_manual_report([anchor], {anchor["pmid"]: [reference]})
        self.assertEqual(report["modo"], "manual")
        self.assertEqual(report["estudos"][0]["referencias"][0]["pmid"], "12345678")
        self.assertEqual(report["estudos"][0]["referencias"][0]["triagem"]["origem"], "curadoria_manual")


if __name__ == "__main__":
    unittest.main()
