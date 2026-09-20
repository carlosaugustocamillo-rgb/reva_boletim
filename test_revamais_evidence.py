import unittest
from unittest.mock import patch

import revamais_editorial as editorial
import revamais_service as service


class RevaMaisEvidenceTest(unittest.TestCase):
    def test_consensus_citation_is_not_promoted_to_abstract(self):
        reference = service._normalizar_referencia_consensus({
            "texto": "Silva et al. Exercise rehabilitation. Journal, 2024.",
            "doi": "10.1000/test",
        }, 0)
        self.assertEqual(reference["material"], "bibliographic_only")
        self.assertEqual(reference["evidence_content"], "")
        packet = editorial.evidence_packet([reference])
        self.assertEqual(packet[0]["material"], "bibliographic_only")
        self.assertEqual(packet[0]["content"], "")

    def test_mixed_recommendation_record_is_quarantined(self):
        report = (
            "Clinical summary with enough content to remain outside the references section.\n\n"
            "References\n"
            "Immediate Response: terminate exertion. Emergency Magnet Protocol: use a magnet. "
            "Karaoguz et al. Exercise participation. https://doi.org/10.5543/tkda.2024.24952"
        )
        references = service._extrair_referencias_da_secao_consensus(report)
        self.assertEqual(references[0]["material"], "invalid_or_mixed")
        self.assertTrue(references[0]["exclude_from_bibliography"])

    def test_incomplete_doi_is_not_published_as_a_link(self):
        reference = service._normalizar_referencia_consensus({
            "texto": "Network meta-analysis",
            "doi": "10.1038/s41598-",
            "link": "https://doi.org/10.1038/s41598-",
        }, 0)
        self.assertEqual(reference["doi"], "")
        self.assertEqual(reference["link"], "")

    @patch("revamais_service.enriquecer_referencias_pubmed", side_effect=lambda refs: refs)
    def test_report_narrative_is_shared_as_identified_evidence(self, _mock_enrich):
        report = (
            "Exercise rehabilitation after device implantation\n\n"
            "A systematic review included randomized studies in clinically stable adults. "
            "The intervention improved functional capacity compared with usual care. "
            "Adverse events were similar between groups, although certainty was limited. "
            "Recommendations should be individualized according to device programming.\n\n"
            "References\nSilva et al. Exercise rehabilitation. https://doi.org/10.1000/test"
        )
        evidence, readiness = service.preparar_evidencias_editoriais([
            {"texto": "Silva et al.", "doi": "10.1000/test", "material": "bibliographic_only"}
        ], report)
        self.assertTrue(readiness["ready"])
        self.assertEqual(readiness["usable_sources"], 1)
        report_source = next(item for item in evidence if item.get("source_id") == "consensus-report:1")
        self.assertNotIn("References", report_source["evidence_content"])

    @patch("revamais_service.enriquecer_referencias_pubmed", side_effect=lambda refs: refs)
    def test_bibliography_alone_stops_generation_readiness(self, _mock_enrich):
        evidence, readiness = service.preparar_evidencias_editoriais([
            {"texto": "Silva et al.", "doi": "10.1000/test", "material": "bibliographic_only"}
        ], "References\nSilva et al. https://doi.org/10.1000/test")
        self.assertFalse(readiness["ready"])
        self.assertEqual(readiness["bibliographic_only"], 1)
        self.assertEqual(len(evidence), 1)

    @patch("revamais_service.enriquecer_referencias_pubmed", side_effect=lambda refs: refs)
    def test_existing_pubmed_abstract_is_ready_without_consensus(self, _mock_enrich):
        evidence, readiness = service.preparar_evidencias_editoriais([{
            "pmid": "123",
            "texto": "Supervised walking trial",
            "resumo": "The trial observed improved mobility after supervised walking.",
            "fonte": "PubMed",
        }])
        self.assertTrue(readiness["ready"])
        self.assertEqual(evidence[0]["material"], "abstract")

    @patch("revamais_service.enriquecer_referencias_pubmed", side_effect=lambda refs: refs)
    def test_legacy_citation_labeled_as_abstract_is_downgraded(self, _mock_enrich):
        citation = "Silva et al. Exercise rehabilitation. Journal, 2024. https://doi.org/10.1000/test"
        evidence, readiness = service.preparar_evidencias_editoriais([{
            "texto": citation,
            "resumo": citation,
            "evidence_content": citation,
            "material": "abstract",
            "doi": "10.1000/test",
            "fonte": "Consensus",
        }])
        self.assertFalse(readiness["ready"])
        self.assertEqual(evidence[0]["material"], "bibliographic_only")
        self.assertEqual(evidence[0]["evidence_content"], "")


if __name__ == "__main__":
    unittest.main()
