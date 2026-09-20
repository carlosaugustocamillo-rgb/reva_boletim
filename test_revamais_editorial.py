import tempfile
import unittest

import revamais_editorial as editorial


def sample_result():
    content = (
        "<h1>Caminhar com segurança</h1>"
        "<p>A caminhada supervisionada pode ajudar a mobilidade.</p>"
        '<img src="https://example.com/science.png" alt="Ciência">'
        "<h2>O que as evidências mostram?</h2>"
        "<p>O estudo observou melhora de mobilidade.</p>"
        '<div class="references"><h3>Referências</h3><ul><li>Study</li></ul></div>'
    )
    return {
        "tema": "Caminhada",
        "titulo": "Caminhar com segurança",
        "calendar_index": 4,
        "email_requested": True,
        "html_content": content,
        "html_full": (
            "<html><body><main><!-- REVAMAIS_CONTENT_START -->"
            + content
            + "<!-- REVAMAIS_CONTENT_END --></main></body></html>"
        ),
        "instagram_assets": [{
            "asset_id": "instagram-slide-1",
            "type": "image",
            "url": "https://example.com/slide.png",
            "name": "Slide 1",
            "texto_base": "Caminhar pode ajudar a mobilidade.",
        }],
        "visual_assets": [
            {
                "id": "newsletter-opening",
                "kind": "newsletter_opening",
                "label": "Abertura",
                "url": "https://example.com/opening.png",
                "prompt": "A patient walking in a park",
            },
            {
                "id": "newsletter-science",
                "kind": "newsletter_science",
                "label": "Ciência",
                "url": "https://example.com/science.png",
                "prompt": "Mobility mechanism",
            },
        ],
        "referencias_utilizadas": [{
            "pmid": "123",
            "texto": "Walking study",
            "resumo": "The study observed improved mobility after supervised walking.",
            "link": "https://pubmed.ncbi.nlm.nih.gov/123/",
            "tipo_estudo": "trial",
        }],
    }


def passing_audit(_payload):
    return {
        "issues": [],
        "claims": [{
            "claim": "O estudo observou melhora de mobilidade.",
            "location": "Seção científica",
            "verdict": "supported",
            "reason": "O resumo descreve o achado.",
            "supports": [{
                "source_id": "pmid:123",
                "quote": "The study observed improved mobility after supervised walking.",
            }],
        }],
        "appraisals": [{
            "source_id": "pmid:123",
            "scope": "Ensaio sobre caminhada supervisionada e mobilidade.",
            "caution": "Aplicar o resultado à população descrita no estudo.",
            "observations": [{
                "kind": "design_scope",
                "explanation": "O material descreve caminhada supervisionada.",
                "quote": "supervised walking",
            }],
        }],
    }


class RevaMaisEditorialTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.base_dir = self.temp.name

    def tearDown(self):
        self.temp.cleanup()

    def create_audited(self):
        draft = editorial.create_draft(self.base_dir, sample_result(), source_task_id="task-1")
        return editorial.audit_draft(self.base_dir, draft["id"], client=None, audit_fn=passing_audit)

    def test_sanitizer_removes_active_content_and_unsafe_urls(self):
        clean = editorial.sanitize_html(
            '<p onclick="steal()">Texto</p><script>alert(1)</script>'
            '<a href="javascript:alert(2)" target="_blank">link</a>'
        )
        self.assertEqual(clean, '<p>Texto</p><a target="_blank" rel="noreferrer noopener">link</a>')

    def test_audit_supports_approval_with_literal_source_quote(self):
        draft = self.create_audited()
        self.assertEqual(draft["status"], "pending_review")
        self.assertTrue(editorial.review_payload(draft)["can_approve"])
        self.assertEqual(draft["audit"]["appraisals"][0]["source_id"], "pmid:123")

        approved = editorial.approve_draft(self.base_dir, draft["id"], draft["sha256"])
        self.assertEqual(approved["status"], "approved")

    def test_invalid_support_quote_becomes_blocking(self):
        draft = editorial.create_draft(self.base_dir, sample_result())

        def invalid_audit(_payload):
            value = passing_audit(_payload)
            value["claims"][0]["supports"][0]["quote"] = "A sentence absent from the source."
            return value

        blocked = editorial.audit_draft(self.base_dir, draft["id"], None, audit_fn=invalid_audit)
        self.assertEqual(blocked["status"], "blocked")
        self.assertTrue(any(issue["severity"] == "blocking" for issue in blocked["audit"]["issues"]))

    def test_incomplete_audit_cannot_skip_claims_or_source_appraisal(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        blocked = editorial.audit_draft(
            self.base_dir,
            draft["id"],
            None,
            audit_fn=lambda _payload: {"issues": [], "claims": [], "appraisals": []},
        )
        reasons = " ".join(issue["reason"] for issue in blocked["audit"]["issues"])
        self.assertEqual(blocked["status"], "blocked")
        self.assertIn("nenhuma afirmação", reasons)
        self.assertIn("não foram conferidos", reasons)

    def test_text_edit_creates_child_revision_and_rebuilds_email_preview(self):
        previous = self.create_audited()
        edited_html = (
            "<h1>Nova versão segura</h1><p>Texto revisado com extensão suficiente para continuar válido e auditável.</p>"
            "<h2>O que as evidências mostram?</h2><p>O estudo observou melhora de mobilidade.</p>"
        )
        edited = editorial.save_edited_draft(
            self.base_dir,
            previous["id"],
            previous["sha256"],
            "Nova versão segura",
            edited_html,
            previous["content"]["instagram_assets"],
        )
        self.assertNotEqual(edited["id"], previous["id"])
        self.assertEqual(edited["parent_id"], previous["id"])
        self.assertEqual(edited["status"], "auditing")
        self.assertIn("Nova versão segura", edited["content"]["html_full"])
        self.assertIn("https://example.com/opening.png", edited["content"]["html_full"])
        self.assertNotIn('<div class="references">', edited["content"]["html_full"])
        untouched = editorial.load_draft(self.base_dir, previous["id"])
        self.assertNotIn("Nova versão segura", untouched["content"]["html_content"])

    def test_single_asset_regeneration_preserves_other_assets_and_updates_html(self):
        previous = self.create_audited()
        updated = editorial.save_regenerated_asset(
            self.base_dir,
            previous["id"],
            previous["sha256"],
            "newsletter-science",
            "https://example.com/science-v2.png",
            "Revised mobility mechanism",
        )
        science = next(item for item in updated["visual_assets"] if item["id"] == "newsletter-science")
        opening = next(item for item in updated["visual_assets"] if item["id"] == "newsletter-opening")
        self.assertEqual(science["url"], "https://example.com/science-v2.png")
        self.assertEqual(len(science["versions"]), 2)
        self.assertEqual(opening["url"], "https://example.com/opening.png")
        self.assertIn("science-v2.png", updated["content"]["html_content"])

    def test_unapproved_revision_cannot_be_published(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        with self.assertRaises(editorial.RevaMaisEditorialError):
            editorial.approved_draft(self.base_dir, draft["id"], draft["sha256"])

    def test_partial_publication_progress_is_saved_without_losing_approval(self):
        audited = self.create_audited()
        approved = editorial.approve_draft(self.base_dir, audited["id"], audited["sha256"])
        progress = editorial.save_publication_progress(
            self.base_dir,
            approved["id"],
            approved["sha256"],
            {"status": "partial", "campaign_id": "campaign-1", "email_content_set": True},
        )
        self.assertEqual(progress["status"], "approved")
        self.assertEqual(progress["sha256"], approved["sha256"])
        self.assertEqual(progress["publication_progress"]["campaign_id"], "campaign-1")


if __name__ == "__main__":
    unittest.main()
