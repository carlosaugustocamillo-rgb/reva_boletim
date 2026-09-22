import tempfile
import unittest
import os
import re

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
        self.previous_remote_persistence = os.environ.get("REVAMAIS_EDITORIAL_PERSIST_REMOTE")
        os.environ["REVAMAIS_EDITORIAL_PERSIST_REMOTE"] = "false"
        self.temp = tempfile.TemporaryDirectory()
        self.base_dir = self.temp.name

    def tearDown(self):
        self.temp.cleanup()
        if self.previous_remote_persistence is None:
            os.environ.pop("REVAMAIS_EDITORIAL_PERSIST_REMOTE", None)
        else:
            os.environ["REVAMAIS_EDITORIAL_PERSIST_REMOTE"] = self.previous_remote_persistence

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

    def test_invalid_support_quote_requires_revision_without_safety_block(self):
        draft = editorial.create_draft(self.base_dir, sample_result())

        def invalid_audit(_payload):
            value = passing_audit(_payload)
            value["claims"][0]["supports"][0]["quote"] = "A sentence absent from the source."
            return value

        blocked = editorial.audit_draft(self.base_dir, draft["id"], None, audit_fn=invalid_audit)
        self.assertEqual(blocked["status"], "needs_revision")
        self.assertTrue(any(issue["severity"] == "needs_revision" for issue in blocked["audit"]["issues"]))

    def test_admin_can_explicitly_override_a_failed_audit_and_approve_it(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        blocked = editorial.audit_draft(
            self.base_dir, draft["id"], None,
            audit_fn=lambda _payload: {"issues": [], "claims": [], "appraisals": []},
        )
        released = editorial.override_audit_block(
            self.base_dir, blocked["id"], blocked["sha256"],
            "As fontes foram revisadas manualmente pela equipe editorial.",
        )
        self.assertEqual(released["status"], "pending_review")
        self.assertTrue(editorial.review_payload(released)["can_approve"])
        self.assertIn("manualmente", released["manual_override"]["reason"])
        approved = editorial.approve_draft(self.base_dir, released["id"], released["sha256"])
        self.assertEqual(approved["status"], "approved")

    def test_failed_audit_can_be_retried_without_regenerating_content(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        blocked = editorial.audit_draft(
            self.base_dir, draft["id"], None,
            audit_fn=lambda _payload: {"issues": [], "claims": [], "appraisals": []},
        )
        restarted = editorial.restart_audit(self.base_dir, blocked["id"], blocked["sha256"])
        self.assertEqual(restarted["status"], "auditing")
        self.assertEqual(restarted["content"]["title"], blocked["content"]["title"])
        self.assertIn("Nova conferência", restarted["audit"]["issues"][0]["reason"])

    def test_incomplete_audit_cannot_skip_claims_or_source_appraisal(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        blocked = editorial.audit_draft(
            self.base_dir,
            draft["id"],
            None,
            audit_fn=lambda _payload: {"issues": [], "claims": [], "appraisals": []},
        )
        reasons = " ".join(issue["reason"] for issue in blocked["audit"]["issues"])
        self.assertEqual(blocked["status"], "audit_error")
        self.assertIn("nenhuma afirmação", reasons)
        self.assertIn("não foram conferidos", reasons)

    def test_dangerous_incompatibility_remains_a_hard_block(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        audit = passing_audit({})
        audit["issues"] = [{
            "severity": "hard_block",
            "location": "Orientação de segurança",
            "reason": "A recomendação contradiz a fonte e pode causar dano.",
        }]
        blocked = editorial.audit_draft(
            self.base_dir, draft["id"], None, audit_fn=lambda _payload: audit,
        )
        self.assertEqual(blocked["status"], "blocked")
        self.assertFalse(editorial.review_payload(blocked)["can_approve"])

    def test_auto_repair_is_opt_in_and_limited_to_one_attempt_per_version(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        audit = passing_audit({})
        audit["issues"] = [{
            "severity": "needs_revision",
            "location": "Abertura",
            "reason": "Reduzir a certeza da afirmação.",
        }]
        revision = editorial.audit_draft(
            self.base_dir, draft["id"], None, audit_fn=lambda _payload: audit,
        )
        repairing = editorial.start_auto_repair(
            self.base_dir, revision["id"], revision["sha256"],
        )
        self.assertEqual(repairing["status"], "repairing")
        self.assertEqual(repairing["repair_attempts"], 1)
        repairing["status"] = "needs_revision"
        editorial.save_draft(self.base_dir, repairing)
        with self.assertRaises(editorial.RevaMaisEditorialError):
            editorial.start_auto_repair(
                self.base_dir, repairing["id"], repairing["sha256"],
            )

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

    def test_email_preview_constrains_every_article_image(self):
        draft = editorial.create_draft(self.base_dir, sample_result(), source_task_id="task-1")
        email_content = draft["content"]["html_full"].split(editorial.CONTENT_START, 1)[1].split(
            editorial.CONTENT_END, 1
        )[0]
        images = re.findall(r"<img\b[^>]*>", email_content, flags=re.I)
        self.assertEqual(len(images), 2)
        for image in images:
            self.assertIn('class="body-img"', image)
            self.assertIn('width="520"', image)
            self.assertIn("width:100%", image)
            self.assertIn("max-width:520px", image)
            self.assertIn("height:auto", image)

    def test_mailchimp_normalization_repairs_legacy_content_without_touching_header(self):
        legacy = (
            '<html><body><img class="header-img" src="https://example.com/header.png">'
            f'{editorial.CONTENT_START}<p><img src="https://example.com/legacy.png"></p>'
            f'{editorial.CONTENT_END}</body></html>'
        )
        normalized = editorial.email_safe_document(legacy)
        header = re.search(r'<img[^>]+header\.png[^>]*>', normalized).group(0)
        article = re.search(r'<img[^>]+legacy\.png[^>]*>', normalized).group(0)
        self.assertNotIn('width="520"', header)
        self.assertIn('width="520"', article)
        self.assertIn("max-width:520px", article)

        draft = editorial.create_draft(self.base_dir, sample_result())
        draft["content"]["html_full"] = legacy
        preview = editorial.review_payload(draft)["content"]["html_full"]
        self.assertIn("max-width:520px", preview)

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
        self.assertIn("science-v2.png", updated["content"]["html_full"])
        science_email_image = re.search(
            r'<img[^>]+src="https://example.com/science-v2\.png"[^>]*>',
            updated["content"]["html_full"],
        ).group(0)
        self.assertIn('width="520"', science_email_image)
        self.assertIn("max-width:520px", science_email_image)

    def test_schedule_edit_creates_approved_child_and_updates_visible_date(self):
        source = sample_result()
        source["data_publicacao"] = "20/09/2026"
        source["data_iso"] = "2026-09-20T12:00:00+00:00"
        source["html_content"] += "<p>Edição de 20/09/2026.</p>"
        source["html_full"] = source["html_full"].replace(
            "<!-- REVAMAIS_CONTENT_END -->", "<p>Edição de 20/09/2026.</p><!-- REVAMAIS_CONTENT_END -->"
        )
        draft = editorial.create_draft(self.base_dir, source)
        audited = editorial.audit_draft(self.base_dir, draft["id"], None, audit_fn=passing_audit)
        approved = editorial.approve_draft(self.base_dir, audited["id"], audited["sha256"])
        updated = editorial.save_schedule_draft(
            self.base_dir, approved["id"], approved["sha256"], "2026-10-01T09:30:00-03:00"
        )
        self.assertEqual(updated["status"], "approved")
        self.assertEqual(updated["metadata"]["data_publicacao"], "01/10/2026")
        self.assertIn("01/10/2026", updated["content"]["html_content"])
        self.assertNotIn("20/09/2026", updated["content"]["html_content"])

    def test_schedule_edit_keeps_documented_manual_override(self):
        draft = editorial.create_draft(self.base_dir, sample_result())
        blocked = editorial.audit_draft(
            self.base_dir, draft["id"], None,
            audit_fn=lambda _payload: {"issues": [], "claims": [], "appraisals": []},
        )
        released = editorial.override_audit_block(
            self.base_dir, blocked["id"], blocked["sha256"], "A equipe conferiu as fontes manualmente antes de publicar.",
        )
        updated = editorial.save_schedule_draft(
            self.base_dir, released["id"], released["sha256"], "2026-10-01T09:30:00-03:00"
        )
        self.assertEqual(updated["status"], "pending_review")
        self.assertTrue(updated.get("manual_override"))
        self.assertTrue(editorial.review_payload(updated)["can_approve"])

    def test_restores_review_payload_after_local_draft_is_lost(self):
        original = self.create_audited()
        payload = editorial.review_payload(original)
        path = editorial.draft_path(self.base_dir, original["id"])
        path.unlink()
        restored = editorial.restore_draft_snapshot(self.base_dir, payload)
        self.assertEqual(restored["id"], original["id"])
        self.assertEqual(restored["content"]["title"], original["content"]["title"])
        self.assertTrue(editorial.draft_path(self.base_dir, original["id"]).exists())

    def test_list_drafts_returns_only_latest_revision_for_each_edition(self):
        previous = self.create_audited()
        editorial.save_schedule_draft(
            self.base_dir, previous["id"], previous["sha256"], "2026-10-01T09:30:00-03:00"
        )
        rows = editorial.list_drafts(self.base_dir)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["metadata"]["data_publicacao"], "01/10/2026")

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
