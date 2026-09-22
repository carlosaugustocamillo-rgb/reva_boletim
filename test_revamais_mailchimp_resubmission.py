import unittest
from unittest.mock import MagicMock, patch

import revamais_service as service


def approved_draft():
    return {
        "content": {
            "title": "Exercício e saúde",
            "html_full": (
                "<html><body><!-- REVAMAIS_CONTENT_START -->"
                '<p><img src="https://example.com/image.png"></p>'
                "<!-- REVAMAIS_CONTENT_END --></body></html>"
            ),
        },
        "metadata": {
            "tema": "Exercício e saúde",
            "calendar_index": 4,
            "calendar_title": "Exercício e saúde",
        },
    }


class RevaMaisMailchimpResubmissionTest(unittest.TestCase):
    def test_email_resubmission_does_not_repeat_calendar_completion(self):
        campaigns = MagicMock()
        campaigns.create.return_value = {"id": "new-campaign"}
        previous = {
            "status": "success",
            "campaign_id": "old-campaign",
            "email_content_set": True,
            "email_scheduled": True,
            "calendar_completed": False,
        }
        with patch.object(service.mc, "campaigns", campaigns), patch.object(
            service, "marcar_tema_revamais_concluido"
        ) as complete_calendar:
            result = service.finalizar_publicacao_revamais(
                approved_draft(),
                publicar_email=True,
                criar_whatsapp=False,
                schedule_time="2026-09-27T10:30:00+00:00",
                previous_publication=previous,
                force_new_email_campaign=True,
            )

        self.assertEqual(result["status"], "success")
        self.assertEqual(result["campaign_id"], "new-campaign")
        self.assertTrue(result["email_scheduled"])
        complete_calendar.assert_not_called()

    def test_scheduled_campaign_is_success_when_schedule_response_is_lost(self):
        campaigns = MagicMock()
        campaigns.create.return_value = {"id": "new-campaign"}
        campaigns.schedule.side_effect = TimeoutError("response timeout")
        campaigns.get.return_value = {"id": "new-campaign", "status": "schedule"}
        with patch.object(service.mc, "campaigns", campaigns):
            result = service.finalizar_publicacao_revamais(
                approved_draft(),
                publicar_email=True,
                criar_whatsapp=False,
                schedule_time="2026-09-27T10:30:00+00:00",
                force_new_email_campaign=True,
            )

        self.assertEqual(result["status"], "success")
        self.assertTrue(result["email_scheduled"])
        self.assertEqual(result["mailchimp_status"], "schedule")

    def test_retry_resumes_partial_campaign_instead_of_creating_a_duplicate(self):
        campaigns = MagicMock()
        campaigns.schedule.side_effect = RuntimeError("already scheduled")
        campaigns.get.return_value = {"id": "new-campaign", "status": "schedule"}
        previous = {
            "status": "partial",
            "campaign_id": "new-campaign",
            "email_content_set": True,
            "email_scheduled": False,
            "calendar_completed": False,
            "mailchimp_submissions": [{"campaign_id": "old-campaign"}],
        }
        with patch.object(service.mc, "campaigns", campaigns):
            result = service.finalizar_publicacao_revamais(
                approved_draft(),
                publicar_email=True,
                criar_whatsapp=False,
                schedule_time="2026-09-27T10:30:00+00:00",
                previous_publication=previous,
                force_new_email_campaign=True,
            )

        self.assertEqual(result["status"], "success")
        self.assertEqual(result["campaign_id"], "new-campaign")
        campaigns.create.assert_not_called()


if __name__ == "__main__":
    unittest.main()
