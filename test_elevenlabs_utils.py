import unittest

from elevenlabs_utils import diagnosticar_erro_elevenlabs


class FakeElevenLabsError(Exception):
    def __init__(self, detail, status_code):
        super().__init__(detail.get("message", "erro"))
        self.body = {"detail": detail}
        self.status_code = status_code


class ElevenLabsUtilsTest(unittest.TestCase):
    def test_identifica_plano_incompativel_com_clone(self):
        error = FakeElevenLabsError(
            {
                "code": "subscription_required",
                "status": "ivc_not_permitted",
                "message": "Instantly cloned voices are not available on your current plan.",
                "request_id": "req-123",
            },
            401,
        )

        result = diagnosticar_erro_elevenlabs(error)

        self.assertEqual(result["code"], "elevenlabs_plan_required")
        self.assertEqual(result["request_id"], "req-123")
        self.assertIn("Creator", result["message"])

    def test_identifica_cota_esgotada(self):
        error = FakeElevenLabsError(
            {"code": "quota_exceeded", "message": "Quota exceeded"},
            401,
        )

        result = diagnosticar_erro_elevenlabs(error)

        self.assertEqual(result["code"], "elevenlabs_quota_exceeded")

    def test_identifica_plano_em_erro_http_textual(self):
        error = RuntimeError(
            "Eleven v3 dialogue falhou (403): "
            "Professional voices require a Creator tier subscription or above."
        )

        result = diagnosticar_erro_elevenlabs(error)

        self.assertEqual(result["code"], "elevenlabs_plan_required")


if __name__ == "__main__":
    unittest.main()
