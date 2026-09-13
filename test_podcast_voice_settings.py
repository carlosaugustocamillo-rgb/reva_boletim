"""Testa o payload real de TTS sem importar clientes ou gastar créditos."""

import ast
import os
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from elevenlabs import VoiceSettings


def load_voice_config():
    source = Path(__file__).with_name("boletim_service.py").read_text(encoding="utf-8")
    functions = {
        "_bool_env", "_float_env", "usar_dialogo_eleven_v3",
        "parametros_tts_podcast", "gerar_fala_com_elevenlabs",
    }
    constants = {
        "ELEVEN_VOICE_ID_HOST", "ELEVEN_VOICE_ID_COHOST", "ELEVEN_AUDIO_MODEL",
        "ELEVEN_AUDIO_DIALOGUE_ENABLED", "ELEVEN_AUDIO_FALLBACK_MODEL",
        "ELEVEN_AUDIO_LANGUAGE_CODE",
    }
    nodes = []
    for node in ast.parse(source).body:
        if isinstance(node, ast.FunctionDef) and node.name in functions:
            nodes.append(node)
        elif isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id in constants
            for target in node.targets
        ):
            nodes.append(node)
    namespace = {"os": os, "VoiceSettings": VoiceSettings, "elevenlabs_client": MagicMock()}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "boletim_service.py", "exec"), namespace)
    return namespace


class PodcastVoiceSettingsTest(unittest.TestCase):
    def setUp(self):
        self.env = patch.dict(os.environ, {}, clear=True)
        self.env.start()
        self.addCleanup(self.env.stop)

    def test_requested_voice_profiles_reach_sdk_payload(self):
        config = load_voice_config()
        for speaker, voice_id, speed in (
            ("HOST", "L0Dsvb3SLTyegXwtm47J", 1.1),
            ("COHOST", "uYXf8XasLslADfZ2MB4u", 1.1),
        ):
            with self.subTest(speaker=speaker):
                config["gerar_fala_com_elevenlabs"]("Texto em português.", speaker)
                actual = config["elevenlabs_client"].text_to_speech.convert.call_args.kwargs
                self.assertEqual(actual["text"], "Texto em português.")
                self.assertEqual(actual["voice_id"], voice_id)
                self.assertEqual(actual["model_id"], "eleven_multilingual_v2")
                self.assertEqual(actual["output_format"], "mp3_44100_128")
                self.assertEqual(actual["voice_settings"].model_dump(), {
                    "speed": speed, "stability": 1.0, "similarity_boost": 1.0,
                    "style": 0.0, "use_speaker_boost": True,
                })
                self.assertNotIn("language_code", actual)

    def test_intro_study_and_goodbye_share_voice_settings(self):
        config = load_voice_config()
        for speaker in ("HOST", "COHOST"):
            params = []
            for text in ("Olá, bem-vindos!", "O estudo descreve a intervenção.", "Até a próxima!"):
                config["gerar_fala_com_elevenlabs"](text, speaker)
                actual = config["elevenlabs_client"].text_to_speech.convert.call_args.kwargs
                params.append({**actual, "text": "", "voice_settings": actual["voice_settings"].model_dump()})
            self.assertEqual(params[0], params[1])
            self.assertEqual(params[1], params[2])

    def test_legacy_v3_flag_does_not_override_multilingual_v2(self):
        os.environ.update({"ELEVEN_AUDIO_DIALOGUE_ENABLED": "true", "ELEVEN_AUDIO_DIALOGUE_MODEL": "eleven_v3"})
        config = load_voice_config()
        self.assertFalse(config["usar_dialogo_eleven_v3"]())
        self.assertEqual(config["parametros_tts_podcast"]("HOST")["model_id"], "eleven_multilingual_v2")

    def test_v3_requires_explicit_model_and_flag(self):
        os.environ["ELEVEN_AUDIO_MODEL"] = "eleven_v3"
        self.assertFalse(load_voice_config()["usar_dialogo_eleven_v3"]())
        os.environ["ELEVEN_AUDIO_DIALOGUE_ENABLED"] = "true"
        self.assertTrue(load_voice_config()["usar_dialogo_eleven_v3"]())

    def test_environment_can_override_ids_and_individual_speeds(self):
        os.environ.update({
            "ELEVEN_VOICE_ID_HOST": "custom-host", "ELEVEN_VOICE_ID_COHOST": "custom-manu",
            "ELEVEN_VOICE_SPEED_HOST": "1.0", "ELEVEN_VOICE_SPEED_COHOST": "0.8",
        })
        config = load_voice_config()
        host = config["parametros_tts_podcast"]("HOST")
        manu = config["parametros_tts_podcast"]("COHOST")
        self.assertEqual((host["voice_id"], host["voice_settings"].speed), ("custom-host", 1.0))
        self.assertEqual((manu["voice_id"], manu["voice_settings"].speed), ("custom-manu", 0.8))

    def test_language_override_only_for_other_models(self):
        os.environ.update({"ELEVEN_AUDIO_MODEL": "eleven_turbo_v2_5", "ELEVEN_AUDIO_LANGUAGE_CODE": "pt"})
        config = load_voice_config()
        self.assertEqual(config["parametros_tts_podcast"]("HOST")["language_code"], "pt")

    def test_v3_fallback_does_not_send_unsupported_v2_language_override(self):
        os.environ.update({"ELEVEN_AUDIO_MODEL": "eleven_v3", "ELEVEN_AUDIO_DIALOGUE_ENABLED": "true"})
        params = load_voice_config()["parametros_tts_podcast"]("COHOST")
        self.assertEqual(params["model_id"], "eleven_multilingual_v2")
        self.assertNotIn("language_code", params)

    def test_invalid_values_stop_before_any_paid_request(self):
        for name, value in (
            ("ELEVEN_VOICE_SPEED_HOST", "0.1"), ("ELEVEN_VOICE_SPEED_COHOST", "1.3"),
            ("ELEVEN_VOICE_STABILITY", "100"), ("ELEVEN_VOICE_SIMILARITY_BOOST", "-1"),
            ("ELEVEN_VOICE_STYLE", "nan"), ("ELEVEN_VOICE_STYLE", "inf"),
        ):
            with self.subTest(name=name, value=value), patch.dict(os.environ, {name: value}):
                config = load_voice_config()
                speaker = "COHOST" if name.endswith("COHOST") else "HOST"
                with self.assertRaises(ValueError):
                    config["gerar_fala_com_elevenlabs"]("Teste.", speaker)
                config["elevenlabs_client"].text_to_speech.convert.assert_not_called()

    def test_unknown_speaker_is_not_silently_mapped_to_another_voice(self):
        config = load_voice_config()
        with self.assertRaises(ValueError):
            config["gerar_fala_com_elevenlabs"]("Teste.", "OUTRO")
        config["elevenlabs_client"].text_to_speech.convert.assert_not_called()


if __name__ == "__main__":
    unittest.main()
