"""Exercise the real pipeline functions with external services replaced.

AST loading avoids the legacy module's credential-dependent clients and import
side effects. The executed functions are compiled directly from the source file.
"""
import ast
import asyncio
import copy
import io
import json
import os
import re
import tempfile
import subprocess
import sys
import unittest
from contextlib import redirect_stdout
from datetime import datetime, timedelta
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytz

import pubmed_related
import podcast_editorial
from test_podcast_editorial import editorial_response
from test_pubmed_related import ANCHOR, candidate, decision, llm


ROOT = Path(__file__).resolve().parent


class FixedDatetime(datetime):
    @classmethod
    def today(cls):
        return cls(2026, 9, 4, 18)

    @classmethod
    def now(cls, tz=None):
        result = cls.today()
        return tz.localize(result) if tz is not None else result


def load_functions(filename, namespace, source=None):
    tree = ast.parse(source if source is not None else (ROOT / filename).read_text(encoding="utf-8"))
    functions = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))]
    for node in functions:
        node.decorator_list = []
    exec(compile(ast.Module(body=functions, type_ignores=[]), filename, "exec"), namespace)
    return namespace


def run_pipeline(base_dir, *, related=True, script=True, source=None, options_override=None):
    anchor = {**ANCHOR, "journal": "Clinical Exercise Journal", "ano": "2026", "volume": "1", "issue": "2", "paginas": "10-20"}
    client = llm([decision()])
    screening_response = client.with_options.return_value.chat.completions.create.return_value
    def respond(**kwargs):
        if kwargs.get('response_format', {}).get('json_schema', {}).get('name', '').startswith('podcast_'):
            return editorial_response(**kwargs)
        return screening_response
    client.with_options.return_value.chat.completions.create.side_effect = respond
    dialogue = [{"speaker": "HOST", "text": "Estudo de Silva A. Contexto: Jones B, 2025."}]
    client.chat.completions.create.return_value = SimpleNamespace(choices=[SimpleNamespace(message=SimpleNamespace(content=json.dumps(dialogue)))])
    mailchimp = MagicMock()
    mailchimp.campaigns.create.return_value = {"id": "test-campaign"}
    namespace = {"__file__": str(ROOT / "boletim_service.py"), "os": os, "json": json, "podcast_editorial": podcast_editorial,
                 "re": re, "copy": copy, "datetime": FixedDatetime, "timedelta": timedelta,
                 "pytz": pytz, "BASE_DIR": str(base_dir), "AUDIO_DIR": str(Path(base_dir) / "audios"),
                 "client": client, "mc": mailchimp, "MC_LIST_ID": "test-list", "MC_FROM_NAME": "Test",
                 "MC_REPLY_TO": "test@example.org", "TEMPLATE_HTML_BASE": "<html>{conteudo_aqui}</html>",
                 "CONSULTAS_PRINCIPAIS": ["Exercício: exercise"], "CONSULTAS_DETALHADAS": {"Exercício": "exercise"},
                 **{name: getattr(pubmed_related, name) for name in ("enabled_from_env", "enrich_episode", "context_for_script", "references_for_notes", "save_json")}}
    load_functions("boletim_service.py", namespace, source=source)
    search = MagicMock(return_value=["100"])
    translate = MagicMock(return_value="Tradução literal, sem alterações.\nSegundo parágrafo.")
    namespace.update({"buscar_ids": search, "buscar_info_estruturada": lambda _: [copy.deepcopy(anchor)],
                      "traduzir_resumo": translate,
                      "gerar_brief_spotify": lambda **_: "Descrição científica do episódio."})
    pubmed = MagicMock()
    pubmed.candidates.return_value = ([candidate()], False)
    options = {"resumos": True, "roteiro": script, "revisao_roteiro": False,
               "brief_spotify": True, "audio": False, "mailchimp": True, "firebase": False,
               "referencias_pubmed": related}
    options.update(options_override or {})
    with patch("pubmed_related.PubMedRelatedClient", return_value=pubmed), patch.dict(os.environ, {"PODCAST_PUBMED_RELATED_ENABLED": "true"}), redirect_stdout(io.StringIO()):
        messages = list(namespace["rodar_boletim"](options))
    return {"result": messages[-1], "messages": messages, "pubmed": pubmed, "client": client,
            "mailchimp": mailchimp, "search": search, "translate": translate, "namespace": namespace}


class PodcastPubmedIntegrationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_email_and_detailed_summaries_identical_with_context_on_or_off(self):
        enabled = run_pipeline(self.root / "on")
        disabled = run_pipeline(self.root / "off", related=False)
        for filename in ("boletim_pubmed_2026-09-04.txt", "boletim_detalhado_2026-09-04.txt", "boletim_para_revisao_2026-09-04.txt"):
            self.assertEqual((self.root / "on" / filename).read_bytes(), (self.root / "off" / filename).read_bytes())
        for method in ("create", "set_content", "schedule"):
            self.assertEqual(getattr(enabled["mailchimp"].campaigns, method).call_args_list,
                             getattr(disabled["mailchimp"].campaigns, method).call_args_list)
        self.assertEqual(enabled["search"].call_args_list, disabled["search"].call_args_list)
        self.assertEqual(enabled["translate"].call_args_list, disabled["translate"].call_args_list)
        enabled["pubmed"].candidates.assert_called_once_with("100")
        disabled["pubmed"].candidates.assert_not_called()

    def test_email_only_does_not_call_related_or_llm_screening(self):
        run = run_pipeline(self.root, script=False)
        run["pubmed"].candidates.assert_not_called()
        run["client"].with_options.assert_not_called()
        self.assertIsNone(run["result"]["referencias_pubmed"])
        self.assertFalse((self.root / "referencias").exists())

    def test_email_is_byte_identical_to_committed_pipeline(self):
        baseline = subprocess.run(['git', 'show', 'HEAD:boletim_service.py'], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout
        before = run_pipeline(self.root / 'before', source=baseline)
        after = run_pipeline(self.root / 'after')
        for filename in ('boletim_pubmed_2026-09-04.txt', 'boletim_detalhado_2026-09-04.txt', 'boletim_para_revisao_2026-09-04.txt'):
            self.assertEqual((self.root / 'before' / filename).read_bytes(), (self.root / 'after' / filename).read_bytes())
        for method in ('create', 'set_content', 'schedule'):
            self.assertEqual(getattr(before['mailchimp'].campaigns, method).call_args_list,
                             getattr(after['mailchimp'].campaigns, method).call_args_list)

    def test_editorial_failure_preserves_email_and_never_calls_tts(self):
        with patch.object(podcast_editorial, 'generate_episode_steps', side_effect=ValueError('Source validation failed')):
            run = run_pipeline(self.root)
        self.assertIn('Source validation failed', run['result']['roteiro_erro'])
        self.assertIsNone(run['result']['roteiro_editorial'])
        run['mailchimp'].campaigns.schedule.assert_called_once()

    def test_pending_or_blocked_script_never_uploads_old_segments_or_creates_whatsapp(self):
        from test_podcast_editorial import fake_client
        for blocked in (False, True):
            with self.subTest(blocked=blocked):
                root = self.root / str(blocked)
                (root / 'audios').mkdir(parents=True)
                (root / 'audios' / 'estudo1_completo.mp3').write_bytes(b'old unrelated audio')
                draft = podcast_editorial.generate_episode([ANCHOR], None, fake_client())
                if blocked:
                    draft['status'] = 'blocked'
                    draft['audit']['issues'] = [{'severity': 'blocking', 'location': '100', 'reason': 'Unresolved caution'}]
                    draft['sha256'] = podcast_editorial.fingerprint(draft)
                def generate(*args, **kwargs):
                    yield 'Test generation'
                    return draft
                firebase = MagicMock()
                whatsapp = MagicMock()
                with patch.object(podcast_editorial, 'generate_episode_steps', generate), patch.dict(
                    sys.modules, {'firebase_service': firebase, 'whatsapp_service': whatsapp}
                ):
                    run = run_pipeline(root, options_override={'audio': True, 'firebase': True})
                firebase.upload_file.assert_not_called()
                firebase.update_podcast_feed.assert_not_called()
                whatsapp.create_draft.assert_not_called()
                run['mailchimp'].campaigns.schedule.assert_called_once()
                self.assertIsNotNone(run['result']['roteiro_editorial'])
                self.assertEqual(run['result']['roteiro_editorial']['can_approve'], not blocked)
                loaded = podcast_editorial.load_draft(root, draft['id'])
                self.assertTrue(podcast_editorial.transcript(loaded))
                self.assertFalse(any('Tentando resgatar' in str(m) for m in run['messages']))
                if blocked:
                    self.assertIsNotNone(run['result']['roteiro_erro'])
                    self.assertFalse((root / 'roteiros').exists())

    def test_audio_requires_approval_and_narrates_exact_saved_text_without_new_llm_calls(self):
        import uuid
        from test_podcast_voice_settings import load_voice_config
        run = run_pipeline(self.root)
        ns = run['namespace']
        review = run['result']['roteiro_editorial']
        ns.update(load_voice_config())
        ns.update({'uuid': uuid, 'DATA_DIR': str(self.root), 'INTRO_PATH': str(self.root / 'no-intro.mp3'),
                   'formatar_erro_elevenlabs': str})
        segment = MagicMock()
        segment.__add__.return_value = segment
        segment.__iadd__.return_value = segment
        segment.export.side_effect = lambda path, **kw: Path(path).write_bytes(b'test-mp3')
        audio = MagicMock()
        audio.empty.return_value = audio.silent.return_value = audio.from_file.return_value = segment
        ns['AudioSegment'] = audio
        synth = ns['elevenlabs_client'].text_to_speech.convert
        synth.return_value = [b'test-mp3']
        options = {'resumos': False, 'roteiro': False, 'audio': True, 'mailchimp': False, 'firebase': False,
                   'roteiro_aprovado_id': review['id'], 'roteiro_aprovado_sha256': review['sha256']}
        with redirect_stdout(io.StringIO()):
            rejected = list(ns['rodar_boletim'](options))[-1]
        self.assertIsNotNone(rejected['audio_erro'])
        synth.assert_not_called()
        approved = podcast_editorial.approve_draft(self.root, review['id'], review['sha256'])
        run['client'].reset_mock()
        run['mailchimp'].reset_mock()
        with redirect_stdout(io.StringIO()):
            result = list(ns['rodar_boletim'](options))[-1]
        expected = [turn['text'] for block in podcast_editorial.dialogues(approved) for turn in block]
        self.assertEqual([call.kwargs['text'] for call in synth.call_args_list], expected)
        self.assertIsNotNone(result['audio_download_url'])
        self.assertIsNone(result['audio_erro'])
        run['client'].with_options.assert_not_called()
        run['mailchimp'].campaigns.create.assert_not_called()
        self.assertFalse(options['resumos'])

        # Failure in one block must not create a second, incomplete final episode.
        synth.reset_mock()
        synth.side_effect = [RuntimeError('Synthetic TTS failure')] + [[b'test-mp3']] * 8
        with redirect_stdout(io.StringIO()):
            partial = list(ns['rodar_boletim'](options))[-1]
        self.assertIsNotNone(partial['audio_erro'])
        self.assertIsNone(partial['audio_download_url'])

    def test_new_draft_does_not_run_legacy_transition_rewriters(self):
        source = (ROOT / 'boletim_service.py').read_text()
        tree = ast.parse(source)
        pipeline = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'rodar_boletim')
        called = {n.func.id for n in ast.walk(pipeline) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
        self.assertTrue({'resumo_para_podcast', 'revisar_roteiro_completo_para_podcast',
                         'suavizar_frases_mecanicas', 'forcar_transicoes_ancoradas_no_proximo_titulo'}.isdisjoint(called))

    def test_metadata_reaches_screening_and_only_script_gets_context(self):
        run = run_pipeline(self.root)
        calls = run["client"].with_options.return_value.chat.completions.create.call_args_list
        screening = next(call.kwargs for call in calls if not call.kwargs.get('response_format', {}).get('json_schema', {}).get('name', '').startswith('podcast_'))
        source = json.loads(screening["messages"][1]["content"])
        self.assertEqual(source["ancora"]["pmid"], "100")
        self.assertEqual(source["ancora"]["doi"], "10.123/new")
        writing = next(call.kwargs for call in calls if call.kwargs.get('response_format', {}).get('json_schema', {}).get('name') == 'podcast_write')
        script_prompt = writing["messages"][1]["content"]
        self.assertIn("Silva A", script_prompt)
        self.assertIn("Earlier exercise trial 200", script_prompt)
        campaign_html = run["mailchimp"].campaigns.set_content.call_args.args[1]["html"]
        self.assertNotIn("Earlier exercise trial 200", campaign_html)
        notes = (self.root / "brief_spotify_2026-09-04.txt").read_text()
        self.assertIn("https://pubmed.ncbi.nlm.nih.gov/200/", notes)
        audit = json.loads(Path(run["result"]["referencias_pubmed_path"]).read_text())
        self.assertEqual(audit["estudos"][0]["referencias"][0]["pmid"], "200")
        self.assertEqual(run["result"]["referencias_pubmed_download_url"], "/baixar-referencias-podcast/2026-09-04")

    def test_disabled_run_does_not_advertise_old_reference_report(self):
        run_pipeline(self.root)
        disabled = run_pipeline(self.root, related=False)
        self.assertIsNone(disabled["result"]["referencias_pubmed_path"])
        self.assertIsNone(disabled["result"]["referencias_pubmed_download_url"])
        self.assertNotIn("Referências do episódio", (self.root / "brief_spotify_2026-09-04.txt").read_text())

    def test_route_forwards_context_setting(self):
        import uuid
        namespace = {"uuid": uuid, "save_task": MagicMock(), "processar_boletim_background": MagicMock(),
                     "BackgroundTasks": object, "Request": object, "Response": object,
                     "CampanhaInput": object, "NewsPayload": object, "RevaMaisInput": object,
                     "RevaMaisConsensusInput": object, "RevaMaisConsensusPdfInput": object,
                     "AgentInput": object}
        # Compile only the route under test to avoid unrelated endpoint annotations.
        source = ast.parse((ROOT / "api.py").read_text())
        route = next(node for node in source.body if isinstance(node, ast.FunctionDef) and node.name == "iniciar_boletim")
        route.decorator_list = []
        exec(compile(ast.Module(body=[route], type_ignores=[]), "api.py", "exec"), namespace)
        tasks = MagicMock()
        namespace["iniciar_boletim"](tasks, referencias_pubmed=False)
        options = tasks.add_task.call_args.args[2]
        self.assertFalse(options["referencias_pubmed"])
        self.assertTrue(options["mailchimp"])

    def test_reference_download_validates_date_and_serves_only_report(self):
        tree = ast.parse((ROOT / "api.py").read_text())
        route = next(node for node in tree.body if isinstance(node, ast.AsyncFunctionDef) and node.name == "baixar_referencias_podcast")
        route.decorator_list = []
        namespace = {"os": os, "re": re, "datetime": datetime, "__file__": str(self.root / "api.py"),
                     "JSONResponse": lambda **values: SimpleNamespace(**values),
                     "FileResponse": lambda **values: SimpleNamespace(**values)}
        exec(compile(ast.Module(body=[route], type_ignores=[]), "api.py", "exec"), namespace)
        download = namespace["baixar_referencias_podcast"]
        for value in ("../../.env", "2026-02-31", "2026-09-04\n", "2026-09-04/../.env"):
            self.assertEqual(asyncio.run(download(value)).status_code, 400)
        self.assertEqual(asyncio.run(download("2026-09-04")).status_code, 404)
        report = self.root / "data" / "referencias" / "contexto_pubmed_2026-09-04.json"
        pubmed_related.save_json(report, {"estudos": []})
        self.assertEqual(asyncio.run(download("2026-09-04")).path, str(report))


if __name__ == "__main__":
    unittest.main()
