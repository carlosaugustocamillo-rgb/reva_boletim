"""Thematic search, canonical metadata and newsletter isolation, without paid APIs."""
import json
import ast
import tempfile
import unittest
import uuid
from pathlib import Path
from unittest.mock import patch, MagicMock

import podcast_editorial as editorial
import podcast_topic as topic
from test_podcast_editorial import fake_client
from test_podcast_pubmed_integration import run_pipeline
from test_pubmed_related import ANCHOR


class TopicPodcastTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.options = {'modo_podcast': 'tema', 'tema_podcast': 'Exercício e fadiga',
                        'orientacoes_podcast': 'Discutir limitações', 'somente_curadoria': True,
                        'resumos': True, 'mailchimp': True, 'firebase': True, 'audio': True}

    def test_api_validates_theme_before_queuing_and_forces_podcast_only_options(self):
        from fastapi import BackgroundTasks
        from fastapi.responses import JSONResponse
        tree = ast.parse(Path(__file__).with_name('api.py').read_text())
        route = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'iniciar_boletim')
        route.decorator_list = []
        save = MagicMock()
        ns = {'BackgroundTasks': BackgroundTasks, 'JSONResponse': JSONResponse,
              'uuid': uuid, 'save_task': save, 'processar_boletim_background': MagicMock()}
        exec(compile(ast.Module(body=[route], type_ignores=[]), 'api.py', 'exec'), ns)
        tasks = MagicMock()
        rejected = ns['iniciar_boletim'](tasks, payload={**self.options, 'tema_podcast': ''})
        self.assertEqual(rejected.status_code, 400)
        tasks.add_task.assert_not_called()
        save.assert_not_called()
        ns['iniciar_boletim'](tasks, payload=self.options)
        options = tasks.add_task.call_args.args[2]
        self.assertEqual(options['tema_podcast'], self.options['tema_podcast'])
        self.assertTrue(options['somente_curadoria'])
        for key in ('resumos', 'mailchimp', 'firebase', 'audio'):
            self.assertFalse(options[key])

    def prepare(self):
        with patch.object(topic, 'search_articles', return_value=('exercise AND fatigue', [ANCHOR])):
            result = list(topic.run_topic(self.options, self.root, fake_client()))[-1]
        return {**self.options, 'somente_curadoria': False,
                'curadoria_tema_id': result['curadoria_tema_id'], 'artigos_podcast_aprovados': ['100'],
                'contexto_pubmed_manual': {'estudos': [{'pmid_ancora': '100', 'referencias': []}]}}

    def test_request_validation_and_hard_isolation(self):
        validated = topic.validate_request(self.options)
        for key in ('resumos', 'mailchimp', 'audio', 'firebase'):
            self.assertFalse(validated[key])
        for body in ({**self.options, 'tema_podcast': ' '},
                     {**self.options, 'orientacoes_podcast': 'a' * 3001},
                     {**self.options, 'modo_podcast': 'unknown'},
                     {**self.options, 'somente_curadoria': False, 'curadoria_tema_id': '../x'}):
            with self.assertRaises(ValueError):
                topic.validate_request(body)

    def test_empty_search_and_stale_selection_do_not_generate(self):
        with patch.object(topic, 'search_articles', return_value=('empty', [])):
            result = list(topic.run_topic(self.options, self.root, fake_client()))[-1]
        self.assertEqual(result['artigos_sugeridos'], [])
        options = self.prepare()
        for override in ({'tema_podcast': 'Outro tema'}, {'artigos_podcast_aprovados': ['999']},
                         {'artigos_podcast_aprovados': ['100', '100']},
                         {'contexto_pubmed_manual': {'estudos': [{'pmid_ancora': '999'}]}}):
            client = fake_client()
            with self.assertRaises(ValueError):
                list(topic.run_topic({**options, **override}, self.root, client))
            client.with_options.assert_not_called()

    def test_theme_reaches_planning_writing_and_survives_edit_approval(self):
        client = fake_client()
        result = list(topic.run_topic(self.prepare(), self.root, client))[-1]
        public = result['roteiro_editorial']
        self.assertEqual(public['status'], 'pending_review')
        self.assertEqual(result['mailchimp']['status'], 'skipped')
        for call in client.with_options.return_value.chat.completions.create.call_args_list[:2]:
            self.assertIn('Exercício e fadiga', json.dumps(call.kwargs, ensure_ascii=False))
        draft = editorial.load_draft(self.root, public['id'])
        old_hash = draft['sha256']
        draft['episode_context']['tema'] = 'Alterado'
        self.assertNotEqual(editorial.fingerprint(draft), old_hash)
        edited = editorial.save_edited_draft(self.root, public['id'], public['sha256'], public['text'])
        self.assertEqual(edited['episode_context']['tema'], self.options['tema_podcast'])
        audited = editorial.audit_edited_draft(self.root, edited['id'], fake_client())
        approved = editorial.approve_draft(self.root, audited['id'], audited['sha256'])
        self.assertEqual(editorial.approved_audio(self.root, approved['id'], approved['sha256'])['episode_context'],
                         edited['episode_context'])

    def test_real_pipeline_thematic_entrypoint_never_searches_newsletter_or_sends_email(self):
        with patch.object(topic, 'search_articles', return_value=('query', [ANCHOR])):
            run = run_pipeline(self.root, options_override=self.options)
        self.assertEqual(run['result']['tipo'], 'curadoria_podcast')
        run['search'].assert_not_called()
        run['translate'].assert_not_called()
        run['mailchimp'].campaigns.create.assert_not_called()
        self.assertFalse(list(self.root.glob('boletim*')))

    def test_approved_theme_recovers_isolation_without_client_mode(self):
        result = list(topic.run_topic(self.prepare(), self.root, fake_client()))[-1]
        draft = result['roteiro_editorial']
        approved = editorial.approve_draft(self.root, draft['id'], draft['sha256'])
        firebase = MagicMock()
        firebase.is_firebase_ready.return_value = False  # Prevent paid TTS in this integration test.
        with patch.dict('sys.modules', {'firebase_service': firebase}):
            run = run_pipeline(self.root, options_override={
                'roteiro_aprovado_id': approved['id'], 'roteiro_aprovado_sha256': approved['sha256'],
                'resumos': True, 'mailchimp': True, 'firebase': True, 'audio': True})
        run['search'].assert_not_called()
        run['mailchimp'].campaigns.create.assert_not_called()
        self.assertIn(approved['id'], run['result']['episodio_path'])
        firebase.upload_file.assert_not_called()


    def test_audio_publishes_exact_approved_dialogue_with_unique_files_and_no_campaign(self):
        result = list(topic.run_topic(self.prepare(), self.root, fake_client()))[-1]
        public = result['roteiro_editorial']
        draft = editorial.approve_draft(self.root, public['id'], public['sha256'])
        narrated = []
        def synthesize(dialogue, path, **kwargs):
            self.assertTrue(kwargs['preservar_texto'])
            narrated.append(dialogue)
            Path(path).write_bytes(b'synthetic-audio')
            return path, 1
        class Audio:
            def __len__(self): return 1000
            def __add__(self, other): return self
            @staticmethod
            def silent(**kwargs): return Audio()
            @staticmethod
            def from_file(*args, **kwargs): return Audio()
            def export(self, path, **kwargs): Path(path).write_bytes(b'synthetic-mix')
        firebase, whatsapp = MagicMock(), MagicMock()
        firebase.is_firebase_ready.return_value = True
        firebase.upload_file.return_value = 'https://example.org/audio.mp3'
        firebase.update_podcast_feed.return_value = 'https://example.org/rss.xml'
        overrides = {'uuid': uuid, 'DATA_DIR': str(self.root), 'AudioSegment': Audio,
                     'INTRO_PATH': str(self.root / 'intro.mp3'), 'elevenlabs_client': object(),
                     'ELEVEN_AUDIO_MODEL': 'test', 'ELEVEN_AUDIO_DIALOGUE_MODEL': 'test',
                     'ELEVEN_AUDIO_DIALOGUE_SEED': None, 'usar_dialogo_eleven': lambda: True,
                     'gerar_dialogo_com_eleven': synthesize, 'formatar_erro_elevenlabs': str}
        with patch.dict('sys.modules', {'firebase_service': firebase, 'whatsapp_service': whatsapp}):
            run = run_pipeline(self.root, namespace_override=overrides, options_override={
                'roteiro_aprovado_id': draft['id'], 'roteiro_aprovado_sha256': draft['sha256'],
                'resumos': True, 'mailchimp': True, 'firebase': True, 'audio': True})
        self.assertEqual(narrated, editorial.dialogues(draft))
        self.assertIsNone(run['result']['audio_erro'])
        self.assertIn(draft['id'], run['result']['audio_download_url'])
        self.assertIn(draft['id'], run['result']['brief_spotify_download_url'])
        firebase.update_podcast_feed.assert_called_once()
        self.assertEqual(firebase.update_podcast_feed.call_args.kwargs['episodio_titulo'], draft['script']['title'])
        run['mailchimp'].campaigns.create.assert_not_called()
        whatsapp.create_draft.assert_not_called()


if __name__ == '__main__':
    unittest.main()
