import ast
import copy
import json
import os
import re
import tempfile
import unittest
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import JSONResponse, FileResponse, Response
from fastapi.testclient import TestClient

import podcast_editorial
from test_podcast_editorial import fake_client
from test_pubmed_related import ANCHOR


class EditorialApiTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        app = FastAPI()
        wanted = {'get_ultimo_roteiro', 'get_podcast_draft', 'approve_podcast_draft',
                  'download_podcast_script', 'download_podcast_preview'}
        tree = ast.parse(Path(__file__).with_name('api.py').read_text())
        routes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name in wanted]
        ns = {'app': app, 'os': os, 're': re, 'json': json, '__file__': str(self.root / 'api.py'),
              'JSONResponse': JSONResponse, 'FileResponse': FileResponse, 'Response': Response,
              'podcast_editorial': podcast_editorial}
        exec(compile(ast.Module(body=routes, type_ignores=[]), 'api.py', 'exec'), ns)
        self.client = TestClient(app)
        self.addCleanup(self.client.close)
        self.draft = podcast_editorial.generate_episode([ANCHOR], None, fake_client())
        podcast_editorial.save_draft(self.root / 'data', self.draft)

    def test_latest_has_readable_text_sources_and_canonical_opening(self):
        data = self.client.get('/ultimo-roteiro').json()
        self.assertEqual(data['editorial']['sha256'], self.draft['sha256'])
        self.assertEqual(data['conteudo'], podcast_editorial.dialogues(self.draft))
        self.assertEqual(data['editorial']['evidence'][0]['sources'][0]['pmid'], '100')
        url = data['editorial']['download_url']
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.text, data['editorial']['text'])

    def test_approval_http_contract_stale_hash_and_missing_file(self):
        url = f'/podcast-roteiro/{self.draft["id"]}/aprovar'
        self.assertEqual(self.client.post(url, json={'sha256': 'stale'}).status_code, 409)
        approved = self.client.post(url, json={'sha256': self.draft['sha256']})
        self.assertEqual(approved.status_code, 200)
        self.assertEqual(approved.json()['status'], 'approved')
        self.assertEqual(self.client.get('/podcast-roteiro/' + '0' * 32).status_code, 404)
        self.assertEqual(self.client.get('/podcast-roteiro/not-valid').status_code, 400)

    def test_audio_download_only_allows_known_episode_names(self):
        self.assertEqual(self.client.get('/baixar-audio-podcast/.env').status_code, 400)
        self.assertEqual(self.client.get('/baixar-audio-podcast/episodio_boletim_2026-09-12.mp3').status_code, 404)
        directory = self.root / 'data' / 'audios'
        directory.mkdir()
        (directory / 'episodio_boletim_2026-09-12.mp3').write_bytes(b'test-only')
        response = self.client.get('/baixar-audio-podcast/episodio_boletim_2026-09-12.mp3')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content, b'test-only')

    def test_blocked_draft_is_readable_but_cannot_be_approved(self):
        self.draft['status'] = 'blocked'
        self.draft['audit']['issues'] = [{'severity': 'blocking', 'location': 'PMID 100', 'reason': 'Ressalva divergente'}]
        self.draft['sha256'] = podcast_editorial.fingerprint(self.draft)
        podcast_editorial.save_draft(self.root / 'data', self.draft)
        result = self.client.get('/ultimo-roteiro').json()['editorial']
        self.assertFalse(result['can_approve'])
        self.assertTrue(result['text'])
        self.assertEqual(self.client.post(f'/podcast-roteiro/{self.draft["id"]}/aprovar',
                         json={'sha256': self.draft['sha256']}).status_code, 409)

    def test_background_completion_reports_partial_or_review_instead_of_success(self):
        tree = ast.parse(Path(__file__).with_name('api.py').read_text())
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'processar_boletim_background')
        for result, expected in (({'roteiro_erro': 'Failed'}, 'partial'),
                                 ({'audio_erro': 'Failed'}, 'partial'),
                                 ({'roteiro_editorial': {'status': 'pending_review'}}, 'pending_review'),
                                 ({'mailchimp': {'status': 'scheduled'}}, 'success')):
            saved = []
            ns = {'rodar_boletim': lambda options: iter([result]), 'load_task': lambda task: {},
                  'save_task': lambda task, value: saved.append(copy.deepcopy(value))}
            exec(compile(ast.Module(body=[node], type_ignores=[]), 'api.py', 'exec'), ns)
            ns['processar_boletim_background']('test', {})
            self.assertEqual(saved[-1]['status'], 'completed')
            self.assertEqual(saved[-1]['outcome'], expected)
            if expected != 'success':
                self.assertNotIn('sucesso', saved[-1]['logs'][-1])


if __name__ == '__main__':
    unittest.main()
