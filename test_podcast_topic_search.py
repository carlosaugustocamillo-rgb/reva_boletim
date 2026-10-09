"""Independent topic routes and date-unrestricted lookup, using synthetic sources."""
import ast
import copy
import tempfile
import unittest
import uuid
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from Bio import Entrez
from fastapi import FastAPI, BackgroundTasks
from fastapi.responses import JSONResponse
from fastapi.testclient import TestClient
from pydantic import BaseModel, Field
import podcast_topic as topic

TITLE = ('The Best Exercise Modality and Dose to Reduce Cancer-Related Fatigue in Breast Cancer: '
         'A Systematic Review with Network and Dose-Response Meta-Analyses.')
ARTICLE = {'pmid': '42816300', 'titulo': TITLE, 'resumo_original': 'Synthetic test abstract.',
           'data_publicacao': '2000-01-01'}


class SearchTest(unittest.TestCase):
    def setUp(self):
        self.builder = MagicMock(return_value='breast cancer AND fatigue')
        self.fetch = MagicMock(return_value=[ARTICLE])
        modules = {'revamais_service': SimpleNamespace(gerar_query_pubmed_tema=self.builder),
                   'boletim_service': SimpleNamespace(buscar_info_estruturada=self.fetch)}
        self.modules = patch.dict('sys.modules', modules)
        self.modules.start()
        self.addCleanup(self.modules.stop)
        self.search = patch.object(Entrez, 'esearch').start()
        self.addCleanup(patch.stopall)
        self.read = patch.object(Entrez, 'read', return_value={'IdList': ['42816300']}).start()

    def assert_unrestricted(self):
        for call in self.search.call_args_list:
            self.assertEqual(set(call.kwargs), {'db', 'term', 'retmax', 'sort'})
            self.assertEqual(call.kwargs['sort'], 'relevance')

    def test_pasted_title_resolves_old_paper_without_general_llm_or_weekly_queries(self):
        query, articles = topic.search_articles(TITLE)
        self.assertEqual(articles, [ARTICLE])
        self.assertIn('fatigue[Title]', query)
        self.assertNotIn('"', query)  # Avoid dependence on PubMed phrase index.
        self.assertEqual(self.search.call_count, 1)
        self.builder.assert_not_called()
        self.fetch.assert_called_once_with(['42816300'])
        self.assert_unrestricted()

    def test_different_title_does_not_masquerade_as_exact_match(self):
        self.fetch.side_effect = [[{**ARTICLE, 'titulo': 'A different article'}], [ARTICLE]]
        _, articles = topic.search_articles(TITLE)
        self.builder.assert_called_once_with(TITLE)
        self.assertEqual(self.search.call_count, 2)
        self.assertEqual(articles, [ARTICLE])
        self.assert_unrestricted()

    def test_general_topic_uses_only_that_topic_and_no_date_range(self):
        self.read.side_effect = [{'IdList': []}, {'IdList': ['42816300']}]
        query, articles = topic.search_articles('Fadiga no câncer de mama')
        self.builder.assert_called_once_with('Fadiga no câncer de mama')
        self.assertEqual(query, self.builder.return_value)
        self.assertEqual(articles, [ARTICLE])
        self.assert_unrestricted()

    def test_direct_identifiers_do_not_expand_to_unrelated_articles(self):
        for value in ('42816300', 'PMID: 42816300', 'https://pubmed.ncbi.nlm.nih.gov/42816300/',
                      '10.1016/j.soncn.2026.152354', 'https://doi.org/10.1016/j.soncn.2026.152354'):
            with self.subTest(value=value):
                _, articles = topic.search_articles(value)
                self.assertEqual(articles, [ARTICLE])
        self.builder.assert_not_called()
        self.assert_unrestricted()
        self.read.return_value = {'IdList': []}
        self.assertEqual(topic.search_articles('999999999')[1], [])
        self.builder.assert_not_called()

    def test_empty_result_does_not_fetch_weekly_candidates(self):
        self.read.return_value = {'IdList': []}
        self.assertEqual(topic.search_articles('tema sem resultados')[1], [])
        self.fetch.assert_not_called()
        self.assert_unrestricted()


class TopicRoutesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.states = {}
        self.weekly = MagicMock(side_effect=AssertionError('Weekly pipeline must not run'))
        app = FastAPI()
        names = {'PodcastTopicSearchInput', 'PodcastTopicScriptInput', 'model_to_dict',
                 'queue_podcast_topic', 'buscar_podcast_tema', 'gerar_podcast_tema',
                 'processar_podcast_tema_background', 'processar_boletim_background'}
        tree = ast.parse(Path(__file__).with_name('api.py').read_text())
        nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and n.name in names]
        ns = {'BaseModel': BaseModel, 'Field': Field, 'BackgroundTasks': BackgroundTasks,
              'JSONResponse': JSONResponse, 'uuid': uuid, 'app': app, 'rodar_boletim': self.weekly,
              'save_task': lambda key, state: self.states.update({key: copy.deepcopy(state)}),
              'load_task': lambda key: self.states.get(key)}
        exec(compile(ast.Module(body=nodes, type_ignores=[]), 'api.py', 'exec'), ns)
        modules = patch.dict('sys.modules', {'boletim_service': SimpleNamespace(BASE_DIR=self.tmp.name, client=MagicMock())})
        modules.start()
        self.addCleanup(modules.stop)
        self.client = TestClient(app)
        self.addCleanup(self.client.close)

    def test_dedicated_search_cannot_run_weekly_even_with_conflicting_flags(self):
        with patch.object(topic, 'search_articles', return_value=('cancer[Title]', [ARTICLE])) as search:
            response = self.client.post('/podcast-tema/buscar', json={
                'tema_podcast': TITLE, 'modo_podcast': 'semanal', 'mailchimp': True,
                'resumos': True, 'somente_curadoria': False})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['modo_podcast'], 'tema')
        search.assert_called_once_with(TITLE)
        self.weekly.assert_not_called()
        state = self.states[response.json()['task_id']]
        self.assertEqual(state['status'], 'completed')
        self.assertEqual(state['result']['artigos_sugeridos'], [ARTICLE])
        self.assertTrue(state['result']['busca']['sem_limite_data'])

    def test_blank_topic_rejected_before_scheduling(self):
        for value in ('', '  '):
            response = self.client.post('/podcast-tema/buscar', json={'tema_podcast': value})
            self.assertIn(response.status_code, (400, 422))
        self.assertEqual(self.states, {})

    def test_script_route_rejects_missing_selection_and_never_runs_weekly(self):
        response = self.client.post('/podcast-tema/gerar-roteiro', json={'tema_podcast': TITLE})
        self.assertEqual(response.status_code, 422)
        self.assertEqual(self.states, {})
        self.weekly.assert_not_called()

    def test_provider_failure_remains_error_without_weekly_fallback(self):
        with patch.object(topic, 'search_articles', side_effect=RuntimeError('PubMed unavailable')):
            response = self.client.post('/podcast-tema/buscar', json={'tema_podcast': TITLE})
        self.assertEqual(self.states[response.json()['task_id']]['status'], 'error')
        self.weekly.assert_not_called()


if __name__ == '__main__':
    unittest.main()
