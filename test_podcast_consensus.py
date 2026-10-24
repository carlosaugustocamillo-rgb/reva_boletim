"""Offline supplemental-source and import API regressions."""
import ast
import base64
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from fastapi import FastAPI
from fastapi.responses import JSONResponse
from fastapi.testclient import TestClient
from pydantic import BaseModel, Field

from connected_papers import build_manual_report
from podcast_consensus import build_consensus_import
from podcast_editorial import evidence_packet
from podcast_pdf import save_pdf
from test_pubmed_related import ANCHOR

REPORT = 'Relatório Consensus: síntese dos estudos, com achados e limitações para discussão clínica. ' * 3


class ConsensusTest(unittest.TestCase):
    def test_complete_report_survives_empty_reference_extraction(self):
        result = build_consensus_import(REPORT, {'referencias': [], 'erro': 'Sem referências identificadas.'})
        self.assertEqual(result['candidatos'][0]['evidence_content'], REPORT.strip())
        self.assertTrue(result['parser_warning'])
        packet = evidence_packet([ANCHOR], build_manual_report([ANCHOR], {'100': result['candidatos']}))
        source = packet[0]['sources'][1]
        self.assertEqual(source['content'], REPORT.strip())
        self.assertEqual(source['material'], 'consensus_report')
        self.assertEqual(source['abstract'], '')

    def test_consensus_extract_and_bibliography_keep_distinct_material(self):
        result = build_consensus_import(REPORT, {'referencias': [
            {'texto': 'Trial A', 'material': 'consensus_extract', 'evidence_content': 'Descrição complementar do ensaio.'},
            {'texto': 'Trial B', 'material': 'bibliographic_only'},
        ]})
        packet = evidence_packet([ANCHOR], build_manual_report([ANCHOR], {'100': result['candidatos']}))
        self.assertEqual(packet[0]['sources'][2]['content'], 'Descrição complementar do ensaio.')
        self.assertEqual(packet[0]['sources'][2]['abstract'], '')
        self.assertEqual(packet[0]['sources'][3]['material'], 'metadata_only')
        self.assertEqual(packet[0]['sources'][3]['content'], '')

    def test_independent_pdf_without_pmid_or_doi_is_read_by_writer(self):
        with tempfile.TemporaryDirectory() as root:
            pdf_dir = Path(root) / 'podcast_pdfs'
            pdf_dir.mkdir()
            source = {'source_id': 'manual-pdf:trial-a', 'titulo': 'Trial A'}
            metadata = {'id': 'a' * 32, 'filename': 'trial.pdf', 'sha256': 'test', 'page_count': 1,
                        'extracted_page_count': 1, 'character_count': 36, 'source': source,
                        'pages': [{'page': 1, 'text': 'Additional full text for the episode.'}]}
            (pdf_dir / ('a' * 32 + '.json')).write_text(json.dumps(metadata))
            reference = {**source, 'pdf_document': {'id': 'a' * 32}}
            context = build_manual_report([ANCHOR], {'100': [reference]})
            packet = evidence_packet([ANCHOR], context, base_dir=root)
            self.assertEqual(packet[0]['sources'][1]['material'], 'full_text_pdf')
            self.assertIn('Additional full text', packet[0]['sources'][1]['content'])

    def test_import_api_text_pdf_and_invalid_input(self):
        tree = ast.parse(Path('api.py').read_text())
        nodes = [n for n in tree.body if getattr(n, 'name', '') in {'PodcastConsensusInput', 'importar_consensus_podcast'}]
        with tempfile.TemporaryDirectory() as root:
            ns = {'app': FastAPI(), 'BaseModel': BaseModel, 'Field': Field, 'base64': base64,
                  'save_pdf': save_pdf, 'BASE_DATA_DIR': root, 'JSONResponse': JSONResponse}
            exec(compile(ast.Module(body=nodes, type_ignores=[]), 'api.py', 'exec'), ns)
            extract = Mock(return_value={'referencias': []})
            with patch.dict('sys.modules', {'revamais_service': SimpleNamespace(extrair_referencias_consensus=extract)}), TestClient(ns['app']) as client:
                response = client.post('/podcast-consensus', json={'relatorio': REPORT})
                self.assertEqual(response.status_code, 200)
                self.assertEqual(response.json()['candidatos'][0]['material'], 'consensus_report')
                self.assertEqual(client.post('/podcast-consensus', json={'relatorio': ''}).status_code, 400)
                self.assertEqual(client.post('/podcast-consensus', json={'pdf_base64': 'bad!'}).status_code, 400)
                # Mock only PDF extraction; exercise base64 decoding and route adaptation.
                with patch('podcast_pdf.load_pdf', return_value={'pages': [{'page': 1, 'text': REPORT}]}), patch.dict(ns, {'save_pdf': Mock(return_value={'id': 'a' * 32})}):
                    response = client.post('/podcast-consensus', json={'pdf_base64': base64.b64encode(b'%PDF-test').decode()})
                self.assertEqual(response.status_code, 200)
                self.assertIn(REPORT.strip(), response.json()['candidatos'][0]['evidence_content'])


if __name__ == '__main__':
    unittest.main()
