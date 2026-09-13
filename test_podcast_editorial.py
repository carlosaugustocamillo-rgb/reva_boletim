"""Source-bound contracts and approval state; no API credentials or paid calls."""
import copy
import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import podcast_editorial as editorial
from test_pubmed_related import ANCHOR, candidate


def editorial_response(**kwargs):
    """Synthetic content, explicitly only for contract/regression tests."""
    stage = kwargs['response_format']['json_schema']['name']
    payload = json.loads(kwargs['messages'][1]['content'])
    packet = payload['evidence']
    if stage == 'podcast_plan':
        value = {'central_question': 'O que o resumo permite concluir?', 'studies': []}
        for item in packet:
            main = item['sources'][0]
            support = {'source_id': main['source_id'], 'quote': main['abstract']}
            value['studies'].append({'pmid': item['pmid'], 'question': main['title'],
                'finding': 'O resumo não mostrou diferença na fadiga.', 'supports': [support],
                'appraisal': [{'kind': 'not_reported', 'explanation': 'O resumo não informa perdas.', 'supports': []}],
                'caution': 'O resumo não informa perdas; não podemos avaliar esse aspecto.',
                'connection_to_next': ''})
    elif stage == 'podcast_write':
        value = {'title': 'Teste sintético de contrato',
            'opening': [{'speaker': 'HOST', 'text': 'Eu sou Ivo. Hoje analisamos apenas os resumos disponíveis.', 'source_ids': []}],
            'studies': [], 'closing': [{'speaker': 'COHOST', 'text': 'Eu sou Manu. Até a próxima.', 'source_ids': []}]}
        for item in packet:
            caution = 'O resumo não informa perdas; não podemos avaliar esse aspecto.'
            value['studies'].append({'pmid': item['pmid'], 'spoken_caution': caution, 'dialogue': [
                {'speaker': 'HOST', 'text': 'Não houve diferença na fadiga.', 'source_ids': [item['sources'][0]['source_id']]},
                {'speaker': 'COHOST', 'text': caution, 'source_ids': [item['sources'][0]['source_id']]},
            ]})
    else:
        value = {'issues': []}
    return SimpleNamespace(choices=[SimpleNamespace(finish_reason='stop', message=SimpleNamespace(content=json.dumps(value), refusal=None))],
                           model=kwargs['model'], usage=SimpleNamespace(prompt_tokens=150, completion_tokens=120))


def fake_client():
    client = MagicMock()
    client.with_options.return_value.chat.completions.create.side_effect = editorial_response
    return client


class EditorialTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.context = {'estudos': [{'pmid_ancora': '100', 'modo': 'curadoria_manual', 'referencias': [candidate()]}]}
        self.client = fake_client()
        self.draft = editorial.generate_episode([ANCHOR], self.context, self.client)

    def test_three_global_calls_with_structured_schema_no_temperature(self):
        calls = self.client.with_options.return_value.chat.completions.create.call_args_list
        self.assertEqual(len(calls), 3)
        for call in calls:
            args = call.kwargs
            self.assertTrue(args['response_format']['json_schema']['strict'])
            self.assertNotIn('temperature', args)
            self.assertEqual(args['model'], editorial.DEFAULT_MODEL)
            self.assertLessEqual(args['max_completion_tokens'], 12000)
        self.assertEqual(len(self.draft['usage']), 3)

    def test_model_override_is_podcast_specific(self):
        with patch.dict(os.environ, {'PODCAST_SCRIPT_MODEL': 'gpt-5.5', 'OPENAI_TEXT_MODEL': 'unchanged-email'}):
            draft = editorial.generate_episode([ANCHOR], self.context, fake_client())
        self.assertEqual(draft['model'], 'gpt-5.5')

    def test_metadata_only_not_promoted_to_source_and_manual_mode_preserved(self):
        self.context['estudos'][0]['referencias'][0]['resumo_original'] = ''
        packet = editorial.evidence_packet([ANCHOR], self.context)
        self.assertEqual(packet[0]['selection_mode'], 'curadoria_manual')
        self.assertEqual(packet[0]['sources'][1]['material'], 'metadata_only')
        broken = copy.deepcopy(self.draft['plan'])
        broken['studies'][0]['supports'] = [{'source_id': '100:ref1', 'quote': 'Invented finding'}]
        with self.assertRaises(editorial.EditorialError):
            editorial.validate_plan(broken, packet)

    def test_unknown_source_invented_quote_or_main_finding_from_reference_rejected(self):
        for support in ({'source_id': 'unknown', 'quote': 'text'},
                        {'source_id': '100:main', 'quote': 'Invented finding'},
                        {'source_id': '100:ref1', 'quote': candidate()['resumo_original']}):
            plan = copy.deepcopy(self.draft['plan'])
            plan['studies'][0]['supports'] = [support]
            with self.assertRaises(editorial.EditorialError):
                editorial.validate_plan(plan, self.draft['evidence'])

    def test_omitted_duplicate_and_reordered_studies_rejected(self):
        packet = editorial.evidence_packet([ANCHOR, {**ANCHOR, 'pmid': '101'}], None)
        for ids in (['100'], ['100', '100'], ['101', '100']):
            with self.assertRaises(editorial.EditorialError):
                editorial.validate_script({'studies': [{'pmid': p} for p in ids]}, packet)

    def test_spoken_caution_and_source_ids_required(self):
        for mutation in ('caution', 'source', 'speaker'):
            script = copy.deepcopy(self.draft['script'])
            study = script['studies'][0]
            if mutation == 'caution':
                study['spoken_caution'] = 'A ressalva desapareceu da fala.'
            elif mutation == 'source':
                study['dialogue'][0]['source_ids'] = ['other:main']
            else:
                study['dialogue'][1]['speaker'] = 'HOST'
            with self.assertRaises(editorial.EditorialError):
                editorial.validate_script(script, self.draft['evidence'])

    def test_caution_typography_and_split_turns_do_not_change_spoken_text(self):
        script = copy.deepcopy(self.draft['script'])
        study = script['studies'][0]
        study['spoken_caution'] = study['spoken_caution'].upper().replace(';', ',')
        before = copy.deepcopy(study['dialogue'])
        editorial.validate_script(script, self.draft['evidence'])
        self.assertEqual(study['dialogue'], before)
        study['dialogue'][1]['text'] = 'O resumo não informa perdas.'
        study['dialogue'].append({'speaker': 'HOST', 'text': 'Não podemos avaliar esse aspecto.', 'source_ids': ['100:main']})
        editorial.validate_script(script, self.draft['evidence'])

    def test_negation_and_decimal_changes_stay_blocked_with_pmid(self):
        for spoken, registered in (
            ('O resumo não informa perdas.', 'O resumo informa perdas.'),
            ('A diferença foi de 0.5 pontos.', 'A diferença foi de 05 pontos.'),
        ):
            script = copy.deepcopy(self.draft['script'])
            script['studies'][0]['dialogue'][1]['text'] = spoken
            script['studies'][0]['spoken_caution'] = registered
            with self.assertRaisesRegex(editorial.CautionMismatch, 'PMID 100'):
                editorial.validate_script(script, self.draft['evidence'])

    def test_one_repair_then_full_validation_and_audit(self):
        writes = []
        def respond(**kwargs):
            response = editorial_response(**kwargs)
            if kwargs['response_format']['json_schema']['name'] == 'podcast_write':
                writes.append(kwargs)
                if len(writes) == 1:
                    content = json.loads(response.choices[0].message.content)
                    content['studies'][0]['spoken_caution'] = 'Paráfrase diferente da fala.'
                    response.choices[0].message.content = json.dumps(content)
            return response
        client = fake_client()
        client.with_options.return_value.chat.completions.create.side_effect = respond
        draft = editorial.generate_episode([ANCHOR], None, client)
        self.assertEqual(len(writes), 2)
        self.assertEqual(draft['status'], 'pending_review')
        self.assertTrue(draft['repair_attempted'])
        self.assertIn('original_script', draft)
        self.assertEqual(draft['usage'][-1]['stage'], 'audit')

    def test_persistent_mismatch_preserves_readable_blocked_draft_and_no_unbounded_retry(self):
        def respond(**kwargs):
            response = editorial_response(**kwargs)
            if kwargs['response_format']['json_schema']['name'] == 'podcast_write':
                content = json.loads(response.choices[0].message.content)
                content['studies'][0]['spoken_caution'] = 'Paráfrase diferente da fala.'
                response.choices[0].message.content = json.dumps(content)
            return response
        client = fake_client()
        client.with_options.return_value.chat.completions.create.side_effect = respond
        draft = editorial.generate_episode([ANCHOR], None, client)
        self.assertEqual(len(client.with_options.return_value.chat.completions.create.call_args_list), 3)
        self.assertEqual(draft['status'], 'blocked')
        self.assertIn('PMID 100', draft['audit']['issues'][0]['reason'])
        editorial.save_draft(self.root, draft)
        loaded = editorial.load_draft(self.root, draft['id'])
        self.assertIn('Manu:', editorial.review_payload(loaded)['text'])
        self.assertFalse(editorial.review_payload(loaded)['can_approve'])
        with self.assertRaises(editorial.EditorialError):
            editorial.approve_draft(self.root, draft['id'], draft['sha256'])

    def test_audit_outage_preserves_script_but_cannot_be_approved(self):
        def respond(**kwargs):
            if kwargs['response_format']['json_schema']['name'] == 'podcast_audit':
                raise TimeoutError('Synthetic timeout')
            return editorial_response(**kwargs)
        client = fake_client()
        client.with_options.return_value.chat.completions.create.side_effect = respond
        draft = editorial.generate_episode([ANCHOR], None, client)
        self.assertEqual(draft['status'], 'blocked')
        self.assertIn('Synthetic timeout', draft['audit']['issues'][0]['reason'])
        self.assertTrue(editorial.transcript(draft))

    def test_refusal_incomplete_or_malformed_json_never_becomes_spoken_text(self):
        for reason, content, refusal in [('length', '{}', None), ('stop', 'bad JSON', None), ('stop', '{}', 'refused')]:
            client = MagicMock()
            client.with_options.return_value.chat.completions.create.return_value = SimpleNamespace(
                choices=[SimpleNamespace(finish_reason=reason, message=SimpleNamespace(content=content, refusal=refusal))],
                model='gpt-6-astra', usage=None)
            with self.assertRaises(editorial.EditorialError):
                editorial.generate_episode([ANCHOR], None, client)

    def test_approval_bound_to_exact_version_and_readable_transcript(self):
        editorial.save_draft(self.root, self.draft)
        for sha in (self.draft['sha256'], 'stale'):
            with self.assertRaises(editorial.EditorialError):
                editorial.approved_audio(self.root, self.draft['id'], sha)
        with self.assertRaises(editorial.EditorialError):
            editorial.approve_draft(self.root, self.draft['id'], 'stale')
        editorial.approve_draft(self.root, self.draft['id'], self.draft['sha256'])
        approved = editorial.approved_audio(self.root, self.draft['id'], self.draft['sha256'])
        self.assertEqual(editorial.dialogues(approved), editorial.dialogues(self.draft))
        self.assertIn('Ivo:', editorial.transcript(approved))
        self.assertIn('Manu:', editorial.transcript(approved))
        self.assertNotIn('source_ids', editorial.transcript(approved))

    def test_blocking_audit_and_tampered_content_fail_closed(self):
        self.draft['audit']['issues'] = [{'severity': 'blocking', 'location': '100', 'reason': 'Unsupported result'}]
        self.draft['sha256'] = editorial.fingerprint(self.draft)
        editorial.save_draft(self.root, self.draft)
        with self.assertRaises(editorial.EditorialError):
            editorial.approve_draft(self.root, self.draft['id'], self.draft['sha256'])
        self.draft['script']['title'] = 'Changed after review'
        editorial.save_draft(self.root, self.draft)
        with self.assertRaises(editorial.EditorialError):
            editorial.load_draft(self.root, self.draft['id'])

    def test_paths_cannot_escape_editorial_directory(self):
        for draft_id in ('../../.env', 'a' * 32 + '\n', None, 'abc'):
            with self.assertRaises(editorial.EditorialError):
                editorial.draft_path(self.root, draft_id)


if __name__ == '__main__':
    unittest.main()
