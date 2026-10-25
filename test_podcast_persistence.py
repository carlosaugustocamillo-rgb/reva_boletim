"""Exercise private backups across loss of the Railway container's local disk."""
import copy
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import podcast_editorial as editorial
from test_podcast_editorial import fake_client
from test_pubmed_related import ANCHOR


class PersistenceTest(unittest.TestCase):
    def test_backup_restores_exact_version_after_local_disk_loss(self):
        draft = editorial.generate_episode([ANCHOR], None, fake_client())
        bucket = MagicMock()
        blob = bucket.blob.return_value
        with tempfile.TemporaryDirectory() as original, tempfile.TemporaryDirectory() as replacement:
            with patch.object(editorial, '_backup_bucket', return_value=bucket):
                editorial.save_draft(original, draft)
                backup = blob.upload_from_string.call_args.args[0]
                self.assertEqual(json.loads(backup), draft)
                self.assertEqual(blob.metadata, {'created_at': draft['created_at']})
                blob.make_public.assert_not_called()
                blob.download_as_bytes.return_value = backup.encode()
                recovered = editorial.load_draft(replacement, draft['id'])
                self.assertEqual(recovered, draft)
                self.assertEqual(editorial.transcript(recovered), editorial.transcript(draft))
                self.assertTrue(editorial.draft_path(replacement, draft['id']).exists())
                approved = editorial.approve_draft(replacement, draft['id'], draft['sha256'])
                self.assertEqual(approved['status'], 'approved')
                self.assertEqual(json.loads(blob.upload_from_string.call_args.args[0])['status'], 'approved')

    def test_remote_hash_and_identity_are_checked_before_caching(self):
        draft = editorial.generate_episode([ANCHOR], None, fake_client())
        for change in ('script', 'id'):
            value = copy.deepcopy(draft)
            if change == 'script':
                value['script']['title'] = 'Alterado depois da revisão'
            else:
                value['id'] = 'a' * 32
            bucket = MagicMock()
            bucket.blob.return_value.download_as_bytes.return_value = json.dumps(value).encode()
            with tempfile.TemporaryDirectory() as root, patch.object(editorial, '_backup_bucket', return_value=bucket):
                with self.assertRaises(editorial.EditorialError):
                    editorial.load_draft(root, draft['id'])
                self.assertFalse(editorial.draft_path(root, draft['id']).exists())

    def test_latest_backup_uses_creation_time_not_last_approval(self):
        bucket = MagicMock()
        older = SimpleNamespace(name='internal/revacast/editorial/' + 'a' * 32 + '.json', metadata={'created_at': '2026-10-01T00:00:00Z'})
        newer = SimpleNamespace(name='internal/revacast/editorial/' + 'b' * 32 + '.json', metadata={'created_at': '2026-10-24T00:00:00Z'})
        irrelevant = SimpleNamespace(name='internal/revacast/editorial/../bad.json', metadata={})
        bucket.list_blobs.return_value = [newer, older, irrelevant]
        with patch.object(editorial, '_backup_bucket', return_value=bucket):
            self.assertEqual(editorial.latest_remote_draft_id(), 'b' * 32)

    def test_local_tests_do_not_connect_to_storage(self):
        with patch.dict('os.environ', {}, clear=True):
            self.assertIsNone(editorial._backup_bucket())
        with patch.dict('os.environ', {'RAILWAY_ENVIRONMENT_ID': 'production', 'PODCAST_DRAFT_REMOTE_BACKUP': 'false'}, clear=True):
            self.assertIsNone(editorial._backup_bucket())


if __name__ == '__main__':
    unittest.main()
