import hashlib
import tempfile
import unittest
from pathlib import Path

from abfilter.assets import verify_assets


class AssetTests(unittest.TestCase):
    def test_valid_missing_and_changed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'weights.bin'
            path.write_bytes(b'test')
            expected = {'weights.bin': {'bytes': 4, 'sha256': hashlib.sha256(b'test').hexdigest()}}
            self.assertTrue(verify_assets(directory, expected)[0]['ok'])
            path.write_bytes(b'edit')
            self.assertEqual(verify_assets(directory, expected)[0]['error'], 'SHA-256 differs')
            path.write_bytes(b'longer')
            self.assertEqual(verify_assets(directory, expected)[0]['error'], 'file size differs')
            path.unlink()
            self.assertEqual(verify_assets(directory, expected)[0]['error'], 'missing file')

    def test_empty_and_escaping_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                verify_assets(directory, {})
            with self.assertRaises(ValueError):
                verify_assets(directory, {'../outside': {}})
