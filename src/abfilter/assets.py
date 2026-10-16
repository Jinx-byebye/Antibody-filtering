"""Check local model files against a published checksum manifest."""
import hashlib
from pathlib import Path


def verify_assets(root, assets):
    root = Path(root).resolve()
    results = []
    for name, expected in assets.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f'Asset path leaves the model directory: {name}')
        result = {'file': name, 'ok': False}
        if not path.is_file():
            result['error'] = 'missing file'
        elif path.stat().st_size != expected['bytes']:
            result['error'] = 'file size differs'
        else:
            checksum = hashlib.sha256()
            with path.open('rb') as handle:
                for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                    checksum.update(block)
            result['sha256'] = checksum.hexdigest()
            result['ok'] = result['sha256'] == expected['sha256']
            if not result['ok']:
                result['error'] = 'SHA-256 differs'
        results.append(result)
    if not results:
        raise ValueError('Asset manifest contains no files')
    return results
