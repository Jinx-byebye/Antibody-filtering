#!/usr/bin/env python3
"""Verify independently obtained model assets without loading them into memory."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from abfilter.assets import verify_assets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictor', required=True, choices=['alphafold3', 'protenix', 'esmfold2'])
    parser.add_argument('--models', required=True, type=Path)
    parser.add_argument('--manifest', type=Path, default=ROOT / 'reference/model_assets.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    results = verify_assets(args.models, manifest['predictors'][args.predictor])
    report = {'predictor': args.predictor, 'ok': all(r['ok'] for r in results), 'files': results}
    rendered = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')
    return 0 if report['ok'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
