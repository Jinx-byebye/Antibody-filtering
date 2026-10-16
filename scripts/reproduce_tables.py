#!/usr/bin/env python3
"""Recalculate Tables 2, 3 and 4 on a CPU."""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from abfilter.tables import reproduce, write_tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', type=Path, required=True)
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--populations', type=Path, required=True)
    parser.add_argument('--config', type=Path, default=ROOT / 'configs/evaluation.json')
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    tables, coverage = reproduce(args.dataset, args.scores, args.populations, config)
    write_tables(tables, config, args.output, args.reference)
    (args.output / 'coverage.json').write_text(json.dumps(coverage, indent=2) + '\n')
    print('Wrote Tables 2, 3 and 4 as CSV, Markdown and LaTeX.')
    if args.reference:
        print('All displayed values match the reference tables.')


if __name__ == '__main__':
    main()
