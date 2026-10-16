#!/usr/bin/env python3
"""Convert explicitly indexed raw predictions to validated ensemble files."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from abfilter.data import read_candidates, read_ensemble
from abfilter.predictors.convert import convert_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--predictor', required=True, choices=['alphafold3', 'protenix', 'esmfold2'])
    parser.add_argument('--seeds', required=True, help='Ordered comma-separated inference seeds')
    parser.add_argument('--samples', type=int, default=5)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    candidates = {c.pair_id: c for c in read_candidates(args.candidates)}
    seeds = [int(x) for x in args.seeds.split(',')]
    if len(set(seeds)) != len(seeds) or len(seeds) < 2 or args.samples < 1:
        raise ValueError('at least two distinct seeds and one sample are required')
    records = [json.loads(line) for line in args.manifest.read_text().splitlines() if line.strip()]
    observed = [(r['pair_id'], r['seed'], r['sample']) for r in records]
    expected = {(p, seed, sample) for p in candidates for seed in seeds for sample in range(args.samples)}
    if len(set(observed)) != len(observed) or set(observed) != expected:
        raise ValueError('raw prediction manifest does not cover the complete requested ensembles')
    for row in records:
        output = args.output / args.predictor / row['pair_id'] / f"run_{seeds.index(row['seed']):03d}" / f"sample_{row['sample']:02d}.npz"
        paths = [Path(row[name]) for name in ('structure', 'confidence', 'summary')]
        paths = [p if p.is_absolute() else args.manifest.parent / p for p in paths]
        convert_sample(candidates[row['pair_id']], args.predictor, *paths, output,
                       row['chain_map'], row['structure_plddt_scale'], row.get('confidence_plddt_scale', 100))
    for candidate in candidates.values():
        read_ensemble(args.output, args.predictor, candidate, len(seeds), args.samples)
    print(f'Validated {len(candidates)} complete ensembles.')


if __name__ == '__main__':
    main()
