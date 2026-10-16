#!/usr/bin/env python3
"""Score validated prediction ensembles for all candidates in a manifest."""

import argparse
import csv
import json
from pathlib import Path

from abfilter.data import read_candidates, read_ensemble
from abfilter.score import score_candidate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--predictor", required=True, choices=['alphafold3', 'protenix', 'esmfold2'])
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--samples", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument('--config', type=Path, default=Path(__file__).resolve().parents[1] / 'configs/evaluation.json')
    args = parser.parse_args()
    candidates = read_candidates(args.candidates)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metrics = json.loads(args.config.read_text())['predictor_metrics'][args.predictor]
    with args.output.open("x", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=['pair_id'] + metrics) if args.output.suffix == '.csv' else None
        if writer:
            writer.writeheader()
        for candidate in candidates:
            ensemble = read_ensemble(args.predictions, args.predictor, candidate, args.runs, args.samples)
            record = {
                "pair_id": candidate.pair_id, "predictor": args.predictor,
                "benchmark": candidate.benchmark, "target_id": candidate.target_id,
                "label": candidate.label,
                **score_candidate(candidate, ensemble, args.predictor),
            }
            missing = set(metrics) - set(record)
            if missing:
                raise ValueError(f'{candidate.pair_id}: missing required metrics {sorted(missing)}')
            if writer:
                writer.writerow({key: record[key] for key in ['pair_id'] + metrics})
            else:
                handle.write(json.dumps(record, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
