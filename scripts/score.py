#!/usr/bin/env python3
"""Score validated prediction ensembles for all candidates in a manifest."""

import argparse
import json
from pathlib import Path

from abfilter.data import read_candidates, read_ensemble
from abfilter.score import score_candidate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--predictor", required=True)
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--samples", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    candidates = read_candidates(args.candidates)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        for candidate in candidates:
            ensemble = read_ensemble(args.predictions, args.predictor, candidate, args.runs, args.samples)
            record = {
                "pair_id": candidate.pair_id, "predictor": args.predictor,
                "benchmark": candidate.benchmark, "target_id": candidate.target_id,
                "label": candidate.label,
                **score_candidate(candidate, ensemble),
            }
            handle.write(json.dumps(record, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
