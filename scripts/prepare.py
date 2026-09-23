#!/usr/bin/env python3
"""Validate a private candidate manifest without changing its contents."""

import argparse
import json
from collections import Counter
from pathlib import Path

from abfilter.data import read_candidates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    candidates = read_candidates(args.manifest)
    print(json.dumps({
        "candidates": len(candidates),
        "benchmarks": dict(Counter(item.benchmark for item in candidates)),
        "labels": dict(Counter(item.label for item in candidates)),
        "targets": len({(item.benchmark, item.target_id) for item in candidates}),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
