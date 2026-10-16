#!/usr/bin/env python3
"""Prepare scoring inputs from the public sequence dataset."""

import argparse
import json
from pathlib import Path
from abfilter.dataset import load_dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = load_dataset(args.dataset)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as f:
        for pair_id in sorted(rows):
            f.write(json.dumps(rows[pair_id]["candidate"].to_record(), allow_nan=False) + "\n")
    print(f"Prepared {len(rows)} antibody-antigen pairs.")


if __name__ == "__main__":
    main()
