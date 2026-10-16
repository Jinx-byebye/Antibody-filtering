#!/usr/bin/env python3
"""Check a downloaded antibody dataset without loading prediction files."""

import argparse
import csv
import json
import math
from pathlib import Path


EXPECTED = {"display_enrichment": (1400, 700, 7), "igdesign_spr": (1171, 249, 7),
            "literature_affinity": (1613, 1510, 165)}
AA = set("ACDEFGHIKLMNPQRSTVWY")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(root):
    report = {}
    all_ids = set()
    for name, (count, positive, antigens) in EXPECTED.items():
        folder = root / name
        with (folder / "sequences.csv").open(encoding="utf-8", newline="") as handle:
            sequences = list(csv.DictReader(handle))
        with (folder / "metadata.jsonl").open(encoding="utf-8") as handle:
            metadata = [json.loads(line) for line in handle if line.strip()]
        sequence_map = {row["pair_id"]: row for row in sequences}
        metadata_map = {row["pair_id"]: row for row in metadata}
        require(len(sequences) == len(sequence_map) == count, f"{name}: sequence count")
        require(len(metadata) == len(metadata_map) == count, f"{name}: metadata count")
        require(sequence_map.keys() == metadata_map.keys(), f"{name}: unmatched pair identifiers")
        require(not all_ids.intersection(sequence_map), f"{name}: duplicate pair across datasets")
        all_ids.update(sequence_map)
        require(all(type(x["label"]) is int and x["label"] in (0, 1) for x in metadata), "Invalid label")
        require(sum(x["label"] for x in metadata) == positive, f"{name}: positive count")
        require(len({x["antigen_id"] for x in sequences}) == antigens, f"{name}: antigen count")
        measurements = 0
        for pair_id, sequence in sequence_map.items():
            item = metadata_map[pair_id]
            chains = sequence["antigen_chain_sequences"].split(";")
            require(len(chains) == len(sequence["antigen_chain_ids"].split(";")), "Antigen chain mismatch")
            for chain in [sequence["heavy_chain_sequence"], *chains]:
                require(bool(chain) and not set(chain) - AA, "Invalid amino-acid sequence")
            require(not set(sequence["light_chain_sequence"]) - AA, "Invalid light chain")
            for chain, positions in item["cdr_sequence_positions"].items():
                length = len(sequence["heavy_chain_sequence" if chain == "H" else "light_chain_sequence"])
                require(positions == sorted(set(positions)) and all(1 <= p <= length for p in positions), "Invalid CDR positions")
            if name == "display_enrichment":
                rounds = item["measurements"]
                require(len(rounds) >= 2, "Insufficient selection rounds")
                require(len({r["round"] for r in rounds}) == len(rounds), "Duplicate selection round")
                for row in rounds:
                    n, total = row["read_count"], row["total_read_count"]
                    require(type(n) is int and type(total) is int and 0 <= n <= total and total > 0, "Invalid read count")
                    require(math.isclose(row["frequency"], n / total), "Frequency mismatch")
                first, last = rounds[0], rounds[-1]
                expected = math.log2(((last["read_count"] + .5) / (last["total_read_count"] + 1)) /
                                     ((first["read_count"] + .5) / (first["total_read_count"] + 1)))
                require(math.isclose(expected, item["enrichment_log2_ratio"]), "Enrichment mismatch")
                require(math.isfinite(item["selection_enrichment_log2_ratio"]), "Missing source enrichment")
                measurements += len(rounds)
            elif name == "igdesign_spr":
                total, binding = item["replicate_count"], item["binding_replicate_count"]
                require(total >= 2 and binding == (total if item["label"] else 0), "SPR label mismatch")
                kd = item["reported_kd_nm"]
                require(kd is None or (math.isfinite(kd) and kd > 0), "Invalid SPR affinity")
                measurements += total
            else:
                values = item["measurements"]
                require(bool(values), "Missing affinity evidence")
                for value in values:
                    require(value["kd_comparator"] in {"=", "<", "<=", ">", ">="}, "Invalid comparator")
                    require(math.isfinite(value["kd_value_nm"]) and value["kd_value_nm"] > 0, "Invalid affinity")
                require(item["affinity_summary"]["measurement_count"] == len(values), "Affinity count mismatch")
                measurements += len(values)
        measurement_key = {"display_enrichment": "round_count_records",
                           "igdesign_spr": "reported_spr_repeats",
                           "literature_affinity": "affinity_records"}[name]
        report[name] = {"pairs": count, "positive": positive, "antigen_constructs": antigens,
                        measurement_key: measurements}
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path, help="Directory containing the three dataset folders")
    args = parser.parse_args()
    print(json.dumps(validate(args.dataset), indent=2))
