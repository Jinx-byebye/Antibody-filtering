#!/usr/bin/env python3
"""Build a validated candidate manifest from benchmark source tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

from abfilter.data import Candidate


CDR_RANGES = ((27, 38), (56, 65), (105, 117))


def records(path: Path, *, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def lines(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def unique(rows: list[dict[str, Any]], field: str) -> dict[str, dict[str, Any]]:
    result = {}
    for row in rows:
        key = str(row[field])
        if key in result:
            raise ValueError(f"duplicate {field}: {key}")
        result[key] = row
    return result


def annotations(paths: list[Path]) -> dict[str, set[int]]:
    result: dict[str, set[int]] = {}
    for path in paths:
        for row in records(path):
            identifier = row["Id"]
            positions = result.setdefault(identifier, set())
            index = int(row["seqstart_index"])
            for name, residue in row.items():
                match = re.fullmatch(r"(\d+)[A-Z]*", name or "")
                if match is None or not residue or residue == "-":
                    continue
                number = int(match.group(1))
                if any(start <= number <= end for start, end in CDR_RANGES):
                    positions.add(index)
                index += 1
    return result


def chain_sequences(job: dict[str, Any], path_maps: list[tuple[str, str]]) -> dict[str, str]:
    source_path = str(job["protenix_json_path"])
    for old, new in path_maps:
        if source_path.startswith(old):
            source_path = new + source_path[len(old):]
            break
    path = Path(source_path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    entries = payload if isinstance(payload, list) else [payload]
    result = {}
    for entry in entries:
        for item in entry.get("sequences", []):
            protein = item.get("proteinChain") or item.get("protein") or {}
            sequence = protein.get("sequence")
            identifiers = item.get("id", protein.get("id"))
            if not isinstance(identifiers, list):
                identifiers = [identifiers]
            for identifier in identifiers:
                if identifier is not None and sequence:
                    result[str(identifier)] = str(sequence).upper()
    return result


def cdr_mask(annotation: dict[str, set[int]], identity: str,
             antibody: dict[str, str]) -> list[bool]:
    result = []
    for chain, sequence in antibody.items():
        key = f"{identity}__chain_{chain}"
        if key not in annotation:
            matches = [value for name, value in annotation.items()
                       if name.startswith(f"{identity}__chain_")]
            if len(antibody) != 1 or len(matches) != 1:
                raise ValueError(f"missing or ambiguous CDR annotation: {identity}")
            positions = matches[0]
        else:
            positions = annotation[key]
        if not positions or max(positions) >= len(sequence):
            raise ValueError(f"invalid CDR positions: {identity}")
        result.extend(index in positions for index in range(len(sequence)))
    if not any(result) or all(result):
        raise ValueError(f"invalid antibody CDR scope: {identity}")
    return result


def optional_number(value: Any) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def build(args: argparse.Namespace) -> list[Candidate]:
    path_maps = []
    for mapping in getattr(args, "path_map", []):
        if "=" not in mapping:
            raise ValueError("path mapping must be OLD=NEW")
        old, new = mapping.split("=", 1)
        if not old or not new:
            raise ValueError("path mapping must have non-empty prefixes")
        path_maps.append((old, new))
    composition = records(args.records, delimiter="\t")
    yeast = unique(lines(args.yeast_panel), "pair_id")
    yeast_jobs = unique(lines(args.yeast_predictions), "pair_id")
    igdesign_jobs = unique(lines(args.igdesign_predictions), "job_id")
    source_jobs = unique(lines(args.source_predictions), "pair_id")
    affinity = unique(records(args.source_affinity), "pair_id")
    cdr = {
        "yeast_display": annotations(args.cdr_yeast),
        "igdesign": annotations(args.cdr_igdesign),
        "literature": annotations(args.cdr_source),
    }
    candidates = []
    for row in composition:
        category = row["benchmark"]
        pair_id = row["pair_id"]
        if category.startswith("Yeast"):
            benchmark = "yeast_display"
            job = yeast_jobs[pair_id]
            panel = yeast[pair_id]
            metadata = panel["metadata"]
            identity = pair_id
            evidence = {
                "type": "selection_enrichment",
                "scale": str(metadata.get("grade_semantics") or "ordinal_enrichment"),
                "score": optional_number(metadata.get("evidence_score")),
                "grade": metadata.get("source_grade_name"),
                "tier": metadata.get("evidence_tier"),
            }
        elif category.startswith("IgDesign"):
            benchmark = "igdesign"
            job = igdesign_jobs[pair_id]
            identity = str(job["record_id"])
            evidence = {
                "type": "spr_replicates",
                "replicates": int(float(job["replicates"])),
                "binding_replicates": int(float(job["binding_replicates"])),
                "method": job.get("source_method"),
            }
            expected_labels = {"binder", "positive"} if row["label"] == "positive" else {"nonbinder", "negative"}
            if str(job["label"]).lower() not in expected_labels:
                raise ValueError(f"SPR label mismatch: {pair_id}")
        elif category.startswith("Literature"):
            benchmark = "literature"
            job = source_jobs[pair_id]
            identity = pair_id
            measurement = affinity[pair_id]
            evidence = {
                "type": "reported_affinity",
                "unit": "nM",
                "kd_min": optional_number(measurement.get("kd_min_nm")),
                "kd_median": optional_number(measurement.get("kd_median_nm")),
                "kd_max": optional_number(measurement.get("kd_max_nm")),
                "measurement_count": int(measurement.get("kd_measurement_count") or 0),
                "document_count": int(measurement.get("document_count") or 0),
                "record_count": int(measurement.get("record_count") or 0),
            }
            if str(measurement["label"]).lower() not in {
                "binder" if row["label"] == "positive" else "nonbinder",
                row["label"].lower(),
            }:
                raise ValueError(f"affinity label mismatch: {pair_id}")
        else:
            raise ValueError(f"unsupported benchmark: {category}")
        chain_map = chain_sequences(job, path_maps)
        antibody_ids = [str(value) for value in job["antibody_chain_ids"]]
        antigen_ids = [str(value) for value in job["antigen_chain_ids"]]
        antibody = {key: chain_map[key] for key in antibody_ids}
        antigen = {key: chain_map[key] for key in antigen_ids}
        if "".join(antibody.values()) != row["heavy"] + row["light"]:
            raise ValueError(f"antibody sequence mismatch: {pair_id}")
        if "".join(antigen.values()) != row["antigen"]:
            raise ValueError(f"antigen sequence mismatch: {pair_id}")
        label = {"positive": 1, "negative": 0}[row["label"]]
        if benchmark == "yeast_display" and str(metadata["label"]) != ("binder" if label else "nonbinder"):
            raise ValueError(f"selection label mismatch: {pair_id}")
        candidates.append(Candidate.from_record({
            "pair_id": pair_id,
            "benchmark": benchmark,
            "target_id": row["target"],
            "antibody_sequences": antibody,
            "antigen_sequences": antigen,
            "cdr_mask": cdr_mask(cdr[benchmark], identity, antibody),
            "cdr_annotation": {"scheme": "IMGT", "ranges": [list(item) for item in CDR_RANGES]},
            "label": label,
            "evidence": evidence,
        }))
    if len(candidates) != args.expected or len({item.pair_id for item in candidates}) != len(candidates):
        raise ValueError("candidate count or identity differs from the declared benchmark")
    return candidates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--yeast-panel", type=Path, required=True)
    parser.add_argument("--yeast-predictions", type=Path, required=True)
    parser.add_argument("--igdesign-predictions", type=Path, required=True)
    parser.add_argument("--source-predictions", type=Path, required=True)
    parser.add_argument("--source-affinity", type=Path, required=True)
    parser.add_argument("--cdr-yeast", type=Path, action="append", required=True)
    parser.add_argument("--cdr-igdesign", type=Path, action="append", required=True)
    parser.add_argument("--cdr-source", type=Path, action="append", required=True)
    parser.add_argument("--path-map", action="append", default=[], metavar="OLD=NEW",
                        help="Replace an input-manifest path prefix when files have moved")
    parser.add_argument("--expected", type=int, default=4184)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    candidates = build(args)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        for candidate in candidates:
            handle.write(json.dumps(candidate.to_record(), sort_keys=True) + "\n")
    print(json.dumps({"candidates": len(candidates),
                      "benchmarks": dict(Counter(item.benchmark for item in candidates))}, sort_keys=True))


if __name__ == "__main__":
    main()
