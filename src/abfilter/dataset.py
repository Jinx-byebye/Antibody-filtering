"""Read the public sequence and experimental-data package."""

import csv
import json
from pathlib import Path

from .data import Candidate


DATASETS = {"display_enrichment": "yeast_display", "igdesign_spr": "igdesign",
            "literature_affinity": "literature"}


def load_dataset(root: Path) -> dict[str, dict]:
    result = {}
    for name, benchmark in DATASETS.items():
        folder = root / name
        with (folder / "sequences.csv").open(newline="", encoding="utf-8") as f:
            sequences = list(csv.DictReader(f))
        with (folder / "metadata.jsonl").open(encoding="utf-8") as f:
            metadata = [json.loads(line) for line in f if line.strip()]
        smap = {x["pair_id"]: x for x in sequences}
        mmap = {x["pair_id"]: x for x in metadata}
        if len(smap) != len(sequences) or len(mmap) != len(metadata) or smap.keys() != mmap.keys():
            raise ValueError(f"{name}: duplicate or unmatched pair identifiers")
        for pair_id, sequence in smap.items():
            if pair_id in result:
                raise ValueError(f"Duplicate pair across datasets: {pair_id}")
            item = mmap[pair_id]
            antibody = {"H": sequence["heavy_chain_sequence"]}
            if sequence["light_chain_sequence"]:
                antibody["L"] = sequence["light_chain_sequence"]
            ids = sequence["antigen_chain_ids"].split(";")
            chains = sequence["antigen_chain_sequences"].split(";")
            if len(ids) != len(chains) or len(set(ids)) != len(ids):
                raise ValueError(f"Invalid antigen chains: {pair_id}")
            positions = item["cdr_sequence_positions"]
            if set(positions) != set(antibody):
                raise ValueError(f"CDR chain mismatch: {pair_id}")
            mask = []
            for chain, seq in antibody.items():
                p = positions[chain]
                if len(p) != len(set(p)) or any(type(i) is not int or not 1 <= i <= len(seq) for i in p):
                    raise ValueError(f"Invalid CDR positions: {pair_id}")
                mask.extend(i + 1 in p for i in range(len(seq)))
            candidate = Candidate.from_record({
                "pair_id": pair_id, "benchmark": benchmark, "target_id": item["target_name"],
                "antibody_sequences": antibody, "antigen_sequences": dict(zip(ids, chains)),
                "cdr_mask": mask, "cdr_annotation": {"scheme": item["cdr_numbering_scheme"]},
                "label": item["label"], "evidence": {"type": name},
            })
            result[pair_id] = {"candidate": candidate, "dataset": name,
                               "target_id": sequence["antigen_id"], "label": candidate.label}
    return result
