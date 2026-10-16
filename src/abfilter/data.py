"""Validated candidate metadata and normalized prediction ensembles."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


BENCHMARKS = {"yeast_display", "igdesign", "literature"}


@dataclass(frozen=True)
class Candidate:
    pair_id: str
    benchmark: str
    target_id: str
    antibody_sequences: dict[str, str]
    antigen_sequences: dict[str, str]
    cdr_mask: tuple[bool, ...]
    cdr_annotation: dict[str, Any]
    label: int
    evidence: dict[str, Any]

    @classmethod
    def from_record(cls, record: dict[str, Any]) -> "Candidate":
        required = {
            "pair_id", "benchmark", "target_id", "antibody_sequences",
            "antigen_sequences", "cdr_mask", "cdr_annotation", "label", "evidence",
        }
        if set(record) != required:
            raise ValueError(f"candidate fields differ from schema: {sorted(set(record) ^ required)}")
        pair_id = str(record["pair_id"])
        if not pair_id or "/" in pair_id or "\\" in pair_id or pair_id in {".", ".."}:
            raise ValueError("pair_id must be a non-empty path-safe identifier")
        benchmark = str(record["benchmark"])
        if benchmark not in BENCHMARKS:
            raise ValueError(f"unsupported benchmark: {benchmark}")
        antibody = _sequences(record["antibody_sequences"], "antibody_sequences")
        antigen = _sequences(record["antigen_sequences"], "antigen_sequences")
        if set(antibody) & set(antigen):
            raise ValueError("antibody and antigen chain identifiers overlap")
        mask = record["cdr_mask"]
        if not isinstance(mask, list) or not mask or any(type(value) is not bool for value in mask):
            raise ValueError("cdr_mask must be a non-empty Boolean list")
        if len(mask) != sum(map(len, antibody.values())) or not any(mask) or all(mask):
            raise ValueError("cdr_mask does not match antibody residues or lacks CDR/non-CDR residues")
        annotation = record["cdr_annotation"]
        if not isinstance(annotation, dict) or annotation.get("scheme") != "IMGT":
            raise ValueError("cdr_annotation must declare the IMGT scheme")
        label = record["label"]
        if type(label) is not int or label not in (0, 1):
            raise ValueError("label must be 0 or 1")
        evidence = record["evidence"]
        if not isinstance(evidence, dict) or not isinstance(evidence.get("type"), str):
            raise ValueError("evidence must include a type")
        target_id = str(record["target_id"])
        if not target_id:
            raise ValueError("target_id is required")
        return cls(pair_id, benchmark, target_id, antibody, antigen, tuple(mask), annotation, label, evidence)

    def to_record(self) -> dict[str, Any]:
        return {
            "pair_id": self.pair_id,
            "benchmark": self.benchmark,
            "target_id": self.target_id,
            "antibody_sequences": self.antibody_sequences,
            "antigen_sequences": self.antigen_sequences,
            "cdr_mask": list(self.cdr_mask),
            "cdr_annotation": self.cdr_annotation,
            "label": self.label,
            "evidence": self.evidence,
        }


def _sequences(value: Any, name: str) -> dict[str, str]:
    if not isinstance(value, dict) or not value:
        raise ValueError(f"{name} must map chain identifiers to sequences")
    result = {str(chain): str(sequence).upper() for chain, sequence in value.items()}
    if any(not chain or not sequence.isalpha() for chain, sequence in result.items()):
        raise ValueError(f"{name} contains an invalid chain or sequence")
    return result


def read_candidates(path: Path) -> list[Candidate]:
    with path.open(encoding="utf-8") as handle:
        candidates = [Candidate.from_record(json.loads(line)) for line in handle if line.strip()]
    ids = [candidate.pair_id for candidate in candidates]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate pair_id in candidate metadata")
    if not candidates:
        raise ValueError("candidate metadata is empty")
    return candidates


@dataclass(frozen=True)
class Sample:
    coordinates: np.ndarray
    plddt: np.ndarray
    pae: np.ndarray
    chain_ids: np.ndarray
    antibody_index: np.ndarray
    antigen_index: np.ndarray
    cdr_mask: np.ndarray
    contact_forward: np.ndarray | None
    contact_reverse: np.ndarray | None
    native: dict[str, float]


def read_sample(path: Path, candidate: Candidate) -> Sample:
    with np.load(path, allow_pickle=False) as archive:
        coordinates = np.asarray(archive["coordinates"], dtype=np.float64)
        plddt = np.asarray(archive["plddt"], dtype=np.float64)
        pae = np.asarray(archive["pae"], dtype=np.float64)
        chains = np.asarray(archive["chain_ids"], dtype=str)
        residue_symbols = np.asarray(archive["residue_symbols"], dtype=str)
        forward = np.asarray(archive["contact_forward"], dtype=np.float64) if "contact_forward" in archive else None
        reverse = np.asarray(archive["contact_reverse"], dtype=np.float64) if "contact_reverse" in archive else None
        native = {
            name: float(np.asarray(archive[name]).item())
            for name in ("pair_iptm", "pair_gpde", "ranking_score", "ptm", "actifptm")
            if name in archive
        }
    n = len(plddt)
    if coordinates.shape != (n, 3) or pae.shape != (n, n) or chains.shape != (n,) or residue_symbols.shape != (n,):
        raise ValueError(f"inconsistent prediction dimensions: {path}")
    if not all(np.all(np.isfinite(value)) for value in (coordinates, plddt, pae)):
        raise ValueError(f"non-finite prediction values: {path}")
    if np.any((plddt < 0) | (plddt > 100)) or np.any(pae < 0):
        raise ValueError(f"pLDDT or PAE outside its declared range: {path}")
    if any(not np.isfinite(value) for value in native.values()):
        raise ValueError(f"non-finite native confidence: {path}")
    antibody = np.flatnonzero(np.isin(chains, list(candidate.antibody_sequences)))
    antigen = np.flatnonzero(np.isin(chains, list(candidate.antigen_sequences)))
    if len(antibody) + len(antigen) != n:
        raise ValueError(f'undeclared chains in prediction: {path}')
    if len(antibody) != len(candidate.cdr_mask) or not len(antigen):
        raise ValueError(f"chain roles or antibody annotation do not match: {path}")
    for role, sequences in ((antibody, candidate.antibody_sequences), (antigen, candidate.antigen_sequences)):
        for chain, sequence in sequences.items():
            observed = "".join(residue_symbols[chains == chain].tolist()).upper()
            if observed != sequence:
                raise ValueError(f"chain sequence differs from candidate metadata: {path}")
    antibody_order = "".join(residue_symbols[antibody].tolist()).upper()
    if antibody_order != "".join(candidate.antibody_sequences.values()):
        raise ValueError(f"antibody chain order differs from CDR annotation order: {path}")
    if (forward is None) != (reverse is None):
        raise ValueError(f"contact directions must both be present or absent: {path}")
    if forward is not None:
        shape = (len(antibody), len(antigen))
        if forward.shape != shape or reverse.shape != shape:
            raise ValueError(f"contact matrix does not match partners: {path}")
        if any(np.any(~np.isfinite(value)) or np.any((value < 0) | (value > 1)) for value in (forward, reverse)):
            raise ValueError(f"contact probability outside [0, 1]: {path}")
    return Sample(coordinates, plddt, pae, chains, antibody, antigen,
                  np.asarray(candidate.cdr_mask, dtype=bool), forward, reverse, native)


def read_ensemble(root: Path, predictor: str, candidate: Candidate, runs: int = 10,
                  samples_per_run: int = 5) -> list[list[Sample]]:
    if not predictor or "/" in predictor or "\\" in predictor:
        raise ValueError("predictor must be a path-safe name")
    if runs < 2 or samples_per_run < 1:
        raise ValueError("RIC requires at least two runs and one sample per run")
    directory = root / predictor / candidate.pair_id
    expected = {
        directory / f"run_{run:03d}" / f"sample_{sample:02d}.npz"
        for run in range(runs) for sample in range(samples_per_run)
    }
    observed = set(directory.glob("run_*/sample_*.npz"))
    if observed != expected:
        raise ValueError(f"incomplete or unexpected ensemble for {candidate.pair_id}: expected {len(expected)}, found {len(observed)}")
    return [[read_sample(directory / f"run_{run:03d}" / f"sample_{sample:02d}.npz", candidate)
             for sample in range(samples_per_run)] for run in range(runs)]
