"""Shared geometry and chain-pair helpers for published baselines."""

from __future__ import annotations

import numpy as np

from abfilter.data import Candidate, Sample


def distances(sample: Sample, left: np.ndarray, right: np.ndarray) -> np.ndarray:
    delta = sample.coordinates[left, None, :] - sample.coordinates[None, right, :]
    return np.sqrt(np.sum(delta * delta, axis=-1))


def chain_pairs(sample: Sample, candidate: Candidate):
    for antibody in candidate.antibody_sequences:
        left = np.flatnonzero(sample.chain_ids == antibody)
        for antigen in candidate.antigen_sequences:
            right = np.flatnonzero(sample.chain_ids == antigen)
            if not len(left) or not len(right):
                raise ValueError("declared chain is absent")
            yield left, right


def native(sample: Sample, name: str) -> float | None:
    return sample.native.get(name)
