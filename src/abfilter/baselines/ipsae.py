"""Shared directional ipSAE calculation for the three published normalizations."""

import numpy as np

from abfilter.data import Candidate, Sample
from .common import chain_pairs


def _direction(sample: Sample, left: np.ndarray, right: np.ndarray, mode: str) -> float:
    error = sample.pae[np.ix_(left, right)]
    valid = error < 10.0
    counts = valid.sum(axis=1).astype(float)
    if mode == "residue":
        d0 = np.maximum(1.0, 1.24 * np.cbrt(np.maximum(26.0, counts) - 15.0) - 1.8)[:, None]
    else:
        length = len(left) + len(right) if mode == "chain" else int(np.any(valid, axis=1).sum() + np.any(valid, axis=0).sum())
        d0 = 1.0 if length <= 27 else max(1.0, 1.24 * (length - 15.0) ** (1.0 / 3.0) - 1.8)
    transformed = 1.0 / (1.0 + (error / d0) ** 2)
    rows = np.divide((transformed * valid).sum(axis=1), counts,
                     out=np.zeros_like(counts), where=counts > 0)
    return float(np.max(rows, initial=0.0))


def score(sample: Sample, candidate: Candidate, mode: str) -> float:
    if mode not in {"residue", "domain", "chain"}:
        raise ValueError("unknown ipSAE normalization")
    return max(max(_direction(sample, left, right, mode), _direction(sample, right, left, mode))
               for left, right in chain_pairs(sample, candidate))
