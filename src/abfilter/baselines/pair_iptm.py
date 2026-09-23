"""Predictor-native pairwise interface pTM."""

from abfilter.data import Candidate, Sample


def score(sample: Sample, candidate: Candidate) -> float | None:
    return sample.native.get("pair_iptm")
