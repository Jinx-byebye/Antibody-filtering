"""Predictor-native pairwise geometric predicted distance error."""

from abfilter.data import Candidate, Sample


def score(sample: Sample, candidate: Candidate) -> float | None:
    return sample.native.get("pair_gpde")
