"""Predictor-native ranking confidence."""

from abfilter.data import Candidate, Sample


def score(sample: Sample, candidate: Candidate) -> float | None:
    return sample.native.get("ranking_score")
