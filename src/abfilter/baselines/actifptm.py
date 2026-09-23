"""Predictor-native actifpTM; expected PAE alone is insufficient to derive it."""

from abfilter.data import Candidate, Sample


def score(sample: Sample, candidate: Candidate) -> float | None:
    return sample.native.get("actifptm")
