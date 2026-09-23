"""Coordinate-contact approximation of regional interaction contrast."""

from .data import Sample
from .ric import score as _score


def score(ensemble: list[list[Sample]]) -> float:
    return _score(ensemble, geometric_contacts=True)
