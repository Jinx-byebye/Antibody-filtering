"""ipSAE with residue-specific distance normalization."""

from abfilter.data import Candidate, Sample
from .ipsae import score as _score


def score(sample: Sample, candidate: Candidate) -> float:
    return _score(sample, candidate, "residue")
