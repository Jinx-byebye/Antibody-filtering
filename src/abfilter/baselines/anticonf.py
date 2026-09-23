"""AntiConf, combining antigen-directed pDockQ2 with predictor pTM."""

from abfilter.data import Candidate, Sample
from .pdockq2 import directional


def score(sample: Sample, candidate: Candidate) -> float | None:
    ptm = sample.native.get("ptm")
    if ptm is None:
        return None
    return 0.3 * directional(sample, sample.antigen_index, sample.antibody_index) + 0.7 * ptm
