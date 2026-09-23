"""Whole-antibody mean pLDDT."""

import numpy as np

from abfilter.data import Candidate, Sample


def score(sample: Sample, candidate: Candidate) -> float:
    return float(np.mean(sample.plddt[sample.antibody_index]))
