"""Mean predicted aligned error at the geometric interface."""

import numpy as np

from abfilter.data import Candidate, Sample
from .common import distances


def score(sample: Sample, candidate: Candidate) -> float:
    left, right = sample.antibody_index, sample.antigen_index
    contact = distances(sample, left, right) <= 8.0
    error = sample.pae[np.ix_(left, right)]
    reverse = sample.pae[np.ix_(right, left)].T
    if not np.any(contact):
        return float(max(np.max(error), np.max(reverse)))
    return float(np.mean(0.5 * (error[contact] + reverse[contact])))
