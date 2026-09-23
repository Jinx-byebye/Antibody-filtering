"""Bidirectional pDockQ2."""

import math
import numpy as np

from abfilter.data import Candidate, Sample
from .common import distances


def directional(sample: Sample, left: np.ndarray, right: np.ndarray) -> float:
    contact = distances(sample, left, right) <= 8.0
    if not np.any(contact):
        return 0.0
    interface = np.concatenate((left[np.any(contact, axis=1)], right[np.any(contact, axis=0)]))
    mean_plddt = float(np.mean(sample.plddt[interface]))
    pae = sample.pae[np.ix_(left, right)][contact]
    x = mean_plddt * float(np.mean(1.0 / (1.0 + (pae / 10.0) ** 2)))
    return 1.31 / (1.0 + math.exp(-0.075 * (x - 84.733))) + 0.005


def score(sample: Sample, candidate: Candidate) -> float:
    return max(directional(sample, sample.antibody_index, sample.antigen_index),
               directional(sample, sample.antigen_index, sample.antibody_index))
