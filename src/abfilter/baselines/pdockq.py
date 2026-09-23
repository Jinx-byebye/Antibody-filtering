"""FoldDock pDockQ."""

import math
import numpy as np

from abfilter.data import Candidate, Sample
from .common import distances


def score(sample: Sample, candidate: Candidate) -> float:
    left, right = sample.antibody_index, sample.antigen_index
    contact = distances(sample, left, right) <= 8.0
    count = int(contact.sum())
    if not count:
        return 0.0
    interface = np.concatenate((left[np.any(contact, axis=1)], right[np.any(contact, axis=0)]))
    x = float(np.mean(sample.plddt[interface])) * math.log10(count)
    return 0.724 / (1.0 + math.exp(-0.052 * (x - 152.611))) + 0.018
