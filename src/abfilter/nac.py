"""Confidence outside annotated antibody CDRs."""

from __future__ import annotations

import numpy as np

from .data import Sample


def score(ensemble: list[list[Sample]]) -> float:
    if not ensemble or any(not run for run in ensemble):
        raise ValueError("a complete ensemble is required")
    run_scores = []
    for run in ensemble:
        values = []
        for sample in run:
            confidence = sample.plddt[sample.antibody_index][~sample.cdr_mask]
            if not len(confidence):
                raise ValueError("non-CDR antibody residues are required")
            values.append(float(np.mean(confidence) / 100.0))
        run_scores.append(float(np.median(values)))
    return float(np.median(run_scores))
