"""Candidate-level scores from complete prediction ensembles."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from . import gric, nac, ric
from .baselines import (actifptm, antibody_plddt, anticonf, interface_pae, ipsae_chn,
                        ipsae_dom, ipsae_res, pair_gpde, pair_iptm, pdockq,
                        pdockq2, ranking)
from .data import Candidate, Sample


BASELINES: dict[str, Callable[[Sample, Candidate], float | None]] = {
    "pair_iptm": pair_iptm.score,
    "pair_gpde": pair_gpde.score,
    "interface_pae": interface_pae.score,
    "antibody_plddt": antibody_plddt.score,
    "ranking_score": ranking.score,
    "pdockq": pdockq.score,
    "pdockq2": pdockq2.score,
    "ipsae_d0res": ipsae_res.score,
    "ipsae_d0dom": ipsae_dom.score,
    "ipsae_d0chn": ipsae_chn.score,
    "anticonf": anticonf.score,
    "actifptm": actifptm.score,
}


def score_candidate(candidate: Candidate, ensemble: list[list[Sample]]) -> dict[str, float]:
    """Compute scores without accessing labels or other candidates."""
    if not ensemble or any(not run or len(run) != len(ensemble[0]) for run in ensemble):
        raise ValueError("ensemble is incomplete")
    scores = {"nac": nac.score(ensemble)}
    contact_available = [sample.contact_forward is not None for run in ensemble for sample in run]
    if any(contact_available) and not all(contact_available):
        raise ValueError("contact probability availability changes within the ensemble")
    scores["ric" if all(contact_available) else "gric"] = (
        ric.score(ensemble) if all(contact_available) else gric.score(ensemble)
    )
    for name, calculate in BASELINES.items():
        runs = [[calculate(sample, candidate) for sample in run] for run in ensemble]
        present = [value is not None for run in runs for value in run]
        if any(present) and not all(present):
            raise ValueError(f"{name} is only partially available")
        if all(present):
            scores[name] = float(np.median([np.median(run) for run in runs]))
    return scores
