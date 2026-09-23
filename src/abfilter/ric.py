"""Regional interaction contrast from a complete prediction ensemble."""

from __future__ import annotations

import numpy as np

from .data import Sample


def _rank_regions(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Assign fractional weights at ties around two adjacent rank regions."""
    size = len(values)
    if size == 1:
        return np.ones(1), np.ones(1)
    width = min(8, size // 2)
    order = np.argsort(-values, kind="stable")
    high = np.zeros(size)
    comparison = np.zeros(size)
    start = 0
    while start < size:
        end = start + 1
        while end < size and values[order[end]] == values[order[start]]:
            end += 1
        group = order[start:end]
        high[group] = max(0, min(end, width) - start) / (width * len(group))
        comparison[group] = max(0, min(end, 2 * width) - max(start, width)) / (width * len(group))
        start = end
    return high, comparison


def score(ensemble: list[list[Sample]], *, geometric_contacts: bool = False) -> float:
    """Score RIC, or its 8-A coordinate-contact approximation when requested."""
    if len(ensemble) < 2 or not ensemble[0] or any(len(run) != len(ensemble[0]) for run in ensemble):
        raise ValueError("a complete ensemble with at least two runs is required")
    first = ensemble[0][0]
    cdr = first.cdr_mask
    shape = (int(cdr.sum()), len(first.antigen_index))
    if not shape[0] or not shape[1]:
        raise ValueError("CDR and antigen residues are required")
    contacts = []
    forward = []
    reverse = []
    for run in ensemble:
        run_contacts = []
        run_forward = []
        run_reverse = []
        for sample in run:
            if not np.array_equal(sample.cdr_mask, cdr) or len(sample.antigen_index) != shape[1]:
                raise ValueError("residue mapping changes across the ensemble")
            if geometric_contacts:
                xyz_a = sample.coordinates[sample.antibody_index[cdr]]
                xyz_g = sample.coordinates[sample.antigen_index]
                contact = (np.linalg.norm(xyz_a[:, None] - xyz_g[None, :], axis=-1) <= 8.0).astype(float)
            else:
                if sample.contact_forward is None or sample.contact_reverse is None:
                    raise ValueError("RIC requires contact probabilities in both directions")
                contact = 0.5 * (sample.contact_forward[cdr] + sample.contact_reverse[cdr])
            if contact.shape != shape:
                raise ValueError("contact matrix changes across the ensemble")
            ab = sample.antibody_index[cdr]
            ag = sample.antigen_index
            run_contacts.append(contact)
            run_forward.append(sample.pae[np.ix_(ab, ag)])
            run_reverse.append(sample.pae[np.ix_(ag, ab)].T)
        contacts.append(run_contacts)
        forward.append(run_forward)
        reverse.append(run_reverse)
    contacts = np.asarray(contacts)
    forward = np.asarray(forward)
    reverse = np.asarray(reverse)
    forecast = contacts.mean(axis=1)
    support = forecast.max(axis=-1)
    regions = [[_rank_regions(row) for row in run] for run in forecast]
    run_scores = []
    for ref in range(len(ensemble)):
        other_scores = []
        for other in range(len(ensemble)):
            if ref == other:
                continue
            fsum = forward[ref, :, None] + forward[other, None, :]
            rsum = reverse[ref, :, None] + reverse[other, None, :]
            agreement = 0.5 * (16.0 / (16.0 + fsum) + 16.0 / (16.0 + rsum))
            delta = np.zeros(agreement.shape[:2] + (shape[0],))
            for residue, (high, comparison) in enumerate(regions[ref]):
                delta[:, :, residue] = np.sum(agreement[:, :, residue] * (high - comparison), axis=-1)
            denominator = support[ref].sum()
            pair_scores = np.sum(delta * support[ref] ** 2, axis=-1) / denominator if denominator else np.zeros(delta.shape[:2])
            other_scores.append(float(np.median(pair_scores)))
        run_scores.append(float(np.median(other_scores)))
    return float(np.median(run_scores))
