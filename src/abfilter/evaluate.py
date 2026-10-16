"""Pooled, within-target, and matched-population filtering evaluation."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

import numpy as np

from .utility import utility


ENDPOINTS = ("recall_fpr_10", "recall_fpr_20", "recall_fpr_30",
             "recall_q_05", "recall_q_20", "utility", "recall_at_utility", "precision_at_utility")
MATCHED_ENDPOINTS = tuple(x for x in ENDPOINTS if x != "recall_at_utility")
LOWER_IS_BETTER = {"pair_gpde", "interface_pae"}


def _arrays(rows: list[dict[str, Any]], metric: str) -> tuple[np.ndarray, np.ndarray]:
    scores = np.asarray([float(row[metric]) for row in rows])
    labels = np.asarray([int(row["label"]) for row in rows])
    if not np.all(np.isfinite(scores)) or not np.all(np.isin(labels, [0, 1])):
        raise ValueError("non-finite score or invalid label")
    return (-scores if metric in LOWER_IS_BETTER else scores), labels


def _boundaries(scores: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(-scores, kind="stable")
    ranked = scores[order]
    sorted_labels = labels[order]
    ends = np.flatnonzero(np.r_[ranked[1:] != ranked[:-1], True])
    return ranked[ends], np.cumsum(sorted_labels == 1)[ends], np.cumsum(sorted_labels == 0)[ends]


def _recall_fpr(tp: np.ndarray, fp: np.ndarray, positive: int, negative: int, limit: float) -> float:
    feasible = fp / negative <= limit
    return float(np.max(tp[feasible], initial=0)) / positive


def _recall_budget(scores: np.ndarray, labels: np.ndarray, fraction: float) -> float:
    budget = math.floor(fraction * len(scores))
    if not budget:
        return 0.0
    boundary = np.partition(scores, len(scores) - budget)[len(scores) - budget]
    higher = scores > boundary
    equal = scores == boundary
    expected = labels[higher].sum() + (budget - higher.sum()) * labels[equal].sum() / equal.sum()
    return float(expected / labels.sum())


def evaluate_group(scores: np.ndarray, labels: np.ndarray, *, f_safe: float = .10,
                   f_severe: float = .30) -> dict[str, float]:
    scores, labels = np.asarray(scores, dtype=float), np.asarray(labels)
    if scores.ndim != 1 or labels.shape != scores.shape or not len(scores):
        raise ValueError("scores and labels must be nonempty vectors of equal length")
    if not np.all(np.isfinite(scores)) or not np.all(np.isin(labels, [0, 1])):
        raise ValueError("non-finite score or invalid label")
    positive = int(labels.sum())
    negative = len(labels) - positive
    if not positive or not negative:
        raise ValueError("both classes are required")
    thresholds, tp, fp = _boundaries(scores, labels)
    options = [(math.inf, 0, 0)] + list(zip(thresholds, tp, fp))
    threshold, best_tp, best_fp = max(options, key=lambda item: (
        utility(int(item[1]), int(item[2]), positive, negative, f_safe=f_safe, f_severe=f_severe),
        -int(item[2]) / negative, int(item[1]) / positive, float(item[0])))
    selected = int(best_tp) + int(best_fp)
    return {
        "n": len(scores), "positive": positive, "negative": negative,
        "recall_fpr_10": _recall_fpr(tp, fp, positive, negative, 0.10),
        "recall_fpr_20": _recall_fpr(tp, fp, positive, negative, 0.20),
        "recall_fpr_30": _recall_fpr(tp, fp, positive, negative, 0.30),
        "recall_q_05": _recall_budget(scores, labels, 0.05),
        "recall_q_20": _recall_budget(scores, labels, 0.20),
        "utility": utility(int(best_tp), int(best_fp), positive, negative, f_safe=f_safe, f_severe=f_severe),
        "recall_at_utility": int(best_tp) / positive,
        "fpr_at_utility": int(best_fp) / negative,
        "true_positive": int(best_tp), "false_positive": int(best_fp),
        "precision_at_utility": int(best_tp) / selected if selected else 0.0,
        "threshold": float(threshold),
    }


def _rank_results(results: dict[str, dict[str, float]], endpoints=ENDPOINTS) -> None:
    ranks = {name: [] for name in results}
    for endpoint in endpoints:
        ordered = sorted(results, key=lambda name: -results[name][endpoint])
        start = 0
        while start < len(ordered):
            end = start + 1
            while end < len(ordered) and math.isclose(results[ordered[start]][endpoint],
                                                       results[ordered[end]][endpoint], rel_tol=1e-12, abs_tol=1e-12):
                end += 1
            for name in ordered[start:end]:
                ranks[name].append((start + 1 + end) / 2.0)
            start = end
    for name in results:
        results[name]["average_rank"] = float(np.mean(ranks[name]))


def _shared_metrics(rows: list[dict[str, Any]]) -> set[str]:
    fields = set.intersection(*(set(row) for row in rows))
    return fields - {"pair_id", "predictor", "benchmark", "target_id", "label"}


def _pooled(rows: list[dict[str, Any]], endpoints=ENDPOINTS) -> dict[str, dict[str, float]]:
    results = {}
    for metric in sorted(_shared_metrics(rows)):
        scores, labels = _arrays(rows, metric)
        results[metric] = evaluate_group(scores, labels)
    _rank_results(results, endpoints)
    return results


def _campaign_groups(rows: list[dict[str, Any]]) -> dict[tuple[str, str], list[dict[str, Any]]]:
    groups = defaultdict(list)
    for row in rows:
        if row["benchmark"] in {"yeast_display", "igdesign"}:
            groups[(row["benchmark"], row["target_id"])].append(row)
    return {key: group for key, group in groups.items()
            if sum(row["label"] == 1 for row in group) >= 10
            and sum(row["label"] == 0 for row in group) >= 40}


def _campaign(groups: dict[tuple[str, str], list[dict[str, Any]]]) -> dict[str, dict[str, float]]:
    if {key[0] for key in groups} != {"yeast_display", "igdesign"}:
        raise ValueError("both campaign benchmarks must have eligible targets")
    information = {}
    for key, group in groups.items():
        positive = sum(row["label"] == 1 for row in group)
        negative = len(group) - positive
        information[key] = math.sqrt(2.0 * positive * negative / len(group))
    weights = {key: 0.5 * value / sum(v for k, v in information.items() if k[0] == key[0])
               for key, value in information.items()}
    metrics = set.intersection(*(_shared_metrics(group) for group in groups.values()))
    result = {}
    for metric in sorted(metrics):
        summaries = {key: evaluate_group(*_arrays(group, metric)) for key, group in groups.items()}
        result[metric] = {endpoint: float(sum(weights[key] * summary[endpoint]
                                             for key, summary in summaries.items())) for endpoint in ENDPOINTS}
        result[metric]["n"] = sum(len(group) for group in groups.values())
        result[metric]["target_count"] = len(groups)
    _rank_results(result)
    return result


def evaluate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Evaluate each predictor on its own fixed population and matched subset."""
    output = {}
    for predictor in sorted({row["predictor"] for row in rows}):
        cohort = [row for row in rows if row["predictor"] == predictor]
        if len({row["pair_id"] for row in cohort}) != len(cohort):
            raise ValueError("duplicate candidate in predictor cohort")
        groups = _campaign_groups(cohort)
        matched = [row for group in groups.values() for row in group]
        output[predictor] = {
            "portfolio": _pooled(cohort),
            "campaign": _campaign(groups),
            "matched_portfolio": _pooled(matched, MATCHED_ENDPOINTS),
            "coverage": {"portfolio": len(cohort), "campaign": len(matched), "targets": len(groups)},
        }
    return output
