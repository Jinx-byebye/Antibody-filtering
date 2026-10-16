"""Wilson-bound filtering utility."""

import math


def wilson(success: int, total: int) -> tuple[float, float]:
    if total <= 0 or not 0 <= success <= total:
        raise ValueError("invalid binomial counts")
    z = 1.959963984540054
    p = success / total
    denominator = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / denominator
    radius = z * math.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)


def utility(tp: int, fp: int, positive: int, negative: int, *, f_safe: float = 0.10,
            f_severe: float = 0.30) -> float:
    if not 0 < f_safe < f_severe <= 1:
        raise ValueError("require 0 < f_safe < f_severe <= 1")
    if positive <= 0 or negative <= 0 or not 0 <= tp <= positive or not 0 <= fp <= negative:
        raise ValueError("invalid filtering counts")
    recall_lower = wilson(tp, positive)[0]
    fpr_upper = wilson(fp, negative)[1]
    excess = max(0.0, fpr_upper - f_safe)
    return recall_lower * math.exp(-(fpr_upper / f_severe) ** 2 - (excess / (f_severe - f_safe)) ** 2)
