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


def utility(tp: int, fp: int, positive: int, negative: int) -> float:
    if positive <= 0 or negative <= 0 or not 0 <= tp <= positive or not 0 <= fp <= negative:
        raise ValueError("invalid filtering counts")
    recall_lower = wilson(tp, positive)[0]
    fpr_upper = wilson(fp, negative)[1]
    excess = max(0.0, fpr_upper - 0.10)
    return recall_lower * math.exp(-(fpr_upper / 0.30) ** 2 - (excess / 0.20) ** 2)
