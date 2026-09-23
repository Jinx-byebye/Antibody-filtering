import unittest

import numpy as np

from abfilter.evaluate import _recall_budget, evaluate, evaluate_group
from abfilter.utility import utility, wilson


class EvaluationTests(unittest.TestCase):
    def test_tie_aware_budget(self):
        scores = np.array([1.0] * 20)
        labels = np.array([1] * 10 + [0] * 10)
        self.assertAlmostEqual(_recall_budget(scores, labels, 0.20), 0.2)

    def test_threshold_includes_empty_and_ties(self):
        scores = np.array([2.0, 2.0, 1.0, 1.0])
        labels = np.array([1, 0, 1, 0])
        result = evaluate_group(scores, labels)
        self.assertEqual(result["recall_fpr_10"], 0.0)
        self.assertIn(result["threshold"], (float("inf"), 1.0, 2.0))

    def test_wilson_utility(self):
        self.assertEqual(wilson(0, 10)[0], 0.0)
        self.assertGreater(utility(8, 0, 10, 90), 0.0)

    def test_three_settings(self):
        rows = []
        for benchmark in ("yeast_display", "igdesign"):
            for index in range(50):
                rows.append({
                    "pair_id": f"{benchmark}-{index}", "predictor": "af3",
                    "benchmark": benchmark, "target_id": "target", "label": int(index < 10),
                    "nac": float(50 - index), "ric": float(50 - index),
                })
        result = evaluate(rows)["af3"]
        self.assertEqual(result["coverage"], {"portfolio": 100, "campaign": 100, "targets": 2})
        self.assertEqual(result["campaign"]["nac"]["recall_fpr_10"], 1.0)
        self.assertEqual(result["matched_portfolio"]["ric"]["n"], 100)


if __name__ == "__main__":
    unittest.main()
