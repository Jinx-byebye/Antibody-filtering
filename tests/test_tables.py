import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from abfilter.evaluate import evaluate_group, _rank_results, ENDPOINTS, MATCHED_ENDPOINTS
from abfilter.tables import load_scores, display_value


class TableTests(unittest.TestCase):
    def test_score_coverage_and_duplicates(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'scores.csv'
            path.write_text('pair_id,nac\na,0.1\nb,0.2\n')
            self.assertEqual(load_scores(path, ['nac'], ['a', 'b'])['a']['nac'], .1)
            with self.assertRaisesRegex(ValueError, 'population mismatch'):
                load_scores(path, ['nac'], ['a', 'b', 'c'])
            path.write_text('pair_id,nac\na,0.1\na,0.2\n')
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                load_scores(path, ['nac'], ['a'])
            path.write_text('pair_id,nac\na,nan\n')
            with self.assertRaisesRegex(ValueError, 'non-finite'):
                load_scores(path, ['nac'], ['a'])

    def test_invalid_labels_and_shapes(self):
        for scores, labels in [([], []), ([1, 2], [1]), ([1, 2], [1, 2]), ([1, float('nan')], [0, 1])]:
            with self.assertRaises(ValueError):
                evaluate_group(scores, labels)

    def test_random_ties_and_permutation(self):
        labels = np.array([1] * 10 + [0] * 90)
        first = evaluate_group(np.ones(100), labels)
        second = evaluate_group(np.ones(100), labels[::-1])
        self.assertEqual(first, second)
        self.assertAlmostEqual(first['recall_q_05'], .05)
        self.assertAlmostEqual(first['recall_q_20'], .2)

    def test_direction_and_selected_counts(self):
        labels = np.array([1] * 20 + [0] * 80)
        result = evaluate_group(np.arange(100)[::-1], labels)
        self.assertEqual(result['recall_at_utility'], 1)
        self.assertEqual(result['precision_at_utility'], 1)
        self.assertEqual(result['fpr_at_utility'], 0)
        self.assertEqual(result['true_positive'], 20)

    def test_rank_endpoint_sets(self):
        rows = {'a': {e: .5 for e in ENDPOINTS}, 'b': {e: .5 for e in ENDPOINTS}}
        rows['a']['recall_at_utility'] = 1
        _rank_results(rows, MATCHED_ENDPOINTS)
        self.assertEqual(rows['a']['average_rank'], 1.5)
        _rank_results(rows, ENDPOINTS)
        self.assertLess(rows['a']['average_rank'], rows['b']['average_rank'])

    def test_display_precision(self):
        self.assertEqual(display_value('utility', .146456), '0.1465')
        self.assertEqual(display_value('recall_q_20', .21271), '21.27')


if __name__ == '__main__':
    unittest.main()
