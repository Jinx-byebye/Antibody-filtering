"""Validate the released HER2 cohort and its reproducible reference table."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from abfilter.data import read_candidates


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'her2'


class Her2Tests(unittest.TestCase):
    def test_fixed_framework_and_distinct_cdrs(self):
        candidates = read_candidates(DATA / 'candidates.jsonl')
        self.assertEqual(len(candidates), 1854)
        self.assertEqual(sum(c.label for c in candidates), 757)
        frameworks, cdrs = set(), set()
        for candidate in candidates:
            sequence = ''.join(candidate.antibody_sequences.values())
            frameworks.add(''.join(a for a, c in zip(sequence, candidate.cdr_mask) if not c))
            cdrs.add(''.join(a for a, c in zip(sequence, candidate.cdr_mask) if c))
            expected = []
            for chain, residues in candidate.antibody_sequences.items():
                positions = candidate.cdr_annotation['sequence_positions_1based'][chain]
                expected.extend(i + 1 in positions for i in range(len(residues)))
            self.assertEqual(tuple(expected), candidate.cdr_mask)
        self.assertEqual(len(frameworks), 1)
        self.assertEqual(len(cdrs), 1854)

    def test_reference_table(self):
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([
                sys.executable, str(ROOT / 'scripts' / 'evaluate_external.py'),
                '--candidates', str(DATA / 'candidates.jsonl'),
                '--scores', str(DATA / 'scores'),
                '--predictors', 'alphafold3', 'protenix', 'opendde',
                '--config', str(DATA / 'evaluation.json'), '--output', directory,
            ], check=True, capture_output=True, text=True)
            with (Path(directory) / 'all.csv').open(newline='') as f:
                actual = list(csv.reader(f))
            with (DATA / 'reference' / 'all.csv').open(newline='') as f:
                self.assertEqual(actual, list(csv.reader(f)))
            coverage = json.loads((Path(directory) / 'coverage.json').read_text())
            self.assertEqual(coverage['all']['candidates'], 1854)
            self.assertEqual(coverage['HCDR3_only']['candidates'], 1607)
            self.assertEqual(coverage['HCDR123_changes']['candidates'], 247)


if __name__ == '__main__':
    unittest.main()
