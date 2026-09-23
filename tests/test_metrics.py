import tempfile
import unittest
from pathlib import Path

import numpy as np

from abfilter import gric, nac, ric
from abfilter.baselines import antibody_plddt, interface_pae, ipsae_chn, pdockq
from abfilter.data import Candidate, Sample, read_ensemble, read_sample
from abfilter.score import score_candidate


def candidate() -> Candidate:
    return Candidate.from_record({
        "pair_id": "example", "benchmark": "yeast_display", "target_id": "target",
        "antibody_sequences": {"A": "AAA"}, "antigen_sequences": {"B": "AA"},
        "cdr_mask": [True, False, False], "cdr_annotation": {"scheme": "IMGT"},
        "label": 1, "evidence": {"type": "enrichment", "value": 2.0},
    })


def sample(contact=True, plddt=80.0) -> Sample:
    coordinates = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=float)
    pae = np.full((5, 5), 5.0)
    forward = np.array([[0.9, 0.1], [0.2, 0.1], [0.1, 0.2]]) if contact else None
    return Sample(coordinates, np.array([60, plddt, plddt, 50, 50], dtype=float), pae,
                  np.array(["A", "A", "A", "B", "B"]), np.arange(3), np.arange(3, 5),
                  np.array([True, False, False]), forward, forward, {"pair_iptm": 0.5})


class MetricTests(unittest.TestCase):
    def test_nac_is_non_cdr_and_hierarchical(self):
        ensemble = [[sample(plddt=70), sample(plddt=90)], [sample(plddt=80), sample(plddt=100)]]
        self.assertAlmostEqual(nac.score(ensemble), 0.85)
        self.assertAlmostEqual(antibody_plddt.score(sample(), candidate()), (60 + 80 + 80) / 3)

    def test_ric_requires_contact_and_constant_agreement_has_zero_contrast(self):
        ensemble = [[sample(), sample()], [sample(), sample()]]
        self.assertAlmostEqual(ric.score(ensemble), 0.0)
        self.assertTrue(np.isfinite(gric.score(ensemble)))
        with self.assertRaisesRegex(ValueError, "contact probabilities"):
            ric.score([[sample(False)], [sample(False)]])

    def test_baseline_boundaries(self):
        self.assertGreater(pdockq.score(sample(), candidate()), 0)
        self.assertAlmostEqual(interface_pae.score(sample(), candidate()), 5)
        self.assertGreater(ipsae_chn.score(sample(), candidate()), 0)

    def test_score_never_uses_label(self):
        first = candidate()
        other = Candidate.from_record({**first.to_record(), "label": 0})
        ensemble = [[sample()], [sample()]]
        self.assertEqual(score_candidate(first, ensemble), score_candidate(other, ensemble))

    def test_complete_ensemble_required(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "incomplete"):
                read_ensemble(Path(directory), "af3", candidate(), 2, 1)

    def test_loader_checks_sequences(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.npz"
            item = sample()
            np.savez(path, coordinates=item.coordinates, plddt=item.plddt,
                     pae=item.pae, chain_ids=item.chain_ids,
                     residue_symbols=np.array(list("AAAAA")))
            self.assertEqual(len(read_sample(path, candidate()).antibody_index), 3)
            np.savez(path, coordinates=item.coordinates, plddt=item.plddt,
                     pae=item.pae, chain_ids=item.chain_ids,
                     residue_symbols=np.array(list("BAAAA")))
            with self.assertRaisesRegex(ValueError, "sequence differs"):
                read_sample(path, candidate())


if __name__ == "__main__":
    unittest.main()
