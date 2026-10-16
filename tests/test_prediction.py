import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from abfilter.data import Candidate, read_sample
from abfilter.msa import balance, validate_a3m
from abfilter.predictors.convert import convert_sample
from abfilter.predictors.inputs import command, prepare, resource_paths


class PredictionTests(unittest.TestCase):
    def test_protenix_omits_paired_search(self):
        candidate = self.candidate()
        inputs = prepare(candidate, 'protenix', {'AA': Path('/msa/aa.a3m'), 'G': Path('/msa/g.a3m')}, [101, 102])
        for item in inputs[0]['sequences']:
            chain = item['proteinChain']
            self.assertIsNone(chain['pairedMsaPath'])
            self.assertTrue(chain['unpairedMsaPath'])

    def test_interpreter_symlink_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root / 'base_python'
            base.touch()
            entry = root / 'venv_python'
            entry.symlink_to(base)
            self.assertEqual(resource_paths({'python': str(entry)})['python'], entry)
            self.assertNotEqual(resource_paths({'python': str(entry)})['python'], base)

    def test_smoke_configs_preserve_model_settings(self):
        root = Path(__file__).resolve().parents[1] / 'configs'
        for predictor in ('alphafold3', 'protenix', 'esmfold2'):
            paper = json.loads((root / 'predictors' / f'{predictor}.json').read_text())
            smoke = json.loads((root / 'smoke' / f'{predictor}.json').read_text())
            self.assertEqual(smoke['seeds'], paper['seeds'][:2])
            self.assertEqual(smoke['samples'], 1)
            for key, value in smoke.items():
                if key not in ('purpose', 'seeds', 'samples'):
                    self.assertEqual(value, paper[key], f'{predictor}: {key}')

    def candidate(self):
        return Candidate.from_record({'pair_id': 'example', 'benchmark': 'igdesign', 'target_id': 'target',
            'antibody_sequences': {'H': 'AA'}, 'antigen_sequences': {'A': 'G'},
            'cdr_mask': [True, False], 'cdr_annotation': {'scheme': 'IMGT'}, 'label': 1,
            'evidence': {'type': 'synthetic'}})

    def test_msa_balancing(self):
        sequences = {'a': 'A'*100, 'b': 'A'*80, 'c': 'A'*20}
        jobs = balance(sequences, 2)
        self.assertEqual(set(sum(jobs, [])), set(sequences))
        self.assertEqual(len(sum(jobs, [])), len(sequences))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'test.a3m'
            path.write_text('>query\nAA\n>hit\nAaA\n')
            validate_a3m(path, 'AA')
            with self.assertRaises(ValueError):
                validate_a3m(path, 'AG')

    def test_no_msa_search_in_af3_command(self):
        config = {'templates': False, 'samples': 5, 'flash_attention_implementation': 'xla'}
        paths = {'python': 'python', 'runner': 'run_alphafold.py', 'models': 'models', 'databases': 'db'}
        result = command('alphafold3', config, paths, 'input.json', 'output')
        self.assertIn('--norun_data_pipeline', result)
        self.assertIn('--num_diffusion_samples=5', result)

    def test_conversion_and_missing_contacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cif = root / 'model.cif'
            columns = ['group_PDB', 'label_atom_id', 'label_comp_id', 'label_asym_id', 'label_seq_id',
                       'Cartn_x', 'Cartn_y', 'Cartn_z', 'B_iso_or_equiv']
            cif.write_text('data_test\nloop_\n' + ''.join('_atom_site.'+c+'\n' for c in columns) +
                           'ATOM CB ALA H 1 0 0 0 80\nATOM CB ALA H 2 1 0 0 90\nATOM CA GLY A 1 2 0 0 70\n#\n')
            full = root / 'full.json'
            confidence = {'token_chain_ids': ['H', 'H', 'A'], 'pae': np.ones((3, 3)).tolist(),
                          'contact_probs': np.full((3, 3), .5).tolist()}
            full.write_text(json.dumps(confidence))
            summary = root / 'summary.json'
            summary.write_text(json.dumps({'ptm': .5, 'ranking_score': .6, 'chain_pair_iptm': [[1, .4], [.3, 1]]}))
            output = root / 'sample.npz'
            convert_sample(self.candidate(), 'alphafold3', cif, full, summary, output, {'H': 'H', 'A': 'A'}, 100)
            sample = read_sample(output, self.candidate())
            self.assertEqual(sample.native['pair_iptm'], .4)
            self.assertEqual(sample.contact_forward.shape, (2, 1))
            esm_confidence = root / 'confidence.npz'
            np.savez(esm_confidence, plddt=np.array([.8, .9, .7]), pae=np.ones((3, 3)))
            esm_output = root / 'esm_sample.npz'
            convert_sample(self.candidate(), 'esmfold2', cif, esm_confidence, summary,
                           esm_output, {'H': 'H', 'A': 'A'}, 100, 1)
            esm_sample = read_sample(esm_output, self.candidate())
            np.testing.assert_allclose(esm_sample.plddt, [80, 90, 70])
            self.assertIsNone(esm_sample.contact_forward)
            with self.assertRaisesRegex(ValueError, 'pLDDT disagree'):
                convert_sample(self.candidate(), 'esmfold2', cif, esm_confidence, summary,
                               root / 'wrong_units.npz', {'H': 'H', 'A': 'A'}, 100, 100)
            del confidence['contact_probs']
            full.write_text(json.dumps(confidence))
            with self.assertRaises(KeyError):
                convert_sample(self.candidate(), 'alphafold3', cif, full, summary, root/'other.npz', {'H': 'H', 'A': 'A'}, 100)


if __name__ == '__main__':
    unittest.main()
