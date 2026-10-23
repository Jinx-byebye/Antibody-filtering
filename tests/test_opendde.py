"""Synthetic tests for lossless OpenDDE output and ensemble identities."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from abfilter.data import Candidate, read_sample
from abfilter.predictors.opendde import read_prediction
from abfilter.predictors.outputs import discover
from abfilter.predictors.convert import convert_sample
from abfilter.predictors.inputs import command


class OpenDDETests(unittest.TestCase):
    def test_lossless_conversion_and_sample_identity(self):
        candidate=Candidate.from_record(dict(pair_id='example',benchmark='igdesign',target_id='test',
            antibody_sequences={'H':'AA'},antigen_sequences={'A':'G'},cdr_mask=[True,False],
            cdr_annotation={'scheme':'IMGT'},label=1,evidence={'type':'synthetic'}))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            full=root/'example_full_data_sample_3.npz'
            summary=root/'example_summary_confidence_sample_3.json'
            structure=root/'example_sample_3.cif'
            structure.touch()
            np.savez(full,structure__chain_id=['H','H','A'],structure__res_id=[1,2,1],
                structure__atom_name=['CB','CB','CA'],structure__res_name=['ALA','ALA','GLY'],
                atom_to_token_idx=[0,1,2],token_asym_id=[0,0,1],token_pair_pae=np.ones((3,3)),
                coordinate=[[0,0,0],[1,0,0],[2,0,0]],atom_plddt=[.8123456,.9,.7],
                contact_probs=np.full((3,3),.5),seed=101,sample_index=0)
            summary.write_text(json.dumps(dict(ptm=.5,ranking_score=.6,
                chain_pair_iptm=[[1,.4],[.3,1]],chain_pair_gpde=[[0,2],[3,0]],gpde=99)))
            sample,seed,index=read_prediction(full,candidate)
            self.assertEqual((seed,index),(101,0))
            self.assertAlmostEqual(sample.plddt[0],81.23456)
            self.assertEqual(sample.native['pair_gpde'],2)
            records=discover(root,'opendde',candidate,[101],1)
            self.assertEqual(records[0]['sample'],0)
            self.assertEqual(Path(records[0]['structure']),structure.resolve())
            output=root/'normalized.npz'
            convert_sample(candidate,'opendde',structure,full,summary,output,{'H':'H','A':'A'},100)
            normalized=read_sample(output,candidate)
            np.testing.assert_array_equal(normalized.plddt,sample.plddt)
            summary.write_text(json.dumps(dict(ptm=.5,ranking_score=.6,
                chain_pair_iptm=[[1,.4],[.3,1]],gpde=99)))
            with self.assertRaises(KeyError):
                read_prediction(full,candidate)

    def test_command_preserves_budget(self):
        root=Path(__file__).resolve().parents[1]
        config=json.loads((root/'configs/predictors/opendde.json').read_text())
        args=command('opendde',config,{'python':Path('/python')},Path('/input'),Path('/output'))
        self.assertEqual(args[args.index('--seeds')+1],','.join(map(str,range(101,111))))
        self.assertEqual(args[args.index('--sample_diffusion.N_sample')+1],'5')
        self.assertEqual(args[args.index('--sample_diffusion.N_step')+1],'200')
        self.assertTrue(Path(args[1]).is_file())


if __name__=='__main__':
    unittest.main()
