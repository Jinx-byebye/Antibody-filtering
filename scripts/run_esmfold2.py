#!/usr/bin/env python3
"""ESMFold2 inference adapter for a local installation and local model weights."""
import argparse
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--models', type=Path, required=True)
    parser.add_argument('--esmc', type=Path, required=True)
    parser.add_argument('--ccd', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    os.environ['HF_HUB_OFFLINE'] = '1'
    os.environ['TRANSFORMERS_OFFLINE'] = '1'
    os.environ['ESMCFOLD_CCD_PATH'] = str(args.ccd.resolve())
    import numpy as np
    import torch
    from esm.models.esmfold2 import ESMFold2InputBuilder, MSA, ProteinInput, StructurePredictionInput
    from transformers.models.esmfold2.configuration_esmfold2 import ESMFold2Config
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    if not torch.cuda.is_available():
        raise RuntimeError('ESMFold2 requires a CUDA GPU')
    parameters = json.loads(args.config.read_text())
    if parameters['precision'] != 'mixed' or parameters['kernel_backend'] != 'none':
        raise ValueError('This adapter supports mixed precision and kernel_backend=none only')
    config = ESMFold2Config.from_pretrained(args.models, local_files_only=True)
    config.esmc_id = str(args.esmc.resolve())
    model = ESMFold2Model.from_pretrained(args.models, config=config, local_files_only=True,
                                        dtype=torch.float32, esmc_precision='bf16').to('cuda').eval()
    task = json.loads(args.input.read_text())
    proteins = []
    for chain in task['chains']:
        msa = MSA.from_a3m(chain['msa'], remove_insertions=True,
                           max_sequences=parameters['msa_max_sequences'])
        if not msa.entries or msa.sequences[0].replace('-', '').upper() != chain['sequence']:
            raise ValueError('ESMFold2 MSA query does not match the input sequence')
        proteins.append(ProteinInput(id=[chain['id']], sequence=chain['sequence'], modifications=None, msa=msa))
    request = StructurePredictionInput(sequences=proteins)
    builder = ESMFold2InputBuilder()
    def array(value):
        return value.detach().cpu().float().numpy() if hasattr(value, 'detach') else np.asarray(value)
    for seed in parameters['seeds']:
        directory = args.output / f'seed_{seed}'
        directory.mkdir(parents=True, exist_ok=True)
        results = builder.fold(model, request, num_loops=parameters['loops'],
                               num_sampling_steps=parameters['diffusion_steps'],
                               num_diffusion_samples=parameters['samples'], seed=seed,
                               lm_mask_pct=parameters['lm_mask_pct'], early_exit=False,
                               lm_dropout=parameters['lm_dropout'], msa_max_depth=parameters['msa_max_depth'],
                               msa_column_mask_rate=parameters['msa_column_mask_rate'], complex_id=task['name'])
        if parameters['samples'] == 1 and not isinstance(results, list):
            results = [results]
        if len(results) != parameters['samples']:
            raise ValueError('ESMFold2 returned an unexpected number of samples')
        for sample, result in enumerate(results):
            plddt, pae = array(result.plddt), array(result.pae)
            if not np.isfinite(plddt).all() or not np.isfinite(pae).all():
                raise ValueError('non-finite ESMFold2 confidence outputs')
            (directory / f'sample_{sample}.cif').write_text(result.complex.to_mmcif())
            np.savez_compressed(directory / f'full_confidence_sample_{sample}.npz', plddt=plddt,
                                pae=pae, pair_chains_iptm=array(result.pair_chains_iptm))
            summary = {'ptm': float(array(result.ptm)), 'iptm': float(array(result.iptm)),
                       'pair_chains_iptm': array(result.pair_chains_iptm).tolist()}
            (directory / f'summary_confidence_sample_{sample}.json').write_text(json.dumps(summary) + '\n')


if __name__ == '__main__':
    main()
