"""Write model inputs using validated local, unpaired MSAs."""
import json
import os
from pathlib import Path


def resource_paths(values):
    # Resolving a venv Python symlink bypasses its pyvenv.cfg and installed packages.
    return {key: Path(os.path.abspath(value)) for key, value in values.items()}


def prepare(candidate, predictor, msa, seeds):
    chains = {**candidate.antibody_sequences, **candidate.antigen_sequences}
    if predictor == 'alphafold3':
        return {'name': candidate.pair_id, 'dialect': 'alphafold3', 'version': 1,
                'modelSeeds': seeds, 'sequences': [
                    {'protein': {'id': chain, 'sequence': sequence,
                                 'unpairedMsa': msa[sequence].read_text(),
                                 'pairedMsa': '', 'templates': []}}
                    for chain, sequence in chains.items()]}
    if predictor in ('protenix', 'opendde'):
        return [{'name': candidate.pair_id, 'sequences': [
            {'proteinChain': {'id': [chain], 'count': 1, 'sequence': sequence,
                              'unpairedMsaPath': str(msa[sequence]), 'pairedMsaPath': None}}
            for chain, sequence in chains.items()]}]
    if predictor == 'esmfold2':
        return {'name': candidate.pair_id, 'chains': [
            {'id': chain, 'sequence': sequence, 'msa': str(msa[sequence])}
            for chain, sequence in chains.items()]}
    raise ValueError(f'Unsupported predictor: {predictor}')


def command(predictor, config, paths, input_path, output):
    if config['templates'] is not False:
        raise ValueError('The paper configuration disables templates')
    if predictor == 'opendde':
        runner = Path(__file__).resolve().parents[3] / 'scripts/run_opendde.py'
        return [str(paths['python']), str(runner), '--model_name', config['model_name'],
                '--input_json_path', str(input_path), '--dump_dir', str(output),
                '--seeds', ','.join(map(str, config['seeds'])),
                '--model.N_cycle', str(config['cycles']),
                '--sample_diffusion.N_sample', str(config['samples']),
                '--sample_diffusion.N_step', str(config['diffusion_steps']),
                '--dtype', config['dtype'], '--device', 'cuda', '--use_msa', 'true',
                '--use_template', 'false', '--use_rna_msa', 'false', '--use_tfg_guidance', 'false',
                '--need_atom_confidence', 'true', '--triangle_attention', 'cuequivariance',
                '--triangle_multiplicative', 'cuequivariance']
    if predictor == 'alphafold3':
        return [str(paths['python']), str(paths['runner']), f'--json_path={input_path}',
                f'--output_dir={output}', f'--model_dir={paths["models"]}',
                f'--db_dir={paths["databases"]}', '--run_inference=true', '--norun_data_pipeline',
                f'--num_diffusion_samples={config["samples"]}',
                f'--flash_attention_implementation={config["flash_attention_implementation"]}']
    if predictor == 'protenix':
        return [str(paths['runner']), 'pred', '-i', str(input_path), '-o', str(output),
                '-s', ','.join(map(str, config['seeds'])), '-e', str(config['samples']),
                '-c', str(config['cycles']), '-p', str(config['diffusion_steps']),
                '-n', config['model_name'], '--use_default_params', 'true', '--use_msa', 'true',
                '--use_template', 'false', '--dtype', config['dtype'], '--enable_cache', 'true',
                '--need_atom_confidence', 'true']
    raise ValueError('ESMFold2 uses the included Python runner')
