"""Discover complete raw output ensembles without accepting partial predictions."""
import re


def discover(root, predictor, candidate, seeds, samples):
    records = []
    chain_map = {c: c for c in (*candidate.antibody_sequences, *candidate.antigen_sequences)}
    if predictor == 'alphafold3':
        for structure in root.rglob('*_model.cif'):
            match = re.fullmatch(r'seed-(\d+)_sample-(\d+)', structure.parent.name)
            if match:
                stem = structure.name.removesuffix('_model.cif')
                records.append((int(match[1]), int(match[2]), structure,
                                structure.with_name(stem + '_confidences.json'),
                                structure.with_name(stem + '_summary_confidences.json')))
    elif predictor == 'protenix':
        for structure in root.rglob('*_sample_*.cif'):
            match = re.fullmatch(r'seed_(\d+)', structure.parent.parent.name)
            sample = re.search(r'_sample_(\d+)\.cif$', structure.name)
            if match and sample:
                stem = structure.name[:sample.start()]
                index = int(sample[1])
                full = structure.with_name(f'{stem}_full_data_sample_{index}.json')
                if not full.exists():
                    full = full.with_suffix('.json.gz')
                records.append((int(match[1]), index, structure, full,
                                structure.with_name(f'{stem}_summary_confidence_sample_{index}.json')))
    elif predictor == 'esmfold2':
        for seed in seeds:
            for sample in range(samples):
                directory = root / f'seed_{seed}'
                records.append((seed, sample, directory / f'sample_{sample}.cif',
                                directory / f'full_confidence_sample_{sample}.npz',
                                directory / f'summary_confidence_sample_{sample}.json'))
    expected = {(seed, sample) for seed in seeds for sample in range(samples)}
    if len(records) != len(expected) or {(r[0], r[1]) for r in records} != expected:
        raise ValueError(f'{candidate.pair_id}: expected {len(expected)} samples; found {len(records)}')
    result = []
    for seed, sample, structure, confidence, summary in sorted(records):
        if not all(path.is_file() for path in (structure, confidence, summary)):
            raise ValueError(f'{candidate.pair_id}: incomplete raw output for {seed}/{sample}')
        result.append({'pair_id': candidate.pair_id, 'seed': seed, 'sample': sample,
                       'structure': str(structure.resolve()), 'confidence': str(confidence.resolve()),
                       'summary': str(summary.resolve()), 'chain_map': chain_map,
                       'structure_plddt_scale': 100, 'confidence_plddt_scale': 1 if predictor == 'esmfold2' else 100})
    return result
