"""Normalize protein-only predictor outputs with explicit chain mappings."""

import gzip
import json
import shlex
from pathlib import Path

import numpy as np

from abfilter.data import read_sample

AMINO_ACIDS = dict(zip(
    'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
    'ARNDCQEGHILKMFPSTWYV'))


def structure_arrays(path, plddt_scale):
    """Read the protein atom loop; preserve structure chain and residue order."""
    residues, columns = {}, []
    reading = False
    for line in path.read_text().splitlines():
        if line.startswith('_atom_site.'):
            columns.append(line.strip())
            reading = True
            continue
        if not reading:
            continue
        if line.startswith('#') or not line.strip():
            break
        fields = shlex.split(line)
        if len(fields) != len(columns):
            raise ValueError(f'{path}: unsupported atom-site record')
        row = dict(zip((c.removeprefix('_atom_site.') for c in columns), fields))
        if row.get('label_alt_id', '.') not in ('.', '?', 'A', '1'):
            continue
        chain, residue = row['label_asym_id'], int(row['label_seq_id'])
        key = chain, residue
        atoms = residues.setdefault(key, {})
        name = row['label_atom_id']
        if name in atoms:
            raise ValueError(f'{path}: duplicate atom {key}/{name}')
        atoms[name] = row
    coordinates, confidence, chains, symbols = [], [], [], []
    for (chain, residue), atoms in residues.items():
        atom = atoms.get('CB') or atoms.get('CA')
        if atom is None:
            raise ValueError(f'{path}: missing representative atom {chain}/{residue}')
        coordinates.append([float(atom[f'Cartn_{axis}']) for axis in 'xyz'])
        confidence.append(float(atom['B_iso_or_equiv']) * (100 if plddt_scale == 1 else 1))
        chains.append(chain)
        symbols.append(AMINO_ACIDS[atom['label_comp_id']])
    for chain in dict.fromkeys(chains):
        indices = [index for (observed, index) in residues if observed == chain]
        if indices != list(range(1, len(indices)+1)):
            raise ValueError(f'{path}: nonconsecutive residue mapping for chain {chain}')
    if not coordinates:
        raise ValueError(f'{path}: no protein residues')
    return {'coordinates': np.array(coordinates), 'plddt': np.array(confidence),
            'chain_ids': np.array(chains), 'residue_symbols': np.array(symbols)}


def read_json(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        return json.load(handle)


def convert_sample(candidate, predictor, structure, confidence, summary, output,
                   chain_map, structure_plddt_scale, confidence_plddt_scale=100):
    if structure_plddt_scale not in (1, 100) or confidence_plddt_scale not in (1, 100):
        raise ValueError('pLDDT scales must be explicitly 1 or 100')
    if predictor == 'opendde':
        from .opendde import read_prediction
        expected = {**candidate.antibody_sequences, **candidate.antigen_sequences}
        if chain_map != {c:c for c in expected}:
            raise ValueError('OpenDDE lossless output must preserve input chain IDs')
        expected_summary = confidence.with_name(confidence.name.replace('_full_data_', '_summary_confidence_').replace('.npz', '.json'))
        if summary.resolve() != expected_summary.resolve():
            raise ValueError('Summary does not match the lossless confidence sample')
        sample, _, _ = read_prediction(confidence, candidate)
        arrays = dict(coordinates=sample.coordinates, plddt=sample.plddt, pae=sample.pae,
                      chain_ids=sample.chain_ids, contact_forward=sample.contact_forward,
                      contact_reverse=sample.contact_reverse, **sample.native)
        symbols = np.empty(len(sample.chain_ids), dtype='U1')
        for chain, sequence in expected.items():
            symbols[sample.chain_ids == chain] = list(sequence)
        arrays['residue_symbols'] = symbols
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists():
            raise FileExistsError(output)
        temporary = output.with_suffix('.partial.npz')
        np.savez_compressed(temporary, **arrays)
        read_sample(temporary, candidate)
        temporary.replace(output)
        return
    arrays = structure_arrays(structure, structure_plddt_scale)
    raw_chains = arrays['chain_ids'].tolist()
    raw_order = list(dict.fromkeys(raw_chains))
    expected_chains = set(candidate.antibody_sequences) | set(candidate.antigen_sequences)
    if set(chain_map) != set(raw_order) or set(chain_map.values()) != expected_chains or len(set(chain_map.values())) != len(chain_map):
        raise ValueError('chain_map must bijectively map every output chain to a declared input chain')
    if predictor == 'esmfold2':
        with np.load(confidence, allow_pickle=False) as archive:
            full = {key: np.asarray(archive[key]) for key in ('plddt', 'pae')}
        plddt = np.asarray(full['plddt'], dtype=float).reshape(-1) * (100 if confidence_plddt_scale == 1 else 1)
        if not np.allclose(plddt, arrays['plddt'], atol=.2, rtol=0):
            raise ValueError('structure and ESMFold2 confidence pLDDT disagree')
        arrays['plddt'] = plddt
        arrays['pae'] = full['pae']
    else:
        full = read_json(confidence)
        if predictor == 'alphafold3':
            token_chains = full['token_chain_ids']
            arrays['pae'] = np.asarray(full['pae'])
        elif predictor == 'protenix':
            asym = full['token_asym_id']
            order = list(dict.fromkeys(asym))
            if len(order) != len(raw_order):
                raise ValueError('token asymmetry groups do not match structure chains')
            mapping = dict(zip(order, raw_order))
            token_chains = [mapping[x] for x in asym]
            arrays['pae'] = np.asarray(full['token_pair_pae'])
        else:
            raise ValueError(f'Unsupported predictor: {predictor}')
        if list(token_chains) != raw_chains:
            raise ValueError('confidence tokens do not match structure residue order')
        contacts = np.asarray(full['contact_probs'])
        if contacts.shape != arrays['pae'].shape:
            raise ValueError('contact probability dimensions differ from PAE')
        antibody = np.flatnonzero([chain_map[c] in candidate.antibody_sequences for c in raw_chains])
        antigen = np.flatnonzero([chain_map[c] in candidate.antigen_sequences for c in raw_chains])
        arrays['contact_forward'] = contacts[np.ix_(antibody, antigen)]
        arrays['contact_reverse'] = contacts[np.ix_(antigen, antibody)].T
    values = read_json(summary)
    pair = values.get('chain_pair_iptm', values.get('pair_chains_iptm'))
    if pair is None:
        raise ValueError('chain-pair iPTM is required; global iPTM is not substituted')
    pair = np.asarray(pair)
    if pair.shape != (len(raw_order), len(raw_order)):
        raise ValueError('chain-pair iPTM dimensions do not match structure chains')
    arrays['pair_iptm'] = max(pair[i, j] for i, left in enumerate(raw_order)
                             for j, right in enumerate(raw_order)
                             if chain_map[left] in candidate.antibody_sequences
                             and chain_map[right] in candidate.antigen_sequences)
    arrays['ptm'] = float(values['ptm'])
    if predictor != 'esmfold2':
        arrays['ranking_score'] = float(values['ranking_score'])
    if predictor == 'protenix':
        arrays['pair_gpde'] = float(values['pair_gpde'] if 'pair_gpde' in values else values['gpde'])
    arrays['chain_ids'] = np.asarray([chain_map[c] for c in raw_chains])
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    temporary = output.with_suffix('.partial.npz')
    np.savez_compressed(temporary, **arrays)
    read_sample(temporary, candidate)
    temporary.replace(output)
