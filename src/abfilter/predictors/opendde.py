"""Read lossless OpenDDE outputs with explicit residue and chain validation."""
import json
import numpy as np
from abfilter.data import Sample
from abfilter.predictors.convert import AMINO_ACIDS


def read_prediction(path, candidate):
    with np.load(path, allow_pickle=False) as z:
        atom_chain = z["structure__chain_id"].astype(str)
        atom_res = z["structure__res_id"]
        atom_name = z["structure__atom_name"].astype(str)
        atom_symbol = z["structure__res_name"].astype(str)
        mapping = z["atom_to_token_idx"].astype(int)
        pae = np.asarray(z["token_pair_pae"], dtype=float)
        n = len(pae)
        if pae.shape != (n, n) or set(mapping.tolist()) != set(range(n)):
            raise ValueError("Invalid residue/token mapping")
        representatives = []
        for token in range(n):
            ix = np.flatnonzero(mapping == token)
            if len(set(zip(atom_chain[ix], atom_res[ix]))) != 1:
                raise ValueError("Token crosses structure residues")
            cb, ca = ix[atom_name[ix] == "CB"], ix[atom_name[ix] == "CA"]
            chosen = cb if len(cb) else ca
            if len(chosen) != 1:
                raise ValueError("Missing or ambiguous representative atom")
            representatives.append(int(chosen[0]))
        ix = np.asarray(representatives)
        chains = atom_chain[ix]
        symbols = np.array([AMINO_ACIDS[x] for x in atom_symbol[ix]])
        sequences = {**candidate.antibody_sequences, **candidate.antigen_sequences}
        if set(chains) != set(sequences):
            raise ValueError("Chain roles differ from candidate")
        for chain, sequence in sequences.items():
            if "".join(symbols[chains == chain]) != sequence:
                raise ValueError(f"Sequence mismatch: {chain}")
            if atom_res[ix][chains == chain].tolist() != list(range(1, len(sequence)+1)):
                raise ValueError("Residue numbering mismatch")
        ab = np.flatnonzero(np.isin(chains, list(candidate.antibody_sequences)))
        ag = np.flatnonzero(np.isin(chains, list(candidate.antigen_sequences)))
        if "".join(symbols[ab]) != "".join(candidate.antibody_sequences.values()):
            raise ValueError("CDR annotation order mismatch")
        asym = z["token_asym_id"]
        order = list(dict.fromkeys(chains.tolist()))
        if len(set(asym.tolist())) != len(order):
            raise ValueError("Confidence chain groups mismatch")
        for chain in order:
            if len(set(asym[chains == chain].tolist())) != 1:
                raise ValueError("Confidence chain order mismatch")
        coords = np.asarray(z["coordinate"][ix], dtype=float)
        plddt = np.asarray(z["atom_plddt"][ix], dtype=float) * 100
        contacts = np.asarray(z["contact_probs"], dtype=float)
        if contacts.shape != pae.shape or np.min(contacts) < -1e-6 or np.max(contacts) > 1+1e-6:
            raise ValueError("Invalid contact probabilities")
        # Float32 softmax summation may exceed one by one ULP.
        contacts = np.clip(contacts, 0, 1)
        if any(not np.isfinite(x).all() for x in (coords, plddt, pae, contacts)):
            raise ValueError("Non-finite prediction values")
        if np.min(plddt) < 0 or np.max(plddt) > 100 or np.min(pae) < 0:
            raise ValueError("Invalid confidence units")
        seed = int(z["seed"])
        sample_index = int(z["sample_index"])
    summary = json.loads(path.with_name(path.name.replace("_full_data_", "_summary_confidence_").replace(".npz", ".json")).read_text())
    pairs = [(i,j) for i,c in enumerate(order) for j,d in enumerate(order)
             if c in candidate.antibody_sequences and d in candidate.antigen_sequences]
    native = {"ptm": float(summary["ptm"]), "ranking_score": float(summary["ranking_score"])}
    for key, field, reduction in (("pair_iptm", "chain_pair_iptm", max),
                                  ("pair_gpde", "chain_pair_gpde", min)):
        matrix = np.asarray(summary[field], dtype=float)
        if matrix.shape != (len(order), len(order)) or not np.isfinite(matrix).all():
            raise ValueError(f"Invalid {field}")
        native[key] = float(reduction(matrix[i,j] for i,j in pairs))
    return Sample(coords, plddt, pae, chains, ab, ag, np.array(candidate.cdr_mask),
                  contacts[np.ix_(ab,ag)], contacts[np.ix_(ag,ab)].T, native), seed, sample_index
