"""Sequence-hash MSA caches and deterministic length-balanced assignments."""
import hashlib
import json
from pathlib import Path


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def sequence_id(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def unique_sequences(candidates):
    return {sequence_id(sequence): sequence for candidate in candidates
            for sequence in [*candidate.antibody_sequences.values(), *candidate.antigen_sequences.values()]}


def balance(sequences, workers):
    if workers < 1:
        raise ValueError('worker count must be positive')
    assignments, costs = [[] for _ in range(workers)], [0] * workers
    for key in sorted(sequences, key=lambda k: (-len(sequences[k])**2, k)):
        worker = min(range(workers), key=lambda i: (costs[i], i))
        assignments[worker].append(key)
        costs[worker] += len(sequences[key])**2
    return assignments


def validate_a3m(path, sequence):
    entries, current = [], []
    for line in path.read_text().splitlines():
        if line.startswith('#') or not line:
            continue
        if line.startswith('>'):
            if current:
                entries.append(''.join(current))
            current = []
        else:
            current.append(''.join(c for c in line if not c.islower() and c != '.'))
    if current:
        entries.append(''.join(current))
    if not entries or entries[0] != sequence or any(len(row) != len(sequence) for row in entries):
        raise ValueError(f'{path}: MSA query or aligned lengths do not match')


def validate_cache(cache, key, sequence, identity):
    path, metadata = cache / f'{key}.a3m', cache / f'{key}.json'
    if not path.exists() and not metadata.exists():
        return False
    if not path.exists() or not metadata.exists():
        raise ValueError(f'incomplete cache entry: {key}')
    record = json.loads(metadata.read_text())
    if record['identity'] != identity or record['sha256'] != digest(path):
        raise ValueError(f'cache configuration or checksum mismatch: {key}')
    validate_a3m(path, sequence)
    return True


def load_manifest(path, candidates):
    record = json.loads(path.read_text())
    sequences = unique_sequences(candidates)
    if record.get('complete') is not True or set(record['sequences']) != set(sequences):
        raise ValueError('MSA manifest is incomplete or belongs to different candidates')
    result = {}
    for key, sequence in sequences.items():
        item = record['sequences'][key]
        msa = Path(item['path'])
        if not msa.is_absolute():
            msa = path.parent / msa
        if digest(msa) != item['sha256']:
            raise ValueError(f'MSA checksum mismatch: {key}')
        validate_a3m(msa, sequence)
        result[sequence] = msa.resolve()
    return result
