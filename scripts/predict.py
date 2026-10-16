#!/usr/bin/env python3
"""Run local prediction ensembles with explicit resources and durable completion records."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from abfilter.data import read_candidates, read_ensemble
from abfilter.msa import digest, load_manifest
from abfilter.predictors.inputs import prepare, command, resource_paths
from abfilter.predictors.outputs import discover
from abfilter.predictors.convert import convert_sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictor', choices=['alphafold3', 'protenix', 'esmfold2'], required=True)
    parser.add_argument('--candidates', type=Path, required=True)
    parser.add_argument('--msa-manifest', type=Path, required=True)
    parser.add_argument('--config', type=Path)
    parser.add_argument('--resources', type=Path, required=True, help='Local executable/model paths and version identifiers')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--gpus', default='0')
    parser.add_argument('--node-index', type=int, default=0)
    parser.add_argument('--nodes', type=int, default=1)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    config_path = args.config or ROOT / 'configs/predictors' / f'{args.predictor}.json'
    config = json.loads(config_path.read_text())
    resources = json.loads(args.resources.read_text())
    if config['predictor'] != args.predictor or not 0 <= args.node_index < args.nodes:
        raise ValueError('predictor configuration or node assignment mismatch')
    seeds, samples = config['seeds'], config['samples']
    if len(seeds) < 2 or len(set(seeds)) != len(seeds) or samples < 1:
        raise ValueError('at least two unique seeds and one sample are required')
    for key in ('software_revision', 'weights_identity'):
        if not resources.get(key):
            raise ValueError(f'resources must declare {key}; installed defaults are not a version record')
    paths = resource_paths(resources['paths'])
    required = {'alphafold3': {'python', 'runner', 'models', 'databases'},
                'protenix': {'runner'}, 'esmfold2': {'python', 'models', 'esmc', 'ccd'}}[args.predictor]
    if required - set(paths):
        raise ValueError(f'Missing resource paths: {sorted(required - set(paths))}')
    if resources.get('structure_plddt_scale') not in (1, 100):
        raise ValueError('resources must specify structure_plddt_scale as 1 or 100')
    for key, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f'Missing {key}: {path}')
    candidates = read_candidates(args.candidates)
    msa = load_manifest(args.msa_manifest, candidates)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    identity = {'config': config, 'resources': resources, 'candidates_sha256': digest(args.candidates),
                'msa_manifest_sha256': digest(args.msa_manifest)}
    signature = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    devices = args.gpus.split(',')
    if not devices or len(set(devices)) != len(devices):
        raise ValueError('GPU identifiers must be distinct')
    selected = [c for c in candidates if int(hashlib.sha256(c.pair_id.encode()).hexdigest(), 16) % args.nodes == args.node_index]
    assignments = [[] for _ in devices]
    costs = [0] * len(devices)
    for candidate in sorted(selected, key=lambda c: -sum(map(len, [*c.antibody_sequences.values(), *c.antigen_sequences.values()]))**2):
        worker = min(range(len(devices)), key=lambda i: costs[i])
        assignments[worker].append(candidate)
        costs[worker] += sum(map(len, [*candidate.antibody_sequences.values(), *candidate.antigen_sequences.values()]))**2
    (output / f'assignment_{args.node_index}.json').write_text(json.dumps({devices[i]: [c.pair_id for c in group] for i, group in enumerate(assignments)}, indent=2) + '\n')

    def worker(index):
        for candidate in assignments[index]:
            work = output / 'raw' / candidate.pair_id
            work.mkdir(parents=True, exist_ok=True)
            done = work / 'complete.json'
            if done.exists():
                previous = json.loads(done.read_text())
                if previous['signature'] != signature:
                    raise ValueError(f'{candidate.pair_id}: previous configuration differs')
                read_ensemble(output / 'predictions', args.predictor, candidate, len(seeds), samples)
                continue
            lock = work / '.running'
            lock.mkdir()
            try:
                previous_run = work / 'run.json'
                if previous_run.exists() and json.loads(previous_run.read_text())['signature'] != signature:
                    raise ValueError(f'{candidate.pair_id}: cannot resume outputs from a different configuration')
                if (work / 'output').exists() and not previous_run.exists():
                    raise ValueError(f'{candidate.pair_id}: existing raw predictions lack run metadata')
                input_path = work / 'input.json'
                input_path.write_text(json.dumps(prepare(candidate, args.predictor, msa, seeds)) + '\n')
                if args.predictor == 'esmfold2':
                    invocation = [str(paths['python']), str(ROOT / 'scripts/run_esmfold2.py'), '--input', str(input_path),
                                  '--config', str(config_path.resolve()), '--models', str(paths['models']),
                                  '--esmc', str(paths['esmc']), '--ccd', str(paths['ccd']), '--output', str(work / 'output')]
                else:
                    invocation = command(args.predictor, config, paths, input_path, work / 'output')
                environment = {**os.environ, **resources.get('environment', {}), 'CUDA_VISIBLE_DEVICES': devices[index]}
                record = {**identity, 'signature': signature, 'command': invocation, 'gpu': devices[index], 'hostname': os.uname().nodename}
                (work / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
                if args.prepare_only:
                    continue
                subprocess.run(['nvidia-smi', '-i', devices[index]], check=True, stdout=subprocess.DEVNULL)
                try:
                    records = discover(work / 'output', args.predictor, candidate, seeds, samples)
                except ValueError:
                    with (work / 'inference.log').open('a') as log:
                        subprocess.run(invocation, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
                    records = discover(work / 'output', args.predictor, candidate, seeds, samples)
                for row in records:
                    row['structure_plddt_scale'] = resources['structure_plddt_scale']
                    row['confidence_plddt_scale'] = resources.get('confidence_plddt_scale', config.get('confidence_plddt_scale', 100))
                    target = output / 'predictions' / args.predictor / candidate.pair_id / f"run_{seeds.index(row['seed']):03d}" / f"sample_{row['sample']:02d}.npz"
                    if not target.exists():
                        convert_sample(candidate, args.predictor, *(Path(row[k]) for k in ('structure', 'confidence', 'summary')),
                                       target, row['chain_map'], row['structure_plddt_scale'], row['confidence_plddt_scale'])
                read_ensemble(output / 'predictions', args.predictor, candidate, len(seeds), samples)
                (work / 'predictions.jsonl').write_text(''.join(json.dumps(row) + '\n' for row in records))
                done.write_text(json.dumps({'signature': signature, 'samples': len(records)}) + '\n')
            finally:
                lock.rmdir()
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(devices)) as executor:
        list(executor.map(worker, range(len(devices))))
    print(f'{"Prepared" if args.prepare_only else "Validated"} {len(selected)} candidate ensembles on node {args.node_index}.')


if __name__ == '__main__':
    main()
