#!/usr/bin/env python3
"""Run local MMseqs2-GPU searches, one worker per GPU, with validated caches."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from abfilter.data import read_candidates
from abfilter.msa import balance, digest, unique_sequences, validate_a3m, validate_cache


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', required=True, type=Path)
    parser.add_argument('--search', required=True, type=Path, help='Installed colabfold_search executable')
    parser.add_argument('--mmseqs', required=True, type=Path)
    parser.add_argument('--mmseqs-sha256', required=True, help='Approved GPU executable checksum')
    parser.add_argument('--databases', required=True, type=Path)
    parser.add_argument('--database-versions', required=True, type=Path, help='JSON mapping database names to version identifiers')
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--work', required=True, type=Path)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--gpus', default='0', help='Comma-separated local GPU identifiers')
    parser.add_argument('--nodes', type=int, default=1)
    parser.add_argument('--node-index', type=int, default=0)
    parser.add_argument('--threads', type=int, default=8, help='CPU threads per GPU worker')
    parser.add_argument('--finalize', action='store_true', help='Validate all cached sequences without launching searches')
    args = parser.parse_args()
    devices = args.gpus.split(',')
    if not 0 <= args.node_index < args.nodes or len(set(devices)) != len(devices) or args.threads < 1:
        raise ValueError('invalid node/GPU/thread assignment')
    for name in ('search', 'mmseqs', 'databases', 'cache', 'work'):
        setattr(args, name, getattr(args, name).resolve())
    if digest(args.mmseqs) != args.mmseqs_sha256:
        raise ValueError('MMseqs executable does not match the approved GPU build')
    version = subprocess.check_output([str(args.mmseqs), 'version'], text=True).strip()
    help_text = subprocess.check_output([str(args.mmseqs), 'search', '-h'], text=True)
    if '--gpu' not in help_text:
        raise ValueError('MMseqs executable has no GPU search option')
    versions = json.loads(args.database_versions.read_text())
    for database in ('uniref30_2302_db', 'colabfold_envdb_202108_db'):
        if not versions.get(database) or not (args.databases / database).exists():
            raise ValueError(f'missing database or version: {database}')
        if not (args.databases / (database + '.index')).exists():
            raise ValueError(f'missing database index: {database}')
    identity = {'mmseqs_sha256': args.mmseqs_sha256, 'mmseqs_version': version,
                'database_versions': versions, 'gpu': 1, 'search_sha256': digest(args.search),
                'templates': False, 'pair_mode': 'unpaired', 'filter': 1}
    args.cache.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    sequences = unique_sequences(read_candidates(args.candidates))
    assignments = balance(sequences, args.nodes * len(devices))

    def worker(local_index):
        index = args.node_index * len(devices) + local_index
        keys = assignments[index]
        missing = [key for key in keys if not validate_cache(args.cache, key, sequences[key], identity)]
        if not missing:
            return
        directory = args.work / f'worker_{index}'
        directory.mkdir(exist_ok=True)
        lock = directory / '.running'
        lock.mkdir()  # Concurrent launches must not duplicate the same assignment.
        try:
            if shutil.disk_usage(directory).free < 10 * 1024**3 or os.statvfs(directory).f_favail < 10000:
                raise ValueError('insufficient free disk space or inodes for MSA search')
            fasta = directory / 'queries.fasta'
            fasta.write_text(''.join(f'>{key}\n{sequences[key]}\n' for key in missing))
            search_dir = directory / hashlib.sha256(fasta.read_bytes()).hexdigest()
            command = [str(args.search), str(fasta), str(args.databases), str(search_dir),
                       '--db1', 'uniref30_2302_db', '--db3', 'colabfold_envdb_202108_db',
                       '--use-env', '1', '--use-templates', '0', '--pair-mode', 'unpaired',
                       '--filter', '1', '--mmseqs', str(args.mmseqs), '--db-load-mode', '2',
                       '--threads', str(args.threads), '--gpu', '1']
            environment = {**os.environ, 'CUDA_VISIBLE_DEVICES': devices[local_index]}
            gpu = subprocess.check_output(['nvidia-smi', '-i', devices[local_index], '--query-gpu=uuid,name', '--format=csv,noheader'], text=True)
            metadata = {**identity, 'command': command, 'node_index': args.node_index,
                        'gpu': devices[local_index], 'gpu_identity': gpu, 'hostname': os.uname().nodename,
                        'query_sha256': digest(fasta), 'queries': missing, 'threads': args.threads}
            (directory / 'run.json').write_text(json.dumps(metadata, indent=2) + '\n')
            observed_gpu_search = False
            with (directory / 'search.log').open('a') as log, (directory / 'gpu.log').open('a') as gpu_log:
                process = subprocess.Popen(command, env=environment, stdout=log, stderr=subprocess.STDOUT)
                while process.poll() is None:
                    snapshot = subprocess.check_output(['nvidia-smi', '-i', devices[local_index]], text=True)
                    gpu_log.write(snapshot + '\n')
                    gpu_log.flush()
                    observed_gpu_search |= 'mmseqs' in snapshot.lower()
                    time.sleep(2)
                if process.returncode:
                    raise RuntimeError(f'MSA search failed; see {directory / "search.log"}')
            if not observed_gpu_search:
                raise RuntimeError('GPU search process was not observed; inspect gpu.log before accepting results')
            for key in missing:
                source = search_dir / f'{key}.a3m'
                validate_a3m(source, sequences[key])
                shutil.copyfile(source, args.cache / f'{key}.a3m')
                (args.cache / f'{key}.json').write_text(json.dumps({'identity': identity, 'sha256': digest(source)}) + '\n')
        finally:
            lock.rmdir()

    if not args.finalize:
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(devices)) as executor:
            list(executor.map(worker, range(len(devices))))
    missing = [key for key, sequence in sequences.items() if not validate_cache(args.cache, key, sequence, identity)]
    if missing:
        if args.finalize or args.nodes == 1:
            raise ValueError(f'MSA manifest incomplete: {len(missing)} sequences missing')
        print(f'Node complete; {len(missing)} sequences await other nodes. Run --finalize after all nodes finish.')
        return
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    record = {'complete': True, 'identity': identity, 'sequences': {
        key: {'path': str(args.cache / f'{key}.a3m'), 'sha256': digest(args.cache / f'{key}.a3m')}
        for key in sorted(sequences)}}
    args.manifest.write_text(json.dumps(record, indent=2) + '\n')
    print(f'Validated {len(sequences)} unique sequence MSAs.')


if __name__ == '__main__':
    main()
