"""Reproduce the paper tables from explicit candidate populations."""

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from .dataset import load_dataset
from .evaluate import evaluate_group, _rank_results

NAMES = {
    'predictor': 'Predictor', 'metric': 'Metric', 'alphafold3': 'AlphaFold 3',
    'protenix': 'Protenix-v2', 'esmfold2': 'ESMFold2', 'pair_iptm': 'Pair-iPTM',
    'pair_gpde': 'Pair-gPDE', 'interface_pae': 'Interface PAE',
    'antibody_plddt': 'Antibody pLDDT', 'ranking_score': 'Ranking score',
    'pdockq': 'pDockQ', 'pdockq2': 'pDockQ2', 'anticonf': 'AntiConf',
    'ipsae_d0res': 'ipSAE-d0res', 'ipsae_d0dom': 'ipSAE-d0dom', 'ipsae_d0chn': 'ipSAE-d0chn',
    'ric': 'RIC', 'gric': 'gRIC', 'nac': 'NAC', 'utility': 'Utility',
    'recall_fpr_10': 'Recall at 10% FPR', 'recall_fpr_20': 'Recall at 20% FPR',
    'recall_fpr_30': 'Recall at 30% FPR', 'recall_q_05': 'Recall at Q=5%',
    'recall_q_20': 'Recall at Q=20%', 'recall_at_utility': 'Recall at maximum utility',
    'precision_at_utility': 'Precision at maximum utility', 'average_rank': 'Average rank',
}


def load_scores(path, metrics, expected):
    rows = {}
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ['pair_id'] + metrics:
            raise ValueError(f'{path}: expected columns pair_id followed by {metrics}')
        for row in reader:
            pair_id = row.pop('pair_id')
            if pair_id in rows:
                raise ValueError(f'{path}: duplicate pair_id {pair_id}')
            values = {key: float(value) for key, value in row.items()}
            if not all(math.isfinite(value) for value in values.values()):
                raise ValueError(f'{path}: non-finite score for {pair_id}')
            rows[pair_id] = values
    if set(rows) != set(expected):
        raise ValueError(f'{path}: population mismatch: {len(set(expected)-set(rows))} missing, '
                         f'{len(set(rows)-set(expected))} unexpected candidates')
    return rows


def reproduce(dataset, scores, populations, config):
    if not math.isclose(sum(config['dataset_weights'].values()), 1.0) or any(v <= 0 for v in config['dataset_weights'].values()):
        raise ValueError('dataset weights must be positive and sum to one')
    candidates = load_dataset(dataset)
    tables = {name: {} for name in ('table2', 'table3', 'table4')}
    coverage = {}
    for predictor, metrics in config['predictor_metrics'].items():
        population = json.loads((populations / f'{predictor}.json').read_text())
        ids = population['included']
        excluded = population['excluded']
        if (len(ids) != len(set(ids)) or set(ids) & set(excluded)
                or set(ids) | set(excluded) != set(candidates)
                or any(not reason for reason in excluded.values())):
            raise ValueError(f'{predictor}: invalid population manifest')
        values = load_scores(scores / f'{predictor}.csv', metrics, ids)
        groups = defaultdict(list)
        for pair_id in ids:
            row = candidates[pair_id]
            if row['dataset'] in config['within_antigen_datasets']:
                groups[(row['dataset'], row['target_id'])].append(pair_id)
        groups = {key: group for key, group in groups.items()
                  if sum(candidates[p]['label'] for p in group) >= config['minimum_positive']
                  and sum(1-candidates[p]['label'] for p in group) >= config['minimum_negative']}
        info = {}
        for key, group in groups.items():
            positive = sum(candidates[p]['label'] for p in group)
            info[key] = math.sqrt(2 * positive * (len(group)-positive) / len(group))
        if set(key[0] for key in groups) != set(config['dataset_weights']):
            raise ValueError(f'{predictor}: missing eligible dataset')
        weights = {key: config['dataset_weights'][key[0]] * value /
                   sum(v for k, v in info.items() if k[0] == key[0]) for key, value in info.items()}
        matched = sorted(p for group in groups.values() for p in group)
        if 'within_antigen' in population and matched != sorted(population['within_antigen']):
            raise ValueError(f'{predictor}: within-antigen population changed')
        coverage[predictor] = {'across_antigen': len(ids), 'within_antigen': len(matched),
                               'antigens': len(groups), 'within_antigen_ids': matched}
        for name in tables:
            tables[name][predictor] = {}
        for metric in metrics:
            def summarize(selected):
                direction = config['score_directions'][metric]
                if direction not in (-1, 1):
                    raise ValueError(f'Invalid score direction: {metric}')
                return evaluate_group(np.array([values[p][metric] * direction for p in selected]),
                                      np.array([candidates[p]['label'] for p in selected]),
                                      f_safe=config['f_safe'], f_severe=config['f_severe'])
            tables['table2'][predictor][metric] = summarize(ids)
            tables['table4'][predictor][metric] = summarize(matched)
            summaries = {key: summarize(group) for key, group in groups.items()}
            tables['table3'][predictor][metric] = {
                endpoint: sum(weights[key] * row[endpoint] for key, row in summaries.items())
                for endpoint in config['endpoints']['table3']}
        for name in tables:
            _rank_results(tables[name][predictor], config['endpoints'][name])
    return tables, coverage


def display_value(endpoint, value):
    if endpoint == 'utility':
        return f'{value:.4f}'
    if endpoint == 'average_rank':
        return f'{value:.2f}'
    return f'{100 * value:.2f}'


def write_tables(tables, config, output, reference=None):
    output.mkdir(parents=True, exist_ok=True)
    mismatches = []
    for name, predictors in tables.items():
        columns = ['predictor', 'metric'] + config['endpoints'][name] + ['average_rank']
        rows = []
        for predictor, metrics in predictors.items():
            for metric, values in metrics.items():
                rows.append({'predictor': predictor, 'metric': metric, **{
                    key: display_value(key, values[key]) for key in columns[2:]}})
        with (output / f'{name}.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        lines = [' | '.join(NAMES[key] for key in columns), ' | '.join(['---'] * len(columns))]
        lines += [' | '.join(NAMES.get(row[key], row[key]) for key in columns) for row in rows]
        (output / f'{name}.md').write_text('\n'.join(lines) + '\n')
        escape = lambda text: text.replace('_', r'\_').replace('%', r'\%')
        tex = [r'\begin{tabular}{ll' + 'r' * (len(columns)-2) + '}', r'\hline']
        tex += [' & '.join(escape(NAMES[key]) for key in columns) + r' \\', r'\hline']
        tex += [' & '.join(escape(NAMES.get(row[key], row[key])) for key in columns) + r' \\' for row in rows]
        tex += [r'\hline', r'\end{tabular}']
        (output / f'{name}.tex').write_text('\n'.join(tex) + '\n')
        if reference is not None:
            with (reference / f'{name}.csv').open(newline='') as handle:
                reader = csv.DictReader(handle)
                expected = list(reader)
                if reader.fieldnames != columns:
                    raise ValueError(f'{name}: reference columns do not match')
            expected = {(r['predictor'], r['metric']): r for r in expected}
            actual = {(r['predictor'], r['metric']): r for r in rows}
            if set(actual) != set(expected):
                raise ValueError(f'{name}: reference metric coverage does not match')
            for key, row in actual.items():
                for column in columns[2:]:
                    if row[column] != expected[key][column]:
                        mismatches.append(f'{name}/{key}/{column}: {row[column]} != {expected[key][column]}')
    if mismatches:
        raise ValueError('Reference differences:\n' + '\n'.join(mismatches))
