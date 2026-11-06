#!/usr/bin/env python3
"""Evaluate an external single-antigen cohort with the main-table endpoints."""
import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import numpy as np
from abfilter.data import read_candidates
from abfilter.evaluate import evaluate_group, _rank_results
from abfilter.tables import load_scores, write_tables


def evaluate_cohort(candidates, values, metrics, config):
    results = {}
    ids = [c.pair_id for c in candidates]
    labels = np.asarray([c.label for c in candidates])
    for metric in metrics:
        scores = np.asarray([values[p][metric] for p in ids]) * config['score_directions'][metric]
        results[metric] = evaluate_group(scores, labels,
            f_safe=config['f_safe'], f_severe=config['f_severe'])
    _rank_results(results, config['endpoints']['table3'])
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates', type=Path, required=True)
    p.add_argument('--scores', type=Path, required=True)
    p.add_argument('--predictors', nargs='+', required=True)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    candidates = read_candidates(a.candidates)
    if len({c.target_id for c in candidates}) != 1 or len({c.benchmark for c in candidates}) != 1:
        raise ValueError('This external evaluation requires exactly one antigen and benchmark')
    if len({c.pair_id for c in candidates}) != len(candidates):
        raise ValueError('Duplicate candidate')
    config = json.loads(a.config.read_text())
    if config['endpoints']['table2'] != config['endpoints']['table3']:
        raise ValueError('Main-table endpoint sets differ')
    groups = {'all': candidates}
    for name in sorted({c.evidence.get('group') for c in candidates} - {None}):
        groups[name] = [c for c in candidates if c.evidence.get('group') == name]
    tables, full, coverage = {}, {}, {}
    for group, subset in groups.items():
        if (sum(c.label for c in subset) < config['minimum_positive']
                or sum(1-c.label for c in subset) < config['minimum_negative']):
            raise ValueError(f'Insufficient class counts in {group}')
        tables[group], full[group] = {}, {}
        coverage[group] = dict(candidates=len(subset), positive=sum(c.label for c in subset),
                               negative=sum(1-c.label for c in subset), antigens=1)
    for predictor in a.predictors:
        metrics = config['predictor_metrics'][predictor]
        values = load_scores(a.scores / predictor / f'{predictor}.csv', metrics,
                             [c.pair_id for c in candidates])
        for group, subset in groups.items():
            full[group][predictor] = evaluate_cohort(subset, values, metrics, config)
            tables[group][predictor] = full[group][predictor]
    a.output.mkdir(parents=True, exist_ok=True)
    table_config = {**config, 'endpoints': {g:config['endpoints']['table3'] for g in groups}}
    write_tables(tables, table_config, a.output)
    (a.output / 'full_results.json').write_text(json.dumps(full, indent=2) + '\n')
    (a.output / 'coverage.json').write_text(json.dumps(coverage, indent=2) + '\n')
    (a.output / 'evaluation.json').write_text(json.dumps(config, indent=2) + '\n')
    print(json.dumps(coverage), flush=True)


if __name__ == '__main__':
    main()
