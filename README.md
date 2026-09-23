# Structure-based antibody filtering

This repository contains the scoring and evaluation code for a study of antibody candidate filtering from predicted antibody-antigen complexes. It compares two selection settings: a shared threshold across antigen targets, and separate ranking within each target. Regional interaction contrast (RIC) summarizes the relative support for predicted CDR contacts across prediction runs. Non-CDR antibody confidence (NAC) summarizes antibody pLDDT outside annotated CDRs. Published confidence baselines are implemented separately in `src/abfilter/baselines/`.

## Installation

Create the environment with `conda env create -f environment.yml`, then run `conda activate abfilter`. Alternatively, install with `python -m pip install -e .` in a Python 3.10+ environment. The package requires NumPy; no training is performed. Run the tests with `python -m unittest discover -s tests`.

## Inputs

`data/candidates.jsonl` contains one record per antibody-antigen candidate. Each record has a stable `pair_id`, `benchmark` (`yeast_display`, `igdesign`, or `literature`), `target_id`, antibody and antigen sequences keyed by chain ID, a Boolean `cdr_mask` in predicted antibody-residue order, `cdr_annotation` with `scheme: "IMGT"`, a binary `label`, and an `evidence` object describing the experimental measurement and units. Labels are used only by the evaluation script. Candidate scoring never reads them.

The `scripts/prepare.py` command builds this manifest from a candidate table, the experimental evidence tables, prediction-input manifests, and IMGT annotations. It checks candidate IDs, labels, chain sequences, and CDR positions before writing. Use `python scripts/prepare.py --help` for its required source paths. The optional `--path-map OLD=NEW` replaces a path prefix when an input manifest points to a moved file. For the yeast-display panel, the available evidence score is an ordinal enrichment grade, not a calibrated affinity measurement or a reconstructed log-fold enrichment value. Literature/patent records retain the reported affinity summaries; IgDesign records retain SPR replicate counts. Paths to source files are supplied at runtime and are not embedded in this repository.

Prediction files use the normalized layout `predictions/<predictor>/<pair_id>/run_000/sample_00.npz`. A run corresponds to an inference seed. Each NPZ contains `coordinates` (C-beta coordinates, or C-alpha for glycine), `plddt` on a 0-100 scale, residue-level `pae`, and per-residue `chain_ids` and one-letter `residue_symbols`. When available, `contact_forward` and `contact_reverse` are antibody-by-antigen matrices of predicted contact probabilities. Optional scalar fields are `pair_iptm`, `pair_gpde`, `ranking_score`, `ptm`, and `actifptm`. The chain sequences and CDR mask must correspond to the same residue order as these arrays; the loader checks this explicitly. All scores require complete ensembles; the paper used ten runs with five samples each.

RIC requires predicted contact probabilities. For a predictor without those outputs, gRIC replaces them with 8-Angstrom geometric contacts; this is an approximation, not the same metric. Predictor-native baselines are omitted when their required scalar is unavailable. The evaluator does not impute missing baselines.

## Usage

```bash
python scripts/prepare.py --help  # Supply the source tables and --output data/candidates.jsonl
python scripts/score.py --candidates data/candidates.jsonl --predictions predictions --predictor af3 --output results/af3.jsonl
python scripts/evaluate.py results/af3.jsonl --output results/evaluation.json
```

The evaluator reports portfolio-wide filtering, within-target filtering on eligible yeast-display and IgDesign targets, and portfolio-wide filtering on that same target-specific population. Target-specific analysis requires at least ten positive and forty negative candidates per target. The two benchmarks have equal total weight; targets within each benchmark are weighted by the square root of their positive-negative information. Scores are evaluated by recall at three false-positive-rate limits, recall at two testing budgets, and a Wilson-bound filtering utility with associated precision. Tied scores are handled as groups.

## Data availability

The candidate metadata, experimental labels, prediction ensembles, and derived score files are not included in this review-stage code repository. They will be released if the paper is accepted, subject to the source datasets' redistribution terms. The repository therefore documents and tests the computation but does not currently reproduce the paper tables without those inputs. Raw predictor outputs also need to be normalized to the NPZ schema described above; this conversion is not yet part of the repository.
