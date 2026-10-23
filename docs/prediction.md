# Recompute from sequences

This workflow requires Linux, CUDA GPUs, local MSA databases, and independently
installed predictor software and weights. Use separate environments for CPU
tools and each predictor. Conversion and scoring have been checked against nine
real complete ensembles. Small GPU inference, conversion, and scoring tests have
also passed for all three models using existing validated MSAs. Fresh local MSA
search validation remains pending.

## Prepare candidates and MSAs

```bash
python scripts/prepare_candidates.py --dataset data/antibody_binding_datasets \
  --output data/candidates.jsonl
python scripts/prepare_msa.py --help
```

The MSA command requires paths to `colabfold_search`, an approved MMseqs2-GPU
executable, and a local database directory. Supply the executable SHA-256 with
`--mmseqs-sha256`. Searches explicitly use `--gpu 1`, one worker per listed GPU,
and must show an MMseqs process in the recorded GPU snapshots. There is no CPU
or remote-server fallback.

The database directory must contain indexed `uniref30_2302_db` and
`colabfold_envdb_202108_db` databases. `--database-versions` supplies a JSON object
mapping those names to installed release identifiers. Searches use environmental
sequences, filtering, unpaired MSAs, and no templates.

Sequences are deduplicated by SHA-256 and balanced by squared length before
assignment. `--gpus 0,1,2,3` uses four workers. `--nodes N --node-index I` partitions
the same manifest across machines; use the same number of GPU workers per node.
Share or consolidate caches before finalization. Once all nodes finish, run the
same command with `--finalize` to validate complete coverage and write the manifest.

Every MSA is checked against its query and aligned length. Cache metadata includes
binary, search, and database identities and output checksums. Logs and assignments
remain under `--work`. After interruption, remove a worker's `.running` directory
only after confirming that worker is no longer active. Do not run keepalive
workloads alongside structure inference.

## Predictor installations

- [AlphaFold 3](https://github.com/google-deepmind/alphafold3): install its Python
  environment, obtain authorized weights, and locate `run_alphafold.py`.
- [Protenix](https://github.com/bytedance/Protenix): install the `protenix` entry
  point and Protenix-v2 checkpoint. The adapter supplies local MSA paths using
  the documented [input format](https://github.com/bytedance/Protenix/blob/main/docs/infer_json_format.md).
- ESMFold2 needs the Biohub ESM and Transformers forks, ESMFold2 and ESMC-6B
  weights, and the CCD file. The inspected source revisions are recorded in
  `reference/predictor_installations.json`; a generic current Transformers
  installation is not a substitute.

### Verify model files

`reference/model_assets.json` records SHA-256 hashes and sizes of the inspected
weights and supporting assets. Check files obtained through each provider's
authorized distribution before inference:

```bash
python scripts/verify_model_assets.py --predictor alphafold3 --models /path/to/af3_weights
python scripts/verify_model_assets.py --predictor protenix --models /path/to/checkpoint
python scripts/verify_model_assets.py --predictor esmfold2 --models /path/to/models
```

The first directory must contain `af3.bin`, the second `protenix-v2.pt`, and the
third the `ESMFold2/` and `ESMC-6B/` subdirectories. The command streams files
from disk, needs no GPU, and exits unsuccessfully for missing or changed files.
Use `--output results/model_check.json` to retain the report. The manifest is an
installation snapshot, not proof that every historical prediction used those
exact files. It does not grant permission to redistribute model weights.

`paper_software_revision` and `paper_weights_sha256` in the predictor configurations
are null where historical verification is incomplete. Null means unknown, not
latest. Do not fill these fields with current defaults and claim a paper-matched run.

## Inference

```bash
python scripts/predict.py --predictor alphafold3 \
  --candidates data/candidates.jsonl --msa-manifest data/msa.json \
  --resources resources.json --gpus 0,1 --output predictions_run --prepare-only
```

`--prepare-only` writes inputs and commands without inference; remove it to run.
Choose `protenix` or `esmfold2` for another predictor. `--nodes` and `--node-index`
partition candidate IDs across machines. Each GPU processes its assigned candidates.
Completion records require all requested outputs to be converted and validated.

`resources.json` is machine-local and must not be committed. Required fields:

| Field | Meaning |
| --- | --- |
| `software_revision` | Exact installed source revision/version |
| `weights_identity` | Checkpoint name and checksum or immutable revision |
| `structure_plddt_scale` | 1 or 100: structure B-factor confidence scale |
| `confidence_plddt_scale` | 1 or 100: NPZ confidence scale; ESMFold2 paper configuration uses 1 |
| `paths` | Local executable and model paths, listed below |
| `environment` | Optional explicit model environment variables |

AlphaFold 3 paths: `python`, `runner`, `models`, `databases`. Protenix paths:
`runner` (installed executable); configure its local model/cache paths through
`environment`. ESMFold2 paths: `python`, `models`, `esmc`, `ccd`. Relative paths
resolve from the working directory.

Configurations retain ten seeds (101-110), five samples each, and no templates.
AlphaFold 3 consumes prepared MSAs without its data pipeline. Protenix uses ten
cycles, 200 diffusion steps, and bf16. ESMFold2 uses 20 loops, 100 sampling steps,
MSA depth 1,024, up to 16,384 sequences, LM dropout 0.3, LM masking 0, and MSA
column masking 0.1; its folding model stays fp32 while ESMC uses bf16. These
parameters come from existing launchers; checkpoint verification is separate.

For a smoke test, use the separate two-seed/one-sample configurations in
`configs/smoke/`, a small candidate manifest, and its matching validated MSA
manifest. Other inference settings are unchanged. These profiles are provided
for installation checks and have passed with existing validated MSAs. They do
not establish full-ensemble reproduction or fresh MSA-search correctness.

```bash
python scripts/predict.py --predictor alphafold3 \
  --config configs/smoke/alphafold3.json \
  --candidates data/smoke_candidates.jsonl --msa-manifest data/smoke_msa.json \
  --resources resources.json --gpus 0 --output smoke_predictions
python scripts/score.py --predictor alphafold3 \
  --candidates data/smoke_candidates.jsonl --predictions smoke_predictions/predictions \
  --runs 2 --samples 1 --output results/smoke_alphafold3.csv
```

Use the corresponding model configuration and resources for Protenix or ESMFold2.
Do not modify the paper configuration or mix smoke-test predictions into paper
ensembles. Smoke-test scores are not expected to match the released scores.

## Conversion and scoring

### OpenDDE

Install the source revision and GPU dependencies in a separate environment:

```bash
git clone https://github.com/aurekaresearch/OpenDDE.git
git -C OpenDDE checkout 9d0d9c975178074c1ce722d2d6eefac8b7057062
python3.11 -m venv opendde-env
opendde-env/bin/python -m pip install './OpenDDE[gpu]'
export OPENDDE_ROOT_DIR=/path/to/opendde_data
bash OpenDDE/scripts/download_opendde_data.sh --skip-search-database
```

Use the general-purpose `opendde.pt`, not the antibody-antigen-specific
checkpoint. `configs/predictors/opendde.json` records the exact source revision,
checkpoint revision, byte count, and SHA-256. The prediction entry point verifies
the source revision and checkpoint hash before starting.

For OpenDDE, `resources.json` uses paths `python` (its environment interpreter),
`source` (the source checkout), and `models` (the runtime data directory containing
`checkpoint/opendde.pt` and `common/`). Set `software_revision` to the configured
revision, `weights_identity` to the checkpoint SHA-256, and
`structure_plddt_scale` to 100. The adapter supplies `OPENDDE_ROOT_DIR` and
the source import path automatically.

```bash
python scripts/predict.py --predictor opendde \
  --candidates data/candidates.jsonl --msa-manifest data/msa.json \
  --resources resources.json --gpus 0,1 --output opendde_predictions
python scripts/score.py --predictor opendde \
  --candidates data/candidates.jsonl --predictions opendde_predictions/predictions \
  --output results/opendde.csv
```

The paper uses ten seeds (101-110), five samples per seed, ten recycling cycles,
200 diffusion steps, bf16, local unpaired MSAs, no templates, and no constraint
guidance. Protein MSAs are supplied by the shared validated MSA manifest; no MSA
search is performed inside inference. The output adapter changes serialization,
not inference: full confidence tensors and atom coordinates are saved losslessly
to NPZ. Sample identities are read from the arrays rather than confidence-ranked
filenames. Pair-gPDE is the minimum antibody-to-antigen entry in the chain-pair
matrix; global gPDE is not substituted. Predicted contact probabilities support
RIC directly. All 4,184 benchmark candidates completed 50-sample ensembles.

### Common output format

Inference writes normalized arrays under
`predictions_run/predictions/<predictor>/<pair_id>/run_000/sample_00.npz`.
For existing outputs, `scripts/convert_predictions.py` takes a JSONL manifest
with one record per sample: `pair_id`, `seed`, `sample`, `structure`, `confidence`,
`summary`, `chain_map`, `structure_plddt_scale`, and optional
`confidence_plddt_scale`. Paths resolve relative to that manifest. `chain_map`
maps every output chain to its candidate chain ID.

The converter validates sequences, residue order, dimensions, chain roles, and
units. Normalized pLDDT is 0-100; PAE is in Angstroms. Coordinates use C-beta, or
C-alpha if C-beta is absent. RIC requires predicted contact probabilities;
ESMFold2 uses geometric-contact gRIC. Raw confidence files remain available,
including any PDE arrays.

```bash
python scripts/score.py --candidates data/candidates.jsonl \
  --predictions predictions_run/predictions --predictor alphafold3 \
  --output results/alphafold3.csv
```

Missing predictions or inconsistent mappings invalidate the candidate; there is
no partial-ensemble averaging. Check recomputed scores against released scores
before using them as an exact paper reproduction.

A `.csv` output contains the exact score columns consumed by
`reproduce_tables.py`. A `.jsonl` output additionally contains candidate metadata
for inspection. The scorer rejects missing required metrics in either format.
