# Data and score files

Two separately distributed archives are required for CPU reproduction. Their
public download addresses and licenses will be added when the release is published.
No download link is available yet.

Archive checksums are in `reference/data.sha256` and `reference/scores.sha256`.
Run `sha256sum -c` with those files from the directory containing the downloaded
archives (on macOS, use `shasum -a 256 -c`). After extraction, each archive also
contains a `SHA256SUMS` file for checking individual files.

The dataset archive contains `display_enrichment/`, `igdesign_spr/`, and
`literature_affinity/`. Each folder has `sequences.csv`, `metadata.jsonl`, and a
README. The archive's `FIELDS.md` defines every experimental metadata field and
unit. Read the assay-specific metadata before treating display enrichment as an
affinity measurement. `pair_id` is the stable join key, not a row number.

The score archive contains:

```text
antibody_scores/
  scores/alphafold3.csv
  scores/protenix.csv
  scores/esmfold2.csv
  populations/alphafold3.json
  populations/protenix.json
  populations/esmfold2.json
  SHA256SUMS
```

Scores retain full floating-point precision and contain no labels. Missing
predictions are recorded in the population manifests, not represented by zero
scores. The dataset and score checksums detect accidental changes to a release.

RIC, gRIC, NAC, pair-iPTM, pDockQ, pDockQ2, and ipSAE are dimensionless; PAE and
gPDE are in Angstroms. Predictor-native ranking scores and AntiConf retain their
original scales. Whole-antibody pLDDT in the released score files follows the
original scoring outputs: 0-100 for AlphaFold 3 and Protenix-v2, and 0-1 for
ESMFold2. This positive scale change does not affect ranks or selected candidates.
NAC is on a 0-1 scale for all predictors. Normalized prediction files always use
0-100 pLDDT; the score exporter must preserve the predictor-specific score scale.

To create candidate metadata for structure scoring:

```bash
python scripts/prepare_candidates.py --dataset data/antibody_binding_datasets \
  --output data/candidates.jsonl
```

This checks sequence/metadata joins and constructs the CDR mask from the supplied
1-based sequence positions. It does not rerun antibody numbering or change labels.
