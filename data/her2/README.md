# AI-designed antibodies targeting HER2

This independent evaluation contains 1,854 AI-designed antibodies targeting one
HER2 construct: 757 author-reported SPR binders and 1,097 non-binders. The parent
antibody control is excluded. The sequences and experimental measurements are
from **Absci Corporation (2023)**,
[Unlocking de novo antibody design with generative artificial intelligence](https://github.com/AbSciBio/unlocking-de-novo-antibody-design).
The original data license is retained in [SOURCE_LICENSE](SOURCE_LICENSE).

All candidates share the same concatenated non-CDR sequence under the supplied
IMGT annotation, while all 1,854 concatenated CDR sequences are distinct.
There are 1,607 candidates in the HCDR3-only design group and 247 with additional
HCDR1/2 changes. NAC excludes CDR residues only from its confidence average:
the prediction itself includes the complete antibody and antigen.

## Contents

| File | Contents |
| --- | --- |
| `candidates.jsonl` | Sequences, chain roles, experimental labels, CDR annotation and source metadata |
| `scores/<predictor>/<predictor>.csv` | Full-precision scores linked by `pair_id` |
| `evaluation.json` | Score directions, utility parameters, eligibility and ranking endpoints |
| `reference/all.csv`, `reference/all.md` | Complete-cohort results for three predictors |
| `SHA256SUMS` | SHA-256 checksums of the dataset, scores, configuration and reference tables |

Each candidate has ten prediction seeds (101-110), with five samples per seed,
for each included predictor. Large structures and residue-level confidence
arrays are not included. This cohort is evaluated separately from the original
4,184-candidate benchmark. ESMFold2 is not included in this complete-cohort release.

### Candidate fields

- `pair_id`: unique candidate identifier shared with the score files.
- `benchmark`: `absci_her2`.
- `target_id`: HER2 construct identifier; `HER2_P04626_23_652` denotes residues 23-652 of UniProt P04626.
- `antibody_sequences`: amino-acid sequences keyed by heavy (`H`) and light (`L`) chain.
- `antigen_sequences`: amino-acid sequence keyed by antigen chain (`A`).
- `label`: original SPR binding call, 1 for binder and 0 for non-binder.
- `cdr_mask`: Boolean mask in antibody-chain order (H then L); true marks a CDR residue.
- `cdr_annotation`: IMGT scheme and one-based CDR positions in each chain sequence.
- `evidence.type`: experimental evidence type (`author_reported_SPR`).
- `evidence.group`: design group, `HCDR3_only` or `HCDR123_changes`.
- `evidence.kd_nm`: reported dissociation constant in nanomolar, or null when unavailable.
- `evidence.source_file`, `evidence.source_csv_line`: original source CSV and one-based line number, including the header.
- `evidence.replicate_consistency_verified`: whether replicate consistency was independently verified; false does not imply inconsistent measurements.

Score definitions and units follow the [data documentation](../../docs/data.md).
All score files have exactly the same candidate set; labels are read from the
candidate file, not inferred from confidence scores.

## Reproduce the results

From the repository root, after installing the CPU evaluation dependencies:

```bash
python scripts/evaluate_external.py \
  --candidates data/her2/candidates.jsonl \
  --scores data/her2/scores \
  --predictors alphafold3 protenix opendde \
  --config data/her2/evaluation.json \
  --output results/her2
```

The command produces CSV, Markdown and LaTeX tables, unrounded results, class
counts and separate design-group analyses. No GPU or prediction files are needed.
Average ranks use the eight endpoints listed in `evaluation.json`, independently
within each predictor. Maximum utility selects a threshold using the labels;
fixed-budget recall uses the top floor(Q*N) candidates with label-independent
boundary-tie handling. At Q=20%, 370 candidates are selected, giving a random
selection expectation of 19.96% recall.

For Protenix-v2, NAC has the best average rank (1.62), compared with 5.75 for
pair-ipTM. NAC recall at Q=20% is 23.91%, versus 20.48% for pair-ipTM and 24.57%
for whole-antibody pLDDT. OpenDDE NAC achieves 32.63% recall, compared with 32.89%
for pair-ipTM. See the [full comparison](reference/all.md) for all endpoints.
