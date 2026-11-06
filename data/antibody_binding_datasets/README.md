# Antibody binding and enrichment datasets

This collection contains the 4,184 antibody-antigen pairs used in the paper.
It combines display-based selection, replicated surface plasmon resonance
(SPR) measurements, and affinity measurements reported in the literature
and patents. The three datasets retain their original experimental meaning:
selection enrichment is not a direct measurement of binding affinity.

| Folder | Antigen constructs | Positive | Negative | Total |
| --- | ---: | ---: | ---: | ---: |
| `display_enrichment` | 7 | 700 | 700 | 1,400 |
| `igdesign_spr` | 7 | 249 | 922 | 1,171 |
| `literature_affinity` | 165 | 1,510 | 103 | 1,613 |
| Total | 179 | 2,459 | 1,725 | 4,184 |

The antigen total is the sum across datasets, not a count of distinct proteins
across the whole collection. Different constructs of one protein can have
different sequences.

## Files

Each dataset contains two data files and a short README:

- `sequences.csv` contains identifiers and amino-acid sequences only.
- `metadata.jsonl` contains one JSON object per antibody-antigen pair, including
  its label, CDR annotation, and experimental evidence. Multiple measurements
  are stored together in the same object.

Join these files using `pair_id`, not the row number. An antibody tested against
two antigens represents two pairs. Identical antibody sequences share an
`antibody_id`; identical ordered antigen chain sequences share an `antigen_id`.
The pair identifiers are retained so that experimental records can be matched
to the paper's prediction and score files.

All fields are defined in [FIELDS.md](FIELDS.md). `summary.json` records dataset
sizes. `SHA256SUMS` provides file checksums when the package is assembled.

## Reading the data

The files use UTF-8. CSV files have a header. In CSV, an empty light-chain field
means that no light chain is included. In JSON, `null` means that a value was
not available in the supplied source, not that the measured value was zero.

```python
import csv
import json
from pathlib import Path

folder = Path("display_enrichment")
with (folder / "sequences.csv").open(newline="", encoding="utf-8") as f:
    sequences = {row["pair_id"]: row for row in csv.DictReader(f)}
with (folder / "metadata.jsonl").open(encoding="utf-8") as f:
    metadata = {row["pair_id"]: row for row in map(json.loads, f)}

assert sequences.keys() == metadata.keys()
pair_id = next(iter(sequences))
print(sequences[pair_id]["heavy_chain_sequence"])
print(metadata[pair_id]["label_description"])
```

## Scope and sources

Sequences are the chains used for structure prediction in the paper. They
should not automatically be interpreted as complete expressed constructs:
signal peptides, purification tags, constant domains, or other experimental
components may not be represented. Chain boundaries and repeated antigen
chains are preserved.

Labels are the labels used in the paper; this export does not redefine them.
The collection is selected for evaluation, so its positive fraction should
not be interpreted as a prospective experimental hit rate.

The display dataset includes sequence-level counts after quality control,
not raw sequencing reads. The IgDesign files currently provide aggregate
SPR binding calls across repeats and the source's reported K_D values where
available, but not sensorgrams or per-repeat affinity fits.
The literature dataset retains reported values and bounds with source
document identifiers; multiple records are not necessarily independent
experimental repeats.

## Source attribution

The IgDesign and
literature/patent portions retain their source attribution; no new blanket
redistribution permission is asserted here.
