# Preparing the sequence dataset

The data package is distributed separately from the scoring code. It contains
three folders: `display_enrichment`, `igdesign_spr`, and `literature_affinity`.
Each has `sequences.csv`, `metadata.jsonl`, and a README. A shared field reference
explains every column, experimental unit, and missing-value convention.

The data files preserve the 4,184 antibody-antigen pairs and labels used in the
paper. They do not contain model predictions or metric scores.

## Checking a downloaded copy

The validator uses only the Python standard library:

```bash
python scripts/validate_dataset.py /path/to/dataset
```

It checks candidate counts, labels, sequence and metadata correspondence,
antigen chain boundaries, CDR sequence positions, count-derived frequencies
and enrichment, SPR repeat counts, and affinity bounds. File checksums can be
verified from inside the dataset directory with `shasum -a 256 -c SHA256SUMS`.

## Experimental measurements

The display files contain both the original selection-data enrichment signal
and a final-to-first log2 frequency ratio recalculated from the included round
counts. Those counts have undergone sequence quality control; they are not
raw sequencing reads.

The IgDesign files include repeat-level binding-count summaries and the source's
reported candidate-level K_D where available. A missing numerical affinity is
not replaced by zero. Individual SPR curves are not included.

The literature/patent files include individual affinity records and document
identifiers. Comparators such as `<` and `>` are retained. Only exact values
enter the minimum, median, and maximum affinity summaries.

## Before publication

The prepared dataset has not yet been uploaded to a public archive. The final
archive link, citation, and approved data license still need to be added.
The two-file dataset layout is distinct from the prediction loader's
`candidates.jsonl` format; do not pass `metadata.jsonl` directly to `score.py`.
