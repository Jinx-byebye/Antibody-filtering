# Literature and patent affinity

This dataset contains 1,613 antibody-antigen pairs spanning 165 antigen
constructs, with 1,510 positive and 103 negative labels from the paper.

`sequences.csv` contains the sequences used for prediction. `metadata.jsonl`
contains the individual affinity records, source document identifiers,
measurement methods, original affinity text, and values converted to nM.

Always interpret `kd_value_nm` together with `kd_comparator`. For example,
`1000` with `>` means K_D > 1000 nM, not K_D = 1000 nM. The minimum, median,
and maximum summaries use exact measurements only and are `null` when no
exact measurement is available. The paper labels are preserved independently
of these descriptive summaries.

Multiple records for one pair can originate from one or several documents.
They need not be independent experimental repeats or share identical assay
conditions. Both source record identifiers and document identifiers are kept
to support inspection of the original evidence.

See [the field reference](../FIELDS.md) for every column and metadata field.
