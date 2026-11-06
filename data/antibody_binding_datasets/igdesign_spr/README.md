# IgDesign SPR

This dataset contains 1,171 designed antibody-antigen pairs across seven
targets, with 249 positive and 922 negative labels. The source is the
[IgDesign dataset](https://github.com/AbSciBio/igdesign).

`sequences.csv` contains the heavy-chain, light-chain, and antigen sequences
used for prediction. `metadata.jsonl` contains the design method, CDR
annotations, and reported numbers of SPR repeats and binding-positive repeats.
It also includes the source's reported K_D values for all 249 positive
candidates. The remaining 922 candidates have no reported numerical K_D;
their `reported_kd_nm` is `null`, not zero.

Positive candidates bind in all reported repeats and negative candidates bind
in none, with at least two repeats per candidate. These are replicated binding
calls. Numerical affinities, where reported, are provided separately and do
not replace the repeat-based labels. The current export does not
include individual sensorgrams, per-repeat response curves, or per-repeat K_D fits;
`individual_measurements` is therefore `null`.

Antigen sequences describe the prediction inputs. Where available, their
construct descriptions and sequence references are retained. Do not infer
unreported assay-construct details from the prediction input alone.

See [the field reference](../FIELDS.md) for every column and metadata field.

The source CSV files were retrieved at revision
`70431eef0afaf0496d7d84e22dfdc1980ec9e70e`. Candidate identities, CDR sequences,
design methods, and binding repeat counts were checked against these files.
The source file name, revision, and checksum are included in each record.
