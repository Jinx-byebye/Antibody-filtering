# Field reference

All three datasets use the same sequence-file columns. Experimental metadata
differs by dataset because the assays measure different outcomes.

## sequences.csv

| Field | Meaning |
| --- | --- |
| `pair_id` | Identifier of the antibody-antigen pair; joins the two data files. |
| `antibody_id` | Sequence-derived identifier for the heavy chain and optional light chain. |
| `antigen_id` | Sequence-derived identifier for the ordered list of antigen chains. |
| `heavy_chain_sequence` | Heavy-chain amino-acid sequence used in prediction, in one-letter notation. |
| `light_chain_sequence` | Light-chain amino-acid sequence; empty when no light chain was included. |
| `antigen_chain_ids` | Antigen chain identifiers separated by semicolons. |
| `antigen_chain_sequences` | Corresponding amino-acid sequences separated by semicolons, in the same order. |

A homodimer has two chain identifiers and two sequences, even when the
sequences are identical. Do not remove the repeated chain or concatenate
chains without retaining their boundaries. Antibody and antigen identifiers
use the first 20 hexadecimal characters of a SHA-256 digest of their ordered
sequence lists, prefixed by `ab_` and `ag_`, respectively. Antibody lists include
an empty second sequence when no light chain is present.

## Shared metadata

| Field | Meaning |
| --- | --- |
| `pair_id` | Identifier matching a row in `sequences.csv`. |
| `target_name` | Target name or short target identifier used in the benchmark. |
| `label` | Paper label: `1` for positive and `0` for negative. Its experimental meaning is given below. |
| `label_description` | Plain-language description of the assigned label. |
| `cdr_numbering_scheme` | `IMGT`, the numbering scheme used to identify CDRs. |
| `cdr_sequence_positions` | Object mapping antibody chain identifiers (`H`, optionally `L`) to CDR residue positions. Positions are one-based indices in the supplied chain sequence, not IMGT residue numbers. |

## Display-based enrichment

| Field | Meaning |
| --- | --- |
| `display_method` | `yeast display` or `phage display`. |
| `selection_enrichment_log2_ratio` | Original final-stage versus first-stage log2 frequency contrast retained from the candidate-selection data. Frequencies were smoothed and processed by the original multi-round sequence-quality procedure. |
| `enrichment_log2_ratio` | Log2 ratio of the final to the first available round's smoothed sequence frequency, calculated from the included count tables. See the formula below. |
| `measurements` | List of round-specific sequence counts, in selection order. |
| `measurements[].round` | Source round name, such as `R0`, `R1`, or `R4`. Missing round numbers are not filled in. |
| `measurements[].read_count` | Number of reads assigned to the exact heavy-chain amino-acid sequence after sequence quality control. |
| `measurements[].total_read_count` | Total accepted read count across all sequences in that round, not just the 200 benchmark candidates. |
| `measurements[].frequency` | `read_count / total_read_count`, without a pseudocount. |

For count `n` and round total `N`, the smoothed frequency is `(n + 0.5)/(N + 1)`.
`enrichment_log2_ratio` is `log2(f_final/f_first)` using these smoothed
frequencies. Positive values indicate increased frequency; zero indicates
equal smoothed frequencies. The first available round need not be an
unselected input library. This descriptive quantity is calculated for the
public export and is not a replacement for the original multi-round selection
and label-assignment procedure. It is not a K_D or a binding probability.

`selection_enrichment_log2_ratio` preserves the existing selection-data value;
`enrichment_log2_ratio` is recalculated from the included round count tables.
They need not be identical because the original selection procedure applies
additional sequence-quality filtering. Neither quantity alone determines the
benchmark label.

A zero count means that the exact sequence was absent from an available,
fully read round table. It does not mean an unavailable file or an unperformed
experiment. The count tables sum duplicate read counts by exact reconstructed
heavy-chain sequence after sequence quality control; they do not report
unique molecules or raw FASTQ read totals.

## IgDesign SPR

| Field | Meaning |
| --- | --- |
| `source_record_id` | Identifier linking the candidate to the reconstructed IgDesign sequence record. |
| `design_method` | Sequence-design method recorded by the source; not the assay method. |
| `designed_cdrs` | Heavy-chain CDRs designated for design in the source, such as `HCDR123`. |
| `assay_method` | `SPR` (surface plasmon resonance). |
| `replicate_count` | Number of reported assay repeats. |
| `binding_replicate_count` | Number of repeats reported to show binding. |
| `reported_kd_nm` | Candidate-level K_D in nM from the source's `KD (nM)` column; `null` when not reported. This is not a list of per-repeat fits. |
| `reference_antibody` | Reference antibody named in the IgDesign source data. |
| `individual_measurements` | Currently `null`: individual repeat measurements and sensorgrams were not supplied by the source used for this export. |
| `antigen_construct_description` | Description of the antigen construct used for prediction, when recorded. |
| `antigen_sequence_reference` | Sequence source or residue range used for that prediction construct, when recorded. |
| `source_url` | Original IgDesign repository. |
| `source_file` | Name of the source CSV file. |
| `source_revision` | Git revision identifying the source-file version. |
| `source_file_sha256` | SHA-256 checksum of that source CSV file. |

Positive candidates show binding in every reported repeat; negative
candidates show binding in none. Each has at least two repeats. A binding
call is not a numerical K_D measurement. Prediction construct descriptions
do not by themselves establish the exact soluble construct used in SPR.

## Literature and patent affinity

| Field | Meaning |
| --- | --- |
| `measurements` | List of affinity records associated with this pair. |
| `measurements[].source_record_id` | Identifier of the extracted source record. |
| `measurements[].document_id` | Source document identifier, such as a patent publication number. |
| `measurements[].antibody_name` | Antibody name as recorded in the source. |
| `measurements[].antigen_name` | Antigen name as recorded in the source. |
| `measurements[].antigen_uniprot_id` | UniProt accession, when available. |
| `measurements[].reported_kd` | Original reported affinity text, retaining its reported unit. |
| `measurements[].kd_value_nm` | Numerical value converted to nanomolar; interpret together with `kd_comparator`. |
| `measurements[].kd_comparator` | `=`, `<`, `<=`, `>`, or `>=`. Bounds are not exact measurements. |
| `measurements[].assay_method` | Reported assay method, such as SPR, BLI, or ITC. |
| `measurements[].assay_context` | Short experimental-context description from the source extraction, when available. |
| `affinity_summary.measurement_count` | Number of included affinity records, including bounds. |
| `affinity_summary.exact_measurement_count` | Number of records with comparator `=`. |
| `affinity_summary.document_count` | Number of distinct source document identifiers. |
| `affinity_summary.exact_kd_min_nm` | Minimum among exact measurements only. |
| `affinity_summary.exact_kd_median_nm` | Arithmetic median among exact measurements only. |
| `affinity_summary.exact_kd_max_nm` | Maximum among exact measurements only. |

The three exact-value summaries are `null` if all records are bounds. Records
may differ in assay conditions and may repeat information across documents;
their count is not an independent replicate count. Inspect individual records
before treating a summary as a single experimental affinity.

## summary.json

Each key is a dataset folder name. `candidates` is the number of pairs,
`positive` and `negative` count the paper labels, and `antigen_constructs`
counts distinct `antigen_id` values within that dataset.
