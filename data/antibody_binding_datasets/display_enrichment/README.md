# Display-based enrichment

This dataset contains 1,400 antibody-antigen pairs: 100 positive and 100
negative candidates for each of seven targets. TfR, spike S2, and GPC3 were
screened by yeast display. CD98hc, EGFR, PLXNA1 ECD, and TNFRSF23 were screened
by phage display.

`sequences.csv` contains sequences and identifiers. `metadata.jsonl` contains
the paper labels, CDR positions, display method, and sequence-level read counts
for each available selection round. Counts are joined by exact heavy-chain
amino-acid sequence within the original selection experiment. Other
experiments for the same target are not pooled into these counts.

Positive labels represent selected enriched candidates. Negative labels
represent baseline or low-enrichment candidates; they are not direct
measurements of absent binding in a soluble assay. Expression, clone growth,
and sampling can also influence sequence frequencies.

For the phage-display targets, successive panning rounds used 5, 4, 2, and
1 microgram of antigen. For the yeast-display targets, the final sorting round
used antigen at 1 nM. Round names follow the source records; not every dataset
starts at `R0`.

The included counts have undergone sequence quality control and are not raw
FASTQ reads. The derived final-to-first enrichment ratio is provided for
convenience; it does not regenerate the paper labels. See [the field
reference](../FIELDS.md) for its definition and the treatment of zero counts.
