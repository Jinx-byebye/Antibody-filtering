# Evaluation populations

The score release contains `populations/alphafold3.json`, `protenix.json`, and
`esmfold2.json`, and `opendde.json`. These manifests accompany the data rather than the source code.
Each lists included candidate identifiers, the within-antigen subset, and excluded
identifiers with reasons. The reproduction script checks the complete partition
against the dataset and requires exactly one score row per included candidate.

AlphaFold 3, Protenix-v2, and OpenDDE use 4,184 candidates for Table 2 and 2,400 for Tables 3
and 4. ESMFold2 uses 3,909 and 2,293, respectively. Successful predictions on new
hardware do not change these paper populations. To evaluate an expanded cohort,
provide a separately named population manifest and do not describe that result as
an exact reproduction of the paper.
