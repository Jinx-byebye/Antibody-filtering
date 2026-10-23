# Validation status

## Checked

- Fresh CPU environment: Python 3.13, NumPy 2.2.6, editable package installation.
- Dataset validation and reproduction directly from the two release archives.
- All 936 original displayed values in Tables 2, 3 and 4 agree with independently extracted
  manuscript references, including average ranks.
- Thirty-one unit tests pass, covering score joins, duplicates, non-finite values,
  ties, ranking, sequences, CDR masks, conversion, and complete-ensemble checks.
- Synthetic single-heavy-chain, heavy/light-chain, and multichain-antigen ensembles
  agree with the original RIC, pDockQ, pDockQ2, ipSAE and AntiConf calculation kernels.
- Nine real ensembles (three configurations for each of the three predictors),
  each containing ten seeds and five samples, have been converted and rescored.
  All 108 compared candidate-metric values match released full-precision scores
  within relative tolerance `1e-8` and absolute tolerance `1e-10`. These cover
  single-heavy-chain antibodies, paired heavy/light chains with one antigen chain,
  and paired heavy/light chains with a multichain antigen. The audit also verified
  that ESMFold2 NPZ pLDDT uses a 0-1 scale while its structure files use 0-100.
- Current model asset sizes and SHA-256 hashes have been collected for all three
  installations, including the six ESMC-6B weight shards. The public asset-check
  command tests missing, changed, and valid files without loading a predictor.

## Still required before claiming sequence-to-table reproduction

The GPU inference smoke test passed for AlphaFold 3, Protenix-v2, and ESMFold2
on 2026-10-16. Each model used one single-heavy-chain complex, two seeds, and
one sample per seed, followed by output conversion and scoring through the
public scripts. The test reused existing MSAs after validating their sequences
and checksums; it did not test a new MSA search or reproduce a full paper ensemble.

- Verify historical software revisions, checkpoint identities, and environment
  dependency locks for all three predictors.
- Validate fresh local MMseqs2-GPU searches and their integration with prediction.

The CPU release uses the original full-precision scores. It does not depend on
the unfinished GPU validation, and synthetic tests do not establish that validation.

## OpenDDE extension

The four-predictor score release reproduces every displayed value in Tables 2
and 3, including average ranks. Adding OpenDDE leaves all original predictor
rows unchanged. OpenDDE contributes 4,184 complete ensembles and 2,400 eligible
within-antigen candidates across 13 antigens.

Public lossless conversion and scoring were checked on three complete ensembles:
a single antibody chain, paired heavy/light chains, and a multichain antigen.
All 13 scores matched the saved evaluation scores exactly for each 50-sample
ensemble. CPU tests cover sample identities, unrounded pLDDT, chain-pair gPDE,
and inference-budget arguments. The public inference wrapper uses the output
adapter and inference flags used for the completed benchmark; it has not been
rerun for a second full benchmark.

Figure 2 uses a separate fixed set of 2,464 candidates. Seven regional scores
are compared by per-antigen AUC, fixed-budget recall, and maximum utility, with
10,000 paired antigen bootstrap resamples. Original three-predictor numerical
results are retained unchanged when OpenDDE is added.
