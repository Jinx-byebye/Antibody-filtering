# Evaluation

`scripts/reproduce_tables.py` reads labels from the dataset and scores from CSV
files joined by `pair_id`. It rejects duplicates, non-finite values, missing
metrics, and discrepancies with the supplied population manifests. Labels are
never inferred from scores or identifier text.

## Tables

- Table 2 pools all eligible candidates for each predictor and uses a shared
  threshold across antigens.
- Table 3 evaluates each eligible antigen independently. A target must have at
  least ten positives and forty negatives among that predictor's available
  candidates. Display enrichment and IgDesign each receive half the total
  weight. Within each dataset, target weights are proportional to
  `sqrt(2 * positives * negatives / (positives + negatives))`.
- Table 4 pools precisely the Table 3 candidates and selects new shared
  thresholds on this subset. It does not reuse Table 2 thresholds.

Within-antigen endpoints are weighted averages of target-level endpoints, not
ratios calculated from pooled selected counts. Each metric is evaluated on the
same population for a given predictor and table.

## Thresholds, budgets, and ties

Higher scores indicate stronger candidates except for PAE and gPDE, whose signs
are reversed for ranking. The threshold sweep includes every distinct score and
the empty selection. A threshold selects all candidates whose score is at least
that value. Equal utilities are resolved by lower empirical FPR, then higher
recall, then the higher threshold. Empty selections have zero recall and precision.

Recall at an FPR limit uses the best feasible complete score group. Fixed-budget
recall uses `floor(Q * N)` candidates, with `Q = 0.05` or `0.20`. A score tie at
the boundary contributes its expected positive count under uniform random
selection within the tied group. Labels do not resolve score ties.

Utility is `R_lower * W(F_upper)`, using 95% Wilson bounds (z = 1.959963984540054).
With `F = F_upper`, the weight is
`exp(-(F/f_severe)^2 - (max(F-f_safe, 0)/(f_severe-f_safe))^2)`.
The paper uses `f_safe = 0.10` and `f_severe = 0.30`. These are dimensionless FPR
fractions, not percentages. The maximizing threshold is selected using known
labels; it is a retrospective evaluation, not a calibrated cutoff for new targets.

Average rank uses eight endpoints for Tables 2 and 3 and seven for Table 4, which
does not include recall at maximum utility. Ranks are computed from unrounded
values, with average ranks for ties (relative and absolute tolerance `1e-12`).
Percentages are displayed to two decimals, utility to four, and rank to two.

## Configuration

`configs/evaluation.json` defines the following fields. All fields are required;
the script does not supply alternative settings for missing values.

| Field | Meaning and paper setting |
| --- | --- |
| `f_safe`, `f_severe` | Utility FPR scales, 0.10 and 0.30 |
| `minimum_positive`, `minimum_negative` | Per-antigen eligibility counts, 10 and 40 |
| `within_antigen_datasets` | Dataset folder names included in Tables 3 and 4 |
| `dataset_weights` | Normalized dataset weights, 0.5 each |
| `score_directions` | +1 for higher-is-better; -1 for lower-is-better |
| `predictor_metrics` | Ordered metric columns required for each predictor |
| `endpoints` | Ordered reported measures used in each table's average rank |

The reference CSVs were extracted independently from the manuscript. Passing
`--reference reference` checks every displayed value and fails on any mismatch.
It does not replace calculated values with reference values.
## OpenDDE and regional confidence

OpenDDE uses the same frozen evaluation rules and 13 metrics as Protenix-v2.
Its populations contain 4,184 candidates for across-antigen evaluation and
2,400 candidates across 13 antigens for within-antigen evaluation. Average
ranks are calculated independently within each predictor; adding OpenDDE
does not change the original predictors' ranks. The reproduction command
also provides the matched-population global evaluation for OpenDDE in its
Table 4 output; the manuscript's Table 4 retains the original three predictors.

Figure 2 uses a different, fixed matched population: 2,464 candidates, including
1,293 display-based and 1,171 IgDesign pairs, with seven antigens in each benchmark.
The seven scores are NAC, N3AC (excluding only CDR3), CDR3+NAC (excluding CDR1/2),
whole-antibody, CDR, CDR3-only, and interface pLDDT. CDR definitions apply to both
heavy and light chains. Interface pLDDT uses unique contact residues in antibody
variable domains and the antigen, with an 8-Angstrom C-beta cutoff and C-alpha
fallback. Samples without contacts receive zero.

Every score uses a per-sample residue mean, then the median of five samples and
the median of ten seeds. AUC, recall at a 20% budget, and maximum utility are
computed separately for each antigen, then averaged with equal antigen weights.
The 95% intervals use 10,000 paired antigen bootstrap resamples, shared across
scores and predictors. They describe antigen variation, not independently
validated threshold performance.

The figure's numerical summaries are in `reference/regional_confidence/`.
With Matplotlib installed, regenerate its PDF, SVG, and PNG using:

```bash
python scripts/plot_regions.py reference/regional_confidence \
  --output results/regional_confidence
python scripts/render_tables.py --reference reference --output results/tables
```
