# Paper illustrations

Result illustrations are rendered from the numerical reference files.
The original binding-information and ensemble-size panels are retained from
the accompanying manuscript.

| File | Source |
| --- | --- |
| `across-antigen-results.png` | Table 2, four predictors |
| `within-antigen-results.png` | Table 3, four predictors |
| `binding-information.png` | Figure 1 |
| `regional-confidence.png` | Figure 2, seven regional scores and four predictors |
| `ensemble-size.png` | Ensemble-size sensitivity analysis |

To regenerate the two tables and regional-confidence figure, install Matplotlib
and run `python scripts/render_paper_assets.py`. Matplotlib is not needed for
CPU table evaluation. The reference CSV files, rather than PDF crop coordinates,
determine every plotted value.
