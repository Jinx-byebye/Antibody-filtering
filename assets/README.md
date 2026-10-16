# Paper illustrations

These images are rendered directly from the manuscript in
`paper/antibody-filtering.pdf`. Values and artwork are unchanged.

| File | Source |
| --- | --- |
| `across-antigen-results.png` | Table 2, pages 6-7, joined into one table with a single column header |
| `within-antigen-results.png` | Table 3, page 8 |
| `binding-information.png` | Figure 1, page 9 |
| `regional-confidence.png` | Table 7, page 16 |
| `ensemble-size.png` | Figure 2, page 17 |

To regenerate them, install PyMuPDF in a documentation environment and run
`python scripts/render_paper_assets.py`. This dependency is not needed for evaluation.
If the PDF layout changes, review the crop coordinates and rendered images.
