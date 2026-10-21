# Figure and computation map

All final PDFs and 600 dpi PNGs are under `figures/`. Frozen cached redrawing goes through `reproduce.py`; computational scripts below are in `code/`.

| Figure | Original-data computation | Current cached inputs |
|---|---|---|
| Fig.1 | `reproduce_fig1a_no_orientation.py`; shared preprocessing in `reproduce_fig1_fig2.py` | `fig1_fig2_results/S_final_pearson.npz` |
| Fig.2 | `reproduce_fig1_fig2.py` | `fig1_fig2_results/T_basis_and_correlations.npz` |
| Fig.3 | `reproduce_fig3.py` | `fig3_fig4_results/fig3_frozen_predictions.npz` |
| Fig.4 | `reproduce_fig4.py` | `fig3_fig4_results/fig4_coefficients_CR1.npz` |
| Fig.5 | `reproduce_fig5.py` | `fig5_results/fig5_scores.npz` |
| Fig.S1 | `reproduce_figS1_figS2.py` reads final Fig.1 | Final unoriented Pearson cache, fixed-lag bounds |
| Fig.S2 | `reproduce_figS1_figS2.py` reads Fig.2 | Original T Pearson cache, fixed-lag bounds |
| Fig.S3 | `reproduce_figS3.py` | `figS3_S9_results/FigS3_scores.npz` |
| Fig.S4–S6 | `reproduce_figS4_S5_S6.py` | Corresponding monthly caches in `figS3_S9_results/` |
| Fig.S7–S8 | `reproduce_figS7_S8.py` | `figS3_S9_results/FigS7_scores.npz` and `FigS8_scores.npz` |
| Fig.S9 | `reproduce_figS9.py` | `figS3_S9_results/FigS9_scores.npz` |
| Table S1 | `reproduce_fig1_fig2.py` | Editable Excel `figures/TableS1.xlsx`; source CSV `code/fig1_fig2_results/table_S1.csv` |

`publication_style.py` and `publication_revision.py` define primary drawing styles; `supplement_style.py` defines SI drawing styles. Existing computational function signatures remain intact. `redraw_final.py` provides the English final-redraw interface.

The oriented S basis remains in `code/reference_inputs/` for recomputation and independent checks. Its legacy correlation CSV has been removed; `code/fig1_fig2_results/S_top100.csv` describes final Fig.1 only. The obsolete sign-flipped cache and custom-statistic implementation are excluded.
