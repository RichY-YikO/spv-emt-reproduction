# Final figure data and verification records

The publication package keeps final drawing results separate from recomputation inputs.

| Location | Role |
|---|---|
| `code/fig1_fig2_results/S_final_pearson.npz` | Final unoriented S basis/scores, negated scores and standard Pearson inference for Fig.1 and Fig.S1 |
| `code/fig1_fig2_results/T_basis_and_correlations.npz` | Final T spatial/lag data for Fig.2 and Fig.S2; also required for Fig.4 inputs |
| `code/fig1_fig2_results/S_top100.csv` | Final S Pearson summary, replacing the older oriented-S summary |
| `code/fig1_fig2_results/T_top100.csv` | T correlation summary |
| `code/fig1_fig2_results/table_S1.csv` | Full-precision machine-readable Table S1 |
| `code/reference_inputs/` | Original oriented S cache and SPV daily CSV needed by recomputation/verification |
| `figures/TableS1.xlsx` | Original editable Excel copy; numerical agreement with the CSV is verified at tolerance 1e-10 |
| `validation/fig1_pearson/` | Fig.1 inference provenance and generated numerical-check report |
| `validation/layout/` | Fourteen figure layout records; these are not plotting data |

The two unused climatology caches and legacy oriented-S summary were removed from the release. The previous standalone Fig.1 cache directory no longer exists. All retained/moved NPZ and CSV bytes, and the 14 final PDF/PNG figures, are preserved. The Excel copy preserves the original workbook bytes. Dated older verification reports remain historical records rather than current path inventories.
