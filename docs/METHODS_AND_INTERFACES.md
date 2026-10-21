# Methods and interface descriptions

Descriptions extracted from the final source before removing code docstrings. Command-line usage is also available through --help. Numerical algorithms are unchanged.

## publication_revision.py / module

Revised Fig.1/2 display; scientific arrays and caches are never modified.

## publication_style.py / module

Publication drawing primitives. Scientific calculations stay in reproduce_*.py.

## replot_publication.py / module

Replot cached results without scientific recomputation.

python replot_publication.py --plot-only [--output-root DIRECTORY]
Default: replace figure PDF/PNG files in the existing result directories.


## reproduce_fig1_fig2.py / module

Reproduce SPV, Fig. 1 and Fig. 2 from the annual polar-cap NetCDF fields.

Run with the Python environment that has numpy, scipy, xarray, netCDF4,
matplotlib and cartopy.  The optional --surrogates 0 mode is a quick numeric
check; use the default 999 for the report's final null envelopes.


## reproduce_fig1_fig2.py / read_annual_fields

Return the continuous common-calendar axis; missing 10 hPa Feb 29 is NaN.

## reproduce_fig1_fig2.py / anomaly_and_climatology

All75 same-month/day grid-point climatology, finite observations only.

## reproduce_fig1_fig2.py / make_surrogates

One SPV surrogate ensemble shared by both layers and all 200 modes.

## reproduce_fig1_fig2.py / plot_figure

Draw using the shared publication style; inputs are unchanged.

## reproduce_fig1a_no_orientation.py / module

Standalone Fig1(a): eigensolver output without either sign-orientation rule.

Final Fig1/FigS1 use ordinary Pearson correlation of negated raw scores.
Each lag is centered using the mean of its own valid, NEGATED paired scores.
Use --from-cache to recompute inference without reading fields or decomposing.


## reproduce_fig1a_no_orientation.py / draw_full

Full final Fig1 from the approved no-orientation cache, no inference.

## reproduce_fig3.py / module

Reproduce report Fig. 3: 75-fold leave-one-winter-out direct SPV forecasts.

This script reads the original annual 10/500 hPa files. It never uses the
241-day merged files or full-sample U as prediction inputs. Run after
reproduce_fig1_fig2.py only for shared calendar constants and QA comparisons.


## reproduce_fig3.py / leaveout_delta

Fold climate minus All75 climate, per month-day and grid point.

## reproduce_fig3.py / plot_fig3

Draw using the shared publication style; inputs are unchanged.

## reproduce_fig4.py / module

Reproduce report Fig. 4 from the All75 source-field EMT/SPV outputs.

First run reproduce_fig1_fig2.py; this script reads its correct 241-day basis,
the signed 7-Apr projections, and standardized SPV. It fits All/SSW/non-SSW
full-sample descriptive regressions, not leave-one-winter predictions.


## reproduce_fig4.py / ols_cr1

Stable centered/scaled OLS and full CR1 winter-cluster covariance.

## reproduce_fig4.py / plot_fig4

Draw using the shared publication style; inputs are unchanged.

## reproduce_fig5.py / module

Reproduce report Fig. 5 from 100 hPa v/T and annual 10/500 hPa fields.

Requires the verified Fig. 3 archive beside this script. Usage:
    python reproduce_fig5.py
    python reproduce_fig5.py --heat-only


## reproduce_fig5.py / daily_heat_flux

Latitude-weighted zonal eddy covariance at each UTC, then daily mean.

## reproduce_fig5.py / fold_h40

All74 same-month/day anomaly, then 40 true consecutive calendar days.

## reproduce_fig5.py / fit_columns

Centered/scaled, unpenalized OLS for an arbitrary column subset.

## reproduce_fig5.py / plot_fig5

Draw using the shared publication style; inputs are unchanged.

## reproduce_figS1_figS2.py / module

Plot S1 from the final no-orientation Fig1 Pearson cache.

S1 uses negated raw scores centered by their own paired mean (Pearson).
S2 retains the original T Pearson correlations. Both gray bands are saved
fixed-lag phase-surrogate null envelopes, not lag-selection envelopes.
No preprocessing, decomposition or surrogate inference is recomputed here.


## reproduce_figS1_figS2.py / draw

Draw using the shared publication style; inputs are unchanged.

## reproduce_figS3.py / module

Reproduce report Fig. S3: B versus separate S and T extensions.

The 75-fold 10/500 hPa EMT, SPV target and pair mask follow reproduce_fig3.py.
Verified Fig. 3 B predictions are reused; BS/BT are independently refit for
each held-out winter, group branch and lead. No 100 hPa heat flux is used.


## reproduce_figS3.py / plot

Draw using the shared publication style; inputs are unchanged.

## reproduce_figS4_S5_S6.py / module

Reproduce report Fig. S4/S5/S6 by target-month scoring of frozen Fig. 3 predictions.

No EMT, SPV preprocessing or regression is refit by month. The 2000 paired
winter bootstrap uses the report's same PCG64 winter weights as Fig. 3.


## reproduce_figS4_S5_S6.py / monthly_winter_sums

Score each held-out winter and target month on the same B/M pairs.

## reproduce_figS4_S5_S6.py / score

Return [month, metric(r/Q), model(B/M), lead] plus percentile bands.

## reproduce_figS4_S5_S6.py / plot

Draw using the shared publication style; inputs are unchanged.

## reproduce_figS7_S8.py / module

Reproduce report Fig. S7 event-aligned and Fig. S8 rolling-winter scores.

Both figures score frozen leave-one-winter-out Fig. 3 predictions. Fig. S7
reads the original SSWevent.xlsx and verifies its hash and Appendix A dates.


## reproduce_figS7_S8.py / plot_s7

Draw using the shared publication style; inputs are unchanged.

## reproduce_figS7_S8.py / plot_s8

Draw using the shared publication style; inputs are unchanged.

## reproduce_figS9.py / module

Reproduce report Fig. S9: default Sep-Apr versus DJF EMT basis window.

The only experimental change is the season used to estimate spatial U. Every
held-out winter has All74 climate/SPV preprocessing and the same dates in both
schemes; all July-April fields, including the S1 history, are reprojected to
their corresponding basis. B/M are independently refit for each scheme.


## reproduce_figS9.py / plot

Draw using the shared publication style; inputs are unchanged.

## supplement_style.py / module

SI-only drawing styles; primary-figure styles and scientific arrays are untouched.

## verify_final_pearson.py / module

Focused numerical audit of final Fig1/FigS1; no plotting or decomposition.

Optional --previous CACHE checks unchanged raw arrays against the backed-up
custom-statistic cache. Uses independently expressed Pearson formulas.


## verify_supplement_figures.py / module

Validate SI drawn data against caches; no science or image files regenerated.

