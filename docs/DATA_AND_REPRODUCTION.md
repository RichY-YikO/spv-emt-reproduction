# Data contract and original-data reproduction

## Included derived data

NPZ files contain frozen basis vectors, scores, lag statistics, prediction results or figure-level summaries. CSV files expose SPV values, lag correlations, model scores, coefficient results, events and Table S1. All original NPZ and CSV bytes selected for release are verified against `validation/source_data_sha256.json`. Raw NetCDF and the SSW event workbook are not included. The verified Table S1 workbook is included as `figures/TableS1.xlsx`.

The study author confirms that the source data are ERA5 reanalysis from the European Centre for Medium-Range Weather Forecasts (ECMWF) Climate Data Store. Source citation: Copernicus Climate Change Service (2018), ERA5 hourly data on pressure levels from 1940 to present, https://doi.org/10.24381/cds.bd0915c6. The citation applies to the ERA5 pressure-level data used for the 10, 500 and 100 hPa analyses. This attribution is author-supplied; the original 10 hPa files contain no provenance attributes. Exact original download requests and preprocessing details are not recreated by the citation alone. The required local processed-file contracts remain documented below.

**SSW event provenance confirmed by the author.** The SSW event dates were calculated by Shuo Yang using the identification method of Charlton and Polvani (2007; https://doi.org/10.1175/JCLI3996.1). This provenance statement is confirmed by the author; this documentation update does not rerun event detection. The existing workbook hash, date checks and derived event tables are retained unchanged. See SSW_EVENT_PROVENANCE.md for the method citation and author-confirmed event provenance. Data and third-party materials retain their respective terms, separately from the original software MIT license.

## Expected external files

Defaults are relative to this repository; each script supports explicit path arguments.

| Location | Input contract |
|---|---|
| `data/raw/10hpa/1948.nc` … `2023.nc` | Variable `z`, shape `(365,11,144)`, coordinates `latitude` and `longitude`. Time is integer 1…365 and omits Feb.29; the loader assigns calendar dates explicitly. |
| `data/raw/500hpa/1948.nc` … `2023.nc` | Variable `z`, coordinates `lat`, `lon`; decoded calendar dates must be unique. Daily fields on the same 11×144 grid. |
| `data/raw/v100/era5_v_100hPa_YEAR.nc` | Variables `v`, `valid_time`, `pressure_level`, `latitude`, `longitude`; pressure level 100 hPa, four UTC observations per day. |
| `data/raw/t100/era5_temperature_100hPa_YEAR.nc` | Variable `t` and matching coordinates/time; v/T grids must agree. |
| `data/raw/SSWevent.xlsx` | Original workbook, `Sheet1!B2:B52`, 51 date rows; exact SHA256 is checked by the event calculation. |

The 10/500 hPa polar grid is 90°N to 65°N in descending 2.5° steps, longitude −180° to 177.5° in 2.5° steps. The 100 hPa fields cover 45–75°N with a 31×360 grid at 1° resolution, latitude order expected by the implementation is 75°N to 45°N; UTC times are 00,06,12,18. Their loader requires existing study-period subsets: dates outside 1948-07-01 through 2023-06-30 are rejected. Merely downloading unrestricted full calendar years may therefore require prior subsetting.

The common study axis is 1948-07-01 through 2023-06-30. The basis mask uses September–April, excludes September 1 and February 29, and has 18,075 days. SPV targets use November–April excluding February 29, giving 13,575 days. Longitude and latitude must match the contracts exactly; the scripts validate them. Units must match the original study inputs; do not silently switch geopotential to geopotential height.

The event workbook must match SHA256 `f7ed7bea677d7d13b63b5a7a8501e5604f9a3bef36cf7732a46a3d6236c0d238`. The study subset has 47 events in 39 SSW winters. Derived event CSVs are included for inspection; they are not silently substituted for the authenticated workbook in a full recomputation.

## Reproduction levels

1. **Redraw only:** `python reproduce.py`; uses included results, no raw inputs, no statistical recomputation.
2. **Verify Pearson:** `python verify_release.py --mode numerical`; directly verifies cached observations and surrogate bounds.
3. **Recompute Fig.1/S1 inference only:**

```bash
python code/reproduce_fig1a_no_orientation.py --from-cache code/fig1_fig2_results/S_final_pearson.npz
python reproduce.py --figures Fig1 FigS1
```

The first command overwrites the final shared inference cache using the stored raw scores. Make a separate working copy before an inference or full raw-data run; the package's checksums intentionally identify changes. It always starts from `scores_raw`, so it cannot repeatedly negate an already negated sequence.

4. **Original-data recomputation:** use a separate copy of the repository, prepare all external files, and execute from the repository root in this order:

```bash
python code/reproduce_fig1_fig2.py --s-dir data/raw/10hpa --t-dir data/raw/500hpa --surrogates 999
python code/reproduce_fig3.py --s-dir data/raw/10hpa --t-dir data/raw/500hpa
python code/reproduce_fig4.py
python code/reproduce_fig5.py --s-dir data/raw/10hpa --t-dir data/raw/500hpa --v-dir data/raw/v100 --heat-t-dir data/raw/t100
python code/reproduce_figS1_figS2.py
python code/reproduce_figS3.py --s-dir data/raw/10hpa --t-dir data/raw/500hpa
python code/reproduce_figS4_S5_S6.py
python code/reproduce_figS7_S8.py --event-file data/raw/SSWevent.xlsx
python code/reproduce_figS9.py --s-dir data/raw/10hpa --t-dir data/raw/500hpa
python reproduce.py
```

Run each script with `--help` for supported overrides. The final S/T drawing data share `code/fig1_fig2_results`; recomputation reference inputs reside in the sibling `reference_inputs` directory. Explicit custom result directories use a sibling `reference_inputs` directory too. Defaults of S3 and S4–S9 now point to the same canonical SI folder. Full raw-data execution is computationally expensive and has not been repeated as part of preparing this release. The full Fig.1/2 entry builds oriented reference/predictor caches and also invokes the final unoriented Fig.1 workflow; the latter retains standard Pearson correlation.

## Core method notes

Climatology is calculated by calendar month-day at each grid point. The anomaly matrix is centered again on the basis window before eigendecomposition. If `A=WΣUᵀ`, spatial vectors are columns of `U`, covariance is `AᵀA/(N−1)`, and a daily score is `(anomaly−basis_mean)·U_j`. For final Fig.1 the raw eigensolver directions are kept, then both score and spatial direction are negated. Other prediction branches retain their original fitted orientation and fold-specific training rules.

SPV preprocessing area-weights the 10 hPa polar anomaly by cos(latitude), fits separate cross-year intercept and slope for each November–April calendar month, subtracts those monthly fitted trends and standardizes the pooled residual with ddof=0. This note describes the existing implementation, not a change in preprocessing.

The fixed-lag surrogate procedure preserves the observed lag-pair masks. It randomizes Fourier phases independently within the original 93 continuous SPV segments, preserves the DC and even-length Nyquist phases, and recenters each surrogate on its own valid paired dates. For Fig.1 the signed correlation with largest absolute value is selected separately for each mode and surrogate over lags 0–60. Strict comparison preserves the smallest lag in exact ties. Fig.S1 does not perform this selection for its envelope.

Prediction and inference details retained from source docstrings are in METHODS_AND_INTERFACES.md; exact reproducible algorithms are in the final scripts. No test here establishes causal identification or adds a Granger interpretation beyond the implemented analyses.

## Current result and reference paths

`S_final_pearson.npz` retains the original raw vectors, raw scores, negated scores and final Pearson inference. `S_top100.csv` is its final 100-mode export. The oriented S reference in `code/reference_inputs/S_basis_and_correlations.npz` is used by Fig.4 recomputation and independent numerical checks; it must not replace the final S cache. `spv_daily.csv` in that reference directory supplies SPV for cached-inference recomputation.

The original-data Fig.1/2 entry writes S reference/SPV data to the sibling reference directory, final S/T data to the selected results directory, PDFs/PNGs to the associated `figures` tree and layout/provenance records to the associated `validation` tree. It no longer saves the unused climatology caches or the legacy oriented-S correlation CSV. A standalone raw-S run requires the reference/SPV files generated by the complete Fig.1/2 entry.

PDF validation accepts `--layout-dir DIRECTORY` on both `verify_release.py` and `code/verify_publication.py`; without it, the shipped `validation/layout`, fresh redraw `outputs/validation/layout`, and legacy PDF-adjacent records are resolved automatically.
