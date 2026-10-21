# SPV EMT reproduction

Local archive candidate for version `v1.0.0`, prepared 2026-10-19. It has not been published or uploaded; this preparation date is not a publication date.

Archive title: **Code and reproduction materials for “State-resolved multimodal contributions to stratospheric polar vortex predictability”**

Software author and copyright holder: **Shuo Yang**.

This repository contains the final analysis and plotting code, frozen numerical results, CSV tables, and 14 publication figures for a stratospheric polar vortex (SPV) eigen microstate reproduction. Software authorship, copyright, archive title and version are confirmed by the author. The SSW event dates were calculated by Shuo Yang following Charlton and Polvani (2007). Repository URL and paper/software DOI are omitted until available.

## Quick start

Use Python 3.10. The original numerical environment used Python 3.10.13.

```bash
python -m venv .venv
```

Activate with `.venv\Scripts\activate` on Windows or `source .venv/bin/activate` on Linux/macOS, then run:

```bash
python -m pip install -r requirements.txt
python verify_release.py
python reproduce.py
```

Alternatively: `conda env create -f environment.yml`, then `conda activate spv-emt`.

The shipped final files are in `figures/pdf/` and `figures/png/`. Layout records are stored in `validation/layout/`, separately from figure files. Redrawing writes to `outputs/pdf/`, `outputs/png/` and `outputs/validation/layout/` and does not replace the shipped final figures or recompute scientific caches.

```bash
python reproduce.py --figures Fig1 FigS1
python reproduce.py --figures Fig3 Fig4 Fig5 --output-dir my-output
python code/redraw_final.py --figures Fig2
python reproduce.py --help
```

`Fig1a` optionally draws only the top Fig.1 panel. Otherwise figure names are `Fig1`–`Fig5` and `FigS1`–`FigS9`. The English wrapper preserves figure selection without requiring a non-ASCII script name.

## Contents

| Location | Content |
|---|---|
| `code/` | Comment-free final computation, drawing and verification modules |
| `code/*_results/` | Required frozen NPZ results and CSV tables |
| `code/fig1_fig2_results/` | Final S/T drawing caches and three CSV exports; final S cache is `S_final_pearson.npz` |
| `code/reference_inputs/` | Oriented S reference and SPV daily values required for recomputation and numerical verification |
| `figures/TableS1.xlsx` | Editable Table S1, checked against `table_S1.csv` |
| `validation/layout/` | Figure layout records used by PDF verification |
| `figures/pdf/`, `figures/png/` | 14 original final vector PDFs and 600 dpi previews |
| `docs/` | Methods, raw-data contract, figure map, release and upload instructions |
| `validation/` | Source-file hashes and executed validation reports |

Historical backups, duplicate result directories, the obsolete custom-statistic cache and scripts that require historical local backups are excluded. The oriented S reference cache is retained only in `code/reference_inputs/` for recomputation and independent verification. Final Fig.1 and Fig.S1 use `S_final_pearson.npz`; `S_top100.csv` now describes that final Pearson result. Unused climatology caches and the legacy oriented-S correlation CSV have been removed.

## Scientific convention

For Fig.1 and Fig.S1, the cached unoriented eigenvectors and raw scores are retained. Each valid lag-paired sample uses `x = -a_raw` and is centered by `mean(x)`:

```text
xc = -a_raw - mean(-a_raw) = -a_raw + mean(a_raw)
yc = SPV - mean(SPV)
r  = sum(xc * yc) / sqrt(sum(xc**2) * sum(yc**2))
```

This is standard Pearson correlation, exactly the negative of Pearson correlation computed from the same original scores and paired SPV. There is no additional plot-time sign reversal. Fig.1 spatial panels use `-U_raw[:, :6]`. All six Fig.1 spatial panels use the same 1.56-inch nominal diameter as Fig.2, in a 3-column by 2-row arrangement.

Fig.1 selects the signed maximum absolute correlation over lags 0–60 days for both the observation and each surrogate. Fig.S1 uses the fixed-lag surrogate envelope over lags −60–60 days. Both use 999 phase surrogates and seed `2026091206`. These envelopes are distinct. The superseded non-Pearson statistic is not used anywhere in the release workflow.

## Verification

```bash
python verify_release.py
python verify_release.py --mode numerical
python verify_release.py --mode si
python -m pip install -r requirements-validation.txt
python verify_release.py --mode pdf
python verify_release.py --mode all
```

The numerical audit independently checks all 12,100 Fig.1 lag correlations and both envelope types; it does not refit the eigenbasis. SI verification checks plotted data and envelope boundaries against the included caches. PDF verification checks embedded fonts, vector paths, absence of raster image objects and panel labels. `verify_release.py --mode pdf --pdf-dir outputs/pdf` audits a fresh redraw. `validation/RELEASE_VALIDATION.json` records the earlier scientific checks plus the current candidate metadata/integrity verification. Full raw-field recalculation and fresh environment installation are not claimed by these checks.

## Raw data and portability

Raw NetCDF files and the original SSW event workbook are not bundled. No raw-data login is needed for cached redrawing. All source-data defaults are relative to `data/raw/`; original CLI flags such as `--s-dir`, `--t-dir`, `--v-dir`, `--heat-t-dir` and `--event-file` can override them. See [raw-data requirements](docs/DATA_AND_REPRODUCTION.md) before attempting an original-data calculation.

Cartopy may download Natural Earth coastline data on its first map draw. Provide network access or prepopulate its data directory. The final figures use Times New Roman, which is not redistributed; install it locally for matching typography. Without it, Matplotlib may choose a fallback. Exact binary output depends on fonts, rendering libraries and metadata. Unoriented eigenvectors also have mathematically arbitrary signs; reproducing a fresh eigensolver output in a different numerical environment may choose different signs. The included frozen cache preserves the confirmed final sign convention.

## License, citation and event provenance

The original software code is licensed under the MIT License: Copyright (c) 2026 Shuo Yang. The software author in `CITATION.cff` is Shuo Yang. Raw data, third-party materials and fonts retain their own terms; the software license does not automatically license included data or figures. ERA5 data source: Copernicus Climate Change Service (2018), ERA5 hourly data on pressure levels from 1940 to present, https://doi.org/10.24381/cds.bd0915c6.

The SSW event dates were calculated by Shuo Yang using the identification method of Charlton and Polvani (2007; https://doi.org/10.1175/JCLI3996.1). This provenance statement is confirmed by the author; this documentation update does not rerun event detection. See [SSW event provenance](docs/SSW_EVENT_PROVENANCE.md) for references and the relationship to the supplied event files. The original software MIT license does not automatically license data and figures; their respective terms should be retained. This is a local candidate, not a completed Zenodo deposit. No DOI or release date has been invented. See [the checklist](docs/RELEASE_CHECKLIST.md) and [Zenodo instructions](docs/ZENODO_UPLOAD_zh-CN.md).

Chinese instructions: [README_zh-CN.md](README_zh-CN.md). Upload guide: [docs/GITHUB_UPLOAD_zh-CN.md](docs/GITHUB_UPLOAD_zh-CN.md).

## Directory cleanup on 2026-10-20

Scientific arrays and all 14 shipped PDF/PNG figures are unchanged. The final Fig.1/2 data directory contains only `S_final_pearson.npz`, `T_basis_and_correlations.npz`, `S_top100.csv`, `T_top100.csv` and `table_S1.csv`. The Fig.1 calculation record is in `validation/fig1_pearson/`. Historical dated validation reports describe their original execution; current paths and checksums are recorded in `validation/source_data_sha256.json` and `CHECKSUMS.sha256`. See `docs/DIRECTORY_LAYOUT.md`.
