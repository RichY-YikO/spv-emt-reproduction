# Changelog

## v1.0.0 - Fig.1 spatial-panel size update, 2026-10-19

- Matched all six Fig.1 spatial panels to the Fig.2 nominal diameter of 1.56 inches, retaining the 3-column by 2-row arrangement.
- Compacted Fig.1 vertical layout; the top correlation values and scientific data are unchanged.
- Updated the shared drawing code, Fig.1 vector PDF/600 dpi PNG and reference checksums. All Fig.2–5 and Fig.S1–S9 files and all NPZ/CSV results are unchanged.

## v1.0.0 - SSW provenance documentation update, 2026-10-19

- Author confirms that Shuo Yang calculated the SSW dates following Charlton and Polvani (2007), DOI 10.1175/JCLI3996.1.
- Resolved current SSW-source placeholders in both READMEs, data documentation, the release checklist and Zenodo description/instructions.
- Preserved all Python source, event-date CSVs, NPZ caches, final figures and archive directory structure; no detection, prediction or plotting was rerun.

## v1.0.0 - local candidate prepared 2026-10-19

- Confirmed sole software author and copyright holder Shuo Yang.
- Added the formal MIT LICENSE and CITATION.cff, with the author-confirmed archive title and version.
- Added the author-supplied ERA5 pressure-level data citation and DOI.
- Added Zenodo copy-and-paste metadata and upload instructions; no account login or upload performed.
- Kept SSW event catalogue provenance explicitly pending. No actual publication date, ORCID, repository URL or software/paper DOI is invented.
- All Python source, numerical caches, CSVs and final figures are unchanged from the validated draft.

## 1.0.0-draft - 2026-10-18

- Prepared an independent local release package from the confirmed final figures and caches.
- Removed explanatory Python comments and module/function/class docstrings while preserving algorithms and command-line help strings.
- Replaced workstation-specific raw-input defaults with repository-relative defaults and retained explicit CLI overrides.
- Added an English final-redraw entry point, integrity verification, dependencies, data/method documentation, publication templates and upload instructions.
- Included the final standard-Pearson Fig.1/Fig.S1 cache, unchanged final 14 PDFs and PNGs, and required reference/forecast/SI caches.
- Excluded obsolete custom-statistic code/cache, duplicate SI outputs, local backups and historical-backup-dependent validators.
- Authorship, copyright ownership, verified upstream data provenance and formal public-release metadata remain pending; no remote repository has been created or populated.

## 2026-10-20 directory cleanup for v1.0.0 candidate

- Selected final S/T drawing data; removed unused climatologies and the legacy S summary.
- Separated required S/SPV reference inputs and moved layout/provenance JSON out of figures/data directories.
- Included original Table S1 Excel; updated computation, redraw and validation paths without changing scientific caches or shipped figures.
