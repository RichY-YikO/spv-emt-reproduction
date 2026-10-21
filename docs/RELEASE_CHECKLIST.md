# v1.0.0 archive candidate checklist

Prepared locally on 2026-10-19. No Zenodo/GitHub login, upload or publication has occurred. Preparation date is not a public release date.

## Confirmed by the author

- [x] MIT copyright holder: Shuo Yang; year: 2026. Formal LICENSE supplied.
- [x] Sole software author: Shuo Yang; family name Yang, given name Shuo. Formal CITATION.cff supplied.
- [x] Archive title: Code and reproduction materials for “State-resolved multimodal contributions to stratospheric polar vortex predictability”
- [x] Version: v1.0.0.
- [x] ERA5 source: Copernicus Climate Change Service (2018), ERA5 hourly data on pressure levels from 1940 to present, https://doi.org/10.24381/cds.bd0915c6.
- [x] Computational Python files, caches, CSVs and 14 final PDF/PNG pairs remain unchanged; source hashes and updated checksums verified.

## Pending before public release

- [x] SSW event date provenance confirmed: calculated by Shuo Yang following Charlton and Polvani (2007), DOI 10.1175/JCLI3996.1. See SSW_EVENT_PROVENANCE.md.
- [ ] Retain applicable source-data and figure redistribution terms. Event-date authorship/provenance is confirmed; the original software MIT license does not automatically license third-party material.
- [ ] Review the final Zenodo draft and confirm the actual public-release date. It is intentionally omitted from CITATION.cff until known.

## Optional or later

- ORCID and affiliation may be added if supplied; neither is invented.
- GitHub repository owner/URL is optional for direct Zenodo upload.
- Paper DOI may be added when available. The ERA5 dataset DOI must not be used as this software's DOI.
- Zenodo software DOI may be reserved later, then copied into CITATION.cff and README before generating the final upload ZIP.
- The documented fresh-environment install and full raw-data calculation have not been rerun; do not claim these tests were completed.

After reviewing intentional metadata changes, run `python verify_release.py --refresh-checksums`, then `python verify_release.py`, and regenerate the ZIP. The source data/figure reference hashes are not rewritten by this command. Keep the current original draft folder/ZIP as a backup.
