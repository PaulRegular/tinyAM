# NOAA West Coast Yelloweye Rockfish — 2025 assessment

## Stock and assessment

NOAA Stock SMART ID 10426 is the U.S. West Coast Yelloweye Rockfish unit. It
is distinct from the Gulf of Alaska unit (ID 10425) and the Canadian inside
and outside Yelloweye units. The register treats it as a provisional new
unit pending cross-inventory reconciliation.

The 2025 update was the adopted assessment used for management advice and
uses data through 2024. Its native Stock Synthesis files and report identify
version 3.30.23.2. This is a single-sex model with two areas, California and
Oregon/Washington, and a January spawning month.

## Evidence and limits

- The native input has seven fishery catch fleets: six use biomass catch in
  metric tonnes and Washington recreational catch is in numbers. The model
  uses fleet-specific length compositions and conditional age-at-length
  observations, plus marginal age compositions for some fisheries and
  surveys. The report describes age samples and data updates through 2024.
  These likelihood data are not a complete annual catch-at-age series. Any
  tinyAM catch-at-age input would need a documented reconstruction from
  catches, available age compositions, and weights; missing composition
  years must remain missing rather than be filled with model estimates.
- The model includes recreational and survey abundance indices, length
  compositions, and survey conditional age-at-length observations, including
  WCGBTS and IPHC observations. It does not directly observe abundance at age.
  A tinyAM age-indexed index would need to be reconstructed from aggregate
  indices and compositions and clearly labeled as approximate.
- The fitted growth curve and the assessment's length-weight relationship
  determine weight-at-age. The native report includes the resulting age-based
  output. Maturity is fixed using a length-based schedule, not estimated as a
  maturity-at-age series. The source materials document the ingredients for
  deriving the corresponding age-specific biological surfaces.
- Natural mortality is fixed at 0.044 per year, as in the 2017 benchmark.
  The report describes the fixed maturity and mortality choices and their
  contribution to uncertainty.
- The native input records observation month for survey data and a January
  spawning month. Preserve each series' recorded month when converting survey
  timing to a fractional year.
- `Report.sso` provides annual numbers-at-age, catch-at-age, biomass-at-age,
  and F-at-age. The report gives spawning output and recruitment summaries
  with uncertainty. It reports 2025 spawning output of 477.63 billion eggs
  (95% interval 384–571 billion eggs).

This is eligible for source curation under the age-input screen because the
official native files expose age-composition inputs, fitted biological
surfaces, and age-specific outputs. That decision does not mean the native
likelihood can be represented exactly in tinyAM: catch-at-age and
age-specific abundance observations require explicit transformations, and
the assessment's two-area, single-sex structure is an additional
approximation.

## Sources and local cache

The report and native model archive were retrieved on 2026-10-29 and are
cached under `source_cache/discovery_noaa_yelloweye_2025/`. Their URLs, sizes,
and SHA-256 checksums are listed in the ignored `sources.csv` manifest.

- [2025 West Coast Yelloweye Rockfish assessment report](https://www.pcouncil.org/documents/2026/01/status-of-yelloweye-rockfish-off-the-u-s-west-coast-in-2025.pdf/)
- [2025 Yelloweye Rockfish native model files](https://www.pcouncil.org/documents/yelloweye-rockfish-update-assessment-2025-model-files/)
- [PFMC groundfish assessment listing](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)
- [Pinned PFMC model repository revision](https://github.com/pfmc-assessments/Sebastes_ruberrimus_2025/commit/9356e8267c23a1a5d33bb76f261fcd5622db062c)
