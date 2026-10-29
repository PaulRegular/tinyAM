# NOAA West Coast Sablefish — 2025 assessment

## Stock and assessment

NOAA Stock SMART ID 10466 is the U.S. West Coast Sablefish management unit.
It is separate from the Alaska and British Columbia units, despite sharing the
species name (*Anoplopoma fimbria*). No canonical U.S. West Coast Sablefish
record was present in the database. The register marks it as a provisional new
unit, subject to later cross-inventory reconciliation.

The selected assessment is the 2025 official assessment, which fits data
through 2024. The Pacific Fishery Management Council provides both the detailed
assessment report and a native Stock Synthesis model archive. The assessment is
a two-sex, age-structured population model with ages tracked from 0 to an
accumulated age of 70. The observation age bins span 0-50, with age 50 used as
a plus group.

The report states Stock Synthesis version 3.30.24. The bundled input echo and
`Report.sso` identify version 3.30.23.1, compiled December 5, 2024. The bundle
is titled for the 2025 assessment, has input data through 2024, and reports
outputs for that run, but the version discrepancy should be resolved or
retained explicitly during import.

## Evidence for age-structured source data

- Catch input is biomass in metric tonnes, divided among trawl, hook-and-line,
  and pot retained fleets and corresponding discard fleets. Retained-catch
  age compositions are sampled by year and sex: trawl ages begin in 1983, and
  hook-and-line and pot ages begin in 1986. Discard age compositions are
  derived in the assessment from length observations using a length-age key.
  These inputs do not provide observed catch numbers-at-age for every year.
  A downstream catch-at-age reconstruction would need to use catch totals,
  observed age compositions, and appropriate weights, and must leave years
  without age samples visibly limited rather than fill them with fitted states.
- Biomass indices are supplied for the early and late Triennial surveys, the
  NWFSC Slope Survey, and the NWFSC West Coast Groundfish Bottom Trawl Survey
  (WCGBTS). Native age-composition observations are also supplied for surveys
  in 0-50 bins. WCGBTS age samples cover 2003-2024 except 2020; the historical
  survey composition coverage is intermittent. The age-specific survey
  abundance index is not directly observed. A tinyAM age-indexed series would
  require an explicit transformation of aggregate index and age-composition
  observations, with the resulting assumptions recorded.
- The model uses empirically estimated sex-specific weight-at-age from WCGBTS
  samples for ages 0-30 during 2003-2024, with age 30 as the plus group.
  Earlier weight-at-age is constant. A fixed spatially integrated maturity
  schedule is used in the model.
- Natural mortality is estimated as a constant shared by sexes. The base
  estimate is 0.088 per year, with a lognormal prior whose median is 0.072 and
  standard error is 0.31. The report identifies mortality, absolute stock
  scale, and productivity as important sources of uncertainty.
- The native data files contain observation-month values for surveys. The
  sampled month is July for the early Triennial, NWFSC Slope, and WCGBTS data,
  and August for the late Triennial period; the model spawning month is
  January. A fractional-year conversion for tinyAM must state how these
  monthly conventions are represented.
- `Report.sso` enables annual numbers-at-age, catch-at-age, F-at-age, and
  biomass-at-age output tables. The report also provides annual spawning
  output, recruitment, and exploitation estimates with uncertainty. It reports
  2025 spawning output of 74,409 mt with a 95% interval of 59,331-89,486 mt.

The native bundle and report make this a source-ready age-structured case.
The distinction between observed catch/index inputs and assessment-derived
age surfaces remains important for any later tinyAM translation.

## Sources and local cache

The report and native model archive were retrieved on 2026-10-29 and are
cached under `source_cache/discovery_noaa_sablefish_2025/`. Their URLs, sizes,
and SHA-256 checksums are in the ignored `sources.csv` manifest.

- [2025 West Coast Sablefish assessment report](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/)
- [2025 West Coast Sablefish native model files](https://www.pcouncil.org/documents/sablefish-2025-model-files/)
- [PFMC groundfish assessment document index](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)
