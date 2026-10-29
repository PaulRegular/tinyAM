# Widow rockfish off the US West Coast

## Stock identity and latest detailed assessment

NOAA Stock SMART ID 10568 is *Sebastes entomelas* off the US West Coast. The 2025 detailed assessment treats this as one coastwide unit, distinct from Canadian and Alaska stocks. No matching stock is present in the canonical database. The NOAA inventory snapshot labels its metadata row 2026; the latest detailed official assessment located is the 2025 update, with model data through 2024. The supplemental review recommended the October 2025 model as best scientific information, and the Council's later assessment listing includes the 2025 report and model files.

The 2025 assessment, official model-file ZIP and supplemental review are cached in `source_cache/noaa_pacific_coast_widow_rockfish_20261029/`. Its `sources.csv` records URLs, retrieval date, sizes and SHA-256 checksums. The model bundle includes `2025widow.dat`, `2025widow.ctl`, `Report.sso` and `ss_summary.sso`; the data file spans 1916-2024 and models two sexes and ages 0-40+.

One run-version discrepancy needs to be preserved in any source import: the assessment report identifies SS 3.30.24.1, the supplemental review's model summary says SS 3.30.23.1, and the bundled `ss_summary.sso` identifies SS 3.30.23.2 with a 2025-12-03 run date. The bundled output's estimated female and male M values (0.135 and 0.147) match the published base-model estimates, but that alone does not resolve the version labels. Keep the bundle's run identity and the discrepancy explicit.

## Inputs and assumptions recovered

The October 2025 model uses three catch sources (bottom trawl, midwater trawl and at-sea hake fishery), six abundance indices (three survey and three fishery-dependent), marginal age compositions for three commercial fleets, and conditional age-at-length data for WCGBTS. The native model data provide these input structures; fishery age compositions are sampled data, not a full annual catch-numbers-at-age matrix. A tinyAM catch-at-age input would need a documented reconstruction using catches, age compositions and biological weights. The WCGBTS conditional age-at-length data are not a direct age-specific abundance index; conversion from survey indices and age information would be approximate and must avoid double-counting derived marginal compositions.

The assessment uses an age-structured, single-area, two-sex Stock Synthesis model. It estimates sex-specific growth and natural mortality, with a lognormal prior for M (median 0.10 per year, log-SD 0.31). M estimates in the base model are 0.135 for females and 0.147 for males. Age at 50% maturity is 5.47 years, carried from the 2015 benchmark, and the model uses externally estimated length-weight and fecundity relationships. Steepness is fixed at 0.72 and recruitment SD at 0.60. Survey, fishery, sex, age-error, and timing distinctions are represented in the native model files and require explicit approximation decisions in a tinyAM translation.

## Outputs and screen decision

The fitted `Report.sso` contains numbers-at-age, catch-at-age, F-at-age, biomass-at-age, and other model tables; the assessment reports recruitment and spawning output. Spawning output is reported in billions of eggs, so it should not be silently labeled or compared as spawning biomass. The published report says the Hessian was positive definite, asymptotic variance estimates were available, the final gradient was low, and 57 of 100 jittered initializations returned to the same objective value.

The candidate is marked `new`, `screened` and `eligible`, with import not started. The 2025 bundle supports source import and a documented age-structured approximation. Eligibility does not imply a like-for-like replication: tinyAM would not reproduce the two-sex, multi-fleet Stock Synthesis likelihood or its conditional age-at-length likelihood. The precise age-composition data, age-specific outputs, official-run version discrepancy and native units should be captured at import.

## Sources

- [Status of widow rockfish off the US West Coast in 2025](https://www.pcouncil.org/documents/2026/05/status-of-widow-rockfish-stock-off-u-s-west-coast-in-2025.pdf)
- [Widow Rockfish 2025 Model Files](https://www.pcouncil.org/documents/widow-rockfish-2025-model-files/)
- [Supplemental Review Panel Report for Widow Rockfish](https://www.pcouncil.org/documents/2025/10/supplemental-review-panel-report-for-widow-rockfish.pdf)
- [PFMC groundfish assessment documents](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)
- [NOAA Stock SMART ID 10568](https://apps-st.fisheries.noaa.gov/stocksmart/?app=browse_by_stock&stockid=10568)