# NOAA Pacific Coast Canary Rockfish

## Stock identity and assessment run

NOAA Stock SMART ID 10415 is the Sebastes pinniger Pacific Coast management stock in U.S. waters. It is distinct from the DFO British Columbia coastwide canary rockfish unit, which covers Canadian waters. Do not crosswalk these records solely because they share a species name.

The Stock SMART inventory row is labeled 2025, but the latest full stock assessment located is the 2023 assessment, published in 2024 and fitted to data through 2022. NOAA's species page identifies 2023 as the most recent Pacific Coast assessment. Later harvest actions use projections associated with the 2023 assessment; they do not establish a newer full age-structured assessment. The run identity is confirmed by the assessment report and the bundled Stock Synthesis files.

The official [NOAA report record](https://repository.library.noaa.gov/view/noaa/72904) identifies Status of Canary Rockfish (Sebastes pinniger) along the U.S. West Coast in 2023. The [PFMC assessment archive](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/) lists the report and its 2023 model-file archive. The model ZIP is cached locally under analysis/comp_assessments/source_cache/discovery_noaa_canary_2023/; its sources.csv records retrieval date, URL and checksum. NOAA's PDF download URL was attempted but returned an access-denied response, so the report PDF itself is not cached.

## Source-input screen

The archived Stock Synthesis 3.30.21 run contains data.ss, control.ss, starter.ss, Report.sso, CompReport.sso, and uncertainty/output files. The data file specifies one annual season, two sexes, a 40-year model maximum age, and 31 fleets. It contains fishery catch totals, biomass survey/index observations and age-composition data. The age-composition section uses age bins 1–35; the fitted population extends to age 40. The input file includes 3,829 age-composition records across fishery and survey fleets, with observed years extending from 1980 to 2022. Some observations are conditional age-at-length records rather than marginal age compositions.

Catch-at-age is not supplied as a single observed series. For fleets and years with marginal age compositions, catch totals and number compositions support a documented reconstruction, with source weights-at-age needed to convert catch biomass to numbers. For conditional age-at-length observations, the reconstruction must also combine the relevant length compositions and age-length information. Do not fill unsampled fleets or years from fitted states.

Survey observations are aggregate indices accompanied by survey age-composition data. Survey age-specific abundance can be approximated by combining the aggregate biomass index, observed age composition and source weights-at-age; this is a derived index, not a direct observation-at-age. Survey and composition observations are recorded at month 7 (July). The one-season model and month-coded observations provide a timing basis; retain the exact source month when representing survey timing downstream.

The control and report files provide source growth, weight-length and maturity specifications and fitted biological-at-age/weight-at-age surfaces. The 2023 report describes separate female and male natural-mortality rates that are constant over age. These inputs and settings are recoverable from the model archive, though a canonical import must distinguish original observations, parameterized biology and fitted derived quantities.

## Output screen and decision

Report.sso advertises and contains official numbers-at-age, catch-at-age, F-at-age, annual time-series, biology and natural-mortality tables. The report records the assessment run's final gradient and Hessian, and the archive includes uncertainty files. This is sufficient for an auditable comparison without treating tinyAM estimates as source data.

Screening decision: eligible for source curation. The archive provides source quantities for documented age-specific catch and survey-index reconstructions, biological and mortality specifications, survey timing and comparable fitted outputs. The reconstructed series must remain labeled as derived, and the assessment's composition likelihood is more detailed than tinyAM's input representation. No tinyAM translation or fit was attempted in this inventory-screening step.

## References

- [NOAA assessment report record](https://repository.library.noaa.gov/view/noaa/72904)
- [NOAA Canary Rockfish species and assessment summary](https://www.fisheries.noaa.gov/species/canary-rockfish)
- [PFMC groundfish assessment documents and 2023 model archive](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)
- [NOAA Stock Synthesis 3.30 User Manual](https://repository.library.noaa.gov/view/noaa/23884/noaa_23884_DS1.pdf), for interpreting the month field in age-composition and index observations.