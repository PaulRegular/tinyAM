# Dover sole off the US West Coast

## Stock identity and latest detailed assessment

NOAA Stock SMART ID 10022 is *Microstomus pacificus* assessed as one US West Coast unit from the US-Canada to US-Mexico borders. It is distinct from Dover sole in Canadian waters. No matching stock is present in the canonical database. The PFMC's latest detailed full assessment located is the 2021 benchmark, with input data through 2020. Recent harvest-specification projections continue to use the 2021 assessment; the PFMC listing reviewed on 2026-10-29 does not identify a later full assessment.

The 2021 assessment report, model-file ZIP and STAR review are cached in `source_cache/noaa_pacific_coast_dover_sole_20261029/`. Its `sources.csv` records original URLs, retrieval date, file sizes and SHA-256 checksums. The model bundle identifies Stock Synthesis 3.30.16.02, the input files `2021_dover.dat` and `2021_dover.ctl`, and fitted outputs including `Report.sso`, `ss_summary.sso` and `wtatage.ss_new`.

## Inputs and assumptions recovered

The native data file spans 1911-2020 and contains annual retained-catch series for California and Oregon/Washington, discard estimates, four survey-index series, and age-composition records. The fishery age compositions are by year, sex, and fleet. The assessment used ages up to 60 and two sexes; the source document notes that new fishery age readings were limited after 2010, including no California otoliths read after 2009. This is age-composition input rather than a complete observed catch-numbers-at-age series. Catch-at-age for a tinyAM translation would therefore require a documented reconstruction using catch totals, age compositions and age-specific weights.

Survey age information is present in the native data as conditional age-at-length records, including the NWFSC West Coast Groundfish Bottom Trawl Survey and NWFSC Slope Survey, alongside aggregate survey indices. The data file also contains derived marginal “ghost” age compositions for the WCGBTS; these are not independent observations and must not be treated as an extra likelihood. The assessment report describes the survey age data and their use. Recovering an age-specific abundance index for tinyAM would require a clearly labeled conversion from the aggregate survey index and age information, with the conditional-age-at-length structure and likelihood differences retained in the audit.

The assessment is size-structured and sex-specific. It uses length-based growth, weight, maturity and fecundity relationships, not a directly observed fixed weight-at-age vector. The model output includes derived weight-at-age (`wtatage.ss_new`) and biology tables. Female natural mortality is fixed at 0.108 per year; male mortality is estimated as an offset from the female value. Steepness is fixed at 0.8 and recruitment standard deviation at 0.35. The model input defines survey timing by month (AFSC Slope 10, Triennial 7, NWFSC Slope 9, WCGBTS 8); these are source settings, not a separate claim that each survey occurs on one exact calendar day.

## Outputs and screening decision

The fitted report contains historical numbers-at-age, catch-at-age, F-at-age and biomass-at-age tables, in addition to spawning biomass, recruitment, survey fits and other derived quantities. These model outputs are not substitutes for missing observed catch-at-age inputs. The fit summary reports final phase 6, 1,185 iterations and a final gradient of 5.49e-05. The STAR review recommended the assessment as the best available science and described the model's overall convergence diagnostics as good, while noting age-data limitations and lack of fit to some recent composition data.

The candidate is marked `new`, `screened` and `eligible`, with import not started. The current evidence supports importing the native inputs and reported outputs. Any tinyAM conversion must be explicitly identified as an approximation: tinyAM does not reproduce the two-sex size-structured Stock Synthesis likelihood, conditional age-at-length likelihood, ageing-error treatment, or the full survey/fishery structure. Eligibility here means the age-composition and abundance information plus official outputs are recoverable for a documented age-structured approximation; it does not imply a like-for-like replication.

## Sources

- [Status of Dover Sole along the US West Coast in 2021](https://www.pcouncil.org/documents/2021/10/status-of-dover-sole-microstomus-pacificus-along-the-u-s-west-coast-in-2021.pdf)
- [Dover Sole 2021 Stock Synthesis model files](https://www.pcouncil.org/documents/doversole-2021-model-files-zip-file/)
- [Dover Sole 2021 STAR panel report](https://www.pcouncil.org/documents/dover-sole-stock-assessment-review-star-panel-report/)
- [PFMC groundfish assessment documents](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)
- [NOAA Stock SMART ID 10022](https://apps-st.fisheries.noaa.gov/stocksmart/?app=browse_by_stock&stockid=10022)