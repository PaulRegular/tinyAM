# Pacific ocean perch in the Gulf of Alaska

## Unit and current assessment

NOAA Stock SMART ID 11056 is the Gulf of Alaska (GOA) Pacific ocean perch unit (*Sebastes alutus*). It is distinct from NOAA's Bering Sea/Aleutian Islands unit (Stock SMART ID 11055) and the British Columbia DFO assessment. No existing canonical GOA stock record or cross-organization match was found.

Stock SMART's inventory year is 2024, but the latest completed detailed fit located is the 2023 GOA assessment update, Model 20.1, using data through 2023. The report calls this an operational update: its configuration matches the model accepted in 2021, with updated inputs. The 2024 SAFE schedule lists 2023 as the last full model and 2024 as a harvest projection. The 2025 catch report says the scheduled full assessment was cancelled, and the June 2026 SSC review plan schedules the next GOA POP assessment for December 2026. No later fitted assessment was available by the review date, 2026-10-29.

The 2023 SAFE, 2024 SAFE introduction, 2025 catch report, 2026 SSC review plan, and selected native model files are cached locally under `source_cache/noaa_goa_pop_20261029/`. The model repository files are pinned to commit `e87bbfc60ab8bfc87d05e02ed1b4b39b76059d7e`; `sources.csv` records URLs, retrieval dates, file sizes and SHA-256 checksums.

## Inputs and source screen

The model data file contains total fishery catch in tonnes for 1961–2023. Catch through 2022 is reported as observed; 2023 is an in-year estimate. The source file notes that early catch values were frozen and that the later series uses AKFIN catch data. Fishery age compositions are number proportions for 17 years: 1990, 1998–2002, 2004–2006, 2008, 2010, 2012, 2014, 2016, 2018, 2020 and 2022. Observed bins are ages 2–24 and 25+, with raw sample counts and hauls in the native input. Fishery weight-at-age is available separately in the assessment repository for ages 2–29, in grams per fish. These quantities support a documented conversion to catch numbers-at-age for observed composition years; unsampled catch years must remain missing rather than be filled from fitted states.

The bottom-trawl survey biomass series is in tonnes, with standard errors, for 16 survey years from 1990 to 2023. Survey age proportions and sampling counts are available in 14 years: 1990, 1993, 1996, 1999, and 2003–2021 on the biennial survey schedule. Survey weights-at-age are available in grams per fish. Together, the aggregate biomass, observed age proportions and survey weights support a documented age-specific abundance-index reconstruction for years with both survey inputs. The 2001 biomass estimate includes substituted Eastern GOA biomass because that area was not sampled.

There is one source-export discrepancy to retain in the audit trail: the final 2023 model `.dat` file and SAFE data summary include 1990 survey age composition, while the repository's auxiliary `goa_bts_age_comp.csv` export begins in 1993. Use the final model input and report as the run-matched sources; verify the 1990 row against the original model `.dat` during any import.

Observed compositions use ages 2–25+, while the assessment population model and reported numbers-at-age span ages 2–29. This is an important age-grid difference for a tinyAM comparison; the terminal ages need a documented aggregation, and the assessment's fitted age states must not be used to create observations. The report's age-composition likelihood uses sample-size weights derived from the square root of the nominal sample count. The assessment also fits fishery length compositions from 1963–1977 and 1991–1997; these are not a tinyAM age input and would be omitted from a translation.

The assessment estimates one constant natural-mortality parameter with a prior (prior mean 0.0614, prior CV 0.1). Maturity is estimated internally with a logistic age curve using two female maturity datasets. The supplied age-specific weight vectors are in grams per fish. The model's survey-biomass prediction uses abundance at the survey year without a separate within-year sampling fraction; the `.dat` value `spawn_fract = 5` is spawning month, not survey timing. The actual survey timing fraction is not established by these model files, so any tinyAM fraction must be clearly identified as an approximation to the source model's annual indexing.

## Outputs and eligibility

The pinned repository contains model-based numbers-at-age for 1961–2023 at ages 2–29, plus annual recruitment, fishing mortality, total biomass and spawning biomass. These are point estimates. The SAFE report and model repository also provide MCMC uncertainty for selected aggregate quantities and parameters, but not a complete set of age-specific abundance intervals in the exported numbers-at-age table.

This unit is provisionally eligible for source curation and a later translation attempt. The source provides catch totals, fishery age proportions and fishery weights; aggregate survey biomass, survey age proportions and survey weights; age-specific biological inputs; estimated M; and comparable fitted outputs. Translation would reconstruct catch numbers and survey abundance indices on the available composition years, change the source observation likelihoods, aggregate the reported 25+ composition to tinyAM's selected age grid, omit length-composition likelihoods, and approximate survey timing. These are explicit translation differences, not exact replications of the official likelihood. No values have yet been imported into the canonical database.

Sources:

- [2023 GOA Pacific ocean perch assessment update](https://apps-afsc.fisheries.noaa.gov/Plan_Team/2023/GOApop.pdf)
- [2024 GOA SAFE introduction and assessment schedule](https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=C2+GOA+Groundfish+SAFE+Introduction.pdf&p=ab2770e7-ed6f-4a2f-a69d-5cbd0bda1896.pdf)
- [2025 GOA catch report](https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=C4+GOA+2025+Catch+Reports.pdf&p=f6b9ff85-9a04-46e3-a7ce-0c2a09b6ccee.pdf)
- [2026 SSC assessment review plan](https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=F+2026+Groundfish+Assessment+SSC+Review.pdf&p=2aa31081-76e1-4025-af01-e17fd99ba0d6.pdf)
- [Pinned GOA POP 2023 model files](https://github.com/noaa-afsc/goa_pop/tree/e87bbfc60ab8bfc87d05e02ed1b4b39b76059d7e/2023/mgmt/2020.1-2023)
- [Pinned 2023 assessment data-generation script](https://github.com/noaa-afsc/goa_pop/blob/e87bbfc60ab8bfc87d05e02ed1b4b39b76059d7e/2023/R/2023_analysis.R)
