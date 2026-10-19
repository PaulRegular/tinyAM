# DFO 3Ps Atlantic cod, 2024 assessment

## Assessment identity and source selection

The metadata-only Charbonneau–Keith record `DFO_morhua_3PS` is DFO Atlantic cod in NAFO Subdivision 3Ps. It is represented as `ck_dfo_morhua_3ps`; the assessment import updates this stock identity to “Atlantic cod” and authority code `3Ps`.

The production assessment used the HYBRID model under the framework adopted in 2019. The current detailed report is DFO Research Document 2026/051, which documents the 2024 assessment. Its Appendix A supplies the model input tables. The associated DFO Science Advisory Report 2025/007 summarizes the same 2024 assessment and provides annual aggregate estimates for 2020–2024 in Table 1. DFO’s assessment-history field identifies 2023 as the previous “last full assessment”; the 2024 record is the latest annual assessment documented here. The later DFO Research Document 2026/047 updates catch and survey information through 2025; it is retained in the local cache as a separate follow-up and none of its values are imported into this 2024 assessment.

The local ignored cache is `source_cache/ck_dfo_morhua_3ps/assessment_2024/`. It contains the detailed 2026/051 report (SHA-256 `759b26fe500dca205b669cf9475707d0b67604ca3c0cf00fb98aae0328d23728`), the 2025/007 summary, the 2026/047 later update (SHA-256 `78fedd1cbc0593013a5066d18878594f73204897084ecc9abe020d42986bb965`), and earlier assessment/data-review documents. `accessible_sources_manifest.csv` records retrieval dates, source URLs, roles and checksums. The 2024/015 document is a comparison-model development report for SAM, not the applied HYBRID run.

## Official model and inputs

The applied model spans calendar years 1959–2024 and ages 2–14+, with age 2 as recruitment. It fits annual landings with a censored likelihood and the age composition using continuation-ratio logits. Fishing mortality follows an age-correlated multivariate random walk for ages 3+, with age 2 separate and a break at the 1993 moratorium. Natural mortality varies with a scaled fish-condition index; the mean of 1978–1982 supplies the index values for 1959–1977.

The import script `scripts/database/073_import_dfo_cod_3ps_2024.R` transcribes Appendix A tables A-1 through A-11 into `inputs.csv`:

- Table A-1: calendar-year landings and the lower and upper bounds used by the censored likelihood, in tonnes, 1959–2024. The report says 2024 landings are set to the 2024/25 TAC.
- Table A-2: commercial catch numbers at age, 1959–2023, ages 2–14+. The current table omits the scale label. DFO Research Documents 2023/025 and 2024/002 label the same catch-at-age series in thousands of fish (000s); this confirms the unit convention, but their revised historical values are not copied into the 2024/051 input series.
- Table A-3: the signed annual fish-condition index used to inform M, 1959–2024; 1959–1977 use the 1978–1982 average.
- Tables A-4 through A-8: age-specific indices from the DFO RV spring survey (1983–2024), French ERHAPS (1978–1991), GEAC (1998–2005), Sentinel Gillnet (1995–2023), and Sentinel Line-trawl (1995–2023). Values and zeros are retained at the printed scales. Because the tables do not identify units or multipliers, these are recorded as survey-index scale with units unresolved. The DFO RV table omits 2006 (incomplete survey), 2020 and 2023 (no survey). The report states that zero survey catches are censored in HYBRID.
- Table A-9: maturity proportions at age, 1959–2024, ages 2–14+.
- Table A-10: model stock weights, kg per fish, 1959–2024, ages 2–14+.
- Table A-11: catch weights, kg per fish, 1959–2024; the source prints the last group as age 14 without a plus sign. This is retained as printed and its mapping to the model’s 14+ group is unresolved.

The tables do not give the exact survey timing fractions or all likelihood weights, random-effect covariance settings, priors, and other run controls. Those source gaps keep assumptions and input status `partial`. Do not supply artificial sampling times or treat the summarized index values as exact count units. Any tinyAM timing values must be labelled as translation approximations.

## Outputs and comparison limits

The detailed report presents historical trajectories mainly as figures and does not provide a complete numerical table of fitted N-at-age, F-at-age, M-at-age, q, or annual SSB and recruitment. Table 1 of the associated 2025/007 Science Advisory Report provides annual aggregate SSB, age-2 recruitment, and mean F/M for ages 5–8 for 2020–2024. The 2024 SSB/recruitment intervals and grouped F/M intervals are also stated in Research Document 2026/051; earlier grouped F/M confidence limits are not tabulated. A separate 2025 SSB projection is reported under a specified catch. No public applied HYBRID executable, full parameter file or fitted result object was found in the examined DFO materials. The database records the five-year aggregate series and 2025 projection, without digitizing model curves:

- SSB for 2020–2024: 27.6, 26.7, 27.4, 27.7, and 29.9 kt; 2024 95% CI 24–37.3 kt.
- Beginning-of-year 2025 projected SSB under a 1,550 t catch: 33.2 kt (95% CI 26.4–41.8 kt).
- Age-2 recruitment for 2020–2024: 17.4, 20.0, 34.3, 36.8, and 36.6 million fish; 2024 95% CI 14.3–93.3 million.
- Mean F for ages 5–8 in 2020–2024: 0.04, 0.02, 0.03, 0.03, and 0.04 per year; 2024 95% CI 0.03–0.06.
- Mean M for ages 5–8 in 2020–2024: 0.37, 0.35, 0.35, 0.35, and 0.32 per year; 2024 95% CI 0.31–0.34.

Outputs remain `partial`; the five-year aggregate table improves recent comparison, but the figures are not digitized as exact time series. Total abundance cannot be derived because numerical N-at-age is unavailable. Age-specific N/F/M and annual outputs before 2020 cannot be compared until additional HYBRID output data are recovered.

## tinyAM translation readiness

The translation recipe is in `scripts/translation/stocks/dfo_cod_3ps_2024.R`. It uses reported catch-at-age directly, fits the five index series with independent age-by-survey q values, approximates seasonal timing by within-season midpoints, uses an age-correlated F random walk without the source age-2 separation or 1993 break, and models log M with age effects and the annual condition covariate. These are explicit approximations: tinyAM does not reproduce censored zero-index observations, censored landings intervals, continuation-ratio catch composition, the exact condition-to-M relationship, or all HYBRID process structures. Survey index units remain unresolved, with survey-specific q absorbing the reported scale. A recipe-wiring test checks the year and age grid, timing arguments, covariate mapping and comparison window using a stubbed observation converter; it does not verify conversion or fitting. The local R environment cannot load `tinyAM`, `RTMB`, or the fitting dependencies, so the actual conversion and fit remain blocked before optimization. No convergence result is claimed.

