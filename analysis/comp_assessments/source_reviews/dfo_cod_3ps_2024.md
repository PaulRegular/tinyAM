# DFO 3Ps Atlantic cod, 2024 assessment

## Assessment identity and source selection

The metadata-only Charbonneau–Keith record `DFO_morhua_3PS` is DFO Atlantic cod in NAFO Subdivision 3Ps. It is represented as `ck_dfo_morhua_3ps`; the assessment import updates this stock identity to “Atlantic cod” and authority code `3Ps`.

The production assessment used the HYBRID model under the framework adopted in 2019. The current detailed report is the DFO Research Document 2026/051, which reports the assessment conducted in 2024. DFO’s 2025 stock summary identifies 2024 as the last full assessment. The detailed report’s Appendix A supplies the model input tables. The later DFO Research Document 2026/047 updates catch and survey information through 2025; it is retained in the local cache as a separate follow-up and none of its values are imported into this 2024 assessment.

The local ignored cache is `source_cache/ck_dfo_morhua_3ps/assessment_2024/`. It contains the detailed 2026/051 report (SHA-256 `759b26fe500dca205b669cf9475707d0b67604ca3c0cf00fb98aae0328d23728`), the 2026/047 later update (SHA-256 `78fedd1cbc0593013a5066d18878594f73204897084ecc9abe020d42986bb965`), and earlier assessment/data-review documents. `accessible_sources_manifest.csv` records retrieval dates, source URLs, roles and checksums. The 2024/015 document is a comparison-model development report for SAM, not the applied HYBRID run.

## Official model and inputs

The applied model spans calendar years 1959–2024 and ages 2–14+, with age 2 as recruitment. It fits annual landings with a censored likelihood and the age composition using continuation-ratio logits. Fishing mortality follows an age-correlated multivariate random walk for ages 3+, with age 2 separate and a break at the 1993 moratorium. Natural mortality varies with a scaled fish-condition index; the mean of 1978–1982 supplies the index values for 1959–1977.

The import script `scripts/database/073_import_dfo_cod_3ps_2024.R` transcribes Appendix A tables A-1 through A-11 into `inputs.csv`:

- Table A-1: calendar-year landings and the lower and upper bounds used by the censored likelihood, in tonnes, 1959–2024. The report says 2024 landings are set to the 2024/25 TAC.
- Table A-2: commercial catch numbers at age, 1959–2023, ages 2–14+. The current table omits the scale label. DFO Research Documents 2023/025 and 2024/002 label the same catch-at-age series in thousands of fish (000s); this confirms the unit convention, but their revised historical values are not copied into the 2024/051 input series.
- Table A-3: the signed condition index for M, 1959–2024, ages 2–14+ as one shared annual covariate.
- Tables A-4 through A-8: age-specific indices from the DFO RV spring survey (1983–2024), French ERHAPS (1978–1991), GEAC (1998–2005), Sentinel Gillnet (1995–2023), and Sentinel Line-trawl (1995–2023). Values and zeros are retained at the printed scales. Because the tables do not identify units or multipliers, these are recorded as survey-index scale with units unresolved. The DFO RV table omits 2006 (incomplete survey), 2020 and 2023 (no survey). The report states that zero survey catches are censored in HYBRID.
- Table A-9: maturity proportions at age, 1959–2024, ages 2–14+.
- Table A-10: model stock weights, kg per fish, 1959–2024, ages 2–14+.
- Table A-11: catch weights, kg per fish, 1959–2024; the source prints the last group as age 14 without a plus sign. This is retained as printed and its mapping to the model’s 14+ group is unresolved.

The tables do not give the exact survey timing fractions or all likelihood weights, random-effect covariance settings, priors, and other run controls. Those source gaps keep assumptions and input status `partial`. Do not supply artificial sampling times or treat the summarized index values as exact count units. Any tinyAM timing values must be labelled as translation approximations.

## Outputs and comparison limits

The detailed report presents historical trajectories mainly as figures and does not provide a complete numerical table of fitted N-at-age, F-at-age, M-at-age, q, or annual SSB and recruitment. No public applied HYBRID executable, full parameter file or fitted result object was found in the examined DFO materials. The database therefore records only the numerical summaries stated in the report, without digitizing curves:

- Beginning-of-year 2024 SSB: 29.9 kt (95% CI 24–37.3 kt).
- Beginning-of-year 2025 projected SSB under a 1,550 t catch: 33.2 kt (95% CI 26.4–41.8 kt).
- 2024 age-2 recruitment: 36.6 million fish (95% CI 14.3–93.3 million).
- 2024 mean F for ages 5–8: 0.04 per year (95% CI 0.03–0.06).
- 2024 mean M for ages 5–8: 0.32 per year (95% CI 0.31–0.34).

Outputs remain `partial`; the figures are useful for visual context but are not treated as exact time-series values. Total abundance cannot be derived because numerical N-at-age is unavailable. A full Official–tinyAM trajectory comparison is therefore blocked unless additional HYBRID output data are recovered.

## tinyAM translation readiness

The source inputs are substantially available, but the production likelihoods and process structures are not reproduced directly by tinyAM. In particular, censored zero-index observations, censored landings intervals, continuation-ratio catch composition, age-correlated F random walks with a moratorium break, and condition-based M require explicit approximation. The source-scale index values and missing exact sampling-time fractions also need review before a fit could be interpreted. The current R environment used for this review could not load `tinyAM`, `RTMB`, or the fitting dependencies, so no fit or convergence result is claimed. The stock is recorded as a data-capture checkpoint, not as a completed runnable translation.

