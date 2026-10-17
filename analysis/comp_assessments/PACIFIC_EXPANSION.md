# Pacific expansion

The expansion adds twelve stock cases, not separate sexes, model variants or
regional components. The eight primary catalogue identities below come from
Charbonneau–Keith's PLOS analysis set. Pacific ocean perch is in the broader
catalogue. Coastal hake and the two Canadian cases are explicit exceptions.

| Stock | Catalogue identity | Detailed production source being reviewed |
|---|---|---|
| Eastern Bering Sea Pacific cod | AFSC_BSAI_Gadus_macrocephalus | 2024 Model 24.1, native SS files and Appendix 2.4 workbook |
| Aleutian Islands Pacific cod | AFSC_AI_Gadus_macrocephalus | 2024 Model 24.1, final SAFE report and public native files |
| Aleutian Islands pollock | AFSC_AI_Gadus_chalcogrammus | 2024 Model 15.1, SAFE tables and AMAK framework |
| BSAI Atka mackerel | AFSC_BSAI_Pleurogrammus_monopterygius | 2024 Model 16.0b, SAFE tables and AMAK framework |
| Eastern Bering Sea yellowfin sole | AFSC_EEBSAI_Limanda_aspera | 2024 Model 23.0, native ADMB files |
| Eastern Bering Sea northern rock sole | AFSC_EEBSAI_Lepidopsetta_polyxystra | 2024 Model 24.2, native ADMB inputs and SAFE outputs |
| BSAI Alaska plaice | AFSC_BSAI_Pleuronectes_quadrituberculatus | 2024 Model 24.1b, native SS input and output files |
| BSAI arrowtooth flounder | AFSC_BSAI_Atheresthes_stomais | 2022 Model 18.9; 2024 document is a projection update |
| BSAI Pacific ocean perch | AFSC_alutus_BSAI | 2024 Model 24, native ADMB input and output files |
| Joint US–Canada coastal hake | Catalogue exception | 2026 catch-only update, 2025 detailed assessment and native SS files |
| Strait of Georgia herring | Canadian catalogue exception | 2025 update, SISCAH ensemble and interim estimation model |
| British Columbia Pacific ocean perch | Canadian catalogue exception | 2023 three-area assessment, Research Document 2025/004 and native SS inputs |

Historical catalogue labels do not establish current eligibility. In particular,
Aleutian Islands cod now has an age-structured Tier 3 assessment. The historical
arrowtooth catalogue identifier retains its original spelling; the species is
*Atheresthes stomias*.

The final December 2024 SSC report identifies the selected Alaska models.
NOAA's 2026 harvest-specifications review explains why the 2024 fitted models
remain the relevant detailed production runs after the 2025 funding lapse.
Those documents are cached in `source_cache/pacific_expansion/`, with manifests.
Current-year development directories are not evidence of an official new fit.

Each stock receives a source review and a database checkpoint before its
translation checkpoint. Cache manifests record retrieval dates, checksums and
pinned repository revisions. Numerical source blockers and non-equivalent
comparisons remain explicit; they do not count as completed translations.

## Delivered records and fitting outcomes

All twelve source records, importers, source reviews and source-specific tests
were added. Source files and provenance manifests are cached locally. The
database addition and translation/readiness checkpoint for each case were
committed separately and pushed to `comp-assessments`.

| Case | Fitting outcome | Main qualification |
|---|---|---|
| EBS Pacific cod | Passed from defaults | Empirical catch-length allocation; proxy maturity and plug-in weights |
| AI Pacific cod | Blocked before fitting | Public native run differs from final production run; final inputs incomplete |
| AI pollock | Neither attempt passed | Failed covariance checks; extreme abundance/q scale trade-off |
| BSAI Atka mackerel | Passed after one source-start retry | Source selectivity and observation precision constraints differ |
| EBS yellowfin sole | Neither attempt passed | Large gradients and failed covariance; female subset |
| EBS northern rock sole | Neither attempt passed | Failed curvature and near-zero M; female subset |
| BSAI Alaska plaice | Passed after one source-start retry | Sex pooling and source selectivity constraints differ |
| BSAI arrowtooth | Blocked before fitting | Likelihood compositions, age-length key and survey timing unverified |
| BSAI Pacific ocean perch | Interrupted after more than seven CPU hours | No completed attempt history, usable fit or dashboard |
| Coastal Pacific hake | Passed after one source-start retry | Source M penalty restored in review; M 0.18335, finite effect SE; core comparisons improve but remain approximate |
| SoG herring | Blocked before fitting | Aggregate spawning-biomass survey has no observed survey-age composition |
| BC Pacific ocean perch | Source-output blocked after historical trial | Numerical checks passed, but matching historical N/F/recruitment/SSB outputs have not been recovered |

Nine cases historically received actual optimization attempts, including the
now-blocked BC POP trial. The current status is four numerical passes, three
non-converged fits, one interrupted run and four source blockers. Seven active
new-case dashboards and the historical BC trial are cached under
`results/cache/<assessment_id>/dashboard.html`. Readiness guards prevent fitting
blocked cases. The twelve-runnable-translation target remains unmet. Numerical
passes are not scientific agreement or evidence of useful uncertainty.

The complete cross-case flags, unavailable/non-equivalent comparisons,
fixed-M checks and numerical outcomes are in [PACIFIC_REVIEW.md](PACIFIC_REVIEW.md).
The review covers the updated library, distinguishing new from existing cases.
It can be refreshed without fitting by sourcing
`scripts/translation/review_pacific_expansion.R` from the repository root.
The expansion initially retained existing model specifications. The subsequent [model review](PACIFIC_REVIEW.md#source-and-model-review) documents the selected mortality, initialization and spawning-time changes and rejected alternatives. No core model mathematics changed.

## Final validation

- Database validation passed: 35 stocks, 41 assessment records, 1,105 assumptions,
  175,397 input rows and 105,332 output rows.
- All original records are unchanged relative to the terminology-cleanup
  checkpoint `311922f`: 23 stocks, 29 assessments, 877 assumptions, 93,105 inputs
  and 44,420 outputs. New records are additions, not replacements.
- New source, translation, default-first/retry, source-blocker and review-flag
  tests passed. Preparation with fitting disabled does not optimize.
- Assessment-analysis tests: 55 of 61 files passed after targeted reruns.
  Six failures also occur at the starting revision: `test_age_subsets.R`,
  `test_dfo_herring_sgsl_translation.R`, `test_f_process_translations.R`,
  `test_herring_translation.R`, `test_ices_norway_pout_translation.R` and
  `test_translation.R`. These concern existing fixtures/settings or a wrapped
  message; they were not hidden by weakening assertions or changing old models.
  The parallel-runner test passed with local worker connections permitted.
- Full package test suite passed, with one interactive-only dashboard test skipped.
  The new-case dashboards were rendered separately.
- `R CMD check --no-manual` completed with zero errors, zero warnings and one
  existing NOTE about an unqualified `head` call in `plot_bubbles`.

Local validation logs are retained outside the tracked repository. The final
cross-case reporting changes were checked separately after the full tests;
no model fitting was repeated to improve agreement.

## Model-review validation — 2026-10-17

The [case-by-case review](PACIFIC_REVIEW.md#source-and-model-review) records the subsequent
source-informed model trials and selected changes. The original expansion
outcomes above remain historical where superseded by that review.

- Canonical numerical inputs and official outputs are unchanged from `86e44bf`.
  The northern rock sole assumption note now distinguishes the verified female
  M penalty from unresolved male-prior activation. Database validation passes.
- All nine targeted test files passed, including shared M density/scale,
  unchanged observations, stock recipes and comparison definitions. The EBS
  alternative's projection mapping also passed a final targeted rerun.
- The full package suite passed with its one interactive-only skip.
- All 62 assessment-analysis files were run: 56 pass after a separate dashboard
  test rerun; the same six pre-existing fixture/settings failures listed above
  remain. The dashboard test initially collided with a concurrent render and
  passes when rendered alone. No assertions were weakened.
- The last complete `R CMD check --no-manual` at the expansion baseline had
  zero errors, zero warnings and the existing `plot_bubbles` NOTE. It was not
  repeated for this analysis-only review; core package sources, tests and
  generated documentation are unchanged.
- Selected dashboards and the agreement heatmap were refreshed. The revised
  hake fit passes after one conditional N-start retry; the historical BC POP
  trial passed from defaults but is now retired from comparison. Rejected,
  failed and time-limited trials remain cached separately.

## BC POP output revisit — 2026-10-17

The matching production repository, report tables and related public PBS
repositories did not yield numerical historical N/F/recruitment/SSB outputs.
BC POP is now blocked from fitting and excluded from the comparison summaries
and plots; canonical inputs and official summaries are unchanged. See the
[source review](source_reviews/dfo_pop_bc.md) for cached search evidence and the
proposed U.S. West Coast petrale sole alternative. No replacement was fitted.

The six targeted source, readiness, prior and reporting test files pass, as
does database validation. Other stocks' diagnostic and comparison rows are
unchanged. The agreement heatmap and four trend plots now contain the same
27 cases; the 7-column, 4-row trend layout and grey ocean labels are retained.
The full package suite passes with its one interactive-only skip. The full
assessment-analysis suite and `R CMD check` were not repeated for this
source-readiness/reporting change; package code and canonical tables are unchanged.
