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

A subsequent Canadian-petrale search recovered regional POP biomass draws for
3CD and 5DE. The 5ABC file ends in 2017, so it cannot complete the matching
coastwide comparison. The [BC petrale preflight](source_reviews/dfo_petrale_bc.md)
identifies a Canadian replacement candidate with 1938–2024 numerical SSB
posterior output; native inputs and ensemble details remain to be verified.

The six targeted source, readiness, prior and reporting test files pass, as
does database validation. Other stocks' diagnostic and comparison rows are
unchanged. The agreement heatmap and four trend plots now contain the same
27 cases; the 7-column, 4-row trend layout and grey ocean labels are retained.
The full package suite passes with its one interactive-only skip. The full
assessment-analysis suite and `R CMD check` were not repeated for this
source-readiness/reporting change; package code and canonical tables are unchanged.

## BC petrale import — 2026-10-17

`dfo_petrale_bc_2024` is an additional partial summary/output record, not a
completed replacement for BC POP. The import retains 87 annual female SSB
medians with 95% credible limits, separate published terminal/reference values,
three reported component M scalars and 33 assumptions. The database now has
36 stocks and 42 assessment records; the 35 current detailed records are
unchanged. This official applied assessment is not marked current detailed
because original observations, biology and configuration remain missing.

The explicit stock entry point stops before model preparation or optimization.
No fitted object, uncertainty estimate, dashboard or agreement statistic was
invented. See the [source review](source_reviews/dfo_petrale_bc.md) for the
public-source search, exact coverage and conditions for resuming translation.

Both new-case tests and database validation pass. All 65 assessment-analysis
files were run: 59 pass and the six pre-existing failures listed above remain.
The full package suite passes with one interactive-only skip. Package checking
with `--no-manual --no-tests` has zero errors, zero warnings and the existing
`plot_bubbles`/`head` NOTE; tests ran separately. The check process used a valid
locale after an initial inherited-locale failure. No core package code or
model specification changed, and all pre-existing canonical values and
comparison/fit diagnostics remain unchanged.

## BC Canary import and translation — 2026-10-17

`dfo_canary_bc_2022` adds the official coastwide Run24 inputs and available
outputs, bringing the database to 37 stocks, 43 assessment records and 36
current detailed cases. Source files, matching fixed-biology evaluation and
SSB posterior draws are cached with provenance. All prior canonical records
are unchanged. Completeness remains partial: official numerical N/F-at-age,
annual recruitment and full posterior draws were not recovered.

The female 1977–2022, ages 1–60+ translation preserves joint-sex denominators
when converting sampled catch/survey biomass to female numbers. Missing ages
are not filled; aggregate-only indices are omitted. Its background and
[source review](source_reviews/dfo_canary_bc.md) describe the simplified
likelihood, recruitment, selectivity, M prior and initialization.

Default starts failed; the one permitted source-informed retry converged
(maximum gradient 0.000333; successful SD report, positive-definite Hessian).
SSB mean absolute difference is 27.91%, terminal difference -59.09% and
correlation 0.295. Female M is 36.51% lower, with no defined temporal
correlation for the constant source series. These are flagged limitations,
not a successful reproduction of the assessment. Other population comparisons
remain missing or non-equivalent. The dashboard and initialization histories
are cached, and compact result tables retain existing stocks' results.

Both new-case tests and database validation pass. All 67 assessment-analysis
files were run: 61 pass, with the same six baseline failures. The full package
suite passes with one interactive-only skip; package checking with
`--no-manual --no-tests` has zero errors, zero warnings and the existing
`plot_bubbles`/`head` NOTE. Tests ran separately. Dashboard and plot exports
were verified, retaining ocean/species order, seven columns, four rows and
the smaller grey ocean labels. No core model mathematics changed.
