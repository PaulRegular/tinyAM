# British Columbia petrale sole — translation notes

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/dfo_petrale_bc.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

BC petrale remains a promising Canadian comparison candidate, with numerical
SSB now preserved in the canonical database. Fitting is **source-blocked**.
These source gaps must be resolved before preparing a model. There is no
catch-only substitute, invented `fit_tam()` model, optimizer result or dashboard.

The stock entry point is now the [curated stock script](../scripts/curated/dfo_petrale_bc_2024.R). It records the source
background and stops with the missing-input inventory before `prepare_tam()`
or `fit_tam()`, for both `fit = FALSE` and `fit = TRUE`. An explicit runner
attempt is recorded as `translation_failed`, with no optimizer or uncertainty
values. It is not an actual optimization attempt or a completed translation.
The non-current summary record is excluded from the default detailed-case
batch; it can still be inspected or explicitly requested by assessment ID.

The explicit readiness attempt used committed database revision `2a66fcf` and
stopped before optimization. Both new-case test files pass. All 65
assessment-analysis files were run: 59 pass, with the same six pre-existing
failures in age-subset, Southern Gulf herring, F-process, herring, Norway pout
and general translation fixtures. The stock-script inventory test now explicitly
checks the documented BC POP and BC petrale blockers; its runnable-case
assertions are unchanged. The full package suite passes with one
interactive-only skip. `R CMD check --no-manual --no-tests` has zero errors,
zero warnings and the existing `plot_bubbles`/`head` NOTE; package tests ran
separately. The initial check's invalid inherited locale was corrected for the
check process without changing package code. Logs and the readiness result
remain local. No agreement or convergence statistics exist for this case.
