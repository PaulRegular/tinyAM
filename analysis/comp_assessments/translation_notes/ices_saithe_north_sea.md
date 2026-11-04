# ices saithe north sea: tinyAM translation notes

Inherited settings and diagnostics; no models were rerun. Check the registered recipe and review tier for current settings.

## Current initial abundance

The retained translation estimates older first-year abundance freely. It does
not use source N/F estimates to construct starting values. Current age pooling,
catchability and other scientific settings are defined in the stock script.

## Remaining limitations for translation

- The final 2026 native SAM data/model object was not located. The report
  provides rounded inputs and outputs, but not age-specific uncertainty or
  fitted observation predictions.
- Exact 2026 within-year sampling fractions are absent from the report. The
  cached 2024 native data object has sample times 0.730 for the Q3–Q4 index
  and 0.525 for commercial CPUE; these are useful context, not verified 2026
  values.
- tinyAM's current observation structure does not represent a single
  aggregate relative CPUE observation tuned to exploitable biomass. The
  translation therefore fits the Q3–Q4 age-specific index and records the
  CPUE series as omitted from the fit. It does not divide aggregate CPUE among
  ages or treat its relative scale as absolute biomass.
- The chosen F process, shared source F state and correlated observation errors
  differ from SAM. tinyAM supports age-correlated F RW increments and process-SD
  sharing, but their availability does not make the retained density identical.

The database status remains partial for assumptions, inputs, and outputs
because the current native run object and exact sampling fractions are
unavailable, and the published age-specific surfaces lack reported
uncertainty. No plot-derived values or fitted estimates are being substituted
for inputs.

[Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_saithe_north_sea.md).

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_saithe_north_sea.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The printed F variance keys are `0 1 2 3 3 3 3 3`: ages 3, 4 and 5 have
separate innovation SDs, with ages 6–10+ sharing one. Separate recruitment and shared older-cohort SDs are
already supported; that simpler N variance grouping is not itself a missing
capability. The initial-state treatment and F density still differ.
