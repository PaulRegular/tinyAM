# Northern Shelf (North Sea) haddock: 2026 source review: tinyAM translation notes

The registered stock script defines the current model. Earlier translation notes were transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_haddock_north_sea.md) and are retained below.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_haddock_north_sea.md).

## Current curated translation

The recipe matches the native SAM likelihood over 1972–2026 at ages 0–8+.
It integrates first-year recruitment and older first-year abundance without
an initial density. Recruitment and age-specific F follow random walks, with
correlated F increments across ages. Older-cohort survival SDs are shared
over ages 1–7, with a separate SD at 8+.

The model uses the native catch and survey observation-SD groups and retains
the survey precision multipliers. Each survey and age has a baseline
catchability. The three q powers apply to Q1 age 1 and Q3+Q4 ages 0 and 1.
Observation errors are independent on the log scale.

The fitted catch and survey observations, fixed M, stock weights, maturity
and beginning-of-year spawning match the native inputs. Separate catch weights
support yield reporting. The 2025 catch weights fill the tinyAM reporting grid
in the survey-only final year. The Fbar comparison uses an abundance-weighted
mean over ages 2–4 for both fits. The published SAM Fbar is an arithmetic mean.
The [native SAM checks](../../sam_equivalence/README.md#northern-shelf-haddock-2026)
document equivalence and uncertainty conventions.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_haddock_north_sea.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The historical translation fits ages 0–8+ over 1972–2026 with random-walk recruitment, IID cohort-survival deviations and grouped survival-process SDs, correlated F random walks, fixed M and both standardized surveys. It preserves the zero pre-spawning mortality fractions. Remaining differences include precise initial-state and covariance conventions, observation variances and the omission of the three estimated density-dependent q powers. The q powers are a potential source of observation-model disagreement. No likelihood or uncertainty equivalence is claimed.
