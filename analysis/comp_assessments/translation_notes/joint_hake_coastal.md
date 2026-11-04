# Joint U.S.–Canada coastal Pacific hake: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/joint_hake_coastal.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/joint_hake_coastal.md).

## Expanded-library review

### Current tinyAM candidate: model structure and scale

The current translation is fitted over 1975–2024 at ages 2–15+. The candidate
uses deterministic older-age cohort survival (`N process = "off"`) with
**exponential initial abundance** (`init = "exp"`), rather than a free initial
age vector. Recruitment enters at age 2 as IID log-recruitment deviations
around a constant estimated mean, with the **fixed source log SD of 1.4**.
This does not reproduce the source Beverton–Holt age-0 recruitment,
pre-1975 cohort initialization, or terminal-recruitment rules.

Fishing mortality has **correlated age-specific temporal random walks**
(`F process = "cor_rw"`), allowing correlated annual increments across ages.
This is a flexible approximation to the source's annual fishing intensity and
nonparametric time-varying fishery selectivity, not the same parameterization.

The M specification preserves one stock-wide Gaussian log-M mean effect with
center -1.60944 and fixed SD 0.1, with M residual process off and no free
M intercept. Source M bounds (0.05–0.4) and Bayesian posterior integration
are not imposed. The posterior median M=0.233454 remains an **official output**,
not a fixed input used to force agreement.

#### Survey catchability and scale

The **key improvement in the candidate fit** was switching survey catchability
to `q_link = "logit"`, constraining effective q at each age block to the open
interval (0, 1). `q_form = ~ q_block` estimates separate effective catchability
for ages 2–5 and one shared level for ages 6–15. The user reports that this
change resolved much of the earlier abundance-scale discrepancy and that the
resulting model converges and approximates the accepted assessment reasonably
well. No new objective, gradient, Hessian, or numerical agreement statistics
were supplied or recalculated for this review.

The catchability estimates are **not especially informative**: at least one
is near one, and the reported uncertainty intervals extend approximately from
zero to one. The constraint helps stabilize population scale, but this wide
uncertainty signals limited separation of survey q from abundance and other
fitted processes. It should not be described as precise recovery of the
accepted model's q. These tinyAM intervals are not the official model's
Bayesian credible intervals.

Importantly, the operational Stock Synthesis assessment analytically profiles
**one scalar survey q** (reported posterior median 0.894) and separately
estimates a normalized, nonparametric age-selectivity curve. TinyAM's
age-block coefficients instead represent **effective catchability at age**:
there is no separately identifiable scalar q and normalized selectivity.
The official median 0.894 therefore should not be compared uncritically with
each tinyAM age-block estimate. The operational model's separate selectivity
and q structure is a possible future refinement, but is not needed to
demonstrate approximate agreement in historical population trajectories.

#### Observation uncertainty and remaining limitations

The candidate retains the supplied year-specific `log_index_sd` inputs in
translation (they are no longer discarded). `database_to_tam_obs()` joins
these as `obs$index$relative_sd`; the fit uses
`index_settings = list(q_form = ~ q_block, q_link = "logit", sd_form = ~ 1,
 sd_supplied = ~ relative_sd)`. Under tinyAM's current rules the ordinary
intercept in `sd_form = ~ 1` is removed when `sd_supplied` is given, so this
fits **no additional global index SD multiplier**. In contrast, the accepted
assessment estimates extra acoustic-index SD combined with sampling SD in
quadrature. The reconstructed index-at-age likelihood also treats age cells
from the same aggregate acoustic biomass observation as independent, without
the original age-composition covariance or ageing-error likelihood.

Catch-at-age is reconstructed from native catch biomass and percentage age
compositions; its lognormal SD is estimated (`catch_settings$sd_form = ~1`).
The original separate catch-biomass and Dirichlet–multinomial age-composition
likelihoods are not reproduced. These differences can influence q, abundance
scale and the width of intervals even when the population trajectories agree.

The working candidate ends with a **single direct `fit_tam()` call** inside
`if (do_fit)`; the earlier source-state initialization retry has been removed.
The unchanged comparison limitations remain: marginal official N medians do
not sum to a joint median, modeled ages start at 2, and native Bayesian F-at-age
posterior draws are unavailable. The near-boundary q estimate is an acknowledged
limitation of this pragmatic translation, not a claim of full operational-model
reproduction.

#### Historical recruitment sensitivity

An earlier sensitivity replaced the constant recruitment mean with a
two-year-lagged Beverton–Holt mean and IID age-2 recruitment deviations. It
converged after one source-based initialization retry (objective 1183.769,
maximum gradient 0.000111, positive-definite Hessian and successful uncertainty
estimation; 25 fixed and 749 random effects). Against the preceding constant-mean
candidate, SSB median absolute percent difference changed from 44.78% to 44.79%,
mean difference from 42.23% to 42.25%, and trend correlation from 0.751 to 0.750.
Recruitment differences remained near 55.2% median absolute percent difference,
with terminal difference changing from +346.1% to +347.2%. This sensitivity is
retained as historical context; it is not the current constant-mean recruitment
specification in the translation script.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/joint_hake_coastal.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The translation uses 1975–2024 and ages 2–15+. Native percentages are
normalized outside the database, as in the source, and converted using
original totals and weights before age restriction. All supplied weights
and fecundity products are identical across ages 15–20+, so the observed
15+ grouping does not require an invented biological tail. The supplied
fecundity product divided by population weight gives maturity; first mature
age 2, female fraction 0.5 and January 1 timing are applied separately.

The local dashboard retains native posterior SSB and explicitly approximate
common-age sums of marginal medians. F-at-age and Fbar comparisons remain
unavailable. Available credible intervals are retained only on their native
quantities; no uncertainty is attached to sums or pooled posterior medians.
