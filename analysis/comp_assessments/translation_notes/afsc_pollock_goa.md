# Gulf of Alaska pollock: official-assessment source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_pollock_goa.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_pollock_goa.md).

## Expanded-library review

### Current tinyAM translation

The current candidate fits 1970–2024 at ages 1–10+, matching the official
stock extent and retaining the fixed age-specific M vector. It uses
deterministic survival of older cohorts (`N_settings$process = "off"`) with
exponential initial abundance. Recruitment is IID around an estimated constant
log mean, with its SD fixed at 1.3. This matches the official annual
recruitment-deviation variance and lack of temporal correlation more closely
than the previous tinyAM recruitment RW. The first recruitment state in
tinyAM remains an unpenalized fixed anchor, and exact official initial-age
construction is not reproduced.

Fishing mortality now uses `F_settings$process = "ar1"` with a rising logistic
mean-at-age (`mu_form = ~ logistic(age)`) and reporting ages 3–10. This
regularizes an age/year log-F surface around a static smooth mean rather than
using independent random walks at each age. It is an approximation to the
official model's annual F intensity multiplied by double-logistic selectivity,
whose ascending-limb parameters evolve under penalties; tinyAM does not
reproduce that parameter evolution or a descending fishery-selectivity limb.

The four survey age patterns are represented parsimoniously by *two* shared
logistic effects and four survey-specific intercepts. A reversed age coordinate
(`sel_age = 11 - age`) gives descending age responses for the Shelikof winter
and summer acoustic surveys; unmodified age gives rising responses for the
NMFS bottom trawl and ADF&G crab/groundfish trawl surveys. The formula
`~ 0 + survey + logistic(sel_age, by = survey_type)` therefore fits one
logistic shape for the acoustic surveys and another for the trawl surveys,
not a separate age curve for every survey. Both enter the **logit-q**
linear predictor; they are not literal multiplicative selectivity curves on
the original q scale. The source uses distinct age-selectivity curves: a
descending logistic for Shelikof, double logistics for NMFS bottom trawl and
summer acoustic, and a rising logistic for ADF&G. This shared two-curve
approximation intentionally does not reproduce all those distinctions.

Annual catchability now includes `rw(year, by = adfg_indicator)`, applied
*only* to the ADF&G survey. This improves on the former independent annual
ADF&G fixed effects by penalizing changes between years. The official model
instead penalizes log-q changes using source-supplied, year-specific increment
SDs; tinyAM's formula uses one estimated increment SD and operates on
logit-q rather than log-q. Shelikof uses the observed environmental covariate
as a logit-q predictor; the source fits a latent environmental AR1 process
with observation error, which is not reconstructed here. Bottom-trawl and
summer acoustic catchability have no annual random walk in the candidate,
consistent with the absence of active q deviations in the official fit.
The bottom-trawl prior is also omitted and all tinyAM catchabilities remain
bounded by the logit link, unlike the official positive log-q formulation.

The candidate fits the published fishery catch numbers at ages 1–15
(1975–2023), summing ages 10–15 into the modeled age-10 plus group. Because
the published table is rounded to 0.01 million fish, printed zeros are treated
as missing by tinyAM's lognormal catch-at-age observation likelihood. This
uses one estimated catch log-SD rather than the source's separate total catch
biomass and Dirichlet-multinomial composition likelihoods, including source
age-bin accumulation and age-reading error. Survey numbers-at-age are
reconstructed from native biomass, proportions-at-age and appropriate survey
weights and timing; aggregate-only survey years and the reconstructed Shelikof
pooled age-1–3 observation (represented as an age-3 row) are omitted. Survey
observation SDs use supplied relative log-SDs multiplied by an estimated
survey-specific factor. The source instead uses distinct biomass/composition
likelihoods and composition weights.

The pinned official C++ implementation expresses population N and catch C
in billions of fish, with weights in kg per fish. Its total catch expression
`1e6 * sum(C * wt_fsh)` gives tonnes; survey biomass is expressed in million
tonnes. For survey reconstruction, source biomass is converted to kg and
allocated across ages by the native composition and survey weights. The
published catch numbers are converted from millions to fish, and official
N/recruitment and SSB/biomass outputs retain their native million-fish and
thousand-tonne reporting scales. The existing staged native data and unit
checks are unchanged by this translation update.

Female maturity is multiplied by 0.5 in tinyAM. **SSB definitions still
differ:** the official calculation is
`sum(N * exp(-0.21 * Z) * wt_spawn * 0.5 * mat)` using spawning weights,
whereas the candidate's tinyAM SSB is
`sum(N * wt_pop * 0.5 * mat)` at the beginning of the year because no
`spawn_time` was supplied. The native preparation sets
`wt_spawn = wt_srv1` and `wt_pop = wt_srv2`. Numerical common-definition
summaries continue to use translated stock weights and maturity, while the
dashboard retains the official native SSB. Yield reporting continues to use
source total catches and predictions multiplied by source catch weights;
these are conditional median predictions without a lognormal mean correction.

The candidate remains a structural approximation: logistic age curves are
grouped by survey type, fitted q effects use a logit rather than log link, and
the source selectivity priors, environmental measurement process, ADF&G
increment-specific penalties, age-reading error, and composition likelihoods
are not fitted directly. The current package-default fit converged with
optimizer code 0 (`relative convergence (4)`), objective 1,358.922, maximum
absolute gradient 0.000196, successful uncertainty estimation, and a
positive-definite Hessian. Residual checks still flag group location, spread,
and outliers.

Over 1970-2024, the current comparison summary reports mean absolute
percentage differences of 29.2% for total abundance, 46.6% for N-at-age,
50.4% for recruitment, 29.7% for SSB under the shared-biology definition, and
26.1% for biomass. Their terminal differences are +83.9%, +104.1%, +79.2%,
+11.1%, and +48.8%; trend correlations are 0.937, 0.950, 0.949, 0.924, and
0.912. This SSB comparison uses official N with the translated stock weights
and maturity at the start of year; it is not the official native SSB
definition. A separate calculation using the candidate's fitted N/F/M with
the source spawning weights and 0.21 spawning fraction gives 38.7% mean
absolute percentage difference, +7.4% terminal difference, and 0.870 trend
correlation against native official SSB. This derived series has no propagated
uncertainty and is not currently part of the exported comparison summary.

The database retains 371 approximate F-at-age cells derived from rounded
N-at-age and fixed M, but they do not form a complete age-year surface; current
F-at-age and Fbar comparisons are unavailable. Fixed M matches by construction.
An otherwise identical fit using `N_settings$init = "free"` and default
initialization failed numerical checks (optimizer code 1, false convergence,
maximum absolute gradient 93,800, unavailable finite fixed-effect covariance,
and non-positive-definite Hessian). The converged exponential initializer is
therefore retained. No source-based starting values or other model changes
were used.

### Catchability-link sensitivity

Because the official model uses positive log-q, a default-start fit was tested
with the same translated observations and formulas but `q_link = "log"` in
place of the current `"logit"`. The optimizer stopped at code 0 with a small
gradient (0.000055), but fixed-effect covariance was unavailable and the
Hessian was not positive definite, so the fit failed tinyAM's numerical
checks. Aggregate trend correlations were mostly high (0.90–0.96), but median
absolute percent differences increased to 59.6% for SSB, 51.5% for biomass,
44.1% for total abundance, 51.3% for N-at-age and 41.2% for recruitment.
This candidate is not retained. The source's positive q scale is a real
structural difference, but this test cannot establish that the bounded q link
is preferable: the reconstructed age-specific observations and independent
lognormal likelihood may leave q and absolute abundance weakly identified.

## Paper-working-set review before migration

Only distinct notes from the earlier main review are retained below. Their settings and diagnostics may differ from the expanded-library recipe.

### Current tinyAM translation

The candidate is intended as a structurally informed approximation rather
than an exact operational replication. In particular, the logistic age curves
are grouped by survey type, the fitted q effects use a logit rather than log
link, and the source's selectivity priors, environmental measurement process,
ADF&G increment-specific penalties, age-reading error and composition
likelihoods are not fitted directly. Official F-at-age is unavailable for
direct comparison. The older custom paired-age q blocks and unpenalized
ADF&G annual-factor construction are no longer used. No new numerical
agreement measures, convergence codes or objective values are recorded here:
those require the revised candidate fit and comparison outputs. Previous
fitting results should not be attributed to this candidate.
