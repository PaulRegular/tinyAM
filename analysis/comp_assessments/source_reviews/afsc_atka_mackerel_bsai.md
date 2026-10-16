# Bering Sea/Aleutian Islands Atka mackerel

## Official run and sources

`afsc_atka_mackerel_bsai_2024` represents the 2024 operational update of
Model 16.0b. The December 2024 SSC report retains that model. The 2026 NOAA
harvest-specifications review confirms that the 2025 funding lapse did not
produce a replacement fitted assessment. It is one combined-sex BSAI stock,
not three Aleutian area cases. The Charbonneau–Keith starting identity is
`AFSC_BSAI_Pleurogrammus_monopterygius`.

The [2024 SAFE report](https://files.npfmc.org/SAFE/2024/BSAIatka.pdf) is cached
with a checksum. Earlier framework inputs and code are pinned to
[`atkawham` revision a65171a](https://github.com/JaneSullivan-NOAA/atkawham/tree/a65171af2c6c90cdc31a5d575e3eb6c95f299455/mod22/base).
They clarify unchanged maturity and timing, not the 2024 observations or fits.
The report refers to the AMAK production implementation as available on request;
the accessible framework run ends in 2022 and is not relabelled as a 2024 fit.

## Inputs and likelihood

The model covers 1977–2024 and ages 1–11+. It fits total catch biomass with
CV 5%, survey biomass with year-specific supplied CVs, and separate
multinomial fishery and survey age compositions. The database keeps totals
and compositions separately. Compositions are recovered by normalizing the
report's expanded numbers, retaining report rounding. Table 17.8 is observed
expanded catch-at-age, not fitted predicted catch. Table 17.10 is design-expanded
survey abundance, not raw age-sample counts; its `n` column is retained separately.
Missing early catch cells (`---`) are not converted to zero or normalized as a
complete composition. Only complete reported rows enter the recovered composition.

Original population, catch and survey weights remain distinct, in kg/fish.
The population vector averages the latest three surveys with ages. The terminal
catch input is the assumed full-year 2024 TAC, 72,987 t, rather than the in-year
64,439 t observation. The supplied female maturity ogive is unchanged from
the framework. AMAK halves this ogive internally for combined-sex N.
Its month 7.5 survey and month 8 spawning become 6.5/12 and 7/12 of the year.

Natural mortality is fixed at 0.3. Survival is deterministic given mortality;
recruitment follows Beverton–Holt with fixed steepness 0.8 and an estimated
lognormal recruitment SD. F is annual intensity times flexible, penalized
age selectivity, with temporal variation; ages 10+ share selectivity.
Survey selectivity is time-invariant and nonparametric. Its scalar q is estimated
under a prior with mean 1 and log variance 0.2². Ages 4–10 define the q/selectivity
normalization. This does not imply fixing q to one.

Exact 2024 supplied composition sizes and the unrounded index CVs are not
available in the recovered report tables. Published rounded index CVs are
retained as aggregate relative standard errors. Neither raw otolith counts
nor fitted effective sample sizes are substituted for supplied likelihood weights.

## Outputs and verification

The report supplies begin-year N in millions, normalized selectivity,
fully-selected F, native female SSB and biomass with approximate 95% confidence
bounds, and age-1 recruitment with natural-scale SE. N and selectivity are
embedded tables, visually transcribed into the importer and checked against
recruitment and biomass. F-at-age is explicitly derived as normalized
selectivity times fully-selected F. It has no invented uncertainty. Values from
the adjacent 2022 comparison columns are excluded.

Completeness remains partial because original unrounded inputs, full composition
weights and state covariance were not recovered. This is distinct from the
confirmed absence of a fitted length likelihood or a stochastic N process.

## tinyAM translation and numerical result

The first translation uses 1990–2023 and ages 2–11+, avoiding incomplete early
catch cells and absent age-1 catch. The 2024 total catch and survey biomass have
no matching age compositions, so cannot enter this age-index approximation.
It retains recorded catch/survey weights, fixed M, female maturity and source
timing. It uses deterministic survival, free initial ages, independent F random
walks and constant-mean IID age-2 recruitment. The native age-1 Beverton–Holt
relationship, F penalties, q prior, ageing error and separate composition
likelihood are omitted. Rounded survey CVs supply log SD; their additional
composition uncertainty is not represented. Source observations are not adjusted.

Default initialization failed (optimizer code 1, false convergence; maximum
gradient 0.934; non-positive-definite Hessian). One retry using aligned official
N and F converged (code 0, objective 511.565, maximum gradient 0.000267,
positive-definite Hessian and successful SD report). Data and assumptions were
identical between attempts. Both diagnostic rows are retained in the cached fit.

The dashboard retains native age-1+ SSB and also compares modeled-age August 1
female SSB calculated from rounded official N/F, fixed M and recorded biology,
without uncertainty. This approximate SSB comparison has MAPD 25.7%, terminal
difference −23.6% and correlation 0.966. N-at-age has MAPD 30.5%, terminal
difference −22.2% and mean within-age correlation 0.874. F-at-age has MAPD 265.6%,
terminal difference +33.8% and correlation 0.647; near-zero official juvenile F
makes relative differences especially large. These are flagged for review.
Age-2 recruitment (16.1%, +8.3%, 0.973) and total abundance (19.1%, −12.8%, 0.985)
do not cross the requested thresholds. Fixed-M agreement is input consistency.
