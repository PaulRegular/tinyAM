# U.S. West Coast sablefish

## Curated translation

The curated fit uses IID log-abundance deviations around expected cohort
survivors, with one estimated log-scale SD (0.178). The official Stock
Synthesis model calculates older-cohort survival from F and M. Its annual
population deviations apply to recruitment. The IID cohort process gives
the pooled tinyAM model flexibility to represent unexplained variation in
apparent survival. Differences in sex composition, age coverage and the
observation likelihood are plausible contributors. [Official assessment](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/) · [Native model files](https://www.pcouncil.org/documents/sablefish-2025-model-files/)

On SAIO `312bc721`, the fit passed numerical checks from default starts.
All 42 years use the same comparison definitions as the fixed-gear fit
below. Natural mortality increased from 0.0344 to 0.0456 per year. The
soft Gaussian prior on log M still has support outside the official bounds.

| Quantity | Median absolute difference, N off → IID | IID trend correlation | IID 2024 difference |
|---|---:|---:|---:|
| SSB | 43.2% → 20.7% | 0.714 | −3.8% |
| Total abundance | 49.0% → 31.4% | 0.959 | −39.1% |
| Average F, ages 3–15 | 27.7% → 15.2% | 0.797 | +5.8% |
| Recruitment | 63.4% → 46.6% | 0.892 | −80.9% |

SSB was 59.0% below the sex-pooled reference in 1983. Its median absolute
difference was 5.4% in 2019–2024. The closer recent SSB and average F support
including the comparison in the curated set. Whole-period SSB correlation
fell from 0.815 to 0.714, and abundance tracks the source trend at a lower
level. The reference comparison remains approximate because of the pooled
biology and the reconstructed catch and survey numbers-at-age.

## Earlier preliminary fits

The recipe fits 2006–2023, ages 0–30+, using the original annual empirical
weight schedules. Marginal joint female/male frequencies are normalized once
over both sexes, excluding inactive source rows and preferring joint samples
over separate unsexed samples. Each retained/discard fleet is reconstructed
separately with its own sex-specific catch weights. Only complete removal
years are summed; 2014 remains wholly missing because pot ages are absent.
No catch or survey observation is derived from fitted population states.

WCGBTS totals are combined with their own age samples and weights. Its
unknown physical multiplier is not invented; q remains relative. Age-bin
tails use the empirical age-30 terminal weight. Combined-sex samples use a
1:1 female/male weight mixture. These are explicit approximations, not
claims of exact SS3 likelihood reconstruction. Source age error, sample
weights, other surveys and the recruitment-deviation index are not reproduced.

The population approximation has deterministic pooled survival, IID
recruitment about a mean, exponential initial abundance, AR1 F around a
logistic age mean, and one estimated M. Native stock–recruitment dependence,
recruitment bias adjustment and the M prior are omitted. Female maturity is
recovered from the supplied mature-weight/mid-season-weight ratio. A female
beginning-weight proxy and a fixed 1:1 sex ratio simplify native two-sex
biology. Source N is summed; source fleet F is added before N-weighting rates
over sexes/terminal ages. Pooled uncertainty is unavailable and not invented.

An initial setup passed numerical checks but had incorrect biological
vectors: scalar-sex selection reduced weight lookups to age 0, recycling that
value over ages. The agreement review exposed this through identical relative
differences for abundance, biomass and mature biomass. Its observations were
correctly reconstructed, but the biological setup and comparison are invalid.
The fit and dashboard are retained only in
`results/preliminary/cache/noaa_sablefish_pacific_coast_2025/experiments/biology_lookup/`.
The lookup now expands the sex label over the age vector, and regression
checks compare every weight and maturity cell directly with the original
source schedules. The corrected fit is evaluated separately.

With correct biology, the default attempt passes numerical checks:
optimizer code 0, objective 1298.976, maximum gradient 0.0000588 and a
positive-definite Hessian. The recipe uses package starting values. Constant
M is about 2.4e−10 per year, with log-M SE about 3023. The estimate is
effectively uninformative after omission of the source prior; numerical
convergence does not make this a defensible estimate of M.

Agreement also needs review. Mean absolute percent differences are 39.4%
for abundance, 225.7% for N-at-age, 49.8% for F-at-age, 102.4% for recruitment
and 29.1% for approximate common-definition SSB. Their respective trend
correlations are 0.990, 0.499, 0.842, 0.920 and 0.835. SSB's terminal difference
is +3.6%, but its mean difference and correlation still trigger review.
Estimated M differs by effectively −100%; constant source M has no meaningful
trend correlation. The corrected dashboard and fit are in the main case
cache, clearly separate from the invalid biological setup.

Source facts and native-file provenance remain in the
[SAIO review](https://github.com/PaulRegular/saio/blob/757b59c/discovery/source_reviews/noaa_west_coast_sablefish_2025.md).

## Curated expansion review on SAIO `cecb45ad`

The estimated-M recipe reproduced the near-zero M estimate. It passed numerical checks and gave a median absolute SSB difference of 29.4% over 18 common years, with trend correlation 0.835 and a terminal difference of +3.6%.

Fixing M at the source prior median of 0.072 per year gave an implausibly large population and virtually zero fishing mortality. Both default starts and a restart from the tinyAM estimated-M fit led to the same pattern. Relative to the reference, the median SSB ratio exceeded 300 million in both fixed-M fits. Optimizer and Hessian checks passed, and several estimated parameters were highly correlated. The two starts indicate that the scale problem persists in this pooled catch and survey specification. Sablefish remains preliminary.

## Survey scale sensitivity

The official WCGBTS selectivity has a sharp peak at age 1 and a nearly flat plateau from age 2. With M held at 0.072 per year, estimating separate q for ages 0, 1 and 2+ still gave almost zero F and an implausibly large stock. A diagnostic then divided the reconstructed index by the official fitted q of 1.0497 times the native relative selectivity (approximately 0.077, 0.995 and 0.801). This assumes the biomass index and sampled age composition can be combined into an absolute numbers index. The catch inputs, biology, and F and M settings were identical across trials. [Official assessment, section 3.4.1](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/)

With unit q on that adjusted index, the fit passed numerical checks and gave a median SSB difference of 13.3% over 18 common years, with trend correlation 0.895 and a terminal difference of −17.8%. Compared directly with the published female SSB, the median difference was 15.9% and the correlation was 0.900. Freeing the overall q multiplier on the adjusted index brought q to about `1.4e-9`, median Fbar to about `8.7e-11` per year and the median SSB ratio to about 611 million. Correcting the relative age shape alone leaves the scale trade-off unresolved. The registered recipe stays preliminary because the plausible scale in the diagnostic depended on official fitted catchability. Extending the 2006–2023 window through the full survey and biology overlap would add four survey years with incomplete catch reconstructions.

## Independent scale trials on SAIO `312bc721`

The sablefish source record matched the previous snapshot exactly. Three trials kept fixed M at 0.072 per year and estimated q for ages 0, 1 and 2+. Each changed one aspect of the short, fixed-M fit with grouped q.

- Separable F, `~ factor(year) + logistic(age)` with no F process deviations, recovered finite Fbar (median 0.0073 per year) and q. The median SSB difference was 481% over 2006–2023.
- Supplied survey log-SDs used the native biomass-index SE and a delta approximation for age-composition expansion with estimated mean weight. Sample sizes used the published WCGBTS Francis multiplier of 0.053. Positive-cell log-SDs ranged from 0.18 to 2.25. The approximate independent likelihood again gave nearly zero F and an implausibly large stock. [Official assessment, Table 27](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/)
- Adding 2003–2005 and 2024 extended the survey window and supplied no new complete catch vectors. Median Fbar was 0.021 per year, and the median SSB difference on the original 18 years was 169% with trend correlation 0.820.

All three passed numerical checks. The results expose sensitivity to both F flexibility and survey coverage. They leave the absolute scale unresolved with freely estimated q. The current translation remains preliminary.

## Three-gear trial

The source supports separating trawl, hook-and-line and pot from 2006–2023. Retained and dead-discard catches were reconstructed separately using the native age samples and sex-specific weights, then summed within gear. The source discard totals already reflect mortality. The 2014 retained pot age sample is missing, so pot-at-age remains missing that year. Trawl and hook-and-line are available. Combined catch weights reproduced all 53 complete gear-year source yields to rounding precision. The recorded total yield remains available for 2014. [Official assessment, section 2.1.2](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/)

The trial used independent AR1 F deviations within each gear, shared process parameters and gear-specific normal mean curves. It kept fixed M at 0.072 and free q at ages 0, 1 and 2+. The fit passed numerical checks, but median SSB was about 269 million times the common-definition reference and median Fbar was `1.90e-10` per year. Terminal fitted yield was 6,083 t versus 6,385 t in the source, illustrating that low F and high N can still reproduce catches. The fleet split is feasible with these inputs and leaves the population scale unresolved. The registered pooled recipe remains preliminary.

## Female-only trial

Females made up 62.7% of reconstructed retained catch numbers in the 17 complete years from 2006 to 2023. They still represented 59.3% of reconstructed removals if all unsexed dead discards were assigned to males. The trial used female catches and WCGBTS observations from the joint sex-and-age expansions, female biological inputs and an equal split by number within age for unsexed discards. The 2014 catch vector remained missing because the retained pot age sample was unavailable. [Official assessment, sections 2.1.2, 2.5.2 and 3.4.1](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/)

The trial retained pooled AR1 F, fixed M at 0.072 per year and free survey q at ages 0, 1 and 2+. From default starts it passed numerical checks, with a positive-definite Hessian and maximum gradient of `7.36e-5`. Median Fbar was `1.61e-9` per year, and fitted SSB was about 44 million times the matched female reference on the median ratio across 18 years. The female-only reconstruction leaves the scale problem unresolved with freely estimated survey q. The registered recipe remains preliminary.

## WCGBTS catchability profile on SAIO `312bc721`

The profile fixes a single multiplier on the source catchability pattern in the pooled 2006–2023 fit. The index is divided by `1.0497 * c(0.0769921, 0.994958, 0.8008)` for ages 0, 1 and 2+, respectively, so a multiplier of 1 gives the source normalization. Catch, survey inputs, fixed M of 0.072, AR1 F and IID recruitment were held to the preceding survey-scale sensitivity. At each of 24 multipliers from `1e-11` to 4, the remaining parameters were refitted under the same Laplace likelihood. A sweep from the conditioned fit and a sweep from the free-q fit both passed optimizer and nuisance-gradient checks at every point. The objectives agreed within `8.5e-8`, and the rebuilt objective reproduced both archived fits. [Official assessment and native model files](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)

| WCGBTS multiplier | Marginal NLL above best grid point | 2023 SSB (t) | Median Fbar (per year) |
|---:|---:|---:|---:|
| 1 | 60.616 | 54,903 | 0.0482 |
| 0.1 | 7.351 | 622,667 | 0.00583 |
| 0.01 | 0.713 | 6,292,968 | 0.000608 |
| 0.001 | 0.071 | 62,994,227 | 0.0000611 |

The marginal NLL improved toward the smallest tested multiplier and varied by less than 0.001 from `1e-5` to `1e-11`. Across that range SSB increased by a factor of about one million and Fbar fell by the same factor. Relative to multiplier 1, catch log-RMSE fell from 0.767 to 0.648 and survey log-RMSE fell from 0.875 to 0.820. At low F and fixed M, catches are approximately proportional to N × F and survey predictions to q × N. Increasing N while decreasing F and q can therefore preserve predictions and leave cohort survival nearly unchanged. The tiny positive multiplier from free optimization lies on this flat part of the profile. This translation fits the observations more poorly at the official normalization and provides little information about absolute scale at small multipliers. The registered recipe remains preliminary.

## Earlier source-conditioned sensitivity on SAIO `312bc721`

The official assessment runs from 1890 through 2024 with two sexes, ages 0–70, six removal streams and four biomass survey indices. It also uses an environmental index for recent recruitment deviations. The expanded translation pools the sexes and fits 1983–2024 at ages 0–30+. The first usable fishery and survey age samples are in 1983. Earlier removals in SAIO are tonnes alone. The translation retains separate trawl, hook-and-line and pot streams for retained catch and dead discards. Each active stream-year composition is expanded to numbers with its source biomass and sex-specific weight-at-age. Unsampled stream-years stay missing. The aggregate yield display retains the reported total removals. [Official assessment, sections 2.1–2.3 and 3.4](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/)

The survey input uses 31 year-series combinations with active age samples. Early Triennial has three years, Late Triennial two, NWFSC Slope five and WCGBTS 21. The source estimates a constant, survey-specific catchability against vulnerable biomass. The four fitted values are 1.02287, 2.77008, 0.509779 and 1.04970, respectively. WCGBTS peaks at age 1 and plateaus around 0.80 from age 2. Early Triennial drops sharply after age 3, Late Triennial has a more gradual decline from an age-2 peak, and NWFSC Slope peaks at age 4 before levelling off around 0.78. The native `Report.sso` contains the fitted `AGE_SELEX` `Asel2` rows, `INDEX_1` `Qbase` and additive log-SDs under `PARAMETERS`. The shapes agree across reported years and both sexes. [Native model archive](https://www.pcouncil.org/documents/sablefish-2025-model-files/)

The fitted age-index proxy divides each expanded index by the official catchability at its original age before pooling age 30+. Using a unit multiplier on that scale conditions stock size on the official survey normalization. It also transfers the official age patterns. The native likelihood fits biomass indices and marginal age samples separately. The translation uses one fitted log-SD per survey and independent age observations, so the native aggregate precision, ageing error and composition sample weights are outside this likelihood.

The conditioned trial fixed common natural mortality at the official estimate of 0.087777 per year and used IID log recruitment. It averaged female and male `Asel2` equally and scaled each curve to a peak of 1. Fourteen fleet-specific time blocks reproduced the reported curves over 1983–2024. The fit estimated a common exponent on those curves and an AR1 annual intensity effect for each fleet. Log reported catch biomass, centred within fleet, informed F intensity where ages were missing. Reported biomass also scaled the age observations. The trial had no separate catch-biomass likelihood.

Stock weights average the female and male January schedules. The maturity proxy gives pooled fish half the official female mature mass at age, approximating SSB under a 1:1 sex ratio. SS3's baseline weights apply before 2003 and the oldest supplied weight applies at older ages. The dashboard compares total and fleet F with native rates weighted by numbers across sexes and terminal ages. [Stock Synthesis empirical weight-at-age documentation](https://nmfs-ost.github.io/ss3-doc/SS330_User_Manual_release.html#weight-at-age-file)

The conditioned fit passed the optimizer, gradient and Hessian checks from default starts. Over 1983–2024, the median absolute SSB difference was 13.8% and the Pearson correlation was 0.977 on the pooled reference definition. The 2024 estimate was 46,556 tonnes, 33.7% below the pooled reference and 24.3% below the published two-sex model's female SSB. Median differences were 12.1% for abundance, 58.4% for recruitment and 21.0% for mean F at ages 3–15. The biomass covariate raised median predicted trawl discard before 2006 from 2.5% to 32.3% of reported totals. Those years had no discard age observations.

Estimating four separate survey multipliers in the conditioned model returned singular convergence. All four multipliers fell below 1.3e-8 and fitted SSB increased by orders of magnitude. Without the biomass covariate, the corresponding free-multiplier sensitivity passed numerical checks at 0.55, 0.49, 0.45 and 0.40 of the official values in survey order, with a 98.7% median absolute SSB difference. An unrestricted age–year AR1 F trial for each fleet reached the 20-minute limit. These local sensitivities remain separate from the scoped recipe.

## Scoped fits with estimated selectivity

SAIO supplies the native biomass totals, active and inactive marginal age samples, sex-specific weights, survey timing and index log-SDs. Catch and survey numbers at age are reconstructed from the totals and active age samples. Some years have biomass data only. In particular, the three discard streams have active age samples only in 2006–2023, and the older surveys have three, two and five years with usable ages. The scoped fits retain the 1983–2024 window, six removal streams and four survey series. No estimated official catchability, selectivity or F enters the fitted data or starting values.

The native control file assigns all four survey series the joined double-normal curve (Stock Synthesis pattern 20). It estimates two shape parameters for Early Triennial, three for Late Triennial, two for NWFSC Slope and only the terminal level for WCGBTS. Fixed age-0 terms imply low positive selectivity. Estimated terminal levels give positive older-age plateaus. The tinyAM scopes use independently fitted age levels or capped parametric curves to allow these patterns, with fewer restrictions on the fitted ordering than the source imposes. [Assessment Figure 69](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/), [Native curve implementation](https://github.com/nmfs-ost/ss3-source-code/blob/v3.30.23.1/SS_selex.tpl#L1394-L1474)

The fishing means use fleet scopes and selected changes to age shape at the native 2002, 2011 and 2019 boundaries. Two trials put age–year AR1 deviations directly on log F, with shared process parameters and some older ages grouped. The first gave retained hook-and-line and pot separate double logistics. The second shared a capped quadratic on the log-F scale between the two retained fleets and used a capped logistic and tail shift for the Slope index. Both reached the 20-minute limit before returning a fit.

The constant natural mortality has one Gaussian effect with log-scale mean −2.631 and SD 0.31, using `iid(M_prior_group, sd = 0.31)` around `exp(-2.631)`. This matches the source log-M penalty up to a constant. tinyAM integrates the effect. Its support extends beyond the source's bounds of 0.04–0.13 per year.

The first annual-intensity AR1 sensitivity used the shared retained-gear shape. The hook-and-line and pot discard logistic slopes approached zero, leaving their parameters weakly identified and the fixed-effect Hessian indefinite. Replacing these with one log-F age ramp per discard fleet and ordering WCGBTS catchability as age 0 ≤ age 2+ ≤ age 1 gave a fit that passed numerical checks from default starts. WCGBTS estimated a tie between age 1 and 2+ at the order constraint, and the Slope index estimated an increase from age 4 to its 6+ plateau. Early Triennial estimated low catchability at both ends and Late Triennial peaked at age 2 with a higher older-age plateau than SS3.

The independently fitted annual-intensity model estimated M at 0.0206 per year, below the official hard lower bound. Across 1983–2024 the median absolute difference in SSB was 52.5%, with trend correlation 0.840. Its 2024 estimate was 43,962 tonnes, 37.4% below the sex-pooled official reference. Median differences were 58.7% for abundance, 66.5% for recruitment and 45.5% for average F at ages 3–15. This remains a preliminary translation.

## Bounded natural-mortality prior sensitivity

Stock Synthesis adds the quadratic penalty \(\tfrac12[(\log M+2.631)/0.31]^2\) and constrains M to 0.04–0.13 per year. On the log-M scale this corresponds to a normal prior conditioned on the permitted interval. The bounds retain 94.27% of the original normal mass. Its log mean shifts from −2.631 to −2.63062 and its log SD drops from 0.310 to 0.26593. Its mean on the M scale drops from 0.07555 to 0.07462 per year, with a central 95% interval of 0.04357–0.11922. The SS3 penalty inside the bounds still has the original curvature. [Stock Synthesis prior code](https://github.com/nmfs-ost/ss3-source-code/blob/v3.30.23.1/SS_objfunc.tpl#L1506-L1517)

The one-level Gaussian IID term in tinyAM integrates log M and has no hard bounds. Matching the bounded log moments leaves 1.35% prior probability below 0.04. A stronger sensitivity keeps the same log centre and chooses a log SD of 0.19036, placing 0.1% below 0.04 and 0.096% above 0.13. It penalizes departures inside the interval more strongly than the source penalty.

### Effect of the prior with annual F

The following comparisons retain the same annual AR1 effects on fleet F, survey formulas, observations and 1983–2024 window. Each fit re-estimates the other model parameters.

| M treatment | Fitted M per year | Median absolute SSB difference | SSB correlation | 2024 SSB (t) |
|---|---:|---:|---:|---:|
| Original unbounded IID prior | 0.0206 | 52.5% | 0.840 | 43,962 |
| IID prior matched to bounded log moments | 0.0238 | 50.9% | 0.842 | 45,003 |
| Tighter IID prior with 0.1% below 0.04 | 0.0312 | 46.8% | 0.848 | 47,768 |
| Diagnostic with exact source penalty and bounds | 0.0400 | 40.5% | 0.848 | 52,404 |
| Diagnostic fixed at the official fitted M | 0.0878 | 40.3% | 0.823 | 112,113 |

The exact-bound diagnostic optimizes log M as an outer parameter of the RTMB objective and retains the Gaussian penalty, matching the source treatment of M. It reached the lower bound with a positive log-M objective gradient of 29.5. The largest projected gradient was 0.0012. Fixing M at the official estimate raised terminal SSB to 59.6% above the pooled reference, with a similar median SSB difference over the full series. M changes stock scale substantially, and agreement over the time series remains limited across these settings.

### Fleet F processes with the tighter prior

The age–year AR1 trial removed the annual AR1 formula term and initialized the 3,108 latent F states from the annual model's fitted fleet F. The starting values came from tinyAM. This trial reached 30 minutes before returning a fit.

The IID trial used the same mean selectivity, survey formulas and tighter M prior. It returned after 25.6 minutes with false convergence (optimizer code 1), a non-positive-definite Hessian and maximum gradient about 9.96e9. Its comparisons are diagnostic only. The annual-AR1 baseline with the tighter M prior converged from default starts in about 103 seconds, matching the warm-start fit. The median SSB difference was 46.8% and the 2024 estimate was 32.0% below the pooled official reference.

## Fixed-gear selectivity

The native retained hook-and-line curve is a joined double-normal with a 2002 peak shift and sex-specific adjustments. Retained pot mirrors it. The preliminary recipe shares a piecewise curve on the log-F scale, allowing a steep rise, broader descent, older-age plateau and a change on the younger limb in 2002. Each fleet still has its own baseline and annual intensity. The fixed-gear discard curves use separate estimated age groups. Hook-and-line groups ages 1, 2–3 and 4+ through 2010, then uses finer ages through 5 and shifts the younger-age pattern in 2019. Pot pools its nearly flat ages 2+ through 2010 and 3+ thereafter. The official age-0 selectivity is small and positive. Both discard streams contain positive age-0 catches, so age 0 has its own estimated level in each. [Official assessment, Figure 68](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/) · [Native controls and output](https://www.pcouncil.org/documents/sablefish-2025-model-files/)

The first grouped discard fit allowed separate old-age levels for pot and estimated too steep a decline after age 6. The retained fit uses the simpler plateau and shares coefficients across the shifted hook-and-line periods. It passed numerical checks from default starts on SAIO `312bc721`. Using years with catch ages, the mean absolute differences in annual F-at-age profiles normalized to a peak of one fell from 0.245 to 0.089 for each retained fleet, from 0.349 to 0.050 for hook-and-line discards, and from 0.131 to 0.074 for pot discards. The official comparison weights F across sexes by abundance. Retained profiles now peak at age 6, down from ages 10–12. Their fitted age-0 selectivity is lower than the native reference. The pot approximation smooths over the native curve's slight decline at older ages.

The hook-and-line discard profile peaks at ages 2–3 through 2010 and at age 5 in both later periods. The native peak is at age 4 in 2011–2018 and age 5 from 2019. SSB's median absolute difference fell from 46.8% to 43.2%, and its trend correlation fell from 0.848 to 0.815. Terminal SSB was 27.4% below the pooled reference and fitted M was 0.0344 per year. The recipe stays preliminary.
