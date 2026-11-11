# Southern Gulf of St. Lawrence 4T American plaice: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ck_dfo_platessoides_4t.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_dfo_platessoides_4t.md).

## Expanded-library review

### tinyAM translation

The first translation uses the source years 1976–2015 and ages 4–16+, with commercial catch, both abundance indices, source weight and maturity surfaces, and estimated-M age groups. tinyAM's population and fishing-mortality processes and its catch likelihood simplify the VPA; it does not reproduce the source's VPA catch equation, M random-walk priors, age-specific q and observation-variance parameters, or MCMC integration. The 2003 RV index remains omitted. The 2003 age-16+ zero weight is preserved as supplied when preparing observations and will be evaluated as a potential source limitation rather than silently changed.

The first fit used package-default initialization, with no explicit starting values. It passed the numerical convergence checks (optimizer code 0, maximum absolute gradient 0.00139, successful standard-error calculation, positive-definite Hessian). Comparing the cached fit with the newly recovered outputs gives mean absolute percent differences of about 661% for N-at-age, 66% for total abundance, 89% for age-4 recruitment, and 90% for the digitized native SSB series; trend correlations are 0.865, 0.947, 0.944, and 0.962, respectively. The grouped mortality comparison uses official-N-weighted means of tinyAM ages. F has about 70% mean absolute percent difference, a -91% terminal difference, and 0.914 trend correlation; M has about 66%, a -76% terminal difference, and 0.453 trend correlation. These figure-based rate comparisons are approximate and F coverage is incomplete. The fit is numerically usable, but its mortality and population scale do not closely reproduce the VPA.

## Mortality sensitivity on SAIO `cecb45ad`

The source fixes the innovation SD at 0.075 on the log scale and anchors its initial mortality with priors. Both tinyAM fits use 0.6 and 0.3 per year for ages 4–9 and 10+ in 1976. The baseline random walks start in 1977, with first states of 0.275 and 0.119 and an estimated common increment SD of 0.697. The sensitivity fixes the 1976 rates throughout. All other settings, the observed inputs and the 40 comparison years were held fixed. Both fits passed numerical checks.

| Quantity | Baseline median absolute difference | Fixed-M median absolute difference | Baseline terminal difference | Fixed-M terminal difference |
|---|---:|---:|---:|---:|
| SSB | 37.8% | 97.4% | −14.7% | +95.9% |
| Recruitment | 89.6% | 73.0% | −87.9% | +85.9% |
| Abundance | 66.3% | 68.8% | −70.8% | +77.7% |

SSB trend correlation was 0.962 for the baseline and 0.981 with fixed M. These comparisons use source median N-at-age with the translated weights and maturity. The resulting 2015 reference is 59,928 t. Against the VPA's published SSB series, digitized from Figure 23 with its 2015 median of 55,023 t, the median differences are 44.4% and 115.6%, respectively. A comparable annual Fbar series remains unavailable.

With fixed M, fitted abundance increased and survey catchability fell. For example, 1976 age-4 abundance rose 7.5-fold while RV catchability at age 4 fell from 0.617 to 0.052. Survey predictions changed by less than 0.5% on the geometric mean scale for either survey. The sensitivity exposes compensation among mortality, abundance and catchability. It leaves the effect of the VPA's initial-M priors unresolved. The registered recipe retains the two random walks.

The random-walk fit estimates 65% fewer fish than the source in 1976, while estimating 57 times as many fish at age 16+. Those oldest fish account for 80% of fitted SSB that year and 5% of the reference. This age composition helps explain why SSB can be high despite low total abundance.

A third trial used two `rw(year, by = M_age_group, sd = 0.075)` effects for M. They start from the supplied 1976 rates and use the VPA's fixed increment SD, including the first transition. The trial prepared successfully with the same observed data. From default starts the optimizer returned false convergence, a maximum gradient of `7.08e8` and a non-positive-definite Hessian. It produced no validated comparisons. The registered recipe remains preliminary.

## Follow-up scale checks

Starting the fixed-SD random walks from the converged constant-M fit passed numerical checks. Median abundance and recruitment differences fell to 22.3% and 14.2%, compared with 66.3% and 89.6% in the baseline. SSB agreement was mixed. The mean absolute difference improved from 90.4% to 60.8%, while the median increased from 37.8% to 61.2% and the terminal difference shifted from −14.7% to +65.7%.

A separate diagnostic supplied the official time-varying M medians from Figure 21. The 10+ series was complete. For ages 4–9, 1993, 1997 and 2005 were linearly interpolated between reported neighboring years. The fit passed numerical checks from defaults. Common-definition median SSB difference fell to 29.7%, with trend correlation 0.997 and a terminal difference of −25.2% (−18.5% against the published median). The initial 16+ abundance fell from 793 million to 5.18 million, versus 13.82 million in the VPA. This supports a role for the mortality trajectory in the scale mismatch. The result conditions on the VPA mortality estimates, and the surviving biomass gap keeps the comparison preliminary.
