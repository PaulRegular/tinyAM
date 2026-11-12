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

## Initial age structure and catch precision on SAIO `312bc721`

Both trials used the same source observations and the fixed-SD M random walks described above. Replacing free first-year abundances with `N_settings$init = "exp"` imposed a mortality-based initial decline. It reduced initial age-16+ abundance from 77.0 to 16.45 million, compared with 13.82 million in the VPA, and reduced the median total abundance difference from 22.3% to 9.2%. The median SSB difference remained 55.4%. Ages 6–15 accounted for most of the +57.1% terminal SSB excess.

Fixing the catch log-SD to 0.05 approximated the VPA's treatment of catches as an input. Relative to the exponential-init fit, median recruitment difference improved from 26.7% to 17.0% and median SSB difference increased to 78.2%. Both warm-started fits passed numerical checks. The results show that regularizing the first-year age composition and tightening catch alone leave large mature-abundance differences.

## Mobile Sentinel timing and residuals

All 169 Sentinel indices matched Table 4. The 2003 vessel exclusion applies to the RV survey. The Sentinel survey ran in August, so the recipe explicitly uses 0.625 in place of the shared 0.75 sampling-time approximation in the selected SAIO rows. The resulting registered specification converged from defaults, with a 38.9% median SSB difference and a −11.5% terminal difference. [2016 assessment](https://waves-vagues.dfo-mpo.gc.ca/Library/365839.pdf) · [Survey timing](https://publications.gc.ca/collections/collection_2016/mpo-dfo/Fs70-7-2016-016-eng.pdf)

The official VPA also underpredicted Sentinel early and overpredicted it late. In the source-SD exponential-init diagnostic, changing timing alone left this pattern. Adding one log-linear catchability trend shared across Sentinel ages reduced Sentinel log-RMSE from 0.690 to 0.389. The 2015 age-16+ prediction fell from 235 to 106 against 83 observed. Estimated q fell about 16% per year. Terminal SSB remained 62.7% above the common-definition reference. The trend is a diagnostic of the difference between survey trajectories, with population scale and survey availability still confounded.

The VPA estimated 26 time-constant catchabilities, one per age for each of the two surveys, with separate lognormal variances by survey and age. The ratio of catches to abundance at ages 15 and 16+ agrees within the rounding of Tables 2 and 7 in all 40 years. With shared M at those ages, that supports a shared oldest-age F state. [2016 assessment, pp. 6–7 and Tables 2 and 7](https://waves-vagues.dfo-mpo.gc.ca/Library/365839.pdf)

A trial combined that F block, time-constant survey catchability by age, mid-August Sentinel timing, free initial abundance, the 0.075 M random-walk SD and a fixed 0.05 catch log-SD to approximate the VPA's direct use of catch. It fixed initial M at the published prior means because the prior widths were unavailable. The trial passed numerical checks. First-year age-16+ abundance fell from 77.0 million in the warm-start fit to 12.0 million, versus 13.8 million in the VPA. Catch log-RMSE fell from 0.220 to 0.020, while Sentinel log-RMSE rose from 0.688 to 0.707. SSB trend correlation was 0.994, with a 54.1% median absolute difference. Fitted 2015 SSB was 27,524 t versus 59,928 t on the common definition and 55,023 t in the published VPA. The joint trial leaves the source of the remaining scale difference unresolved. The registered recipe remains preliminary.

## Fixed effects for F by year and age on SAIO `312bc721`

The VPA estimates 12 terminal abundance parameters and reconstructs cohorts and historical F from catch and survival. The random walks in its objective apply to M. A forward `fit_tam()` approximation turns off the F process and estimates fixed year-by-age effects, with common F at ages 15 and 16+. A fixed catch log-SD of 0.05 keeps positive catch predictions close to the inputs. The VPA instead uses catch directly in the reconstruction. tinyAM integrates the Gaussian M effects, while the assessment uses M penalties in its ADMB objective and reports MCMC summaries. The two published initial-M prior means are available, and the prior widths are unspecified. [Assessment, pp. 6–7](https://waves-vagues.dfo-mpo.gc.ca/Library/365839.pdf#page=15)

The first trial estimated all 480 F coefficients with recruitment and M held to the preceding trial's specifications. It hit the evaluation limit and had a non-positive-definite Hessian. Six year-by-age blocks have no positive catch, and their F coefficients were weakly identified. Sharing each of those six coefficients with the nearest supported year at the same age gave 474 coefficients and passed numerical checks. A seventh printed zero is in a shared oldest-age block with a positive catch at age 15. All catch and index observations were preserved.

| Trial | Median absolute SSB difference | SSB trend correlation | 2015 SSB |
|---|---:|---:|---:|
| Previous F random-walk trial | 54.0% | 0.994 | 27,524 t |
| Fixed year-by-age F with shared coefficients for the six unsupported cells | 62.6% | 0.989 | 22,215 t |
| Fixed F and free annual recruitment | 63.3% | 0.987 | 21,314 t |
| Fixed F and free annual recruitment, estimating both initial M rates freely | 75.2% | 0.983 | 14,833 t |
| Fixed F and free annual recruitment, supplying official annual M | 9.9% | 0.997 | 55,243 t |

The comparison uses the same 40 years and SSB calculated from Table 7 abundance, with 59,928 t in 2015. The annual-M diagnostic supplies the 1976–2015 rates for both age groups from Figure 4 of the rebuilding plan. Its median absolute SSB difference against the rounded Figure 2 trajectory was 1.8%, and its 2015 estimate was close to the published posterior median of 55,023 t. It passed numerical checks and reduced the median absolute differences in abundance and recruitment to 15.4% and 19.4%. Against the fixed F and recruitment control, Sentinel log-RMSE fell from 0.728 to 0.690 and RV log-RMSE rose from 0.402 to 0.456. This demonstrates the scale sensitivity to M when F and recruitment are free. Agreement in the conditional fit depends on using official estimates of M. The registered recipe remains preliminary. [Rebuilding plan, Figs. 2 and 4](https://www.dfo-mpo.gc.ca/fisheries-peches/ifmp-gmp/groundfish-poisson-fond/american-plaice-plie-canadienne-eng.html)
