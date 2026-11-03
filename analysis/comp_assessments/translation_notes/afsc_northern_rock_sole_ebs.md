# Eastern Bering Sea northern rock sole: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_northern_rock_sole_ebs.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_northern_rock_sole_ebs.md).

## Expanded-library review

### tinyAM translation and numerical result

The recipe represents female ages 5–20+ over 1982–2023. It converts original
biomass and joint compositions using both sexes' empirical weights before
retaining female observations and restricting ages. Female biology and
February 1 spawning are retained. Female M is estimated as one constant;
the official estimate is not converted to a fixed input. Native male states
remain canonical but are outside this fitted female population.

The independent F random walks, constant-mean age-5 recruitment, logistic
survey q and independent age observations approximate the source curves,
stock–recruit relationship and joint composition likelihood. Source M/q
priors are omitted, and the different source composition/biomass timings
remain an approximation. The original observations are not tuned to the fit.

Neither permitted initialization attempt passed numerical checks. Default
starts returned code 1 (singular convergence), objective 3866.355, maximum
gradient 1.203 and a non-positive-definite Hessian after 193 seconds. One
aligned female N/F retry returned code 0, objective 3963.453 and gradient
0.0788 after 148 seconds, but still failed the gradient tolerance and Hessian
checks. Fixed-parameter SEs were non-finite; an SD-report object alone does
not mean uncertainty estimation succeeded. No further retry or model change
was made.

The returned retry drives M close to zero, the scalar q baseline extremely
high and initial young-age N upward, while some old initial states approach
zero. These are failed-fit estimates, not evidence of stock status or model
agreement. The source priors and joint observation structure differ, but the
two attempts do not isolate the cause. An exploratory dashboard and numerical
comparisons are retained with the failure diagnostics; female-only comparisons
and rounded-derived F/SSB are labelled explicitly.

An exponential initial-abundance sensitivity was also fit with package-default
starts and the same remaining settings. It did not resolve the instability:
optimizer code 1, gradient 0.217, non-positive-definite Hessian, and
non-finite fixed-effect covariance. N-at-age and SSB median percent differences
rose to about 100% and 864%, respectively, and abundance and recruitment
differences became orders of magnitude larger. The free initial-age treatment
is therefore retained as the closer source representation, although neither
initialization has produced a numerically reliable fit.
