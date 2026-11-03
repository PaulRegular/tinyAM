# Aleutian Islands walleye pollock: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_pollock_ai.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_pollock_ai.md).

## Expanded-library review

### tinyAM translation and numerical result

The translation uses 1978–2018 and ages 1–15+, ending at the latest matching
age observations. Aggregate-only survey/catch observations after 2018 cannot
enter this approximation. It retains original catch weights and uses the
recorded annual survey weights as an explicit population-weight approximation.
It fits deterministic survival, free initial abundance, a Beverton–Holt median
recruitment curve with estimated parameters and IID SD, independent F random
walks and one estimated M. The source fixed steepness, M prior, q normalization,
F penalties, composition likelihood and ageing error are not reproduced.
Survey time 0.5 approximates summer; March spawning and GOA maturity are retained.

Default initialization failed (code 1, objective 642.037, gradient 2.53e8 and
a non-positive-definite Hessian). The one permitted retry changed only starting N
and M. It reduced the objective to 625.738 and the gradient to 2.18e-5, but
still failed with code 1, a non-positive-definite Hessian and unavailable
fixed-parameter covariance. The SD-report object does not provide usable
uncertainty in this case. Both attempt records are retained; no further model
tuning or initialization retry was conducted.

A source-prior sensitivity approximated the reported M prior (mean 0.2, CV 0.2)
with a moment-matched lognormal penalty on log M, holding the observations and
other model settings fixed. The package-start attempt failed (code 1, objective
635.205, gradient 4,976 and non-positive-definite Hessian). A single aligned
source-start retry reached objective 628.863 and gradient 3.31e-5, but still
failed the fixed-parameter covariance check (code 1; no usable uncertainty).
The resulting N scale was orders of magnitude beyond the official estimate, so
its point comparisons are not interpretable. This prior approximation alone
does not resolve the q/N/M identifiability problem; it was not retained.

Usable predictions permit an **exploratory** dashboard with the numerical failure
visible. N rises to roughly 10^12–10^16 fish, F falls to 10^-12–10^-8 and survey
q to 10^-10–10^-8. Unit checks confirm this is an extreme fitted scale trade-off,
not confusion between the report's billions, recruitment millions and tinyAM's
individual fish. Removing the source's fixed-q normalization and M prior is an
important identifiability limitation of this approximation. It must not be
interpreted as a successful reproduction. Native SSB is labelled approximate
because population weights differ; historical F-at-age is unavailable and
arithmetic Fbar is non-equivalent. All available abundance/SSB comparisons and
estimated M are flagged for review; missing F comparisons are not passes.

## Paper-working-set review before migration

Only distinct notes from the earlier main review are retained below. Their settings and diagnostics may differ from the expanded-library recipe.

### tinyAM translation and numerical result

Default initialization failed with singular convergence (code 1), objective
625.738, gradient 0.000419 and weak stock–recruit shape support. The one permitted
retry changed only starting N and M. It reached the same objective and also
failed (code 1, gradient 0.0000218, non-positive-definite Hessian and non-finite
fixed-parameter covariance). A returned SD-report object is not usable uncertainty
in that situation. Both attempt records are retained; no further model tuning
or initialization retry was conducted.
