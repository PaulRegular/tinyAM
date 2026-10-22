# Formula-driven F validation

`deterministic_f_recovery.R` simulates and refits logistic mean F with the F
residual process off. It uses 40 years, seven ages, deterministic cohort
survival, RW recruitment, a survey with constant q, and independent lognormal
catch/index errors (log SD 0.15). No rescue fits or alternative starts are used.
All attempts and warnings are saved under ignored `results/`; the function
returns the final usable fit as an `example_fit` attribute for local inspection.

Seed 1129, 20 assessments:

| Quantity | Truth | Median estimate |
|---|---:|---:|
| Asymptotic F | 0.4 | 0.4021 |
| Logistic midpoint age | 3 | 3.0198 |
| Logistic slope | 0.8 | 0.7908 |

All 20 passed numerical checks, including the gradient and positive-definite
Hessian. This validates a relatively well-informed generating model, not the
usefulness of constant selectivity in every real assessment. No settled stock
recipe was changed.

Historical F equals its formula mean. Formula coefficients still have fitted
uncertainty; no absolute F states, residual SD or residual correlations are
created. Existing Gaussian formula terms, when used, retain their own random
effects: `process = "off"` turns off residual F variation, not formula effects.
Projections continue from terminal F using the existing multiplier/catch control.

Tests cover ordinary/mono/logistic means, exact equations, finite derivatives,
simulation, shared stock-recruit initialization, warm starts, projected catch,
uncertainty, information-criterion counts and a hindcast. A zero-catch test
exposed a diagnostic issue: deliberately zero projected catches have log
prediction `-Inf`, which is valid rather than numerical failure. This exception
is restricted to those declared catch targets; nonfinite historical predictions
still fail. Derived projected F intervals remain exactly zero in that case.

Default models retain identical parameter lists, objectives, gradients and
reported quantities against the pre-expansion snapshot. A local dashboard was
rendered from the same simulated example, with no browser opening in tests.
