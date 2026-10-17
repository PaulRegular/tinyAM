# Initial abundance: source-checked SAM trials

The six cached SAM runs below all have `initState = 0` and empty `initN`.
Their older initial states have no initial density. tinyAM's new `latent`
option integrates older initial states without estimating an initial SD.
It retains **r0 as an estimated fixed effect**, so it is closer to this SAM
boundary, not identical.

`scripts/translation/trial_initialization.R` regenerates both baseline and
latent trials, comparisons, warnings and local dashboards. The source branch
was `e18c18f`; trials used initialization implementation `b18873c` plus the
unit-only migration. All observations and other scientific settings were held
unchanged. The candidate starts at the baseline's fitted population; no
optimizer retries or model simplification were used. Raw objectives were not
used to compare fixed and integrated initial states.

| Stock | Baseline → candidate | Candidate adjusted max gradient | Optimizer | pdHess | Retained |
|---|---|---:|---:|---|---|
| Blue whiting | free → latent | 0.000391 | 0 | yes | latent |
| Northeast Arctic cod | free → latent, fitted start | 0.000560 | 0 | yes | free; see below |
| Northern Shelf haddock | free → latent | 0.000692 | 0 | yes | latent |
| Western Baltic herring | free → latent | 0.00301 | 1, false convergence | yes | free |
| North Sea saithe | free → latent | 0.00129 | 1, false convergence | yes | free |
| North Sea whiting | exp → latent | 0.000816 | 0 | yes | latent |

Saithe has a large raw gradient (~51.8) at a catchability increment bound;
the table uses the bound-adjusted gradient. Both rejected candidates fail
optimizer checks despite positive-definite Hessians. No rescue fits were made.

The retained scripts were then fitted as written. Blue whiting, haddock and
whiting passed (adjusted gradients 0.000541, 0.00119 and 0.00149). The cod
candidate failed from ordinary package starts: adjusted gradient about
1.72 billion and non-positive-definite Hessian. Its existing **free** recipe
is therefore retained despite the successful fitted-start experiment. Whiting's
existing source-based starting values now also include the source first-year
older abundances; these are starting values, not constraints on the estimates.

## Agreement and trade-offs

Median absolute percent differences and trajectory correlations are assessed
separately. Fixed M is excluded as evidence of improvement.

- **Northeast Arctic cod, experimental fitted-start candidate:** SSB difference falls from 7.15% to 6.30%; SSB trend
  correlation rises from 0.963 to 0.992. N and F comparisons also improve, but
  starting-value sensitivity prevents retaining this candidate.
- **Northern Shelf haddock:** agreement is essentially unchanged; median N/F
  differences fall slightly. The initialization better represents the source boundary.
- **Blue whiting:** N/F/recruitment differences fall slightly and their trend
  correlations rise. SSB difference rises from 4.40% to 4.71%, while its trend
  correlation rises from 0.985 to 0.986. Deterministic subsequent survival remains
  a material difference from SAM.
- **North Sea whiting:** N and SSB differences fall from 15.31% to 13.89% and
  9.21% to 7.96%. Recruitment difference is essentially unchanged, but its trend
  correlation falls from 0.868 to 0.843. This trade-off remains visible;
  initialization does not improve every comparison.

Other SAM recipes are unchanged: the first-state configuration was not verified
from a matching cached native fit in these trials. The database is unchanged.
Detailed differences, every warning and failure, fitted objects and dashboards
remain ignored under `results/cache/initialization_trials/`. Aggregate result
CSVs are not refreshed from these targeted runs.
