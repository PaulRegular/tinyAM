# Compare a fitted SAM assessment with tinyAM

This workflow starts from `stockassessment::fitfromweb("WKCOD_combined_99")`.
It does not reconstruct or refit SAM, compare likelihood equality, or change
any tinyAM model equations.

```r
sam_fit <- stockassessment::fitfromweb("WKCOD_combined_99")
tam_obs <- sam_to_tam_obs(sam_fit)
settings <- sam_to_tam_settings(sam_fit)
audit <- sam_to_tam_audit(sam_fit, settings)
# Explicitly select a complete biological input period before fitting.
tam_fit <- do.call(fit_tam, c(list(obs = tam_obs), settings))
sam_comparison <- sam_to_tam_comparison(sam_fit)
vis_tam(model_list = list(SAM = sam_comparison, tinyAM = tam_fit))
```

Install/load this branch and run `source("analysis/comp_sam/001_compare.R")`
from the repository root for the complete worked case. The script selects
**1983–2022**, because original maturity is absent for 1963–1982. It subsets
all input tables explicitly. Neither missing observations nor maturity are
replaced using SAM fitted parameters.

The reference is cached in `cache/`, with retrieval provenance and an RDS
checksum in CSV files. To retrieve a new reference, explicitly remove the
cached reference first. Reference and fitted RDS files and generated HTML
assets are ignored by Git; scripts and exported review tables are retained.
SAM estimates initialize the tinyAM optimizer without fixing any coefficients
or states. Settings are unchanged after audit: IID survival errors, free
initial abundance, independent temporal F random-walk increments, supplied
original M, and exact q/observation-SD sharing within observation tables.

## Read the outputs

- `results/SAM_tinyAM_dashboard.html`: native population trends, age-specific
  states, q, catch/survey predictions, parameters and available intervals.
- `results/common_definitions.html`: shared-definition population trends.
- `results/settings.csv` and `audit.csv`: applied settings and all reviewed
  assumptions, including unsupported features and unresolved legacy fields.
- `results/diagnostics.csv`: optimizer status, gradient, Hessian and uncertainty
  diagnostics, elapsed time, and parameter counts.
- `results/trend_agreement.csv`: descriptive relative differences, trend
  correlation and terminal-year differences. No agreement threshold or ranking.
- Other CSVs retain original inputs, native reports with uncertainty, predictions,
  q comparisons, shared-definition trends, and N/F comparisons by age.

Common SSB uses original maturity, weight and spawning timing for both fitted
N/F surfaces. Common Fbar is arithmetic over ages 2–4. Common catch biomass
uses original catch weights (unavailable for 2022). No SEs or confidence
intervals are invented for recomputed common summaries. Native confidence
intervals retain each model's own definitions and uncertainty calculation.
`se_scale = "log"` means the estimate and interval are natural-scale, while
`se` is a log-scale standard error; it approximates a CV only for small errors.

## Interpretation of the first case

The cached SAM reference records optimizer code 1 ("false convergence (8)")
despite a stored maximum fixed-effect gradient of 8.4e-9, a positive-definite
Hessian and finite reported SEs. These native diagnostics are retained; SAM is
not refitted or silently certified as converged.

The translated tinyAM fit converged (code 0, maximum absolute gradient 0.000403,
positive-definite Hessian and successful uncertainty estimation). Shared-definition
SSB and recruitment follow SAM closely in trend (correlations 0.997 and 0.998).
Their mean absolute relative differences are 3.7% and 4.5%; 2022 tinyAM estimates
are 13.1% and 2.3% lower, respectively. Arithmetic Fbar has correlation 0.970,
mean absolute relative difference 9.4%, and is **51.7% higher in 2022**.
Predicted catch biomass has correlation 0.997, but differs by +34.8% in its last
available year, 2021. Agreement in long-term trends does not establish agreement
in recent fishing mortality or projections.

SAM uses correlated F increments and estimated maturity in this case; the
approximation uses independent increments and original maturity. SAM's initial
states are integrated, while tinyAM estimates free first-year abundance as fixed
parameters. The SAM reference also conditions on 1963–1982 data, whereas the
explicit tinyAM fit begins in 1983. Missing legacy `initState` and
`logNMeanAssumption` settings remain unresolved in the audit. Native Fbar differs
in weighting: SAM is arithmetic, tinyAM abundance-weighted. These differences
are documented explanations to investigate, not demonstrated causal allocations
of the observed discrepancies. No likelihood values or AICs are used to rank
models with different inputs and likelihoods.

Next scientific questions are whether correlated F increments account for the
terminal-age F differences, how the longer SAM history influences 1983 boundary
states, and how fitted versus original maturity affects native SSB. Resolving
these questions requires explicitly controlled sensitivity analyses; this bridge
does not add new model capabilities or undocumented rescue adjustments.
