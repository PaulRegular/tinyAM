# Simulated-cod example exploration

The primary selection criterion is absolute mean annual relative SSB error below
10% over all displayed years, with successful numerical checks. This is a
single-realization teaching criterion, not an estimate of repeated-sampling bias.
Mean absolute relative error and terminal error are reported alongside the signed
mean. Abundance and recruitment errors are secondary checks. No annual years are
removed to improve the criterion.

The original 80 x 80 simulation, seed 438, packaged inputs and local raw caches
are retained. Trials use separate ignored files under
`data-raw/cache/sim_cod_exploration/`; each attempt records settings, diagnostics,
warnings/failures and numerical results. Run one fitting/simulation job at a time.

1. On the original data, replace IID N deviations and constant M with exact
   cohort survival and age/year AR1 M, retaining the supplied M baseline of 0.3.
   Test the existing separable and AR1 F specifications separately with default
   starts and fitting controls.

   Both numerical fits passed. SSB mean relative errors were -63.41% (separable
   F) and -56.59% (AR1 F); terminal errors were -83.38% and -80.79%. Thus shifting
   process variation to M alone did not solve the recovery problem.
2. Retain the closed N and AR1 M models, but replace the broad q age groups with
   RV logistic and Sentinel double-logistic curves. The previous q group 1--3
   pooled ages with very different RV length-selectivity; it could distort
   inferred recruitment and apparent survival. This is a separate q hypothesis,
   evaluated from default starts, not a retry of the earlier model.

Further trials and their outcomes will be recorded here before selecting the
final example. No exploratory trial changes package likelihoods or settled
assessment specifications.
