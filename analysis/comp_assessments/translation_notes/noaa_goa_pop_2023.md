# Gulf of Alaska ocean perch: preliminary setup review

The recipe reconstructs numbers only in overlapping source biomass and age-
composition years, beginning in 1990. It retains the observed 25+ group and
older population weights/maturity through age 29+. The effective survey
timing of zero follows the native beginning-year prediction equation, not a
claimed calendar sampling date. Reported fitted maturity is an explicit
fixed biological plug-in, with the source 0.5 female fraction. Age error,
length likelihoods, time-varying fishery selectivity and the M prior are
omitted, with each difference identified in the script background.

An initial setup mistakenly specified `M_settings = list(process = "off",
mu_form = ~1)` without clearing the inherited supplied-M default. The actual
model fixed M at 0.2 instead of estimating a constant. Default and one source-
state attempt both failed uncertainty checks (gradients 0.00119 and 0.000266).
These are setup trials, not results for the intended model. Their fit,
attempt table and warnings are preserved under the ignored case cache in
`experiments/fixed_M_setup/`. The current recipe explicitly sets
`mu_supplied = NULL`; a regression check examines the prepared settings,
not only the settings list in the script. Evaluation of that corrected
model is recorded separately.

With M actually estimated, the default attempt returned optimizer code 0
and gradient 0.000421, but an invalid covariance/non-positive Hessian. The
single source-N retry returned code 1 (singular convergence), gradient
0.00265 and the same uncertainty failure. Both reached objective 651.544.
These outcomes are retained in the main case cache. No successful comparison
is claimed and no additional initialization or model rescue was attempted.

Source facts and remaining gaps are in the
[SAIO review](https://github.com/PaulRegular/saio/blob/757b59c/discovery/source_reviews/noaa_goa_pacific_ocean_perch.md).
