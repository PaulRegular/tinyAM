# 4T winter flounder: preliminary translation

The recipe retains 1973–2016, ages 2–12+, the original catch and trawlable
survey count surfaces, RV weights and maturity. Derived compositions and the
unresolved-scale aggregate biomass index are excluded to avoid double
counting. It approximates the report's Model 2, with deterministic survival,
AR1 recruitment at phi = plogis(−4), static logistic F, two survey curves
changing after 1984, and M estimated in two age groups and five-year blocks.
Source M and q priors are not reproduced. The source explicitly identifies
weak q information without its prior; that limitation matters here.

During setup an intended relative catch-SD multiplier was instead treated as
an absolute supplied SD (1 ordinarily, 20 in 1984), because supplying SDs
removes the formula intercept. That setup converged but is not the intended
source-scaled model. Its fit and dashboard are retained separately under
`results/preliminary/cache/dfo_sgsl_winter_flounder_2017/experiments/sd_setup/`.
The current recipe explicitly uses fixed log-SDs 0.1 and 2, with no free
intercept. Applying aggregate catch SDs independently to age counts is an
approximation, stated in the background and tested on the prepared data.

Only terminal 2016 source population estimates have been recovered. A usable
fit cannot establish trend agreement from that single point, and source
uncertainty is unavailable. The main prose reports Model 2 while Model 7 has
the lowest AIC; the translation does not conceal this run-selection ambiguity.
No aligned official age-state starts are available for a numerical retry.

The source-scaled default attempt did not converge: optimizer code 1
(singular convergence), maximum gradient 0.0574, invalid covariance and
non-positive Hessian. No source-state retry was possible. The current result
and full diagnostics are retained in the main case cache; its earlier setup
dashboard is not presented as a successful result for this model.
