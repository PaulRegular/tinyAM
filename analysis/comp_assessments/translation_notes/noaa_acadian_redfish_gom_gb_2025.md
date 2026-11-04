# Acadian redfish: preliminary 2025 translation

The recipe uses SAIO preliminary revision
`6f08c06c49e4b73591427a1495bfb4568caf49f0`, restricting the fit to 1969–2024
when catch age samples are available. Catch numbers are reconstructed only in
sampled years from biomass, number proportions and the original catch-weight
schedule. Unsampled age shares remain missing. Survey totals and their own
age shares produce relative age indices without an invented physical
multiplier. These derived ages share totals; independent lognormal age errors
approximate the source aggregate/composition likelihoods.

The recipe follows deterministic survival, AR1 recruitment, fixed M = 0.05,
static fishery age effects and annual F effects. Fishery ages 9+, fall survey
ages 4+ and spring survey ages 7+ share effects. Initial exponential abundance
uses first-year mortality rather than the source's separate equilibrium F.
SSB comparisons are labelled approximate because January and spawning weight
schedules differ in the official model.

Default and one source-state initialization retry both returned optimizer
code 0, but neither produced a valid fixed-effect covariance. Maximum
absolute gradients were 0.000988 and 0.000666. The case remains numerically
non-converged. Both attempts, warnings and the fitted object are retained in
`results/preliminary/cache/noaa_acadian_redfish_gom_gb_2025/`; no successful
comparison or dashboard is claimed.

Source equations, data and gaps remain in the
[SAIO review](https://github.com/PaulRegular/saio/blob/e2ef4ae682bf4e80fd4b85edd0c19f1c62780c08/discovery/source_reviews/noaa_acadian_redfish_2025.md).
