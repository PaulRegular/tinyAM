# West Coast sablefish: preliminary translation

The recipe fits 2006–2023, ages 0–30+, using the original annual empirical
weight schedules. Marginal joint female/male frequencies are normalized once
over both sexes, excluding inactive source rows and preferring joint samples
over separate unsexed samples. Each retained/discard fleet is reconstructed
separately with its own sex-specific catch weights. Only complete removal
years are summed; 2014 remains wholly missing because pot ages are absent.
No catch or survey observation is derived from fitted population states.

WCGBTS totals are combined with their own age samples and weights. Its
unknown physical multiplier is not invented; q remains relative. Age-bin
tails use the empirical age-30 terminal weight. Combined-sex samples use a
1:1 female/male weight mixture. These are explicit approximations, not
claims of exact SS3 likelihood reconstruction. Source age error, sample
weights, other surveys and the recruitment-deviation index are not reproduced.

The population approximation has deterministic pooled survival, IID
recruitment about a mean, exponential initial abundance, AR1 F around a
logistic age mean, and one estimated M. Native stock–recruitment dependence,
recruitment bias adjustment and the M prior are omitted. Female maturity is
recovered from the supplied mature-weight/mid-season-weight ratio. A female
beginning-weight proxy and a fixed 1:1 sex ratio simplify native two-sex
biology. Source N is summed; source fleet F is added before N-weighting rates
over sexes/terminal ages. Pooled uncertainty is unavailable and not invented.

An initial setup passed numerical checks but had incorrect biological
vectors: scalar-sex selection reduced weight lookups to age 0, recycling that
value over ages. The agreement review exposed this through identical relative
differences for abundance, biomass and mature biomass. Its observations were
correctly reconstructed, but the biological setup and comparison are invalid.
The fit and dashboard are retained only in
`results/preliminary/cache/noaa_sablefish_pacific_coast_2025/experiments/biology_lookup/`.
The lookup now expands the sex label over the age vector, and regression
checks compare every weight and maturity cell directly with the original
source schedules. The corrected fit is evaluated separately.

With correct biology, the default attempt passes numerical checks:
optimizer code 0, objective 1298.976, maximum gradient 0.0000588 and a
positive-definite Hessian. The recipe uses package starting values. Constant
M is about 2.4e−10 per year, with log-M SE about 3023. The estimate is
effectively uninformative after omission of the source prior; numerical
convergence does not make this a defensible estimate of M.

Agreement also needs review. Mean absolute percent differences are 39.4%
for abundance, 225.7% for N-at-age, 49.8% for F-at-age, 102.4% for recruitment
and 29.1% for approximate common-definition SSB. Their respective trend
correlations are 0.990, 0.499, 0.842, 0.920 and 0.835. SSB's terminal difference
is +3.6%, but its mean difference and correlation still trigger review.
Estimated M differs by effectively −100%; constant source M has no meaningful
trend correlation. The corrected dashboard and fit are in the main case
cache, clearly separate from the invalid biological setup.

Source facts and native-file provenance remain in the
[SAIO review](https://github.com/PaulRegular/saio/blob/757b59c/discovery/source_reviews/noaa_west_coast_sablefish_2025.md).

## Curated expansion review on SAIO `cecb45ad`

The estimated-M recipe reproduced the near-zero M estimate. It passed numerical checks and gave a median absolute SSB difference of 29.4% over 18 common years, with trend correlation 0.835 and a terminal difference of +3.6%.

Fixing M at the source prior median of 0.072 per year gave an implausibly large population and virtually zero fishing mortality. Both default starts and a restart from the tinyAM estimated-M fit led to the same pattern. Relative to the reference, the median SSB ratio exceeded 300 million in both fixed-M fits. Optimizer and Hessian checks passed, and several estimated parameters were highly correlated. The two starts indicate that the scale problem persists in this pooled catch and survey specification. Sablefish remains preliminary.

## Survey scale sensitivity

The official WCGBTS selectivity has a sharp peak at age 1 and a nearly flat plateau from age 2. With M held at 0.072 per year, estimating separate q for ages 0, 1 and 2+ still gave almost zero F and an implausibly large stock. A diagnostic then divided the reconstructed index by the official fitted q of 1.0497 times the native relative selectivity (approximately 0.077, 0.995 and 0.801). This assumes the biomass index and sampled age composition can be combined into an absolute numbers index. The catch inputs, biology, and F and M settings were identical across trials. [Official assessment, section 3.4.1](https://www.pcouncil.org/documents/2026/02/status-of-the-sablefish-off-the-u-s-west-coast-in-2025.pdf/)

With unit q on that adjusted index, the fit passed numerical checks and gave a median SSB difference of 13.3% over 18 common years, with trend correlation 0.895 and a terminal difference of −17.8%. Compared directly with the published female SSB, the median difference was 15.9% and the correlation was 0.900. Freeing the overall q multiplier on the adjusted index brought q to about `1.4e-9`, median Fbar to about `8.7e-11` per year and the median SSB ratio to about 611 million. Correcting the relative age shape alone leaves the scale trade-off unresolved. The registered recipe stays preliminary because the plausible scale in the diagnostic depended on official fitted catchability. Extending the 2006–2023 window through the full survey and biology overlap would add four survey years with incomplete catch reconstructions.
