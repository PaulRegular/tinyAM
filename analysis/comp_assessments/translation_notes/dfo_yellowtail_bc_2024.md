# British Columbia Yellowtail Rockfish

The 2024 DFO base model has one coastwide, predominantly trawl fishery,
separate female and male populations, and six fishery-independent biomass
series. SAIO revision `c36cc22` preserves the native total catches, survey
indices, and marginal age compositions as separate observations. The fixed
source biology includes the native sex-specific weights and the female
reproductive weights used for January spawning. [DFO advisory report](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41294804.pdf)

The reference reproduces the pinned B1.R02v2/02.01 point fit with Stock
Synthesis 3.30.22.1. The rerun reports all four series from 1935 to 2024.
Female numbers multiplied by native reproductive weight reproduce reported
SSB. Summing both sexes gives total abundance and age-0 recruitment. Average F
uses instantaneous F weighted by numbers across sexes and ages 9–20. The
native report's exploitation rate uses a different definition. Public SSB
posterior summaries are preserved separately. Over the 1979–2018 window, the
point-fit SSB differs by a median 10.4% from the public posterior median
(`r = 0.965`) and lies inside the public 90% interval in all 40 years.
[Native base inputs](https://github.com/pbs-software/pbs-synth/tree/16f1d4368ad7d3e854219ca1c2d9d6353252c861/models/input/2024/YTR/BC/Base/B1.R02v2/02.01)

The translation fits ages 0–45+ from 1979 to 2018, with catch-at-age in 38
years and 20 paired biomass and age samples from three surveys. It converts
each joint-sex age sample using its own biomass total and source weights, then
sums numbers across sexes. The lognormal likelihood uses those reconstructed
age observations. Native aggregate-only survey observations, ageing error,
and composition weights stay outside the tinyAM likelihood. The translation
uses one pooled M and a 1:1 sex ratio at age for SSB.

## Single-fleet trials

The curated recipe uses a separate logistic curve for each survey and
deterministic cohort survival. The WCVI midpoint was near age 200, with very
large uncertainty on the midpoint and scale. The native QCS and WCVI curves
are approximately flat by age 12. A sensitivity fit used monotone
catchability over ages 1–5, 6–7, 8–9, 10–11 and 12+ for those surveys. A
second sensitivity also freed the 1979 initial abundance at older ages.
The IID trial adds cohort deviations to the original logistic fit with
`update(fit, N_settings = list(process = "iid"), start_par = NULL)`.

| Fit | Median absolute SSB difference | SSB correlation | 2018 SSB difference | Average F correlation | Numerical checks |
|---|---:|---:|---:|---:|---|
| Survey-specific logistics | 7.7% | 0.719 | +19.1% | 0.914 | Passed |
| Synoptic survey plateaus | 10.7% | 0.239 | +5.9% | −0.125 | Passed |
| Synoptic plateaus and free initialization | 16.4% | 0.671 | −2.8% | 0.911 | Passed |
| Survey-specific logistics and IID cohort deviations | 19.1% | 0.785 | +1.1% | 0.916 | Failed |

Scores use all 40 common years and the native female SSB calculation.
Median absolute differences in average F, abundance and recruitment rose
from 12.0%, 32.0% and 37.4% in the logistic baseline to 30.6%, 49.1% and
59.2% with IID deviations. Estimated M fell from 0.094 to 0.067 per year.
The IID fit returned singular convergence and a maximum gradient of 0.082
against a tolerance of 0.01. The logistic fit with deterministic survival
was retained for the curated comparison. All trials fell short of the
0.90 SSB trend target. The curated comparison uses the shorter time series
and pooled sexes. The figures show the remaining trend discrepancy.

## Why SSB tracks the reference less closely than F

The native fit estimates a 13% decline in SSB from 1979 to 1983. The logistic
fit estimates a 5% increase, primarily because of fish already alive in 1979.
Later, the native fit estimates 72.8 million recruits in 1990, compared
with 21.7 million in tinyAM. The 1997 and 2000 cohorts show similar
attenuation. As these cohorts mature, the fitted SSB peaks differ in timing
and magnitude. Substituting the native F trajectory into the tinyAM cohort
equations raises SSB correlation from 0.719 to 0.758. Substituting both
initial abundance and recruitment raises it to 0.866, with a much larger
magnitude difference. These substitutions hold the other fitted inputs
fixed and illustrate sensitivity to each trajectory.

The native proportion female is often below one-half at mature ages. Using
one-half on the native abundance increases median SSB by 24% and preserves
the broad trend (`r = 0.977` with native SSB). That helps explain how the
pooled fit can have fewer fish and still approximate female spawning
biomass. The sparse surveys limit what the fit can learn about the cohorts.
Across ages 31–45, 70–89% of reconstructed survey observations are zeros
excluded by the lognormal likelihood. The few positive records are
consistent with the high fitted old-age catchability. The native
composition likelihood retains information from zero proportions.

The WCVI observations use the correct May survival timing. The expansion
uses midyear weights, while the native survey prediction uses early-season
weights. Recomputing the six reconstructed WCVI samples with the native
early-season weights increases numbers by 2.6–6.0%. Both fits above retain
the same weights so the IID comparison changes only the N process.
