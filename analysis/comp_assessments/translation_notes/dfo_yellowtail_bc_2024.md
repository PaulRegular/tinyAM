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

All three fits passed tinyAM's numerical checks. The first used a logistic
curve for each survey. The WCVI midpoint was near age 200, with very large
uncertainty on the curve and scale. The native QCS and WCVI curves are
approximately flat by age 12. The preliminary recipe therefore estimates
monotone catchability over ages 1–5, 6–7, 8–9, 10–11, and 12+ for those
surveys. The final sensitivity kept those groups and freed the 1979 initial
abundance at older ages. It estimated one coefficient per older age and
retained the annual recruitment model.

| Fit | SSB median absolute difference | SSB trend correlation | 2018 SSB difference | Average F trend correlation |
|---|---:|---:|---:|---:|
| Survey-specific logistics | 7.7% | 0.719 | +19.1% | 0.914 |
| Synoptic survey plateaus | 10.7% | 0.239 | +5.9% | −0.125 |
| Synoptic plateaus and free initialization | 16.4% | 0.671 | −2.8% | 0.911 |

Scores use all 40 common years and the native female SSB calculation. The
logistic fit had median absolute differences of 12.0% in average F, 32.0% in
total abundance, and 37.4% in recruitment. The survey plateaus resolve the
poorly identified logistic shape but weaken trend agreement. Free
initialization improves the trend correlation relative to that fit, with
larger median differences in all four quantities and one fitted parameter
at a bound. The SSB trend correlation fell short of the 0.90 curation
target in every trial. The recipe remains preliminary.
