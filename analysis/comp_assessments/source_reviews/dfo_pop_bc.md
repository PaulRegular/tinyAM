# British Columbia Pacific ocean perch

## Production source

The detailed 2023 assessment was published as Research Document 2025/004.
This Canadian-managed assessment is a documented catalogue exception, not
three new stock cases. The production coastwide run is **R21v3 / 21.01**,
Stock Synthesis 3.30.20, with 2,000 retained posterior draws. The report and
native files are cached under `source_cache/dfo_pop_bc/`; the pbs-synth
revision is `16f1d4368ad7d3e854219ca1c2d9d6353252c861`.
Use the validated 310-page Waves PDF. Two earlier DFO/publications downloads
were HTML responses despite their PDF filenames and are not assessment evidence.

Preserve areas **5ABC, 3CD and 5DE**, with no movement between them, female/male
states, three trawl catches and six survey series. Historical years are
1935–2023 and ages 0–60+; original observed age bins are 1–60+. The 2023 catch
input repeats 2022 provisionally. It is not a final 2023 observed landing.

## Inputs and fixed biology

Native catches and survey indices are tonnes. Joint-sex age compositions are
number proportions: female and male cells sum jointly to one. Original supplied
N, partition, age-error keys and sample years are retained. Francis multipliers
are 2.7923, 3.1274, 2.8821, 0.5744, 1.0766, 1.2071, 0.7348 and 0.4874 for
fleets 1–8. These are not raw otolith counts. Native survey month 6 is time
5/12, not 0.5. All length-selectivity patterns are zero; age selectivity has
double-normal and linked shapes with sex/area differences.

Growth, length CV and weight–length parameters are **fixed**, regardless of
the generic control-file header saying growth is estimated. No empirical
weight matrix was supplied. Recover beginning/midyear population weights by
evaluating these original fixed relationships with the matching official
SS3 3.30.20 release, `-maxfn 0 -nohess -nox`. The small source-only recovery
script records provenance and checksums. The resulting SS parameter file
retains initial values; no fitted population states are used. Integrating the
length distribution avoids substituting weight at mean length for mean weight.
Weights are time-invariant, distinct by sex/timing. All-fleet midyear weights
are identical because there is no length selectivity. Only these fixed biology
calculations are imported; the generated N, F, M, objective and SSB are **not
official estimates**.

The original female maturity vector is zero at ages 0–4, then 0.118685573,
0.174868158, 0.247826408, 0.337837982, 0.442989478, 0.558730517,
0.677852837, 0.791029151, 0.88791952, 0.95869125, 0.995652911, then one.
Female and male M are estimated constants shared across areas, with normal
prior mean 0.06, SD 0.018 and bounds 0.02–0.2. Recruitment uses estimated
Beverton–Holt steepness, fixed log SD **0.9**, annual deviations/bias adjustment
and time-varying allocation among the three areas.

## Official outputs and unresolved definition

Table F.4 gives female M median 0.05229 (0.04365–0.06146) and male M 0.05939
(0.05050–0.06902). These are **90% credible limits**, not SEs or 95% intervals.
Table F.5 reports January 2024 female spawning biomass 61,965 t with 90%
limits 44,390–90,825 t. Native subarea CSV medians and limits are retained
separately; marginal component medians/limits are not summed into coastwide
quantiles. January 2024 is a projection boundary, not historical 2023.
Historical posterior N/F/recruitment/SSB trajectories were not recovered
numerically; figures and initial-state evaluation output do not fill that gap.

A definition issue needs confirmation from final native posterior files or
the assessment authors. The public R21v3 control specifies fecundity option 2
with coefficients 1/0 and first mature age 9.49256964. Matching SS3 reads the
latter as integer 9 and calculates reproductive output separately from mature
biomass. Its code explicitly stores both `SSB_yr` and `SSB_B_yr`. The public
reporting code also handles `SpawnBio`, while the paper tables label B in
tonnes. Without the final posterior extraction, do not silently equate those
quantities or "repair" the official assessment. Published values are retained
as reported and SSB comparison is marked non-equivalent.

## Translation

The recipe uses female ages 5–60+, 2003–2022. Each source area is converted
from biomass and joint-sex proportions using both sexes and all observed ages
in the weight denominator. Retain females after conversion, then pool the
three catch areas only when **all three** have compositions. A missing area
does not become zero catch or a partial coastwide observation. Modern synoptic
surveys retain their identity and time 5/12. A pooled female population and
survey-specific q approximate fixed area shares; original spatial recruitment
variation cannot be reproduced by q alone.

Use deterministic N survival, free initial ages, IID constant-mean recruitment
with log SD 0.9, independent F random walks and estimated constant female M
without its source prior/bounds. Survey-specific four-df natural cubic age
curves allow domes without inventing a new model feature. Original female
beginning-year weights/maturity define tinyAM January SSB, without another
half-sex factor. Source zero age proportions remain zero and are treated as
missing by tinyAM's lognormal observation likelihood.

Normal package starts are tried first. A failed numerical check permits one
retry using only the official female M median; unavailable N/F starts retain
defaults. Data, formulas and assumptions are identical between attempts.
M is the only recovered historical comparable core output. Native 90% limits
remain canonical and are omitted from the common 95%-interval dashboard path.
Unavailable series and unresolved SSB are not agreement passes.
