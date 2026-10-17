# British Columbia Pacific ocean perch

## Production source

**Current comparison status: blocked (2026-10-17).** A complete matching
historical population comparison has not been recovered. A subsequent search
found regional biomass draws for 3CD and 5DE, but not the matching 5ABC series
(see the follow-up below).
The source record stays in the database; it is excluded from comparison
tables and figures. Earlier fitting results below are historical experiments.

The detailed 2023 assessment was published as
[Research Document 2025/004](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41282565.pdf).
This Canadian-managed assessment is a documented catalogue exception, not
three new stock cases. The production coastwide run is **R21v3 / 21.01**,
Stock Synthesis 3.30.20, with 2,000 retained posterior draws. The
[public production directory](https://github.com/pbs-software/pbs-synth/tree/16f1d4368ad7d3e854219ca1c2d9d6353252c861/models/input/2023/POP/CST/Base/R21v3/21.01)
contains inputs, rather than those fitted draws. The report and
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

## Translation — initial expansion (historical)

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

The default-start fit passed numerical checks: optimizer code 0, objective
2199.5445, maximum gradient 0.0014052, positive-definite Hessian and finite
covariance, 12.78 seconds. The dashboard renders. No initialization retry was
allowed or needed by those checks. However, female M is about **2.87e-9** versus
the official median 0.05229, almost −100%. Its log-M coefficient is −19.6677
with SE **1924.9**, a practically unidentified mortality estimate. A numerical
pass does not make this biologically usable. Omission of the informative M
prior and pooling nonmixing areas are important limitations; this is flagged
separately from numerical convergence, without tuning or another fit.

## Source and model review — 2026-10-17

Replaced the omitted M prior with a positive lognormal approximation matching the native natural-scale normal mean 0.06 and SD 0.018. One shared log-M effect is estimated, not fixed to an official value; its density contributes once. Default starts pass (code 0, gradient 0.00164, PD Hessian); M is 0.04053 with effect SE 0.22890, versus official median 0.05229. The M difference is 22.5%, still flagged. Native bounds, exact normal density, posterior integration and sex/area structure differ. Historical N/F/R/SSB remain unavailable or non-equivalent; no agreement is claimed for them. Full M-surface intervals are not automatically propagated. Canonical observations and native 90% limits are unchanged.

See [the complete twelve-case review](../PACIFIC_REVIEW.md#source-and-model-review) for paired median differences, correlations, diagnostics and rejected trials. Canonical numerical inputs and outputs are unchanged.

## Public-output revisit — 2026-10-17

Rechecked the complete, non-truncated repository inventory at the same
`16f1d4368ad7d3e854219ca1c2d9d6353252c861` revision: 1,447 entries and 196
files under the 2023 POP paths. Only the master branch exists; there are no
tags or releases. The matching production directory contains inputs, while
the reporting material contains scripts and terminal/reference-point summaries.
No matching fitted `Report.sso`, parameter/posterior file, R object or numerical
historical N/F/recruitment/SSB table was recovered. The cached research report's
Appendix F presents historical trajectories graphically; its parameter and
reference-point tables do not provide those annual series.

Also searched the complete inventories of `pbs-assess/pacea`,
`pbs-software/pbs-data` and `pbs-software/pbs-tools`. The POP datasets in
PBSdata are biological age/maturity observations. PBStools' `pop.rep` is a
legacy example with a different age/year structure, not this SS3 production
assessment. Neither is a substitute for the 2023 fitted outputs. Inventories,
branch/release responses and supporting files are cached in
`source_cache/dfo_pop_bc/output_search/`, with URLs, revisions, retrieval dates
and checksums in the existing manifest.

The available constant sex-specific M, 2024 boundary biomass and terminal
subarea/reference-point summaries cannot support the intended historical
comparison. Recovering new fitted states from an unverified rerun, mixing
older assessments or treating graph coordinates as exact tables would not
resolve that gap. Reopening this case needs matching numerical historical
outputs from the assessment authors or a public release, plus confirmation of
the native spawning-output definition. No claim is made that the authors do
not hold those outputs. No correspondence was sent.

### Follow-up through the Canadian petrale search — 2026-10-17

The [gftrends repository](https://github.com/pbs-assess/gftrends/tree/7f1c59b8cc3895a3477204ac16da6cf15e2466c6)
provides regional POP biomass posterior objects. The 3CD and 5DE files each
contain 1,965 draws for 1935–2024. Their processing code identifies exports
from the 2023 SS3 assessment, so this is a meaningful recovery of historical
output beyond the original production-directory search. The 5ABC object,
however, contains 1,000 draws for 1940–2017, despite metadata naming the newer
assessment. It must not be combined with the other areas as the 2023 result.

Thus the earlier finding of no numerical historical population outputs was
too broad. Partial regional biomass is available; a complete matching
coastwide comparison remains blocked. N/F/recruitment, the native biomass
definition and joint area-pooling conventions remain unresolved. The files,
processing code and revisions are cached with the
[BC petrale preflight](dfo_petrale_bc.md); no canonical outputs or fitted model
have been changed in this follow-up.

## Suggested alternative: U.S. West Coast petrale sole

The [council's assessment list](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)
provides the [detailed 2023 assessment](https://www.pcouncil.org/documents/status-of-petrale-sole-eopsetta-jordanialong-the-u-s-west-coast-in-2023/)
and its [official model archive](https://www.pcouncil.org/documents/petrale-sole-2023-model-files-zip-file/). The
archive was downloaded and inspected, not inferred from a link label. It
contains `petrale_data.ss`, `petrale_control.ss`, `starter.ss`, `forecast.ss`,
`Report.sso`, `CompReport.sso` and covariance output. `Report.sso` contains
`TIME_SERIES`, `NUMBERS_AT_AGE` and `F_AT_AGE` sections, including historical
recruitment and fishing rates. This is a different stock from Canadian
petrale sole, and is a proposed replacement, not a new completed case.

The report uses data through 2022, female/male states and ages through 40+.
Native spawning output is in **trillions of eggs**; it must remain distinct
from mature biomass. Official N, weights and maturity could support a separate
common-definition biomass comparison, subject to curation. Fleet/season
aggregation, likelihoods, ageing error and uncertainty still need review
before fitting. No new translation or model fit was attempted. No matching
identifier was found in the locally pinned Charbonneau catalogue, so any
addition needs an explicit catalogue exception, not an invented identifier.

The report and 25 MB model archive are cached with checksums under
`source_cache/pfmc_petrale_2023_preflight/`. The public-output availability is
verified, but this is a preflight, not a claim of replication readiness.
