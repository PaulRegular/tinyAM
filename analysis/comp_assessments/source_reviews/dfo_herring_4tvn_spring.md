# Southern Gulf spring-spawning herring

## Assessment record

The canonical database represents the accepted 2024 assessment to 2023, the
most recent assessment with detailed numerical tables available. DFO Research
Document 2024/058 is a data and support document; it does not provide a full
description of the model run. The 2026 accepted assessment to 2025 is recorded
separately as summary-only because its detailed research document is still in
preparation. No 2026-run data or estimates are mixed into the 2024 record.

The source files are cached locally under
`source_cache/dfo_sgsl_herring_spring_2023/`. The cache includes the 2024
support document, the 2022 methods report used for method context, the 2026
science advisory report, the DFO acoustic biomass CSV and data dictionary,
the 2016/060 acoustic table used to clarify units, and a file manifest.

## Inputs and assumptions recovered

The 2024 support document reports fixed- and mobile-gear spring catch-at-age
and weight-at-age, spring fixed-gear CPUE by age, and spring-spawner acoustic
age-index values for ages 2-10 in 1994-2023. It also reports numerical January
1 abundance, biomass-at-age, and fishing-mortality-at-age estimates. The
fishery age range is 2-11+, with age 11 as the plus group. The 2023
mobile-gear catch-at-age is represented by zeros because the annual landings
table reports no mobile-gear spring catch that year; this expansion is labeled
as a reconstruction in the input records.

Table 15 is not a table of raw fish processed from trawl samples. The report
describes the acoustic trawl samples as weighted by acoustic density and uses
them to allocate the survey signal by spawning component and age. Figure 19
labels the values as abundance-at-age in numbers, and the 2022/068 methods call
the series an age-disaggregated acoustic abundance index. The table does not
state a scale multiplier. DFO Research Document 2016/060, Table 16 (p. 36),
explicitly labels the same historical acoustic series as thousands of fish.
For example, its 1994 age-4 value is 100,087, compared with 100,062 in the
2024 table; its 2006 age-4 value is 24,601 in both reports. These small
historical revisions do not imply a change of units. The canonical database
therefore retains the 2024 values in `thousand fish`, with the unit source and
resolution recorded. The shared converter multiplies them by 1,000 to obtain
fish. No manual multiplier is needed in the stock recipe.

The biomass provides an independent scale check: the 2023 Table 15 spring
values sum to 55,619 thousand fish. The later official CSV's 7,358 t of
spring biomass implies a mean weight of about 0.132 kg per fish on that scale,
whereas treating the table as individual fish would imply about 132 kg per
fish. This is a scale cross-check only; the later CSV is not added to the
2024 fit. The previous importer incorrectly called the observations weighted
sample counts and later treated the multiplier as unresolved. Both descriptions
have now been corrected using the historical unit documentation.

The accepted 2024 model's acoustic likelihood is only partly recoverable. The
2022/068 methods describe a multivariate-logistic age-composition likelihood,
with small age proportions grouped with adjacent ages until they exceed 0.01,
and a separate lognormal acoustic biomass likelihood. Spring acoustic biomass
uses ages 4-8 and its likelihood weight is 3; CPUE biomass weight is 1. The
2024/058 support document says the 2020 SCA model was updated, but does not
repeat the 2022 likelihood details or confirm that the weights and age groups
were unchanged. The tinyAM translation uses the published age-index values at
ages 2-10 as direct age-specific observations; this differs from the source's
composition-plus-aggregate-biomass likelihood, whose biomass component uses
ages 4-8. This uncertainty remains recorded rather than presented as a
verified 2024 setting.

The Open Government biomass CSV and dictionary are cached with provenance.
The dictionary states that spring acoustic biomass is in tonnes, with values
for 1994-2025. The 1994-2023 rows are kept on the separate 2026 summary record;
2024-2025 are not attached to the 2024 assessment. The 2024 translation uses
Table 15 values alone, converted from thousand fish to fish, to avoid mixing the later CSV
release into the accepted 2024 record. The latest file is useful as a
cross-check, but it is not confirmed to reproduce the 2024 model's acoustic
biomass inputs: its 2022 total is 27,209 t versus 27,268.7 t in the 2024
support report, while its 2023 total (19,363 t) agrees with the reported
19,363.8 t after rounding. The 2021 total in the current file (23,147.68 t)
differs from the 2022/068 methods report's 37,953.1 t. The reason for that
revision and the exact age-aggregated biomass series used by the 2024 run
remain unresolved.

The most recent detailed methods description is DFO Research Document
2022/068. It describes estimated age-2 recruitment, logistic selectivity in
three periods, separate log-M random walks for ages 2-6 and 7-11+, a 0.2
initial-M prior mean, fixed M-increment SD 0.075, and a random walk in the
spring CPUE catchability. These settings are useful context but may not capture
changes made for the 2022-2023 assessment.

The maturity schedule is knife-edge between ages 3 and 4: ages 2-3 are
immature, and ages 4-11+ are mature. The methods report says beginning-year
weights combine fixed- and mobile-gear weights and then use a geometric mean
across adjacent age-year cells. It does not provide the complete fitted
stock-weight surface or specify exactly how the two gear series are combined.

## Derived source summaries

The database now includes total January 1 abundance and biomass across ages
2-11+, calculated by summing the published age-specific MLEs. It also includes
January 1 mature biomass at age and its sum, calculated from published biomass
and the knife-edge maturity schedule. This uses the source biomass surface,
not the approximate weights chosen for the tinyAM fit. Native reported 4+
totals are retained separately. Rounding in the age tables can make a derived
sum differ slightly from a reported aggregate.

These rows are labeled `derived_source_output`, with formulas and sources in
their notes and no SEs or intervals. The mature-biomass sum is stored under
`SSB` for the common January 1 comparison used by tinyAM. It is not the
accepted assessment's April 1 SSB: applying pre-spawning mortality would
require the source's estimated M surface, which is unavailable numerically.

## Translation limits

The tinyAM recipe sums fixed- and mobile-gear catches into one fleet and uses
published age-specific spring CPUE at a seasonal midpoint of 0.25. It now also
uses the 1994-2023 Table 15 acoustic age series at a late-season midpoint of
0.75. The source model used age composition and aggregate biomass likelihoods;
tinyAM instead treats the reported age series as age-specific lognormal index
observations in fish after resolving the table units. This is not the accepted
model's likelihood. The two surveys have separate non-decreasing age-q curves
and separate observation SDs. CPUE catchability varies in three blocks (1990-1999, 2000-2009 and
2010-2021), rather than the source's annual random walk. Acoustic catchability
is constant through time. A log link allows q to act as an index scaling
coefficient without imposing an upper bound. The user's original multiplication
of acoustic values by 1,000 was correct; the correction is now made through
the canonical unit rather than a manual fit-only adjustment. Resolving the unit
raises fitted acoustic q by a factor of 1,000 (maximum about 0.612) while
leaving population states essentially unchanged under the log link. A plausible
q supports the scale check but is not, by itself, evidence for a unit conversion.

The recipe builds a fit-only stock-weight approximation from the published
gear weights. It combines the two gear values with an arithmetic mean weighted
by catch numbers at that age; if both gear catches are zero, it uses the mean
of available positive weights. It interpolates missing values within age and constructs
beginning-year weights from the geometric mean of age `a - 1` in year `t - 1`
and age `a` in year `t`. Age 2 and the first modeled year use same-year
weights because the earlier age-year cells are unavailable. These values are
not added to the canonical source-input table. The gear-weighting rule is an
inference, not a fully specified source instruction. Its median absolute
percentage difference from the weights implied by the reported biomass/N
surface is 0.073%, compared with 2.63% for the former equal geometric gear
mean. Source biomass/N is used only to check this rule, never as an input to
the fit. Boundary cells and rounding account for some remaining differences.

The source publishes no numerical April 1 spawning biomass or age-specific M
surface in the detailed output tables. Its tabulated biomass is January 1
biomass; the derived mature-biomass comparison is explicitly labeled January
1 rather than accepted April 1 SSB. The tinyAM M process is centered on the
source's 0.2 initial-M prior mean and split at ages 2-6 and 7-11+, but it uses
AR1 rather than the source random walk because tinyAM's random-walk initial
state is unpenalized and the source prior cannot be applied as a matching
penalty. MLE M estimates are not used as fixed inputs or comparison values.

The translation retains the revised model's exponential initial abundance,
IID N process for older ages, AR1 F, and two AR1 M blocks. F mean age blocks
are 2-3, 4-5, 6-7 and 8-11; Fbar uses ages 6-8. The source model does not
include the same older-age N process and estimates initial cohorts differently.
Combined catch-at-age is converted from thousand fish to fish. Both catch and
index settings exclude zeros and missing values instead of estimating filled
observations. Source N, F and recruitment provide starting values only;
catchability starts use the corresponding, sorted observation rows.

## Refinement trials (October 2026)

Each trial changed one assumption from the revised model. The previously
unresolved-scale, log-q version was used for the trials below; converged states from that fit
were also used as starting values where needed. The source-like choices are
approximations using existing tinyAM settings, not reconstructions of the
accepted likelihood.

| Trial | Converged / positive Hessian | Median absolute % difference: total N | January 1 mature biomass | Decision |
|---|---|---:|---:|---|
| Revised input model | Yes / yes | 42.2 | 17.0 | Baseline |
| Previously unresolved acoustic scale and log q | Yes / yes | 42.2 | 17.0 | Log link retained; acoustic units subsequently corrected to thousand fish |
| Restrict CPUE to ages 4-10 and acoustic to 4-8 | Yes / yes | 32.8 | 22.6 | Do not retain: mixed agreement; keep the more detailed published indices |
| M random walk | Yes / yes | 49.7 | 27.5 | Do not retain: worse agreement; no matching initial-M prior |
| AR1 rather than IID older-age N process | Yes / yes | 53.6 | 27.1 | Do not retain: closer terminal biomass but worse historical agreement |
| Separate survey baselines with a three-degree-of-freedom natural spline for CPUE year | Yes / yes | 50.5 | 31.8 | Do not retain: worse agreement than CPUE blocks |
| No older-age N process | No / no | — | — | Do not retain: false convergence from source and fitted-state starts |
| Free initial abundance | No / no | — | — | Do not retain: false convergence/evaluation limit |
| Random initial abundance | No / no | — | — | Do not retain: false convergence |
| Final, with catch-number-weighted gear weights | Yes / yes | 42.2 | 19.0 | Retain: better supported biological inputs; overall output agreement is mixed |

The final fit has optimizer code 0, objective 1173.766 and a positive-definite
reported Hessian. Its raw maximum gradient is about 18.24 at a zero monotone-q
increment; a positive derivative there satisfies the lower-bound optimality
condition. The package's convergence check projects that component to zero.
After correcting the acoustic unit, the maximum projected gradient is
0.0000538, below the 0.01 tolerance. Acoustic q ranges from 0.220 to 0.612.
The rerun changes N, F, M and derived population estimates by less than
0.0011% relative to the preceding fit, while acoustic predictions scale by
1,000 to numerical tolerance. The dashboard and comparison summaries use
this corrected database revision.
The log-q fit reproduces the original population-state estimates to
numerical tolerance. Gear weights enter derived biomass, not the number-based
observation likelihood, so their correction does not change fitted N or F.

The final mean absolute percentage difference in January 1 mature biomass is
24.9%, compared with 25.3% before the gear-weight correction. Its terminal-year
difference improves from -36.5% to -33.8%, while its median difference worsens
from 17.0% to 19.0%. Total biomass's terminal difference improves from -57.0%
to -55.4%. This is a modest biological correction, not a general improvement
in replication: terminal abundance remains 64.5% lower and recruitment 88.5%
lower than the source MLEs. Matching source aggregate/composition likelihoods,
acoustic weighting, M priors and fixed process SDs remains outside this
translation. Failed trial fits are not evidence that these alternatives are
scientifically inappropriate; only the tested starts/settings are ruled out
for this recipe.

The AR1 N sensitivity brings terminal January 1 mature biomass to 5.8% below
the source, but its historical median difference is 27.1% and mean difference
35.0% (versus 17.0% and 25.3% for the baseline with the same weights). Its
median total-abundance difference also rises to 53.6%. It is therefore not
selected merely for a closer terminal estimate; this trade-off warrants
inspection if a correlated N process is explored further.

The stock does not yet have a Charbonneau identifier in the local crosswalk.
The identifier is left blank until a source crosswalk can confirm it.

## Sources

- [DFO Research Document 2024/058](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41256384.pdf)
- [DFO Research Document 2022/068](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41091589.pdf)
- [DFO Research Document 2016/060, Table 16: acoustic units](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/365860.pdf)
- [DFO Science Advisory Report 2026/028](https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-6/Fs70-6-2026-028-eng.pdf)
