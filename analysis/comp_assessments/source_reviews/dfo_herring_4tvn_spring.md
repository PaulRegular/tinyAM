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
and a file manifest.

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
state a scale multiplier, so the values are kept unchanged as
`number (index scale; multiplier not stated)`. The previous importer
incorrectly called them weighted sample counts and excluded them from tinyAM.

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
Table 15 values alone on their native scale to avoid mixing the later CSV
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

## Translation limits

The tinyAM recipe sums fixed- and mobile-gear catches into one fleet and uses
published age-specific spring CPUE at a seasonal midpoint of 0.25. It now also
uses the 1994-2023 Table 15 acoustic age series at a late-season midpoint of
0.75. The source model used age composition and aggregate biomass likelihoods;
tinyAM instead treats the reported age series as age-specific lognormal index
observations on the source scale. This preserves the published values without
inventing a multiplier, but it is not the accepted model's likelihood. The
existing q formula is unchanged and therefore shares catchability effects
between the two surveys.

The recipe builds a fit-only stock-weight approximation from the published
gear weights. It combines the two gear values with a geometric mean when both
are available, interpolates missing values within age, and constructs
beginning-year weights from the geometric mean of age `a - 1` in year `t - 1`
and age `a` in year `t`. Age 2 and the first modeled year use same-year
weights because the earlier age-year cells are unavailable. These values are
not added to the canonical source-input table.

The source publishes no numerical April 1 spawning biomass or age-specific M
surface in the detailed output tables. Its tabulated biomass is January 1
biomass, so it is not labeled as SSB. The tinyAM M process is centered on the
source's 0.2 initial-M prior mean and split at ages 2-6 and 7-11+, but it uses
AR1 rather than the source random walk because tinyAM's random-walk initial
state is unpenalized and the source prior cannot be applied as a matching
penalty. MLE M estimates are not used as fixed inputs or comparison values.

The stock does not yet have a Charbonneau identifier in the local crosswalk.
The identifier is left blank until a source crosswalk can confirm it.

## Sources

- [DFO Research Document 2024/058](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41256384.pdf)
- [DFO Research Document 2022/068](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41091589.pdf)
- [DFO Science Advisory Report 2026/028](https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-6/Fs70-6-2026-028-eng.pdf)
