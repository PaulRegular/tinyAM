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
science advisory report, and a file manifest.

## Inputs and assumptions recovered

The 2024 support document reports fixed- and mobile-gear spring catch-at-age
and weight-at-age, spring fixed-gear CPUE by age, acoustic weighted-trawl
sample counts, and numerical January 1 abundance, biomass-at-age, and
fishing-mortality-at-age estimates. The reported age range is 2-11+, with age
11 as the plus group. The 2023 mobile-gear catch-at-age is represented by zeros
because the annual landings table reports no mobile-gear spring catch that
year; this expansion is labeled as a reconstruction in the input records.

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
published age-specific spring CPUE at a seasonal midpoint of 0.25. That CPUE
series is not the same observation likelihood as the source model's aggregate
index plus age composition. The acoustic table contains weighted-trawl sample
counts used to form age composition, not an abundance index compatible with
tinyAM, so it is not fitted.

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
