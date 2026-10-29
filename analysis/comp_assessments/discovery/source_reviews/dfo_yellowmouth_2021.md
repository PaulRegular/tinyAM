# Yellowmouth Rockfish

Reviewed 2026-10-29. DFO's current inventory lists Yellowmouth Rockfish as a Pacific stock. The latest detailed assessment located is the coastwide 2021 assessment, published as Research Document 2022/010 and Science Advisory Report 2022/001. It covers PMFC areas 3CD and 5ABCDE and excludes area 4B. No canonical Yellowmouth stock record was found.

## Assessment and source evidence

The 2021 production assessment uses Stock Synthesis 3 (version 3.30.17.01) as an annual, two-sex, age-structured model. It tracks ages 1-60, with age 60 as an accumulator, and combines bottom and midwater trawl catch into one commercial fishery. Commercial catch is estimated annually from 1935 through 2021. Commercial age-frequency samples cover 28 years from 1979-2019; quarterly samples were weighted by sampled-trip catch and expanded with quarterly landings. Four fishery-independent trawl biomass-index series are used: QCS, WCVI, WCHG and GIG. Survey index years span 1967-2020, and survey age compositions cover 25 years. The report downweights survey age compositions to 0.25 because they did not consistently track cohorts; commercial age compositions provide the main recruitment information.

The report describes sex-specific growth and length-weight relationships, a time-invariant female maturity ogive, and the biological data used to derive them. Its composite base case pools five component assessments with fixed scalar M values of 0.04, 0.045, 0.05, 0.055 and 0.06. An estimated-M sensitivity had poor MCMC diagnostics and is not the adopted central run. The assessment reports annual female spawning-biomass and age-0 recruitment trajectories with posterior intervals through 2022, and tabulates terminal B2022 quantiles. The reported fishing-pressure measure is exploitation; numeric N-at-age and instantaneous F-at-age arrays were not recovered in this screen.

The published report contains age-composition figures and detailed input-processing descriptions, but a run-matched native Stock Synthesis bundle or separate numeric age-composition matrices were not located. The quarterly catch expansion is described, so catch-at-age can be reconstructed from catch totals and compositions if the plotted observations can be recovered with adequate precision. The exact within-year survey timing was not verified. These are explicit source-recovery and translation limitations, not reasons to substitute fitted states for observations.

The 2025/26 Groundfish IFMP anticipated an updated assessment in 2024. The 2026/27 Pacific plan summary lists other new advice priorities but no later Yellowmouth assessment. No more recent detailed stock assessment was found in DFO's public CSAS search; the 2021 assessment remains the latest detailed run located for this review.

## Screening decision

Identity is a provisional new coastwide unit, screened eligible for source curation. The official report documents the age-composition inputs, annual catches, four aggregate survey indices, age-specific biological relationships, mortality scenarios, and official SSB/recruitment outputs. Source curation must preserve native aggregate observations, composition weights and the composite fixed-M run structure. Numerical age arrays, fractional survey timing and age-specific fitted N/F outputs remain unavailable until further recovery. No canonical values were imported.

## Sources

- [DFO Research Document 2022/010](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41091449.pdf), detailed 2021 assessment.
- [DFO Science Advisory Report 2022/001](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/4103336x.pdf), accepted advice and reported stock status.
- [DFO CSAS publication record](https://csas-sccs.dfo-mpo.gc.ca/publications/96b903ca-fd67-4275-8236-84df0597b400).
- [2025/26 Groundfish IFMP](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/4127698x.pdf), which anticipated an update in 2024.
- [2026/27 Pacific Groundfish plan summary](https://www.pac.dfo-mpo.gc.ca/fm-gp/mplans/ground-fond-ifmp-pgip-sm-eng.html), checked for later advice priorities.
- [DFO peer-review proceedings, 2022/003](https://publications.gc.ca/collections/collection_2022/mpo-dfo/fs70-4/Fs70-4-2022-003-eng.pdf), indexed in the Government of Canada Publications catalogue; its download endpoint returned an archived-information page during this review.

The detailed report, advisory report, IFMP, plan page and proceedings download response are cached under the gitignored `source_cache/discovery_dfo_yellowmouth/`; `sources.csv` records retrieval dates, URLs, sizes and SHA-256 checksums.
