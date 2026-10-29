# DFO Western Component Pollock

## Unit and assessment run

The Maritimes inventory entry refers to Western Component Pollock, *Pollachius virens*. DFO's science unit is NAFO 4Xopqrs plus the Canadian portions 5Yb and 5Zc. This differs from the broader 4X5YbZc Resource Management spatial area. The component structure is also distinct from the eastern component (4VW and 4Xmn).

Do not coalesce this record with NOAA Stock SMART 10515 (Atlantic pollock, 5YZ). The assessment areas overlap in parts of 5Y/5Z but are not the same one-to-one stock definition.

The 2024 HCR update identifies 2022 as the last full assessment and reports no later interim model update at that time. It extends the survey and catch information used for HCR advice through 2024, but does not refit the population model; its 2024 SSB is an index-based three-year estimate, not an extension of the accepted model's SSB series. DFO's 2026 Research Documents 2026/025 and 2026/029 publish detailed data and model review from the 2022 MSE after the advice update. The accepted population model is ASAP Model J_M4B, with Model I treated as a constant-M sensitivity. Its population series ends in 2020.

## Inputs and model assumptions

- Commercial catch-at-age: ages 2-13, 1982-2020. The 2026 input report identifies corrected age assignments for 2017-2020. Source charts provide the series, but a run-matched machine-readable input array was not found.
- Summer RV trawl survey: abundance-at-age, ages 2-13, survey years from 1984 through 2020 with gaps. The index is relative abundance; the 2022 model uses mid-year timing (0.5).
- Summer acoustic survey: biomass-at-age, ages 2-13, 2012 and 2016-2020. The age-specific series is derived from total biomass using length frequencies by tow; it is not a directly aged abundance index. Model timing is mid-year (0.5).
- Weight-at-age: annual survey and beginning-of-year surfaces for model ages 2-13 through 2020; these are shown in the model/data review, but numeric arrays were not located.
- Maturity-at-age: age-year matrix for ages 1-13 in the population model. The data review describes winter-survey maturity samples and early/recent periods; exact model-input values are not printed as a numeric table.
- Natural mortality: the accepted ASAP model uses the four-block M matrix estimated by VPA Model F. Ages 2-4 are fixed at 0.2 from 1994; ages 5+ have three period blocks (1994-2000, 2001-2010, 2011-2020). ASAP receives the matrix as input; it does not jointly estimate M.
- The accepted model applies equal objective-function weights to commercial catch-at-age and the two survey indices; the report specifies CV 0.2 for catch and CV 0.4 for survey observations. SSB is calculated at start of year with beginning-of-year weights and summed over ages 4+.

## Outputs and source limits

The accepted model report shows numbers-at-age and fishing mortality-at-age, age-2 recruitment, and SSB (ages 4+) through 2020. The time-series surfaces appear in figures rather than downloadable numeric arrays. The reports provide model-fit and retrospective summaries, but the fitted trajectory uncertainty is not available as a full numeric series.

This is an eligible source-curation candidate because official documents verify all required age-structured input/output families and support source-based reconstruction. Any recovered values from figures must be marked as digitized estimates and retain the plot scale and reading uncertainty. No canonical assessment values have been imported in this screening step.

The 2026 research PDFs are available in the official publication viewer, but attempts to retrieve their direct file URLs locally returned the site's archived-page notice. The official 2024 HCR update PDF is cached in the ignored local directory `source_cache/discovery_dfo_western_pollock_2022/`; its `sources.csv` records the URL, retrieval time, byte count, and SHA-256.

## Sources

- [2026 Research Document 2026/025, Data Inputs](https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-5/Fs70-5-2026-025-eng.pdf)
- [2026 Research Document 2026/029, Population Modeling](https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-5/Fs70-5-2026-029-eng.pdf)
- [2024 HCR update, published as Research Document 2025/012](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41291050.pdf)
- [DFO Pollock profile and component boundaries](https://www.dfo-mpo.gc.ca/species-especes/profiles-profils/pollock-goberge-eng.html)
