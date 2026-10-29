# NAFO 4T white hake: source screen

## Identity and current assessment

The latest detailed assessment located is the DFO 2024 assessment to 2022,
documented by Science Advisory Report 2024/055 and Research Document 2025/008.
It is a statistical catch-at-age model for southern Gulf of St. Lawrence white
hake (*Urophycis tenuis*) in NAFO Division 4T, covering 1978-2022 and ages
2-10+. This current DFO management unit is distinct from the NOAA Gulf of
Maine/Georges Bank white hake record. The report notes uncertainty about stock
structure in the St. Lawrence Estuary; the assessed 4T unit should not be
silently expanded to include other regions.

## Inputs and outputs identified

The model uses annual fishery catch in tonnes and fishery age compositions.
Research Document 2025/008 Appendix Table A4 reports commercial catch-at-age
in thousands for 1978-2022, with age bins 0-2, 3-12 and 13+. Those bins do not
match the model's ages 2-10+ directly, so the age mapping needs to be checked
before import.

The native survey likelihood uses age-aggregated biomass indices for ages 2-7
from the September Research Vessel (RV) survey and August Mobile Sentinel (MS)
survey, together with proportions-at-age. Appendix Table A5 reports RV
stratified mean catch rates by age, in fish per tow, for ages 2-14 in
1978-2022. This provides a source for a possible age-specific tinyAM index;
that derived series must remain clearly distinguished from the assessment's
aggregate-index likelihood. The RV index uses strata 415-439. MS observations
cover 2003-2019, and the report notes that the MS survey has not operated since
the last assessment.

The model uses age-specific weights from the fishery, RV and MS data. Appendix
Table A6 gives RV mean weight-at-age in kg for 1978-2022, but the fishery and
MS vectors still need to be recovered or verified. Maturity is knife-edge:
ages 2-3 are immature and age 4+ is mature. The model estimates log-scale
random-walk M series for ages 2-3, 4-5 and 6+.

The official advisory report gives 2022 SSB of 6.1 kt. It reports annual SSB,
age-2 recruitment, age-6+ fishing mortality, and age-group natural-mortality
patterns. The detailed report describes MCMC intervals, but numeric annual
uncertainty and N/F-at-age output arrays were not recovered during this screen.

## Screening decision

This distinct current unit is screened eligible for source curation: its
official assessment provides age-disaggregated fishery and survey information,
age-specific biology, explicit mortality assumptions and comparable official
population outputs. No canonical data have been imported. Before import, resolve
the commercial age-bin mapping; recover or verify fishery and MS weight vectors;
extract the native aggregate indices, age compositions, timing and output
series; and record uncertainty limitations. The source-reported RV age-specific
catch rates should not be substituted for the native aggregate index without
labeling that transformation.

The detailed Research Document PDF and its catalogue record are cached locally
with retrieval metadata and checksums. The official advisory-report catalogue
lists Science Advisory Report 2024/055, but its PDF download redirects to a
Government of Canada archived-web notice, so that PDF could not be cached.

## Sources

- [DFO Science Advisory Report 2024/055](https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs70-6/Fs70-6-2024-055-eng.pdf)
- [DFO Research Document 2025/008](https://open-science.canada.ca/server/api/core/bitstreams/28963929-7d4b-47c6-a56d-c0d705be4621/content)
- [Government of Canada catalogue record for Research Document 2025/008](https://publications.gc.ca/site/eng/9.947511/publication.html)
