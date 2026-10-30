# Pacific Coast black rockfish: 2023 assessment screen

## Current assessment and stock identity

NOAA Stock SMART's assessment identifier 166727 is associated with the 2023
black rockfish assessment suite. The suite contains four regional model
components: Washington, Oregon, northern California, and central California.
For inventory counting, these are grouped as one assessment case; this does not
assert that they are one biological population, and the regions remain separate
in the source files and any later translation.

The inventory's older NOAA:10402 and NOAA:11302 entries identify broad northern
and southern Pacific Coast units with an assessment year of 2007. Those labels
do not map one-to-one to the four 2023 model areas. They are retained as
ambiguous historical identifiers and are not counted as separate current
cases. The 2023 suite is registered once as NOAA:166727.

## Source evidence

The four 2023 components use Stock Synthesis 3 age-structured models. PFMC's
report and model-file links include the model data and control files, fitted
reports, and numeric age-specific population tables. The downloaded model
archives contain sex-specific numbers-at-age and biomass-at-age tables. The
available age grids are ages 0-40 for Washington and Oregon and ages 0-50 for
northern and central California. The reported population histories begin in
1940 for Washington, 1892 for Oregon, and 1875 for both California components.

Inputs include fleet catches, fishery and survey age/length compositions,
abundance indices, and sex-specific biological settings. Catch-at-age and
age-specific survey indices are not all supplied as direct observation
matrices: age compositions, conditional age-at-length samples, aggregate
catches, and survey indices require documented reconstruction for tinyAM.
The native model files preserve the source likelihood inputs and component-
specific fleet and survey settings. Natural mortality and other biological
assumptions differ among regional models and sexes and must remain component-
specific during translation.

The official output archives provide sex-specific abundance-at-age and
biomass-at-age. Spawning output is reported in millions of eggs, not as
spawning biomass in tonnes. This is an eligible source-curation case because
the native inputs support documented age-based reconstruction and the official
files provide numeric age-specific population outputs; eligibility does not
mean the tinyAM observation likelihood is identical to Stock Synthesis.

## Screening decision

Screen NOAA:166727 as one eligible current assessment case with four preserved
regional model components. Keep both legacy 2007 inventory identifiers
ambiguous and ineligible as separate current cases. Do not pool regional
observations or biological assumptions without an explicit translation
decision. No canonical assessment, input, or output records were added and no
tinyAM translation was attempted.

## Sources

- [NOAA-EDAB Stock SMART changelog](https://noaa-edab.github.io/stocksmart/news/index.html)
- [PFMC Groundfish assessment documents and model files](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/)
- [2023 Black Rockfish STAR panel report](https://www.pcouncil.org/documents/2023/11/stock-assessment-review-star-panel-report-for-black-rockfish.pdf/)
- Washington, Oregon, and California assessment reports and their model-file archives are cached locally. URLs, revisions, retrieval date, byte sizes, and SHA-256 checksums are listed in `source_cache/discovery_noaa_black_rockfish_2023/sources.csv`.
