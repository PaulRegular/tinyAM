# NOAA West Coast Chinook candidates (2026)

## Scope and identity

This review covers 12 NOAA Stock SMART rows: Washington Coast Hoko Summer/Fall
(10483) and 11 Puget Sound/eastern Strait of Juan de Fuca management units
(10484–10494). Each has a distinct NOAA stock ID and matching scientific and
common names in the cached 2026 Stock SMART summary. The canonical database has
no *Oncorhynchus tshawytscha* stock, so these are recorded as distinct new
candidate units (
oaa_chinook_<Stock SMART ID>`). This is an inventory
crosswalk only; no canonical numerical records were added.

The cached Stock SMART summary marks all 12 latest entries as 2026 operational
assessments using “Cohort Analysis.” Composition Input Data is level 4 for all
12. The NOAA data dictionary says this category indicates a time series of age
composition data sufficient for assessment methods that use fishery age
compositions. Model Category is 4 (VPA) for 11 units, but 2 (index-based) for
Skokomish Summer/Fall (10485). This difference is material: the metadata are
strong evidence that age-composition information exists, but do not establish
that all 12 use the same assessment structure or expose the same annual input
series. The catch, abundance, and biological input-data scores also vary among
the 12. Stock SMART is summary metadata and directs users to regional assessment
records for details.

## Source screening

The latest reviewed public documents are the 2026 Preseason Report I and the
Review of 2025 Ocean Salmon Fisheries. Preseason Report I gives run-size and
adult-return forecasts. For Puget Sound it says forecasts are described in
annual management reports by management unit, not by individual species, and
that forecasts are generally dominated by age-4 adults. The reviewed reports
do not expose the annual catch numbers-at-age vectors, an age-specific marine
survey index, stock weights-at-age, maturity-at-age, or age-specific natural
mortality needed to verify a tinyAM-ready input set. The reports also do not
provide a complete N-at-age/F-at-age/M-at-age output surface for these units.

Stock SMART does publish aggregate indicators for some units. In its 2026
summary, an estimated exploitation rate is reported for 9 of the 12 candidates
and a spawning-stock point estimate for 6. The abundance basis varies (natural,
hatchery, or both), and these point estimates do not establish a consistent
population output series across the 12 units. The candidate register records
the available point-estimate units and years row by row. The age-composition
metadata make it inappropriate to mark the units ineligible solely because the
summary reports aggregates; instead, eligibility remains **unresolved** until
the annual regional management reports and underlying assessment inputs are
checked.

| NOAA ID | Management unit | Model category | Composition input | Decision |
|---:|---|---:|---:|---|
| 10483 | Washington Coast: Hoko Summer/Fall | 4 | 4 | Unresolved |
| 10484 | Puget Sound: Eastern Strait of Juan de Fuca Summer/Fall | 4 | 4 | Unresolved |
| 10485 | Puget Sound: Skokomish Summer/Fall | 2 | 4 | Unresolved |
| 10486 | Puget Sound: Nooksack Spring Early | 4 | 4 | Unresolved |
| 10487 | Puget Sound: Skagit Summer/Fall | 4 | 4 | Unresolved |
| 10488 | Puget Sound: Skagit Spring | 4 | 4 | Unresolved |
| 10489 | Puget Sound: Stillaguamish Summer/Fall | 4 | 4 | Unresolved |
| 10490 | Puget Sound: Snohomish Summer/Fall | 4 | 4 | Unresolved |
| 10491 | Puget Sound: Cedar River Summer/Fall | 4 | 4 | Unresolved |
| 10492 | Puget Sound: White River Spring | 4 | 4 | Unresolved |
| 10493 | Puget Sound: Green River Summer/Fall | 4 | 4 | Unresolved |
| 10494 | Puget Sound: Nisqually River Summer/Fall | 4 | 4 | Unresolved |

Source and output scores are discovery evidence, not reconstructed likelihood
vectors. No age data, model states, or numeric observations were copied into
the canonical database, and no tinyAM translation was attempted.

## Sources and cache

- [NOAA Stock SMART](https://apps-st.fisheries.noaa.gov/stocksmart/), with a
  [stock-specific view for Hoko](https://apps-st.fisheries.noaa.gov/stocksmart/?app=browse_by_stock&stockid=10483).
- [NOAA Stock SMART Data Dictionary](https://apps-st.fisheries.noaa.gov/stocksmart/StockSMART_DataDictionary.pdf), for definitions of composition-input and model-category scores.
- [PFMC 2026 Preseason Report I](https://www.pcouncil.org/documents/2026/03/preseason-report-1-stock-abundance-analysis-and-environmental-assessment-part-1-for-2026-ocean-salmon-fishery-regulations.pdf/).
- [PFMC Review of 2025 Ocean Salmon Fisheries](https://www.pcouncil.org/documents/review-of-2025-ocean-salmon-fisheries/).

The three PDFs are cached under the gitignored
`source_cache/discovery_noaa_puget_sound_chinook_2026/`; their retrieval date,
URLs, byte sizes, and SHA-256 checksums are in `sources.csv`. The cached NOAA
Stock SMART summary and its pinned repository revision are documented in
`source_cache/discovery_noaa_stocksmart/sources.csv`.
