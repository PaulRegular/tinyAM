# NOAA Gulf of Alaska southern rock sole screen

## Identity and latest official run

NOAA Stock SMART entry 10040 is Gulf of Alaska rock sole, scientific name *Lepidopsetta bilineata* (southern rock sole). The database has no matching Gulf of Alaska record. The existing `afsc_northern_rock_sole_ebs` record is the Bering Sea/Eastern Bering Sea northern rock sole (*L. polyxystra*), so it is not the same unit.

The latest full assessment located is the 2021 GOA northern and southern rock sole SAFE chapter. Southern rock sole's preferred assessment was Model 21.1, which represented central and western GOA separately. The 2023 report provides harvest projections based on the 2021 fit. The 2025 SAFE chapter explicitly reports that no assessment was conducted that year because of the federal shutdown and says the 2021 fit remains the latest full assessment. Thus no later model run should be inferred from the projection or catch-report updates.

## Inputs and outputs

The 2021 assessment covers ages 0-30 and uses total annual rock sole catch split 50:50 between northern and southern rock sole. The report states that fishery input compositions are lengths from 1997-2021, not catch-at-age. Survey inputs include total biomass with uncertainty, length compositions, and conditional age-at-length data. Southern rock sole survey CAAL samples span 1996-2019; the 2021 survey age samples were not available in the assessment report's input table. The assessment therefore is age-structured internally, but its documented fishery observations do not include the catch-at-age series required for the current tinyAM translation screen. Survey CAAL does not replace fishery catch-at-age.

The assessment uses fixed age-based maturity, a length-weight relationship, female M fixed at 0.2, and estimated male M. The preferred 2021 SRS model reports male M of 0.253 in central GOA and 0.271 in western GOA. Exact numerical maturity and weight-at-age vectors were not recovered from these sources. The report contains age-structured model methods, survey biomass, biomass and SSB summaries, mortality and management outputs, but not annual full N-at-age and F-at-age output tables.

The 2023 projection report states the last full assessment was in 2021 and gives projected aggregate biomass, SSB, M, and ABC/OFL by area. The 2025 catch-report SAFE chapter states no assessment was conducted that year and carries forward the 2021 assessment for harvest specifications. Neither projection/reporting update supplies the missing fishery catch-at-age input.

## Screening decision

This is a distinct, current NOAA stock unit and is provisionally `new` relative to the canonical database. It is screened as ineligible for a tinyAM age-structured translation based on the recovered official inputs because the fishery provides total catch and length compositions rather than catch-at-age. The documented input sources also do not provide the required numeric age-specific biology or annual N/F-at-age output series. It remains a valid stock assessment record for future expansion if run-matched age data or an approved age-observation translation becomes available.

## Sources and cached files

- NOAA, 2021, *Assessment of Northern and Southern Rock Sole Stocks in the Gulf of Alaska*, preferred southern rock sole Model 21.1: <https://apps-afsc.fisheries.noaa.gov/Plan_Team/2021/GOAnsrocksole.pdf>.
- NOAA, 2023, *Assessment of the Northern and Southern Rock Sole Stock in the Gulf of Alaska*, off-cycle harvest projection: <https://apps-afsc.fisheries.noaa.gov/Plan_Team/2023/GOAnsrocksole.pdf>.
- NPFMC, December 2025, *GOA Shallow Water Flatfish SAFE chapter*, one-page catch report stating no assessment was conducted in 2025: <https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=4.+Shallow-Water+Flatfish.pdf&p=969301d9-957a-48ed-8d78-094a02ffce04.pdf>.
- NOAA, *2023 Assessment of Northern and Southern Rock Sole Stocks in the Gulf of Alaska*, official assessment page and report link: <https://www.fisheries.noaa.gov/resource/data/2023-assessment-northern-and-southern-rock-sole-stocks-gulf-alaska>.

The 2021 assessment and 2025 SAFE catch-report PDFs are cached under `source_cache/discovery_noaa_goa_rock_sole_2026/`; `sources.csv` records retrieval date, source URL, SHA-256, and role.
