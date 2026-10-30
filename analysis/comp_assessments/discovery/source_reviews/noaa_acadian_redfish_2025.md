# NOAA Acadian redfish source screen

**Screened:** 2026-10-30  
**Candidate:** NOAA Stock SMART 10455, Acadian redfish, Gulf of Maine / Georges Bank  
**Scope:** Discovery screen only; no canonical stock, assessment, input, assumption, or output rows are added.

## Current assessment and source run

The latest inventory identifies assessment 14070 as the 2025 Acadian redfish Management Track assessment, with data through 2024 and full acceptance. The NOAA SASINF archive contains the run-specific report, `BaseWHAMCode.R`, `BaseWHAMInput.DAT`, and the fitted figures-and-tables bundle. The code identifies Run 5 and WHAM 2.1.0.9002 at commit `47506172e9a55376137501608e645ddae8a574e9`. The final September 2025 peer-review panel concluded that the assessment met its terms of reference and represented best scientific information available.

The cached Stock SMART summary still labels the model as ASAP. The assessment-specific SASINF code and the peer-review report document the 2025 transition from ASAP to an ASAP-like WHAM model, so the run-specific sources are used for model identification. The input file retains a stale comment calling it the “Updated 2023 redfish final model”; its 62-year span beginning in 1963, the 2025 Run 5 code, the September 2025 report, and the accepted terminal year 2024 identify the current production run.

## Recoverable inputs and outputs

The fitted model covers 1963-2024 and uses ages 1-25 plus a 26+ group. It has one commercial trawl fleet and two aggregate abundance indices, spring and fall. Catch input includes annual total catch and proportions-at-age. The surveys also include age compositions; their native month values are 4.5 for spring and 10.5 for fall. The ASAP file uses index-unit code 2 for both series; its exact physical-unit definition should be checked against the reader/model convention during import rather than guessed. Some age-composition years are incomplete. The final WHAM setup uses a logistic-normal AR1 age-composition likelihood that treats zeros as missing.

The input supplies three weight-at-age matrices, maturity-at-age schedules, and fixed natural mortality of 0.05 per year across ages. The peer review confirms the same weight matrix is used for the catch fleet, SSB calculation, and projections. Survey timing must be translated from the native month convention without substituting the spawning date or assuming a different survey season.

The fitted figures-and-tables bundle reports abundance-at-age in thousands of fish and fishing mortality-at-age by fleet, as well as recruitment and biomass time series. Parameter estimates include standard errors and 95% intervals. The bundle also contains projection results beyond the terminal historical year; those years must not be mixed with the fitted 1963-2024 series. Uncertainty for every age-year cell in the reported N/F surfaces was not verified in this screen.

## Identity and screening decision

No canonical NOAA Acadian redfish stock was found. Keep this as the distinct Gulf of Maine / Georges Bank NOAA unit. DFO and NAFO also manage redfish units containing the same species, but shared taxonomy alone does not establish that their assessment boundaries are identical. This candidate is **eligible for source curation**: direct catch-age compositions, aggregate survey indices with survey age compositions, run-specific biological schedules and fixed M, native timing, and comparable fitted age-specific and aggregate outputs are available. Resolve the exact index-unit code and retain missing composition years during the import.

The official source files, retrieval date, URLs, sizes, and SHA-256 values are recorded in the gitignored cache at `source_cache/discovery_noaa_acadian_redfish_2025/sources.csv`. The discovery register is updated; the canonical database remains unchanged.
