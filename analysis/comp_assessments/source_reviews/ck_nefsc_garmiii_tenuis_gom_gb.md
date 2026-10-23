# White Hake (Gulf of Maine–Georges Bank)

## Assessment and source selection

The stock record is ck_nefsc_garmiii_tenuis_gom_gb (Charbonneau–Keith identifier NEFSC-GARMIII_tenuis_GOM, GB), White Hake (*Urophycis tenuis*) in the Gulf of Maine–Georges Bank area. The current detailed record is NOAA NEFSC's 2025 Management Track update, ASAP 3.0 Run_005, with data through 2024. The final September 2025 peer-review panel report concluded the assessment met its terms of reference and represented best scientific information available for management. The run archive includes the native RUN_005.DAT input and a 25-page table supplement. The SASI working paper is marked draft for peer review; its status is made clear here and it is used for detailed input/output tables alongside the final review decision.

The 2022 Management Track assessment updated the 2013 SAW56 benchmark and provides background on model history. Attempts to cache those NOAA Institutional Repository PDFs returned Access Denied. Their NOAA links and the failed retrievals are recorded in the local source-cache manifest; no 2022 or 2013 numerical values are mixed into the 2025 record.

## Model and inputs

The ASAP input covers 1963–2024, ages 1–9+ and a sex-combined stock. The single combined catch fleet has catch biomass in every model year and catch numbers-at-age with effective sample size 58 for 1989–2024; the source catch-age ESS is zero before 1989. The source discard and release arrays are zero. The age-specific catch numbers are stored in thousands of fish, with a total-catch biomass observation in metric tons and a separate annual CV. Multiplying catch numbers by the matrix-1 age weights gives an internal cross-check against observed total catch. The supporting report's Table 5 labels its catch-at-age table in metric tons while repeating the native age values; the native ASAP field structure and the weight-based cross-check support the thousands-of-fish interpretation used here. This table-label inconsistency remains a source caveat.

The model fixes natural mortality at 0.20 yr-1 for ages 1–9+. Maturity is fixed by age. Three annual weight matrices are retained: matrix 1 is assigned to catch, matrix 2 to SSB and all five surveys, and matrix 3 to start-of-year biomass. SSB is calculated at 0.25 of the year. The two fishery selectivity blocks are 1963–1997 and 1998–2024 and use age-specific selectivity; fully selected F is the unweighted mean over ages 6–9+.

All five model indices are included: NEFSC spring and fall bottom trawl, ASMFC shrimp trawl, and NEFSC spring and fall bottom longline. The native index and age-composition units are numbers. Survey timing is April (0.25), July (0.50), or October (0.75), as specified by the ASAP input month codes. Aggregate indices have separate CVs; age-composition effective sample sizes are retained on the age rows. Annual q deviations are inactive in Run_005. Exact index observations, age rows, weights, timing, and sample sizes are retained in the database.

The input file also records estimated recruitment deviations, estimated initial numbers-at-age and a stock-recruit relationship. tinyAM's available initial-state and process options do not reproduce ASAP's exact integration and likelihood setup; those choices are described in the translation background rather than represented as equivalent.

## Reported outputs and remaining gaps

Table 14 reports annual SSB, fully selected F, age-1 recruitment, and their CVs for 1963–2024. The database retains those estimates and converts each reported CV to an approximate natural-scale SE as estimate multiplied by CV; the source CV is repeated in the row notes. This is a derived uncertainty scale, not a separately reported SE.

The public run archive did not include a numerical ASAP .rep output. No full fitted N-at-age or F-at-age surface was found in the report tables. Those outputs and their uncertainty remain unavailable; figures are not digitized and estimates are not reconstructed from SSB or recruitment. The output status is therefore partial.

## tinyAM translation

The translation uses the complete 1963–2024 year and 1–9+ age grid, source M, maturity, spawning timing, weights, all five indices, index timing, and available catch-age data. It keeps the original annual catch-biomass and index-CV inputs in the canonical database. Since tinyAM's catch observations are age-specific, its recipe uses native catch numbers-at-age from 1989–2024; total catch biomass from 1963–1988 cannot be fitted as an independent aggregate observation. Source zero catch and survey-age values are retained in the database and omitted from tinyAM's lognormal observation rows, which cannot represent exact zeros.

The source combines total-index and age-composition likelihoods; tinyAM represents the available age-specific numbers as direct age-index observations and applies the aggregate CV as a shared log-scale SD within each survey-year. This does not reproduce ASAP's multinomial composition likelihood or its effective-sample-size weighting. tinyAM also uses its own N, F, catch and q process structures. Comparisons are available for SSB, recruitment and Fbar only. Fbar is an age-aggregated source output and must not be interpreted as a comparison of F-at-age.

## Sources

- [Final September 2025 Management Track Peer Review Panel Report](https://www.fisheries.noaa.gov/s3/2025-11/September-2025-Management-Track-Peer-Review-Panel-Report_Final-508.pdf)
- [2025 White Hake assessment working paper and report tables](https://apps-nefsc.fisheries.noaa.gov/saw/sasi_files.php?year=2025&species_id=8&stock_id=6&review_type_id=3&info_type_id=-1&map_type_id=&filename=2025_HKW_UNIT_RPT.pdf)
- [2025 Run_005 ASAP input](https://apps-nefsc.fisheries.noaa.gov/saw/sasi_files.php?year=2025&species_id=8&stock_id=6&review_type_id=3&info_type_id=3&map_type_id=&filename=RUN_005.DAT)
- [2025 White Hake run archive](https://apps-nefsc.fisheries.noaa.gov/saw/sasi_files.php?year=2025&species_id=8&stock_id=6&review_type_id=3&filename=2025_HKW_UNIT.zip)
- [2022 Management Track assessment](https://repository.library.noaa.gov/view/noaa/55264/noaa_55264_DS1.pdf)
- [2013 SAW56 White Hake assessment summary](https://repository.library.noaa.gov/view/noaa/4374/noaa_4374_DS1.pdf)

The local, gitignored source cache and checksums are recorded in source_cache/ck_nefsc_garmiii_tenuis_gom_gb/manifest.csv.
