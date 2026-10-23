# Beaked redfish in ICES subareas 1 and 2 (`Sebastes mentella`)

## Source review status

The 2026 JRN-AFWG assessment is the latest detailed official assessment located. This review records its published coverage and unresolved source gaps; it does not yet create an assessment event or model translation. The stock remains a catalogue-only database entry until the official input surfaces and their likelihood roles can be resolved without inference.

## Sources

- [JRN-AFWG Report 2026, chapter 6](https://www.hi.no/hi/nettrapporter/imr-vniro-2026-5), assessment and table source. The PDF and structured report page are cached locally with retrieval dates and SHA-256 checksums in `source_cache/ck_ices_afwg_mentella_deep_1_2_2026/manifest.txt`.
- [WKREDFISH benchmark report 2018](https://ices-library.figshare.com/articles/report/Report_of_the_Benchmark_Workshop_on_Redfish_Stocks_WKREDFISH_/19291142), model framework source.
- [2026 ICES advice](https://www.hi.no/en/hi/nettrapporter/imr-vniro-en-2026-3), confirms the current SCAA assessment and advice basis.

## Verified assessment structure

The 2026 report describes an analytical statistical catch-at-age assessment for 1992–2025, configured as the 2018 benchmark baseline. The population matrix and outputs use ages 2–19+, with age 2 recruitment and a 19+ group. The assessment has demersal and pelagic fishing mortality components, fixed natural mortality of 0.05 per year, and annual demersal selectivity. The 2018 benchmark documents autoregressive processes for age-2 recruitment, annual F in both fleets, and annual demersal selectivity. The 2026 report says the fitted assessment was configured as that baseline with 53 estimated parameters.

## Published input and output coverage

The official report contains machine-readable tables for the following source series:

- Total catch-at-age, ages 6–19+, in thousands, with total landings in tonnes (Table 6.8; 1992–2025). The report identifies estimated age compositions in specified years when samples were unavailable.
- Pelagic catch-at-age outside the EEZ, ages 7–19+, in thousands (Table 6.9; 2006–2025). Several years are estimated or lack direct age data, with footnotes describing the method.
- Stock weight-at-age and maturity-at-age (Tables 6.12 and 6.14). Their printed age columns begin at age 6, although the report's data inventory describes these biological series as ages 2–19+.
- Russian autumn survey catch rates by age, ages 0–11, in fish per trawling hour (Table 6.15; survey ended after 2017 because coverage was insufficient).
- Ecosystem survey abundance by age, ages 2–15 and 16+, in millions, for 1996–2009 and 2011–2025 (Table 6.17; August–September).
- Winter survey abundance by age, ages 2–15, in millions, for 1992–2011 and selected later years (Table 6.18b; February).
- Deep pelagic survey proportions-at-age, ages 2–75, in selected survey years (Table 6.19). These are compositions, not absolute abundance indices.
- Fitted population numbers-at-age, ages 2–19+, in thousands, and fleet-year F coefficients (Table 6.20a); pelagic and demersal selectivity by age (Tables 6.20b–c); and recruitment, total biomass, SSB and F summaries (Table 6.21). The report tables do not provide uncertainty for these fitted outputs.

## Unresolved items before import or fitting

1. The source model spans ages 2–19+, but the published stock-weight and maturity tables expose only ages 6–19+. The required ages 2–5 values are not recoverable from these tables. They must not be set to zero or extrapolated without an authoritative source.
2. Table 6.8 is labelled total catch-at-age, while Table 6.9 separately reports pelagic catch-at-age outside the EEZ. The exact model input relationship between the overlapping series and the two reported F fleets has not been verified. Treating both as independent catches could double-count removals.
3. The report documents some catch-age estimates and survey gaps, but does not fully specify all likelihood weights and variance-sharing choices needed to reproduce the 53-parameter run.
4. The native 2026 SAM fit object was not located. The public `stockassessment::fitfromweb("ple.27.7d_2026")` endpoint resolves to a missing `run/model.RData` file (HTTP 404); the available public material includes the report tables but not the fitted object.

These gaps prevent a source-faithful tinyAM fit at ages 2–19+ at this stage. Do not mark inputs complete, infer missing young-age biology, or fit both overlapping catch series. Continue with another stock while searching for the missing original model inputs or a documented age-6+ model specification.
