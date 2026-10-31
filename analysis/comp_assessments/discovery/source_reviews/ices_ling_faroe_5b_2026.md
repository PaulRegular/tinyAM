# ICES Faroes ling (lin.27.5b), 2026

**Identity:** Ling (*Molva molva*) in Division 5.b, Faroese grounds. This is a distinct management unit from the Icelandic ling unit lin.27.5a. The canonical database contains no matching Faroes ling stock record.

**Current assessment:** ICES advice record 7205 and SAG assessment key 23735 identify the 2026 Category 1 assessment. WGDEEP 2026 identifies the final SAM run as lin.27.5b_wgdeep2026_final. The full report is ICES Scientific Reports 8:30 (DOI 10.17895/ices.pub.32051352.v1); the 2021 WKBARFAR benchmark report and its Faroes ling stock annex document the SAM benchmark and age-length-key (ALK) procedures.

## Published inputs

- **Catch numbers-at-age:** Table 3.2.7 gives 1996–2025, ages 3–12+, in thousands of fish. Catch-at-age is built from Faroese commercial landings; the report discusses low sampling and the use of weight-grader data for longliners in 2025.
- **Survey indices-at-age:** Tables 3.2.9 and 3.2.10 give numeric spring-survey indices for ages 4–11 in 1998–2025 and summer-survey indices for ages 3–11 in 1996–2025. Each table includes an effort column. WGDEEP describes spring coverage as February/March and summer coverage as August. The 2026 report says summer indices in 2021–2023 and spring indices in 2022–2024 were downweighted after a change in rockhopper gear; it does not publish numeric weights. The report tables also do not explicitly label the absolute index scale, so verify the index scale, effort treatment and downweight factors against the final SAM data before import.
- **Weights-at-age:** Table 3.2.8 publishes weighted mean catch weights at age for ages 3–12+ in 1996–2025. The table does not label the weight unit. The benchmark stock annex explains that catch weights were selected for the assessment's stock weights rather than spring-survey weights; confirm the native unit before import.
- **Maturity-at-age:** Fixed, time-invariant Faroese survey maturity proportions are published for ages 3–12+: 0.00, 0.04, 0.19, 0.50, 0.79, 0.93, 0.98, 1.00, 1.00, and 1.00.
- **Natural mortality:** Fixed M = 0.15 for every modeled age and year.
- **Survey timing:** Spring survey in February/March; summer survey in August. Month-based fractional-year times would be approximations, not recovered model parameters.

## Published outputs

WGDEEP 2026 Tables 3.2.11–3.2.13 report annual recruitment at age 3, SSB, Fbar for ages 6–10, total biomass, F-at-age, and stock numbers-at-age for 1996–2025. The annual recruitment, SSB, Fbar, and total-biomass series include 95% limits. The age-specific F and stock-number tables are point estimates without age-specific uncertainty. ICES Standard Graphs identifies SSB in tonnes, recruitment in thousands, and F for ages 6–10.

## Screening decision

This is a **new, eligible unit for source curation**. The latest detailed official report publishes numeric catch-at-age, two numeric age-disaggregated survey series, weights, fixed maturity and fixed M, and comparable SSB, recruitment, and mortality outputs. Full N-at-age and F-at-age estimates are also tabulated. Before importing values, verify the survey index scale/effort, downweight factors and explicit weight units against the native SAM run; do not infer those labels from numerical magnitudes. Survey fractional-year times will remain documented approximations if the exact model values are not recoverable.

No canonical stock, assessment, input, or output rows have been added, and no tinyAM translation has been attempted. The Faroes unit remains distinct from ICES lin.27.5a (Icelandic waters). Source files are cached in the gitignored discovery cache; the 2026 WGDEEP PDF is reused from the shared WGDEEP cache rather than copied.

## Sources

- WGDEEP 2026 report, ICES Scientific Reports 8:30: https://doi.org/10.17895/ices.pub.32051352.v1
- ICES 2026 advice record 7205: https://asd.ices.dk/viewAdvice/7205
- ICES Standard Graphs source data, assessment key 23735: https://standardgraphs.ices.dk/ViewSourceData.aspx?key=23735
- ICES SAG XML, assessment key 23735: https://sg.ices.dk/SAG_API/StockXmlDownload/23735
- 2026 final SAM-run directory referenced by WGDEEP: https://stockassessment.org/datadisk/stockassessment/userdirs/user3/lin.27.5b_wgdeep2026_final/
- WKBARFAR 2021 report and stock annex: https://doi.org/10.17895/ices.pub.7920
- The report links the stock annex at https://www.ices.dk/sites/pub/Publication%20Reports/Stock%20Annexes/2021/lin.27.5b_SA.pdf; the complete annex material is also present in the cached WKBARFAR report. The direct legacy ICES URL returned HTTP 401 during retrieval.

Cached source URLs, retrieval dates, file sizes, revisions, and SHA-256 checksums are recorded in source_cache/discovery_ices_ling_5b_2026/sources.csv. The reused 2026 WGDEEP PDF and its checksum are also recorded in source_cache/discovery_ices_gss_5b6a_2026/sources.csv.
