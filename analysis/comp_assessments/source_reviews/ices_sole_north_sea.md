# North Sea sole source review

## Assessment selected

The stock is ICES `sol.27.4`, Solea solea in Subarea 4 (North Sea), matched to the existing Charbonneau–Keith catalogue record `ck_ices_wgnssk_solea_ns_4`. The latest detailed production assessment recovered is WGNSSK 2025. Its model uses data through 2024 and reports a 2025 forecast. The 2024 WKBFLATFISH benchmark introduced the Stock Synthesis 3 reference model.

ICES also publishes a 2026 assessment-year SAG series. The retrieved 2026 source provides annual aggregate recruitment, SSB and fishing-pressure summaries, but no age-specific inputs or outputs. It is recorded separately as a summary-only update and does not replace the 2025 detailed assessment as the current detailed record.

## Sources and local cache

- 2025 stock chapter: [WGNSSK 2025, North Sea sole chapter](https://ndownloader.figshare.com/files/59358092), from [ICES WGNSSK report, version 4](https://doi.org/10.17895/ices.pub.29085995.v4). Tables 16.2–16.11 provide catch, biological inputs, surveys, N-at-age and F-at-age.
- Framework: [WKBFLATFISH 2024](https://doi.org/10.17895/ices.pub.25471987.v1). The sole chapter documents survey standardization, benchmark choices and SS3 structure.
- Annual summaries: [ICES SAG 2025, key 21217](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=21217) and [ICES SAG 2026, key 23777](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=23777), with the corresponding XML exports.

The PDFs, extracted text and SAG HTML/XML are cached locally under the gitignored `source_cache/ices_sole_north_sea_2025/`. `SOURCES.md` there records retrieval date, source versions, URLs and SHA-256 hashes. Public ICES TAF/GitHub searches did not identify a matching native 2025 SS3 run archive; this import therefore uses the detailed official report tables and the ICES SAG machine-readable outputs.

## Recovered inputs

- Landings-at-age: 1957–2024, ages 1–10+, in thousand fish. Discards-at-age: 2002–2024, ages 1–10+, in thousand fish. Age 0 is not reported in these catch tables and is left absent.
- Landings weights-at-age: ages 1–10+, 1957–2024. The report explicitly gives age-0 landings weight as a constant 0.01 kg; that value is recorded for each year. Discard weights and stock weights cover ages 0–10+ and 1957–2024. All weights are in kg per fish.
- Maturity and natural mortality are fixed age vectors for ages 0–10+. Maturity is knife-edged at age 3. M is not estimated; the reported rates follow the size-based Lorenzen formulation, with exponent −0.288 and scaling based on maximum validated age 26.
- The age-specific survey series are the Q3 BTS index (ages 0–10+, 1990–2024) and the Q3 COAST index combining SNS and DYFS (ages 0–5, 1990–2024). The report’s index tables do not state a physical multiplier, so original values are retained as relative-index units. TinyAM represents Q3 with the quarter midpoint, 0.625.
- Aggregate BTS and COAST index series from Tables 16.8b and 16.9b are retained separately. The report says aggregate survey-biomass indices also inform model scale; their reported values closely match the rounded sums of the age-specific tables. Their dependence is unresolved, so tinyAM must not fit both representations without a deliberate choice.

## Recovered outputs and assumptions

The 2025 report provides point estimates of N-at-age (1957–2024, ages 0–10+) and F-at-age (1957–2024, ages 1–10+); it does not provide age-specific uncertainty. The ICES SAG 2025 XML supplies recruitment, SSB and Fbar summaries and lower/upper bounds. Its bounds are defined as mean ± 2 SD; those bounds are retained directly, and no standard errors are inferred. Recruitment is in thousand fish, SSB in tonnes, and Fbar is the mean F over ages 2–6. The ICES SAG 2026 XML supplies the comparable aggregate summary series through 2026; the 2026 row contains recruitment and SSB but no Fbar.

The production model is single-area and single-sex, with separate landings and discards fleets, age-0 recruitment at the beginning of the year, a Beverton–Holt stock-recruitment relationship with fixed steepness 0.6646, time-varying landings selectivity, time-invariant discards and survey selectivity, and selectivity/catchability constraints from age 5 onward. SSB is calculated at the start of the year. The 2016 age-reading-error matrix is retained in the assessment, but its numerical values were not recovered from the cached report tables.

## Completeness and translation limits

The 2025 assessment has partial input, assumption and output status because native SS3 files, some likelihood weighting details, and the numeric age-reading-error matrix were not recovered. The 2026 summary-only record has no detailed inputs or assumptions and partial aggregate outputs. The tinyAM translation must keep landings and discards distinct in the database, document any catch aggregation, and choose between correlated age-specific and aggregate survey observations so it does not double count information. Fixed M is a source input, not an estimated M output.

## tinyAM translation and fit

The translation fits 1957–2024 at ages 1–10+, because source catch and F-at-age begin at age 1. It sums landings and discards by year and age where both are reported; before the discard series begins in 2002, the model input contains landings only. It fits the two age-specific Q3 indices at timing 0.625, groups q by survey and age 1–4 versus 5+, and omits the closely related aggregate indices. The index values remain on their reported relative scales. Fixed M, annual stock weights, maturity, and start-of-year SSB timing are retained.

tinyAM estimates an IID cohort-survival process, a random walk for log recruitment, and a random-walk F process with mean ages 2–6. The source's fixed-steepness Beverton–Holt curve cannot be represented with steepness fixed in tinyAM. An initial attempt to estimate the tinyAM Beverton–Holt parameters had a non-positive-definite Hessian and undefined standard errors; switching to the simpler log-recruitment random walk resolved that numerical failure without changing observations or other model settings. This omits the source stock-recruitment curve and should be treated as an approximation.

The final fit passed `check_tam()`: optimizer code 0, maximum adjusted gradient about 0.001, positive-definite Hessian, and successful uncertainty estimation. Residual diagnostics still flag some standardized residual spreads outside 0.5–1.5 and some residuals above 3 SD. The cached dashboard is local at `results/cache/ices_sole_north_sea_2025/dashboard.html`.

Common-definition agreement over 1957–2024 is mixed. SSB has 11.7% mean absolute percent difference, −41.8% terminal difference, and 0.957 trend correlation. Recruitment is 29.1%, +5.9%, and 0.662; N-at-age is 27.2%, −41.1%, and 0.913; F-at-age is 31.7%, +106.5%, and 0.731. Total abundance has 16.2% mean absolute percent difference and 0.693 trend correlation. M is identical because it is fixed to the same source input, not independently estimated. These numbers describe this simplified fit and are not evidence of equivalence to the SS3 likelihood.
