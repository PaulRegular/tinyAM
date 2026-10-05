# North Sea saithe source review

## Assessment identified

ICES assessed saithe (*Pollachius virens*) in subareas 4 and 6 and Division 3.a
under stock code `pok.27.3a46`. The 2026 advice is linked to assessment key
22504. Its published summary table reports recruitment at age 3, fishing
pressure for ages 4–7, and annual stock-summary values. For 2025 it reports
recruitment of 214,587 thousand, spawning-stock biomass of 116,722 tonnes,
total biomass of 259,593 tonnes, and F for ages 4–7 of 0.414. The table also
provides lower and upper uncertainty bounds for those quantities.

The [2026 WGNSSK report](https://agris.fao.org/search/es/records/6a96d767dd5645257f1c9d53)
is a detailed 1,183-page report published as ICES Scientific Reports 8:43.
Its [PDF](https://archimer.ifremer.fr/doc/01075/118685/133209.pdf) is publicly
listed. The 2026 ICES advice and source table are available through the
[Advice and Scenarios Database](https://doi.org/10.17895/ices.advice.30932291)
and [Stock Assessment Graphs, assessment key 22504](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22504).

## Assessment structure recovered

The [2024 WKBGAD benchmark report](https://ices-library.figshare.com/articles/report/Benchmark_workshop_on_selected_haddock_and_saithe_stocks_WKBGAD_/25002470)
describes North Sea saithe as an existing SAM assessment. It records the
benchmark's review of life-history inputs, survey and commercial CPUE index
models, and stock weights. The assessment used catch numbers-at-age and
catch weights-at-age for catch fractions raised through InterCatch, a
design-based North Sea IBTS Q3 index for ages 3–8, and a combined commercial
CPUE index scaled to exploitable biomass. The 2026 Stock Assessment Graphs
page confirms that the summary outputs report recruitment at age 3 and F for
ages 4–7.

The 2025 WGNSSK report is available as [ICES Scientific Reports 7:57](https://doi.org/10.17895/ices.pub.29085995).
The 2026 WGNSSK report supersedes it for the current assessment year. The
2024 benchmark is useful for model context, but neither an older stock run nor
benchmark settings can stand in for the accepted 2026 run's actual inputs.

## Remaining gap

The full 2026 report and model inputs have not been recovered into the local
source cache. Direct downloads of ICES/figshare and the linked report PDF fail
with TLS/access errors in this environment. A read of the historical
`NS_saithe_2024_benchmark_final` object through `stockassessment::fitfromweb()`
also failed at the TLS connection. The directory is a historical benchmark
run and is not being used as a substitute for the 2026 accepted assessment.

The accessible 2026 summary table does not provide the accepted run's
numerical catch-at-age, processed survey indices, annual stock and catch
weights, maturity, natural mortality, or age-specific N and F surfaces. The
2026 report is publicly listed, so these are retrieval/extraction gaps rather
than evidence that the inputs do not exist. Until the detailed report and
model files can be inspected and cached, this stock is not ready for a
canonical numerical input record or tinyAM fit. No output values are being
reused as model inputs.

## Sources

- [ICES advice 2026, assessment key 22504](https://doi.org/10.17895/ices.advice.30932291)
- [ICES Stock Assessment Graphs source data, key 22504](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22504)
- [WGNSSK 2026 report record and PDF link](https://agris.fao.org/search/es/records/6a96d767dd5645257f1c9d53)
- [WGNSSK 2025 report](https://doi.org/10.17895/ices.pub.29085995)
- [WKBGAD 2024 benchmark report](https://doi.org/10.17895/ices.pub.25002470)
