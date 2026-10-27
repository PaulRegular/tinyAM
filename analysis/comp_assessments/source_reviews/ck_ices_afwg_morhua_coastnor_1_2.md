# Historical Norwegian coastal cod catalogue unit

## Source identity

The metadata-only Charbonneau–Keith record is `ICES-AFWG_morhua_COASTNOR 1,2`
for Atlantic cod in Norwegian coastal waters. Its Master entry cites the AFWG
2017 report. Before the 2021 benchmark, the official ICES stock code was
`cod.27.1-2coast`.

ICES advice for 2018 classified the former combined unit as category 3.8. The
2018 Standard Graph record describes recruitment at age 2, fishing pressure
over ages 4–7, and stock size as a biomass index. Its annual XML contains
recruitment, catch, fishing pressure, and recreational catch, but no SSB,
N-at-age, F-at-age, or age-structured likelihood inputs. The XML is cached at
`source_cache/ck_ices_afwg_morhua_coastnor_1_2/assessment_2018/`; its manifest
records the source URL and checksum.

## Boundary change and database decision

The 2021 WKBARFAR benchmark split the former combined unit at 67°N. The north
became the age-structured SAM assessment `cod.27.1-2.coastN`; the south remained
data-limited. The current northern assessment is already represented by
`ices_cod_northern_coastal_2026`, but it covers only the northern component.
It cannot supply the former combined stock's missing inputs or outputs.

No separate age-structured production assessment with recoverable inputs and
population surfaces was found for the full historical 1–2 coastal unit. Keep
this catalogue record metadata-only. Do not copy the northern SAM data or
aggregate advice indicators into the combined record; no tinyAM translation
or comparison dashboard is supported.

## Sources

- [ICES 2017 advice for the former combined stock](https://asd.ices.dk/viewAdvice/2582)
- [ICES 2018 Standard Graph source data, assessment key 9915](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=9915)
- [ICES 2018 assessment XML, key 9915](https://sag.ices.dk/SAG_API/StockXmlDownload/9915)
- [WKBARFAR 2021 benchmark](https://doi.org/10.17895/ices.pub.7920)
- [Current northern coastal cod assessment, 2026](https://standardgraphs.ices.dk/ViewCharts.aspx?key=22637)
