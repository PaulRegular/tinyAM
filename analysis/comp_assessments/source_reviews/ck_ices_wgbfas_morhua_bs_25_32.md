# Eastern Baltic cod in subdivisions 25–32

## Assessment selected for the partial database record

The latest recoverable official age-specific output tables are from the 2024
WGBFAS Stock Synthesis assessment. It is a one-area, quarterly,
length-based Stock Synthesis 3.30 model for both sexes combined, with age
classes 0–15+ and spawning time set to month 5. The model includes BITS Q1
survey data through 2024; catch for 2024 was assumed equal to 2023 catch.
Because the Charbonneau–Keith stock row identifies subdivisions 25–32, the
stock identity remains unchanged; the model also includes eastern Baltic cod
from the eastern part of subdivision 24.

The 2024 report tabulates:

- stock numbers at the beginning of the year for ages 1–8, 1946–2024, in
  thousands of fish;
- fishing mortality at ages 1–14, 1946–2023, per year;
- recruitment and SSB relative to the time-series average, plus relative
  fishing mortality for ages 4–6, in Table 2.1.11.

The report does not tabulate age-specific uncertainty. Its age tables do not
include modeled age 0 or ages 9–15+, and age 8 is not identified as a plus
group. For that reason no total abundance is derived from the partial N-at-age
surface. The values are stored as partial official outputs with these limits
recorded on each row.

One F-at-age cell needs source verification: Table 2.1.13 prints age-11 F in
2023 as `0.4`, while the neighboring age cells are around `0.03–0.04` and the
2022 age-11 value is `0.03`. The 2025 WGBFAS report repeats the cell in its
table, marked as not updated in 2024. The database preserves the printed value
and flags it; it is not silently changed to `0.04`.

## Input availability and translation readiness

WGBFAS 2024 Table 2.1.8 describes the input categories and coverage, and the
report discusses the model assumptions. It does not provide the complete
numerical series used for catch, length and age compositions, survey indices,
or biological surfaces in a reusable form. The source review therefore records
the known model structure and assumptions but leaves `inputs.csv` empty and
marks input recovery as not started. It does not treat estimated catch-at-age
or fitted outputs as observed model inputs. The record is not ready for a
tinyAM fit, and no translation is created.

The WGBFAS 2025 report says no new advice was requested for this stock and
reprints its age-specific N and F tables as not updated in 2024. It provides
no newer numerical age-specific population surfaces. Its newer relative
summary values are not treated as a separate at-age assessment output here.

The latest advice record is from 2026 and uses a survey-index trend assessment
(ICES category 3.2) based on standardized BITS Q1/Q4 indices, not an
age-structured population model. It has no age-specific N/F/M estimates. The current index-based assessment is recorded separately as
`ices_cod_eastern_baltic_2026`; the 2024 age-structured outputs remain
available as a distinct historical record. The exploratory SPiCT
analysis in the 2026 report is not the method used for current advice.

Sources:

- [WGBFAS 2024 report](https://doi.org/10.17895/ices.pub.25764978),
  sections 2.1.4–2.1.5 and Tables 2.1.8–2.1.13.
- [WGBFAS 2025 report](https://doi.org/10.17895/ices.pub.29099786),
  section 2.1 and Tables 2.1.11–2.1.13.
- [WGBFAS 2026 report](https://doi.org/10.17895/ices.pub.32455056), section 2.1.
- [ICES SAG source data, 2026 assessment key 22502](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22502)
  and [XML export](https://sag.ices.dk/SAG_API/StockXmlDownload/22502).
- [2026 advice record 6135](https://asd.ices.dk/viewAdvice/6135).
- [WKBALTCOD2 2019 benchmark](https://doi.org/10.17895/ices.pub.4984).

The 2023, 2024, 2025, and 2026 reports or machine-readable files reviewed for
this stock are cached locally under the gitignored
`analysis/comp_assessments/source_cache/ck_ices_wgbfas_morhua_bs_25_32/`
directory. Per-assessment manifests record source URLs and checksums.

## Related catalogue entry

The Master catalogue also contains
`ck_ices_wgbfas_morhua_wbs_25_32`, with source label “BS gadus morhua.pdf” and
area text `WBS,25-32`. That label does not identify a distinct current ICES
assessment code or resolve its apparent mismatch between “western” and the
25–32 area. The western Baltic stock is assessed in subdivisions 22–24; no
independent assessment identity was verified for this catalogue row. It
remains metadata-only, and no data are copied between the two entries.

## Current 2026 index assessment

The current assessment uses a GAM-standardized combined BITS Q1/Q4 biomass
index in kg/hour from 1991–2025, including the eastern part of subdivision 24.
The index combines numbers-at-length with annual weight-at-length. The SAG
stock-size field is this survey index, not SSB. It is recorded as an input and
as the same reported indicator; these are not independent series. The ICES SAG
export also provides 1946–2025 landings and 1966–2025 discards, in tonnes.
They are contextual catch series, with recent data incomplete: the working
group reports Russian landings missing for 2022–2025 and limited discard
sampling. No age-specific population states or separate seasonal survey series
are available from the current method. No tinyAM fit or dashboard is supported
by the current aggregate index alone.

The 2026 record is now the selected current assessment for this stock; the
2024 Stock Synthesis assessment remains a separate historical record. The
2026 source files and their checksums are listed in the existing cache
manifest. See [ICES source data](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22502)
and [WGBFAS 2026](https://doi.org/10.17895/ices.pub.32455056), section 2.1.
