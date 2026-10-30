# ICES Bay of Biscay sea bass (`bss.27.8ab`), 2026 assessment

## Identity and latest assessment

The ICES 2026 stock list and advice record 6153 identify the current unit as
sea bass (*Dicentrarchus labrax*) in divisions 8.a–b, northern and central Bay
of Biscay. The assessment year is 2026 and the advice is for 2027. The latest
production assessment is an update using data through 2025 and the model
structure accepted at the 2025 WKBSEABASS benchmark. Seasonal movement is
accounted for by reallocating part of catches from the northern sea bass unit;
the Bay of Biscay unit remains distinct.

## Inputs and assumptions

- The assessment is an age- and length-based Stock Synthesis model. Commercial
  catches are represented by landings and discards, recreational removals, and
  length compositions. The 2026 report says annual conditional age-at-length
  observations (age-length keys) are combined into a commercial fleet to
  provide catch-at-age estimates. Marginal age compositions by fleet are not
  included in the likelihood. This is a reconstructible age input, but it must
  not be described as a directly observed catch-at-age likelihood.
- Two research surveys, Nourdem Loire (from 2016) and Nourdem Gironde (from
  2019), provide juvenile abundance indices. An adult OTB LPUE index covers
  2000–2025, with 2000 and 2008 excluded from the fitted series. The 2025
  benchmark describes the juvenile indices as covering ages 1–4; the 2026
  update reports aggregate juvenile-index series and does not publish
  age-by-age index values. Survey age composition and timing need to be
  verified before an age-specific tinyAM translation.
- The 2026 update estimates most von Bertalanffy growth parameters and fixes
  length-at-age 30 at 87.7611 cm. Female maturity is length-based, with
  L50 = 42.14 cm and logistic slope -0.37809. These quantities support
  deriving maturity-at-age, but the report does not tabulate a maturity-at-age
  vector. Exact weight-at-age inputs and their conversion remain to be
  extracted from the model input/control files.
- Natural mortality is specified at ages 0–5 as
  0.840918, 0.463849, 0.365109, 0.306347, 0.267492, and 0.24 per year;
  ages 6 and older use 0.24. Recruitment occurs at the beginning of the year
  and follows a Beverton–Holt stock–recruit relationship with fixed steepness
  0.8 and recruitment-deviation SD 0.6.
- Seasonal catch reallocations use exchange rates of 0.28803 for Q3 and
  0.40665 for Q4. This is a stock-mixing treatment, not survey timing.
  Exact fractional-year timing for the three abundance indices was not
  recovered from the reviewed 2026 section.

## Outputs and eligibility

ICES Standard Graphs assessment key 22538 reports age-0 recruitment in
thousands, combined-sex SSB in tonnes, and F averaged over ages 4–15, with a
90% interval. These are annual assessment outputs; the report's data cutoff is
2025, so 2026 is the forecast year. It does not publish historical N-at-age,
F-at-age, or M-at-age arrays in the standard graph.

Screened **eligible for source curation** because the official assessment
provides catch-at-length and conditional age-at-length source quantities for
reconstructing catch-at-age, juvenile abundance indices with age-composition
information, age-based mortality, growth and maturity assumptions, and
comparable annual recruitment, SSB and F outputs. Eligibility does not mean
that numerical inputs are ready for import: the age-length data, exact
weight-at-age derivation, survey timing, and their units still require
source-level extraction and validation. No canonical database records were
added.

## Sources

- [2026 advice record 6153](https://asd.ices.dk/viewAdvice/6153) and [Standard Graphs source data, key 22538](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22538).
- [WGBIE 2026 report, section 11](https://doi.org/10.17895/ices.pub.32790321.v1), especially sections 11.2.1–11.3.
- [WKBSEABASS 2025 benchmark report](https://doi.org/10.17895/ices.pub.25002488.v1), especially the Bay of Biscay stock model annex and its Stock Synthesis data/control files.

The ignored local source cache contains the reviewed report files, advice,
graph page, source data, XML, retrieval dates, and SHA-256 checksums.
