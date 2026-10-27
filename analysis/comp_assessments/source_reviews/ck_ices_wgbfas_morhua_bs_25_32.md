# Eastern Baltic cod in subdivisions 25–32

## Source screening: 2026 assessment

The latest detailed source reviewed is the 2026 ICES Baltic Fisheries
Assessment Working Group (WGBFAS) report, with assessment key 22502 and stock
code `cod.27.24-32`. The Charbonneau–Keith row records the management area
25–32; the current assessment also includes Eastern Baltic cod from part of
Subdivision 24, east of 13°E, where it overlaps with the western stock.

The report states that the analytical category-1 assessment was downgraded to
category 3 in 2024 because of recruitment-related model problems, then further
downgraded in 2026 to a survey-index assessment because biological data,
landing data and sampling were insufficient or unreliable. The 2026 advice
record is valid and category 3.2. The current method uses a standardized BITS
Q1/Q4 survey index from 1991 onward, including the eastern-stock portion of
Subdivision 24. The report's exploratory SPiCT model is not the method used
for current advice.

Sources:

- [WGBFAS 2026 report](https://doi.org/10.17895/ices.pub.32455056), section 2.1.
- [ICES SAG source data, assessment key 22502](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22502) and [XML export](https://sag.ices.dk/SAG_API/StockXmlDownload/22502).
- [2026 advice record 6135](https://asd.ices.dk/viewAdvice/6135).
- [2024 Eastern Baltic cod stock annex](https://doi.org/10.17895/ices.pub.25930969), cited by the 2026 report for stock boundaries and historical assessment details.

The WGBFAS report, SAG XML and advice page are cached locally under the
gitignored `analysis/comp_assessments/source_cache/` directory. The local
manifest records their source URLs and SHA-256 checksums.

## Why no assessment data were imported

The database records age-structured assessments. The current official method
for Eastern Baltic cod is an index-based trend assessment, not an age-structured
population model. Catch and biological tables in the working-group report do
not make them inputs to the current survey-index assessment. The report also
contains an exploratory SPiCT model, but its output is not the official method
used for the 2026 advice. Therefore no age-specific population states,
instantaneous fishing-mortality surface, or tinyAM-compatible full input set
was identified for the current assessment.

This Charbonneau–Keith entry remains metadata-only. No partial catch or survey
series was added to the canonical assessment tables, and no tinyAM fit or
dashboard is produced. Reconsider if ICES returns to an age-structured method
or if the comparison database is expanded to include survey-index assessments.
