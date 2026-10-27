# Herring in divisions 6.a South and 7.b–c

## Source screening: 2026 assessment

The latest detailed public assessment located is the ICES Herring Assessment
Working Group (HAWG) report for 2026, using data through 2025. It reports the
current method for `her.27.6aS7bc` as ICES category 3.2, WKLIFE method 2.2: a
constant harvest-rate rule using the genetically split Malin Shelf Herring
Acoustic Survey index. This is not an age-structured population assessment.

The 2026 advice record (SAG key 21289; advice record 4061) is marked **Invalid**.
The 2025 advice record (SAG key 19342; advice record 3843) is marked **Active**.
The 2026 report is retained as the latest detailed account of the method, but
it is not treated as applied advice.

Sources:

- [HAWG 2026 full report](https://backend.orbit.dtu.dk/ws/portalfiles/portal/440308463/HAWG_2026_-_Full_Report.pdf), sections 5.3 and 5.6 and Tables 5.3.1.1 and 5.6.2.5a.
- [ICES SAG source data, assessment key 21289](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=21289) and [XML export](https://sag.ices.dk/SAG_API/StockXmlDownload/21289).
- [2026 advice record 4061](https://asd.ices.dk/viewAdvice/4061) and [2025 advice record 3843](https://asd.ices.dk/viewAdvice/3843).
- [HAWG 2023 report](https://ices-library.figshare.com/articles/report/Herring_Assessment_Working_Group_for_the_Area_South_of_62_N_HAWG_/22182034), documenting the earlier assessment-method history.

The HAWG 2026 report, SAG source page and XML, and both advice pages are
cached locally under the gitignored `analysis/comp_assessments/source_cache/`
directory. A download of the 2023 report PDF was blocked by the publication
host's challenge page; its online report page remains available at the link
above. The local cache manifest records retrieved files, URLs and checksums.

## Why no assessment data were imported

The assessment database is scoped to age-structured assessments. The current
official production method is a survey-index harvest rule, so its catch and
survey series are not a complete input set for an age-structured fit, and its
reported `Frel` is not instantaneous F or an estimated fishing-mortality
surface. The 2026 report includes age-specific information from exploratory
model work, but those values are not outputs of the adopted production rule.
They are not relabeled as official N-at-age, F-at-age, recruitment or M-at-age.

The 2014–2025 catch and split survey-index series were checked against the
report and SAG XML. For 2025, the series reports 1,378 tonnes of catch, a
7,774-tonne survey biomass index and `Frel` of about 0.1773. The index is on a
reported tonne scale and is not absolute population SSB. These indicators were
not added to the canonical assessment tables because doing so would mix a
category-3 harvest-rule assessment into the age-structured database and would
not support a faithful tinyAM translation or Official–tinyAM dashboard.

The database row for this Charbonneau–Keith entry therefore remains
metadata-only. A separate, already populated row records the historical 2020
SAM assessment for the combined 6.a and 7.b–c area; that older assessment is
not a substitute for the current split-stock method.

## Translation status

No tinyAM fit or dashboard is appropriate for this record. Reconsider only if
the official assessment adopts an age-structured population model or the
analysis scope is expanded to include category-3 harvest-rule assessments.
