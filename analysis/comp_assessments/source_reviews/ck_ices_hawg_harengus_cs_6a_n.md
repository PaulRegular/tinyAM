# North of 56°N herring (6.a.N): source review

## Stock and assessment history

The Charbonneau–Keith entry `ICES-HAWG_harengus_CS 6a (N)` corresponds to
North Sea autumn-spawning herring in ICES Division 6.a North, stock code
`her.27.6aN`. The stock was assessed separately through 2014, jointly with
6.a South and 7.b–c from 2015 to 2021, and separately again from 2022 using
a genetically split acoustic survey. The combined historical assessment and
the two current units must not be treated as the same stock or merged.

HAWG 2026 presents a final 6.aN survey-index assessment based on data through
2025. The published 2026 advice record 4063, assessment key 21308, is marked
**Invalid** in the ICES Advice and Scenarios Database. Do not assume that its
advice was applied, or infer the reason for invalidation. The 2025 ICES ACOM
advice for this stock uses assessment key **19370** and DOI
`10.17895/ices.advice.27202632`; an earlier version of this review incorrectly
identified key 19342 and advice 3843, which refer to the distinct **southern**
stock. The 2025 advice recommends catches of 1,501 tonnes in 2026.

## Assessment method and recovered sources

The 2025 advice describes a survey-biomass trend assessment using a constant
harvest-rate rule. It combines commercial catches, catch-length frequencies,
and the genetically split Malin Shelf Herring Acoustic Survey (MSHAS). Its
2024 input index is 12,557 tonnes; the reported target harvest rate is 0.28,
Itrigger 14,711 tonnes, biomass safeguard 0.854, and precautionary multiplier
0.5. The advice rounds displayed values and calculates with unrounded inputs.
This is not an age-structured population assessment and its index is not an
independently estimated absolute SSB surface. Relative catch/index is a
harvest-rate indicator, not instantaneous fishing mortality.

The detailed 2025 and 2026 HAWG reports provide historical catch numbers and
catch weights by age, along with MSHAS numbers-at-age and biological summaries.
These age tables provide source information, but the adopted harvest-rate rule
is based on aggregate survey biomass and is not an age-specific likelihood.
Commercial catch is apportioned geographically and includes an unquantified
proportion of other herring populations; the survey is genetically split.

Source interpretation requires care:

- HAWG 2026 Table 4.2.2 prints zero catch numbers at each age in 2022, although
  the aggregate series reports 51 tonnes; its Table 4.2.3 calls the 2024
  zero-weight row “no catches” despite 2,157 tonnes of aggregate catch. HAWG
  2025 instead marks the corresponding age data unavailable. Preserve source
  values and flag these contradictions; do not assert that they establish zero
  removals or manufacture an age allocation.
- HAWG 2026 Table 4.4.1.1 labels weight in kilograms, but its 2025 row matches
  grams in Table 4.3.1. Do not silently apply the kilogram label to that row.
- For 2023, HAWG 2025 and 2026 Table 4.3.2 show an index of 22.463 thousand
  tonnes, whereas the ICES SAG export for both assessment keys gives 23,117
  tonnes, as does the final 2025 advice Table 8. These are not rounding
  differences. The adopted 2025 advice is the source for its published index;
  the conflicting working-group age table is not silently treated as the
  decomposition of that value. The cause of the discrepancy is unresolved.

The 2026 HAWG PDF, 2025 HAWG PDF, 2025 advice, SAG XML and associated pages are
cached locally with file hashes and URLs. The report tables and source export
are being checked before import; this review alone does not claim that they are
already represented in the canonical database.

## Translation status

The survey-index rule does not estimate population N/F/M-at-age, recruitment,
or an age-structured population model. A comparison that presents those
quantities as official estimates would be misleading. A tinyAM translation
and comparison dashboard are not supported by this assessment alone.

## Sources

- [HAWG 2026 report](https://doi.org/10.17895/ices.pub.32625759), section 4
  and Tables 4.2.2–4.4.3.1 and 4.6.1.
- [HAWG 2025 report, version 5](https://ices-library.figshare.com/articles/report/Herring_Assessment_Working_Group_for_the_Area_South_of_62_N_HAWG_/28389008), section 4 and associated tables.
- [ICES 2025 6.aN advice](https://doi.org/10.17895/ices.advice.27202632).
- [ICES 2025 6.aN SAG source data, key 19370](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=19370).
- [ICES 2026 6.aN SAG source data, key 21308](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=21308).
- [ICES 2026 advice record 4063](https://asd.ices.dk/viewAdvice/4063), marked Invalid.
