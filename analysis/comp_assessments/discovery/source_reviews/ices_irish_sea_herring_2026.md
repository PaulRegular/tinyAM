# Irish Sea herring: 2026 source screening

ICES stock `her.27.nirs` is herring in Division 7.a north of 52°30′N. The 2026 HAWG report describes the final production run as a SAM assessment with terminal year 2025. It uses catch data from 1980 onward; older catch estimates extend to 1961 but have documented quality concerns. This is a distinct current assessment unit from Celtic Sea herring (`her.27.irls`, south of 52°30′N) and the West of Ireland `her.27.6aS7bc` unit.

## Inputs and assumptions

The 2026 HAWG report provides detailed numeric assessment tables:

- Table 8.2.1.2 reports catch-at-age in thousands for winter-ring ages 1–8+ over 1980–2025; 2009 has no estimate. Sampling reports show the 2025 catch-age samples and seasonal coverage.
- Table 8.3.1 reports ISAS age-specific acoustic abundance in thousands for winter-ring ages 1–8+ over 1994–2025. Table 8.3.2 gives a second acoustic series, the Irish Sea Commercial Acoustic Survey (ISSS), as an aggregate SSB index from 2007–2025.
- Tables 8.5.1–8.5.2 provide age-specific weights in kg. Catch weights use all catch seasons; spawning-stock weights use Q3 catches. The 2009 values are replaced by the mean of 2008 and 2010.
- Table 8.5.3 provides annual maturity ogives at ages 1–8+ from Q3 catch samples, with 2009 replaced by the mean of 2008 and 2010.
- Table 8.5.4 gives fixed natural mortality for all years, borrowed from North Sea herring: 0.787, 0.380, 0.353, 0.335, 0.315, 0.311, 0.304 and 0.304 for ages 1–8+.
- ISAS is conducted in late August and September. The 2025 survey ran 26 August–19 September, with remaining transects completed 7–9 October. ISSS is normally conducted in late September or early October. Exact annual SAM observation-time fractions are not reported; a translation would need to document a date-based approximation.

The report refers to the 2026 stock annex for model settings. The latest `ices-advice` code snapshot located is the 2024 repository, revision `808b6c896474f7210c7f95ee4c7aef5046aaf491`, which predates the 2024 and 2025 configuration reviews described by HAWG. It is not used to supply the current numerical inputs. No 2026 `ices-taf` repository was found in the ICES TAF GitHub organization search.

## Official outputs and limitations

HAWG Table 8.7.1.1 reports annual age-1 recruitment, SSB, Fbar at ages 2–4 and catch with 95% intervals through 2025. The 2026 ICES SAG source data (key 21319) corroborate the assessment table; for 2025 they report recruitment 157,581 thousand, SSB 11,298 t and Fbar 0.342, with intervals. These outputs support a meaningful aggregate comparison. Numeric fitted N-at-age and F-at-age arrays were not recovered from the report or SAG source table.

HAWG documents substantial uncertainty in acoustic scaling: estimated survey catchabilities range from 2.7 to 4.9, and fitted biomass is below direct survey biomass estimates. The assessment treats the fishery and surveys as mixed-stock observations, including Celtic Sea herring; the mixing extent is not quantified and occurs across ages. Natural mortality is borrowed from North Sea herring rather than estimated specifically for the Irish Sea. These are major interpretive limitations for a comparison.

**Screened: eligible for source curation.** Numeric catch-at-age, an age-specific acoustic abundance series, biological inputs, mortality, survey seasons and comparable official SSB/recruitment/F outputs are available. Record the missing raw run bundle, age-specific fitted output arrays, approximate timing requirements, mixed-stock treatment and acoustic-scaling concern during any later import or translation.

## Sources

- [HAWG 2026 full report, section 8 and Tables 8.2.1.2, 8.3.1–8.3.2, 8.5.1–8.5.6, 8.7.1.1](https://doi.org/10.17895/ices.pub.31424315)
- [ICES 2026 SAG graphs, key 21319](https://standardgraphs.ices.dk/ViewCharts.aspx?key=21319)
- [ICES 2026 SAG source data, key 21319](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=21319)
- [ICES 2027 advice](https://doi.org/10.17895/ices.advice.32133760.v1)
- [Older 2024 ICES-advice code snapshot](https://github.com/ices-advice/2024_her.27.nirs/tree/808b6c896474f7210c7f95ee4c7aef5046aaf491/2024_her.27.nirs_assessment)

The report PDF is shared with the Celtic Sea screen and cached under `source_cache/discovery_ices_celtic_sea_herring_2026/`. The 2026 ICES SAG pages are cached in `source_cache/discovery_ices_irish_sea_herring_2026/sources.csv`, with retrieval date and SHA-256 checksums. This review is based on sources retrieved 2026-10-29.
