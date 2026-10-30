# Golden redfish, Iceland and adjacent areas (ICES reg.27.561214)

Review date: 2026-10-30

## Unit and latest assessment

The current ICES unit is Golden redfish (*Sebastes norvegicus*) in subareas 5,
6, 12, and 14. The latest assessment is the 2026 advice assessment (ICES
assessment key 22441; advice record 5102), based on data through 2025. This is a
different management unit from Northeast Arctic golden redfish, `reg.27.1-2`,
which is already represented in the database. The shared species name does not
make the units interchangeable.

The MFRI 2026 production report identifies the model as SAM and gives a
simulation period of 1966–2026 with ages 6–25+, where age 25 is the plus group.
The 2023 MFRI assessment report and the ICES WKBNORTH 2023 benchmark report are listed as background references. Current-run inputs, outputs and model settings below come from the 2026 production report and tables; no historical values are substituted.

## Published input and biology coverage

| Series | Published coverage and scale | Assessment interpretation |
|---|---|---|
| Catch-at-age | Ages 6–25+, 1995–2025 used by the age-based assessment; numbers in thousands. Catch-at-age and total landings are available from 1966, but the report says age data before 1995 are not used apart from pooled age-length keys applied in selected historical years. | Icelandic commercial length samples and age-length keys provide age allocation; estimated age catches are scaled to total landings across the stock areas. These are model catch-at-age inputs, not raw fish counts. |
| Autumn survey | Published IS-SMH age-disaggregated index in numbers, ages 6–25+, 1996–2025. The report says the autumn survey was not conducted in 2011 and that the missing autumn index was replaced by the mean of 2010 and 2012; the published 2011 age-table row should be treated as imputed, not a direct survey observation. | The production report describes autumn length and age data as model inputs. The table values are not reported as a physical abundance unit. |
| Spring survey | Published IS-SMB age-disaggregated index in numbers, ages 6–25+, 1985–2025. | The production report says the spring survey enters SAM as one aggregate biomass series because age data are limited. The at-age table is retained as supporting source material, not treated as a confirmed SAM input. |
| Catch weights | Ages 6–25+, 1966–2025; grams. | Commercial catch weights are reported separately from stock weights. |
| Stock weights | Ages 6–25+, 1966–2025; grams. | The report describes filling occasional missing ages from the alternative source and smoothing stock weights with a two-year moving average. |
| Maturity | Ages 6–25+, 1986–2025; proportion mature. | The report derives age maturity from a fixed length ogive and age-length keys, then smooths over two years below age 15 and four years for ages 15+. |
| Natural mortality | M = 0.05 for ages 6–24 and M = 0.1 for age 25+. | Fixed in the current assessment; M is not an estimated surface. |

The spring and autumn surveys are timed in March and October, respectively.
The reviewed sources do not give the exact within-year fractions or all
observation-likelihood weights.

## Official outputs and screening decision

The MFRI source tables report N-at-age in thousands for ages 6–25+ from
1966–2026 and F-at-age through 2025. The 2026 N row is a forecast/state beyond
the terminal observation year, not another historical observation. The ICES
2026 source table reports age-6 recruitment in thousands, SSB in tonnes, and F
for ages 9–19; aggregate intervals are available. The reviewed at-age tables
do not provide N/F standard errors.

This current unit is screened **eligible** for source curation: the detailed
2026 sources provide catch-at-age, survey age information, weights, maturity,
fixed M, and comparable official outputs. Eligibility does not mean the entire
SAM likelihood is reconstructable. A run-matched native SAM object or data
archive was not located in the reviewed sources. Exact likelihood weights,
observation-error settings, and survey timing fractions remain unresolved.

The MFRI tables mix observed historical series and forecast/state rows. Any
later import must preserve their source type, units, coverage, and forecast
status and must not treat the published spring at-age table as a SAM input
without further confirmation. No canonical stock, assessment, input, or output
rows have been added, and no tinyAM translation has been attempted.

## Sources

- ICES, [2026 advice assessment key 22441](https://asd.ices.dk/viewAdvice/5102).
- ICES, [2026 source data](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22441) and [assessment graphs](https://standardgraphs.ices.dk/ViewCharts.aspx?key=22441).
- MFRI, [2026 production assessment report](https://www.hafogvatn.is/static/extras/images/5_red_2026_1_techreport_en.html), [tables](https://www.hafogvatn.is/static/extras/images/5_red_2026_1_tables_en_v2.html), and [advice](https://www.hafogvatn.is/static/extras/images/5_red_2026_1_advice_en.html).
- MFRI, [2023 golden-redfish assessment report](https://www.hafogvatn.is/static/extras/images/04-goldenredfish_tr1388153.pdf).
- ICES, [WKBNORTH 2023 benchmark report](https://doi.org/10.17895/ices.pub.22304638).
- ICES, [Stock Information Database](https://sid.ices.dk/Default.aspx).

The reviewed sources and local checksums are recorded in the ignored cache at
`source_cache/discovery_ices_golden_redfish_561214_2026/sources.csv`. Candidate
screening details are in `candidates.csv`; the canonical database remains
unchanged.
