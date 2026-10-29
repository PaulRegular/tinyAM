# Celtic Sea herring: 2026 source screening

ICES stock `her.27.irls` is the Celtic Sea herring assessment unit covering Division 7.a south of 52°30′N and divisions 7.g, 7.h, 7.j and 7.k. The 2026 HAWG update is an official age-based ASAP assessment with terminal catch year 2025. Its assessment year runs from 1 April through 31 March, and the report cautions that winter-ring ages differ by one year from calendar age for this autumn/winter spawning stock.

## Identity and boundaries

No canonical stock has authority code `her.27.irls`. Existing records for `her.27.6aN`, `her.27.6aS7bc` and the historical combined `her.27.6a7bc` cover different assessment units and areas; they are not aliases for this stock. The 2026 HAWG report notes that Celtic Sea herring mix with Irish Sea herring (`her.27.nirs`), with the extent unquantified. This is a biological caveat, not evidence that the management assessment units are identical. The report also notes uncertainty about whether area 7.h belongs to the stock area.

## Inputs and assumptions

The 2026 HAWG report provides numerical inputs in the ASAP tables:

- Table 6.6.1.1 reports catch numbers-at-age (winter-ring ages 1–9+) for 1958–2025 and catch and spawning-stock weights-at-age in kg. The assessment assumes no discards.
- Table 6.6.1.2 reports the acoustic survey index in numbers, age-specific numbers for ages 2–7, annual CV 0.5 and sample sizes for 2002–2025. The assessment excludes 2004 (no survey), 2017 (insufficient biological data) and 2024 (no mature fish). Survey tuning uses ages 2–7. Table 6.3.1.1 presents full survey abundance by age in millions; the ASAP input table is in thousands of fish.
- The model uses a time-invariant maturity ogive: 0.5 at winter-ring age 1 and 1 at ages 2–9+. Natural mortality is fixed at 0.767, 0.385, 0.356, 0.339, 0.319, 0.314, 0.307, 0.307 and 0.307 for ages 1–9+. The report says this vector follows the HAWG 2015 final procedure using an SMS run from 2011.
- The acoustic survey index month is October. The 2025 survey ran from 9 to 29 October. The ASAP report identifies the index as numbers, all-year CV 0.5, and ages 2–7 for tuning.
- Catch composition is fitted with year-specific sample sizes from the Irish sampling programme, down-weighted to five for 2015–2025. The assessment fixes fishery selectivity at age 3 and survey selectivity at ages 3–5, with other stated selectivity parameters estimated.

The assessment report contains the detailed numeric tables needed to document the principal catch, survey and biological inputs. The pinned public TAF repository snapshot contains its data-preparation, ASAP-run and output scripts. Its scripts identify the underlying Lowestoft input files and the 2026 ASAP DAT/executable, but those bootstrap data and model files are not included in the public snapshot. Do not treat the scripts alone as the missing raw bundle.

## Official outputs and screening decision

Table 6.6.1.4 reports annual catch, SSB, total biomass, Fbar at ages 2–5 and recruitment at age 1 through 2025. The TAF output script also writes N-at-age and F-at-age tables, but the generated numeric age-specific output files are not present in the public repository snapshot. The published annual series is sufficient to verify comparable population and mortality outputs for the candidate screen; it does not substitute for those absent age-specific files.

**Screened: eligible for source curation.** The report supplies catch-at-age, an age-specific acoustic index, weights, maturity, mortality, timing and comparable official annual outputs. This is an assessment eligibility decision only; it does not import or assert canonical database records. Record the absent raw bootstrap bundle, absent generated age-specific output files, age convention and documented mixing as limits to resolve during curation.

## Sources

- [HAWG 2026 full report, section 6 and Tables 6.6.1.1–6.6.1.4](https://doi.org/10.17895/ices.pub.31424315)
- [ICES SAG assessment graphs, key 21307](https://standardgraphs.ices.dk/ViewCharts.aspx?key=21307)
- [ICES Advice and Scenarios Database, Celtic Sea herring 2026](https://asd.ices.dk/viewAdvice/4064?handler=SubComm)
- [Pinned ICES-advice TAF snapshot](https://github.com/ices-advice/2026_her.27.irls/tree/cf29535caedd3b97d1912a1fd5b37be350e3d916/2026_her.27.irls_assessment)

Downloaded source files and SHA-256 checksums are recorded in `source_cache/discovery_ices_celtic_sea_herring_2026/sources.csv` (gitignored). This review is based on sources retrieved 2026-10-29.
