# Greenland cod identity and source curation

**Review date:** 2026-11-03

**Inventory:** ICES Stock Information Database snapshot retrieved 2026-10-28
**Scope:** Three distinct Greenland cod advice units; N-WISC source curation completed as a partial database import.

## Identity and assessment boundaries

ICES lists West Greenland offshore spawning cod (cod.21.1.osc), West Greenland inshore spawning cod (cod.21.1.isc), and East Greenland-Iceland offshore spawning cod (cod.21.27.1.14) as separate advice units. The West Greenland inshore advice unit contains two separately assessed components: northern inshore cod (cod.21.1a-c, N-WISC) and southern inshore cod (cod.21.1d-f, S-WISC). Keep the inshore advice unit as one stock record and retain component identity in the assessment scope and output region.

The Greenland cod units are distinct from the ICES cod stocks already represented in the canonical database, which cover the North Sea/West of Scotland and Northeast Arctic.

## Screening and curation status

- **West Greenland offshore (cod.21.1.osc): eligible, source curation pending.** The 2026 report describes an analytical SAM assessment over 2000-2025, commercial catch-at-age, age-specific bottom-trawl survey indices, catch and stock weights-at-age, maturity fixed between years, and natural mortality fixed at 0.2 for all ages. It identifies NWWG26_WOSC_v1 and reports recruitment, SSB, catch, and F ages 5-8 with uncertainty intervals. Exact input arrays, age coverage, and survey timing still need source checks.
- **West Greenland inshore (cod.21.1.isc): partially curated.** The 2026 report retains the analytical SAM assessment for N-WISC (cod.21.1a-c); it says the updated S-WISC model (cod.21.1d-f) was not used for advice and the stock receives Category 5 catch-based advice. This import records only N-WISC summary outputs. It does not treat S-WISC as a separate stock or add model results that were not used for official advice.
- **East Greenland-Iceland offshore (cod.21.27.1.14): not eligible for the age-based comparison batch at present.** The 2026 report says the SAM model failed to fit observed catches for some ages (up to 92%) and was not a robust basis for advice. ICES used Category 5 catch-based advice; no comparable official population or F outputs were identified.

## West Greenland inshore run and data checks

The inshore chapter reports fitted years 2000-2025 and names the SAM run NWWG25_N-WISC. Table 15.3.1 reports age-2 recruitment, SSB, total catch, and Fbar for ages 5-8. The historical recruitment, SSB, and Fbar rows include 95% confidence limits. Catch is blank for 2001 because the report says no commercial samples were available to split that year's catch by stock. The report states that M is fixed at 0.2 per year for all ages, maturity is fixed between years, and catch-at-age is calculated annually from commercial catches. It describes annual multimesh gillnet CPUE data for N-WISC in area 1B for 2002-2024, but the chapter does not provide the numeric age-specific input arrays, exact model age range, or fractional-year survey timing. No values were imputed for these missing inputs.

The public ICES SAG source-data page for assessment key 20951 is labelled assessment year 2025. Its latest reported values differ from Table 15.3.1 in the 2026 report (2025 recruitment 6795 and SSB 9242 on the SAG page, versus recruitment 7294 and SSB 10112 in the 2026 report). Key 20951 is therefore retained as a separate earlier published summary and was not mixed into the 2026 record.

A stockassessment::fitfromweb() retrieval was attempted for NWWG25_N-WISC, using its public run/model.RData URL. The host reset the connection from this environment, so the model object and its data arrays could not be cached or inspected. The source record is partial: it contains only the report's verified historical summary outputs and documented assumptions. There are no canonical age-specific input rows, N-at-age outputs, F-at-age outputs, or estimated-M outputs.

## Cached sources

The official [ICES 2026 NWWG report](https://doi.org/10.17895/ices.pub.32414577), chapter 15, is cached at source_cache/discovery_ices_greenland_cod_2026/nwwg_2026_northwest_greenland_inshore_cod.pdf. The local sources.csv manifest records the URL, retrieval date, file size, and SHA-256. The derived assessment_summary_2026.csv transcribes only Table 15.3.1, pages 204-205, and is also recorded in the manifest.

The public [ICES SAG summary for key 20951](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=20951) is cited above only to explain why the older assessment values were not used. It is not the source of the imported 2026 values.

The 2026 advice pages and other Greenland cod chapter PDFs remain cached as described in the source manifest. No tinyAM translation has been started.
