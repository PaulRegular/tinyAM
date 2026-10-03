# North Sea autumn-spawning herring: 2026 source review

Status: canonical record ices_herring_north_sea_2026 imported; inputs, outputs and assumptions remain partial.
Charbonneau identifier: ICES-HAWG_ NS-IV 3a,7d_Clupea_harengus.
Stock: her.27.3a47d, North Sea, Skagerrak/Kattegat and eastern English Channel.

## Current source trail

Official detailed report: ICES 2026 HAWG, 8:16, DOI 10.17895/ices.pub.31424315. The complete 892-page final report is cached from its DTU author repository:
https://backend.orbit.dtu.dk/ws/portalfiles/portal/440308463/HAWG_2026_-_Full_Report.pdf

The report identifies 2025 SSB near 1.68 million tonnes and F ages 2–6 near 0.191. These provide accepted-run checks, not substitutes for native time series.

Current official reproducible repository:
https://github.com/ices-advice/2026_her.27.3a47d
Pinned inspected revision: d2690dac28a0b44f1735193cf8c7c87fe8fcbbbe.
The tree and key scripts are cached in source_cache/ices_herring_north_sea_2026/.
The older ICES-dk/wg_HAWG repository contains historical NSAS artifacts and is not used as the current-run source.

The TAF model.R executes separate single-fleet and multifleet scripts. model_sf.R explicitly names NSAS_HAWG2026_sf. A direct lookup under stockassessment.org user3 returned 404; this does not establish that the fit is unavailable under another owner or via TAF output archives. The repository contains scripts but no native RData objects in its tracked tree. Its DATA.bib still labels some inputs 2021/1947–2020, whereas the scripts are current; verify actual native files and terminal years rather than treating the bibliography as a numerical source.

Next: establish which single-/multifleet outputs define accepted historical assessment versus fleet-wise forecast; search TAF output/data archives and other SAM owners; inspect the report inventory and framework; recover complete inputs and outputs without mixing historical objects or advice projections.

## Accepted-model inventory from report and pinned TAF code

The historical assessment is the single-fleet FLSAM model, with one combined catch fleet. The multifleet model supports fleet-wise forecast calculations and must not be substituted for the reported single-fleet historical outputs. Report Table 2.6.2.7 gives years 1947–2026, ages 0–8 winter rings, plus group 8+, Fbar ages 2–6, and eight fitted index streams: HERAS, IBTS-Q1, IBTS0, IBTS-Q3, LAI-ORSH, LAI-BUN, LAI-CNS and LAI-SNS. Ages are winter-ring classes and must be labeled as such, rather than silently interpreted as chronological ages.

The code trims HERAS to ages 1–8+, IBTS-Q1 to age 1 and IBTS-Q3 to ages 0–5. This resolves a mismatch with abbreviated prose listing only IBTS-Q3 ages 2–5. LAI indices have type partial (FLSAM control fleet type 6), with separate regional spawning-component information, shared q and logP process parameters; they are not abundance-at-age observations. All four LAI sampling fractions are explicitly 0.67 in data_construct_input.R. Other survey timings must be retrieved from fleet.txt or verified configuration evidence.

Table 2.6.1.1 gives fitted coverage: LAI 1973–2025; Q1 1984–2026; Q3 1998–2025; HERAS age 1 1997–2025 and older ages 1989–2025; IBTS0 1992–2026. The report contains numerical catch, M, maturity, weights, all survey inputs, N/F outputs and summary uncertainty in tables 2.6.1.2–14 and 2.6.2.1–7. Some adjacent narrative paragraphs retain prior-year values; use explicitly labeled final numerical tables and crosscheck them against the current accepted run.

The TAF model_sf.R adds 0.02 to the constructed SMS-2023 M surface before fitting. Verify whether the reported input M matrix includes this addition before importing it. Do not add it twice or silently omit it.

## Report-table staging and configuration checks

The source-specific staging extraction now retains 6,131 report-native values in source_cache/ices_herring_north_sea_2026/report_tables_raw.csv. Tables 2.6.1.2–6 cover 1947–2025; the N and F surfaces in Tables 2.6.2.3–4 cover 1947–2026. These are staging values, not canonical imports. Missing catch entries in 1978–1979 are printed as dots and remain missing; survey -1 sentinels are retained in staging for explicit interpretation during import. Table 2.6.1.6 is captioned catch weight but its values and the accompanying narrative identify catch numbers; do not apply the weight variable based on that caption alone.

The printed current control object confirms eight unique catch F states: ages 7 and 8 share a state. N process variance is separate for recruitment age 0 and shared across ages 1–8. F innovation variance is shared over ages 0–1, 2–5 and 6–8. The pinned single-fleet configuration explicitly specifies cor.F = 2 and an AR observation structure for IBTS-Q3; other fitted index observation flags are ID. No power-law catchability parameters are configured. HERAS q is shared across ages 1–2 and across ages 3–8; the four larval streams share one q parameter. Three logP variance parameters are printed. Their statistical interpretation still requires FLSAM implementation or framework documentation, rather than inference from those integer keys.

The 2026 F surface must be retained as reported but distinguished from a year with complete catch data: catch and annual biological input tables end in 2025, while Q1 and IBTS0 observations extend into 2026. Final biological surfaces used for 2026 and their terminal-year construction are not yet verified. The report-derived M surface must also be checked against the explicit +0.02 adjustment in model_sf.R before its input status can be established.

Public TAF access check (2026-10-02): https://taf.ices.dk redirects to the general ICES homepage; https://taf.ices.dk/app/stock returns 404. The official ices-tools-prod/icesTAF package still constructs https://taf.ices.dk/api/Artifacts (R/taf_api.R and R/get.artifacts.R, master inspected on this date). A request for year=2026 and stock=her.27.3a47d also returns 404 without authentication. This is a limitation of this source trail, not proof that the native assessment object is unavailable everywhere. No credentials were requested or used. The portal landing response and source-specific table extraction script are cached locally. Report tables remain a valid fallback; earlier-run TAF artifacts must not be used to fill current numerical gaps.

The wrapped summary Table 2.6.2.6 has now been extracted for all 80 years, 1947–2026, including published lower/upper intervals, to report_summary_raw.csv. Recruitment equals the reported age-0 N exactly in every year. Mean F over ages 2–6 agrees with reported Fbar within 0.0006, allowing rounding of four-significant-digit table entries. The checks use the report values directly and do not overwrite them with recomputed quantities.

The source data_construct_input.R explicitly deletes catch observations for 1978–1979 (fishery closure), agreeing with the report dots. Preserve those missing cells; zero catch in other years is a distinct published value, although its likelihood treatment needs implementation verification. Larval tables include 1972, while Table 2.6.1.1 states fitted coverage starts in 1973; the unresolved inclusion rule must be checked rather than silently importing 1972 as a fitted observation.

LAI time-window labels are explicit: SNS columns 0/1/2 mean 16–31 December, 1–15 January and 16–31 January; CNS 0/1/2/3 mean 1–15 September, 16–30 September, 1–15 October and 16–31 October; BUN and ORSH 0/1 mean 1–15 September and 16–30 September. These belong in the season dimension, with no fish age. Their configured common model timing is 0.67 and is distinct from actual sampling-window dates. Native LAI numerical units still require clarification.

Report section 2.4.3 confirms the SMS-2023 update and additive +0.02 offset. It describes three-year moving averages outside the latest SMS period, while preceding historical-method prose describes five-year means. Final numerical inputs, rather than either generic prose description, must determine the accepted M surface and terminal-year extension.
