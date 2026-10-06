# Norway pout source review

## Assessment selected

The selected record is the February 2026 WKBWNG benchmark assessment for
Norway pout in ICES Subarea 4 and Division 3.a (stock code nop.27.3a4). The
report names the accepted final SESAM run NP_Sep2025_bench_final; the cached
fitted summary is under the related archive name NP_Sep2025_bench_final_Mfor.
It matches the report's final settings and the reported 2012 Q4 B_lim: the
cached SSB at 2012.75 is 30,766.56 tonnes, compared with 30,766.76 tonnes in
the report. Its quarterly fitted states span 1984 Q1–2025 Q4; catch and survey
observations extend through 2025 Q3. The fitted summary reports optimizer
convergence code 0, objective 1547.746, and the message “relative convergence
(4)”.

The benchmark is distinct from the 2025 production-advice assessment. WGNSSK
2026 states that no spring advice was issued for Norway pout. Accordingly the
2026 benchmark is the current detailed assessment represented here, but
is_applied is FALSE.

The authoritative benchmark report identifies the selected M-0 scenario:
annual-varying natural mortality at age and quarter from the North Sea SMS
model for 1984–2022, with 2022 values carried through 2025. The report says
this scenario was selected after comparison of fit, retrospective patterns,
and diagnostics. The native fitted object carries the quarterly M values used
over its historical input intervals. The separate M_forecast.csv file contains
forecast values, including 2026, and is not substituted for the fitted-history
input.

## Data recovered

The cached native archive files are in
source_cache/ices_norway_pout_north_sea_2026_benchmark/native/:

- sum.RData, the final fitted SESAM summary and the exact processed data
  retained by that run;
- obs.dat and aux_input.dat, broader source files that include rows not used in
  the selected fitted run;
- xsa.dat, the model configuration data;
- M_forecast.csv, kept for reference but not treated as fitted-history input.

The processed model object contains 1,092 observations: quarterly
catch-at-age and five survey index fleets. The importer uses those processed
rows rather than combining them with the broader candidate rows in obs.dat.
Catch-at-age is reported in millions of fish for ages
0–3+, with Q4 2025 absent. The five survey series are:

| Fleet | Series and coverage | Within-year timing |
|---|---|---:|
| 2 | IBTS Q1, ages 1–3+, 1984–2025 | 0.125 |
| 3 | EGFS Q3, ages 0–3+, 1998–2025 | 0.625 |
| 4 | SGFS Q3, age 0 in 1998–2012 and ages 1–3+ in 1998–2025 | 0.625 |
| 5 | IBTS other countries Q3, ages 0–1 and one grouped age 2–3+ index, 1999–2024 in the fitted object | 0.625 |
| 6 | SGFS Q3 age 0, 2013–2025 | 0.625 |

Timing is the midpoint of the model quarter. The native run uses a quarter
length of 0.25 year and observations halfway through each survey quarter.
Age 3 is the plus group. Maturity is constant at 0, 0.2, 1, and 1 for ages
0–3+. Stock weights are age- and quarter-specific but constant across years;
catch weights vary by age, year, and quarter. The model inputs treat discards
as equal to catch weights, consistent with the report's conclusion that
discards and bycatch are negligible.

Natural mortality is supplied as a time-varying input, not estimated by SESAM.
The input data object stores 668 age-quarter values over intervals from 1984
Q1 through 2025 Q3. The separate forecast file is retained in the local cache
but is not combined with those fitted-history inputs.

There is one unresolved source discrepancy: WKBWNG Section 2.5 lists the
IBTS other-countries Q3 index through 2025, while the selected fitted object
and cached obs.dat both end that series in 2024. The database retains only
the observations present in the fitted object and does not invent a 2025
index value. This series is flagged in the assumptions, and the input
completeness status remains partial until the difference is resolved.

The accepted run has a combined age 2–3+ observation in that IBTS series.
The canonical record identifies it as an age-group index rather than an
age-2-only index. tinyAM does not have a grouped-age index likelihood, so the
translation keeps the row in the database and excludes it from the fitted
tinyAM observation list.

## Model structure and available outputs

Recruitment is at age 0 and enters in Q3. N process variance is shared for
ages 0–2, with a separate variance for age 3+. F has a distinct age state at
each age and quarterly random-walk evolution; its process variance is shared
for ages 0–1 and 2–3, with AR(1) correlation across ages. Fbar is ages 1–2.
Catch observation variance is shared for ages 1–2, with separate age-0 and
age-3 groups. Survey variance and q-sharing follow Table 2.5 of the benchmark
report. EGFS and SGFS each share q between ages 2 and 3+; the other age keys
are separate within each survey fleet. The post-2013 SGFS age-0 index has a
separate q key. Selected F-process increments with large changes in
log(catch + 1) have their standard deviation multiplied by 100 when the change
exceeds 2.9.

The cached fit provides quarterly N-at-age and F-at-age, quarterly SSB and
Fbar with 95% intervals, and annual age-0 recruitment with 95% intervals.
Fitted age-specific N and F point estimates are available from the native
object; age-specific intervals are not recovered here. M-at-age is an applied
input surface and is therefore not described as an estimated output.

The native object is a compact, processed representation of the accepted
model inputs. The observation and auxiliary text files are also cached for
review. The benchmark PDF is publicly available at the source link below, but
the local Windows download client could not retrieve it because its HTTPS
credential negotiation failed. This does not prevent importing the native
model records; the report remains the authority for assessment selection and
model assumptions.

## Sources

- ICES WKBWNG 2026 report, especially Sections 2.3–2.6, Tables 2.1–2.5, and
  Figures 2.1–2.2:
  https://doi.org/10.17895/ices.pub.32019630
- Report PDF:
  https://ices-library.figshare.com/ndownloader/files/64519302
- WKBWNG working documents:
  https://ices-library.figshare.com/articles/report/Benchmark_Workshop_on_Whiting_Norway_Pout_and_Golden_Redfish_WKBWNG_/32019630
- Accepted final fitted run:
  https://stockassessment.org/datadisk/stockassessment/userdirs/user3/NP_Sep2025_bench_final_Mfor/run/sum.RData
- Accepted run input directory:
  https://stockassessment.org/datadisk/stockassessment/userdirs/user3/NP_Sep2025_bench_final_Mfor/data/
- 2025 production advice and Stock Assessment Graphs remain separate records:
  https://doi.org/10.17895/ices.advice.27202770
  https://standardgraphs.ices.dk/ViewSourceData.aspx?key=21230
