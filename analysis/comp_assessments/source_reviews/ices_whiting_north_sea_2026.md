# North Sea whiting source review

## Assessment and source

The represented assessment is the accepted 2026 SAM assessment for whiting
(*Merlangius merlangus*) in ICES Subarea 4 and Division 7.d, stock code
`whg.27.47d`. The assessment is described in the [WGNSSK 2026 report](https://doi.org/10.17895/ices.pub.32676345),
with outputs also available through [ICES Stock Assessment Graphs, key 22488](https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22488).
The [WKBW 2026 benchmark](https://doi.org/10.17895/ices.pub.32019630) documents
recent model changes. No Charbonneau–Keith identifier was confirmed in the
current curated catalogue; the stock is linked here by its ICES code.

The accepted `NSwhiting_2026n` model and data objects, native input files, and
the report are cached under
`analysis/comp_assessments/source_cache/ices_whiting_north_sea_2026/`. The
cache is gitignored. The report SHA-256 is
`9d33846cf100dd8544ba30dbeffc7435a857a71854d318dbc7b383d668b03857`; the
model object SHA-256 is
`6bd2d93ad87b6cbbE6f479346d478796183a24befb6c2d87da97216b4ef5f6b8`; the
native data object SHA-256 is
`689003079c3fede21befdf735d0af68b20a26e631d9ed9214ba06994e1b3914c`.
The import script is
`scripts/database/049_import_ices_whiting_north_sea_2026.R`.

The April `NSwhiting_2026n` run is treated as accepted because its input data
and fitted object reproduce the published recruitment, SSB, Fbar and TSB
values for 1978–2025 in Table 22.22 to the table's rounding precision; its
calculated SSB also reproduces the reported 2026 intermediate-year value. A later June
`NSwhiting_2026n_updateIBTSQ1` variant is cached for reference, but its annual
fitted series do not match the published 2026 assessment table and it is not
used as the canonical run. The selected model reports optimizer convergence
code 0 (`relative convergence (4)`), objective 51.29096.

## Model structure and source data

The model spans ages 0–8+, with recruitment at age 0. Total catch numbers at
age cover 1978–2025. The model also uses total catch weight at age, annual
stock weight, smoothed maturity, natural-mortality observations, and two
survey indices. The native fitted object stores survey indices on an unscaled
relative index scale; the reports do not give a calibrated physical unit.

The IBTS Q1 survey covers ages 1–6+ in 1983–2026 and has sampling time 0.125.
The IBTS Q3 survey covers ages 0–6+ in 1991–2025 and has sampling time 0.625.
Catchability is separate by survey and age, with no density-dependent q power.
The two survey observation processes use lognormal errors with AR(1) age
correlation. Catch observations use lognormal errors with independent age
residuals. The fitted model's native relative precision weights and the
corresponding relative log-scale SD factors are retained in `inputs.csv`.

Stock and catch weights and maturity are used as known inputs. Stock weights
and maturity use a 6+ group, expanded over the model's ages 6–8. Maturity is
estimated from Q1 information for 1991–2021; the 1991 values are held for
earlier years and the 2021 values for later years. The WGSAM natural-mortality
surface is supplied through 2022, with age 6+ grouped; SAM estimates missing
later years using its mortality GMRF. The accepted model has random-walk N and
F processes, with age-specific variance sharing for N and AR(1) correlation
between F process increments across ages. Fbar is ages 2–5.

The native SAM data input contains `propF = 0` and `propM = 0`, which define
the source model's spawning-biomass timing. These are recorded separately from
the Q1 and Q3 survey timing fractions.

## Survey precision detail

The native Q1 and Q3 CV files contain a leading `1` followed by the six or
seven age-labelled CV values shown in the report. The accepted data script
uses `stockassessment::read.ices()` and takes the first declared number of
columns from those files. Consequently, the actual precision weights in
`fit$data$weight` include the leading `1` as the first modeled age and omit
the last age-labelled CV column. This matches the accepted fit's stored weights
and is how the fit is represented in the database. The age-labelled CV values
from the files are retained separately so the distinction is visible; they do
not replace the fitted weights in the tinyAM translation.

## Comparison limits

The report publishes point estimates and 95% intervals for recruitment, SSB,
Fbar, and TSB through 2025, with SSB alone for 2026. The native fit also
contains 2026 state surfaces. Those internal states are retained as
`native_model` outputs, while the summary table leaves unreported 2026
recruitment, Fbar, and TSB unavailable. The age-state output surfaces do not
have report-published uncertainty intervals.

The tinyAM translation is an illustrative approximation. It cannot reproduce
SAM's exact age-correlated F increments, N-process variance sharing, or M GMRF.
It will retain the reported data, q sharing, timing and observation precision
as far as tinyAM's existing model allows. These are model differences, not
corrections to the accepted SAM assessment.
