# Icelandic haddock (had.27.5a)

## Assessment record

The canonical database represents the 2025 benchmark assessment as the current
detailed record. MFRI published a newer 2026 assessment report, so that update
is recorded separately as summary-only. Its age-specific matrices are not
copied from the 2025 assessment.

The 2025 MFRI tables provide catch numbers-at-age, catch weights, stock weights,
maturity, spring and autumn survey indices-at-age, and estimated numbers and
fishing mortality-at-age. The assessment covers ages 1–12+, with age 12 as the
plus group. Catch and biological input tables end in 2024; population estimates
extend through 2025. The database therefore records 1979–2024 as fitted input
years and marks 2025 estimates as beyond the detailed terminal input year.

The 2026 MFRI report says the model uses the 2025 benchmark framework, starts in
1979, tracks ages 1–12+, fixes M at 0.2, and allows selectivity to vary over
time. The 2026 table endpoint was unavailable during this review. Its report
contains summary figures, but age-specific numerical surfaces were not
recovered. No values were digitized from figures and no 2025 values were
carried forward as 2026 data.

## Inputs and assumptions

- Catch-at-age is reported directly in thousands of fish, ages 1–12+.
- Catch weights are from commercial samples. Stock weights and maturity are
  from the March survey; the source reports that pre-1985 stock-weight and
  maturity vectors use the 1985 values.
- The two survey series are IS-SMB (spring, March; 1985–2025) and IS-SMH
  (autumn, October; 1995–2024). The official tables label the indices as
  numbers but do not state a physical unit. Their values are retained on the
  native scale with an unresolved-unit label.
- The database records 0.20 for March and 0.80 for October as within-year
  timing approximations based on the documented survey months. These are not
  claimed to be the exact timing fractions used by SAM.
- M is fixed at 0.2 per year for every age. The report gives pre-spawning
  mortality fractions of 0.4 for F and 0.3 for M.
- The summary tables provide aggregate SSB, reference biomass and age-1
  recruitment, along with low/high values. Age-specific output uncertainty was
  not recovered.

## tinyAM translation

The translation fits 1979–2024, ages 1–12+, with one catch stream, the two
native survey-index series, supplied M, an IID N process and an AR1 F process.
Survey q and observation error are estimated by series and age or survey,
respectively. These choices approximate SAM's state-space model; tinyAM does
not reproduce the source model's exact process covariance, residual
correlation, or all parameter sharing.

The source's native SSB includes partial mortality before spawning. tinyAM
does not apply those fractions. The comparison summary instead labels its
common-definition SSB calculation: accepted N is converted to mature biomass
using the shared translated weight and maturity inputs. The accepted
length-based 45+ cm reference biomass is retained as a native output and is
not compared to tinyAM's total biomass.

## Sources reviewed

- [MFRI 2025 assessment tables](https://dt.hafogvatn.is/astand/2025/2_HAD_en.html)
- [MFRI 2025 technical report](https://www.hafogvatn.is/static/extras/images/02-had_2025_techreport_en.html)
- [MFRI 2026 technical report](https://www.hafogvatn.is/static/extras/images/2_had_2026_1_techreport_en.html)
- [ICES WKICEGAD 2025 benchmark report](https://doi.org/10.17895/ices.pub.28444499.v1)

The numerical MFRI table exports used for the 2025 record are cached locally
under the gitignored source_cache/iceland_haddock_2025/ directory. Their
filenames and access notes are listed in that directory's README.
