# Western Baltic cod, ICES 2025

## Source selection

The WGBFAS 2025 report documents the assessment run `WBCod2025`. The report and
the linked SAM input/configuration files provide detailed observations,
biological inputs, fitted population outputs and settings. The source stock is
ICES cod in subdivisions 22–24 (`cod.27.22-24`). The selected run has catch
observations through 2024 and fitted population estimates through 2025. The
2025 age-0 population is not reported; age-1 recruitment is reported for 2025.

Sources:

- [ICES WGBFAS 2025 report, DOI 10.17895/ices.pub.29099786](https://doi.org/10.17895/ices.pub.29099786), western Baltic cod section 2.3, especially Tables 2.3.15 and 2.3.17–2.3.26.
- [ICES WKBALTCOD2 benchmark, DOI 10.17895/ices.pub.4984](https://doi.org/10.17895/ices.pub.4984).
- [ICES advice, DOI 10.17895/ices.advice.27202560](https://doi.org/10.17895/ices.advice.27202560).
- [Official ICES TAF repository](https://github.com/ices-advice/2025_cod.27.22-24), pinned locally at commit `093cc29bf0e4d47b89adfb209348adc946041b3f`.
- [SAM data and configuration for `WBCod2025`](https://stockassessment.org/datadisk/stockassessment/userdirs/user3/WBCod2025/), cached locally with the fitted object returned by `stockassessment::fitfromweb("WBCod2025")`.

The report identifies the assessment as uncertain because mortality data are
incomplete. It states that indices are used for SSB and recruitment and that
harvest rate is used instead of F for advice. The fitted SAM object still
contains its estimated F surface; it is retained as a model output but must not
be read as an equally reliable advice quantity.

## Inputs and settings

The fitted SAM object covers 1985–2025 and ages 0–7. Total catch-at-age covers
1985–2024, ages 1–7+, in thousand fish; the source model omits unavailable
catch-age cells. The observation likelihood weight is 1 through 2019 and 0.1
for catch observations in 2020–2024. The public report does not explain this
weight change, so the native weight values are preserved without converting
them to an assumed SD or precision multiplier.

The three survey series are BITS Q4 (ages 0–4+, 1999–2024), BITS Q1 (ages
1–4+, 1996–2025) and FEJUCS pound net (age 0, 2011–2024; 2012 is missing).
The report describes the BITS observations as CPUE indices. It does not give a
physical scale multiplier, so native SAM index units are retained. The fitted
sampling times are 0.8477816, 0.1653494 and 0.85, respectively. BITS Q1 has
correlated age residuals in the SAM configuration; tinyAM uses independent
survey residuals.

Natural mortality is fixed at 1985–2024 values from Table 2.3.21: ages 0–7+
are 1.318, 0.598, 0.411, 0.324, 0.274, 0.241, 0.218 and 0.201 per year.
The SAM input surface also supplies 2025 values. Stock, catch, landings and
discard weights, maturity, landings fractions and pre-spawning F/M fractions
are extracted from the native fit data and retained without interpolation.
The report defines SSB at the start of the year; pre-spawning F and M fractions
are zero. Maturity is a fixed ogive in the model, documented in Table 2.3.20.

SAM uses an age-1 recruitment state, an age-0 to age-7 population surface,
age-specific log-F states with ages 5–7 sharing one state, and a shared F
process variance. The SAM configuration has no F residual correlation. The
N-process variance is separate for age 0 and ages 1–7. Catch observation
variance is separate for age 1 versus ages 2–7. BITS Q4 and Q1 have the
survey-age variance-sharing groups in the model configuration; Q1 age residuals
are correlated. Catchability is estimated by survey and age where the SAM
configuration includes a parameter.

## Outputs and checks

The canonical outputs retain native SAM N-at-age, the fitted F states,
recruitment at age 1, SSB, total stock biomass, Fbar for ages 3–5, catch/index
predictions and q where available. Conditional random-effect uncertainty for
N-at-age and F is stored as natural-scale delta-method SEs with 95% log-Wald
intervals; the source report does not publish these age-specific intervals.
The report's Tables 2.3.25–2.3.26 provide grouped total mortality and harvest
rate, respectively; those rounded table values are retained as separate native
outputs. Total abundance is derived from N-at-age, with uncertainty left
blank.

The cached SAM fit converges with optimizer code 0, positive-definite Hessian,
and maximum absolute fixed-effect gradient below `4e-12`. Its 2024–2025 age-1
recruitment, SSB, Fbar and TSB reproduce the rounded values in Table 2.3.23.
The fitted object includes a 2025 age-0 state that is not reported in Table
2.3.24; that unavailable table cell is not imported as a source estimate.

The tinyAM translation uses an IID recruitment process and cohort deviations,
IID F states, fixed source M, source survey timings and survey-by-age q groups.
It cannot reproduce SAM's joint survey residual correlation, native catch
observation weights, or the exact initial-state and process parameterization.
tinyAM estimates F at age 0, whereas the source has no F state at age 0. These
are recorded as approximations; no observations are adjusted to improve
agreement.

The source PDF, extracted text, raw `.dat` inputs, configuration, selected TAF
revision and fitted SAM object are cached under the gitignored
`source_cache/ices_cod_western_baltic_2025/` directory. The local manifest
records retrieval dates, checksums and source URLs.
