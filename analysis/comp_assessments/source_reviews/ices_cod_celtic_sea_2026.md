# Celtic Sea cod, 2026

## Official assessment and sources

This record represents the 2026 WGCSE production assessment of Celtic Sea cod
(`cod.27.7e-k`, ICES divisions 7.e-k), used for 2027 advice. The fitted SAM
object is `Cod_CS_WG2026`, covering 1980-2025 at ages 1-7+; age 7 is the plus
group. The run uses stockassessment 0.12.0 (RemoteSha `1cc464b80f6f`), has
optimizer convergence code 0, objective 279.6445, and a positive-definite
Hessian.

- [WGCSE 2026 report and stock chapter](https://doi.org/10.17895/ices.pub.32625759.v1), especially cod section 7 and Tables 7.2-7.9.
- [WKCELTIC 2020 benchmark report](https://doi.org/10.17895/ices.pub.5983), which the 2026 report identifies as the source of the current model settings.
- [ICES 2027 advice](https://asd.ices.dk/viewAdvice/6150) and [standard assessment graphs](https://standardgraphs.ices.dk/ViewCharts.aspx?key=22676).
- [Official fitted SAM object](https://stockassessment.org/datadisk/stockassessment/userdirs/user3/Cod_CS_WG2026/run/model.RData).

The two reports and native SAM object are cached under
`analysis/comp_assessments/source_cache/ck_ices_wgcse_morhua_cs7e_k/assessment_2026/`.
The ignored `source_manifest.csv` records each file's URL, retrieval date and
SHA-256 checksum. The 2020 cod stock-annex PDF was requested from the ICES
publication archive but returned HTTP 401; no replacement copy was recovered.
The native fitted object and benchmark report still provide the numerical
inputs and fitted configuration used here.

## Inputs and assumptions

- The SAM object contains total catch-at-age, including landings and discards,
  for ages 1-7+ from 1980-2025. Catch at ages 1 and 2 is missing through 2003;
  those cells remain missing in the database.
- The model fits two tuning series. The French OTDEF Q2+3+4 commercial LPUE
  index is aggregate, with one observation per year from 2002-2025. It is
  retained as a relative biomass index on its native kg/hour scale, but is not
  translated to tinyAM because tinyAM currently fits age-specific indices.
- The combined Irish IGFS and French EVHOE Q4 VAST series is age-specific for
  ages 1-3 over 2003-2025. Its model timing is 0.855. Table 7.6 also displays
  age 4, but the final SAM object uses only ages 1-3. The report shows 0.000
  for age 3 in 2025; SAM represents that log-scale observation as missing, so
  it is retained as `NA` and omitted from the tinyAM fit.
- The SAM object supplies annual stock weights, catch weights, landings
  weights, discard weights and maturity. The report identifies the numerical
  input tables as Tables 7.2-7.4. Missing discard weights remain missing.
- Natural mortality is fixed at 0.501, 0.330, 0.264, 0.233, 0.211, 0.211 and
  0.211 per year at ages 1-7+. The assessment does not estimate M. Maturity is
  0, 0.54, 0.93 and then 1 at ages 1-7+, and is constant across years. The
  model's F and M fractions before spawning are both zero.
- Recruitment is at age 1, with age 1 in its own N-process variance group and
  ages 2-7 sharing another group. F has six age states, with age 7+ sharing
  the age-6 state. Its process increments have AR(1) correlation across ages
  and one shared process-variance group. Fbar is ages 2-5.
- The commercial index has one catchability group. VAST has separate q
  parameters at ages 1, 2 and 3. Catch likelihood SD groups are age 1, age 2,
  and ages 3-7; VAST uses one SD group and an age-correlated AR residual
  structure. Both fitted indices and catch use lognormal likelihoods.

## Outputs and tinyAM translation

The database includes SAM N-at-age, F-at-age, SSB, total stock biomass,
age-1 recruitment, Fbar and fitted catch/index predictions. The total-abundance
series is derived as the sum of N-at-age over ages 1-7+; no uncertainty is
assigned to that derived total. For N and F, natural-scale SEs are delta-method
approximations from the conditional log-state uncertainty in the SAM fit
object. SSB, total biomass, recruitment and Fbar retain the source 95%
confidence limits.

The tinyAM recipe fits 1980-2025 at ages 1-7+, using catch and VAST. It uses
source fixed M, source stock weights and maturity, a zero spawning time, three
catch SD groups, age-specific q for VAST, one VAST SD group, an IID N process,
random-walk recruitment and a correlated-random-walk F process. The aggregate
commercial LPUE likelihood and VAST's age-correlated observation errors are
not represented. The source's exact age-state sharing, F covariance and N
initial-state treatment are also simplified. The initial fit uses package
defaults and no source-derived parameter starts.

The final SAM output in 2025 is SSB 98.94 tonnes, total stock biomass 251.39
tonnes, age-1 recruitment 183.20 thousand fish and Fbar ages 2-5 of 1.4348
per year. These values agree with the rounded values in WGCSE Table 7.9.
Agreement with tinyAM is reviewed through the dashboard and comparison
summary; matching definitions do not imply matching fitted models.
