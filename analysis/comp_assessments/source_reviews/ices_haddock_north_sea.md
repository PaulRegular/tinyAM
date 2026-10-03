# Northern Shelf (North Sea) haddock: 2026 source review

Status: canonical record ices_haddock_north_sea_2026 imported. Inputs, outputs and assumptions remain partial pending benchmark/review, unit and initial-state checks.

Charbonneau–Keith identifier: ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus.
The source stock covers North Sea, West of Scotland and Skagerrak (had.27.46a20).

## Accepted-run evidence

Detailed 2026 chapter: https://ndownloader.figshare.com/files/66104291
Working-group report: https://doi.org/10.17895/ices.pub.32676345
Native final run: https://stockassessment.org/datadisk/stockassessment/userdirs/user3/NShaddock_WGNSSK2026_Run1/run/model.RData

The baserun/model.RData object in this directory is historical (ends 2021), so it is not used. The run/model.RData object spans 1972–2026. All 648 annual recruitment, SSB, Fbar and TSB estimates and interval endpoints for 1972–2025 match the final-model Table 8.3.6 within published rounding. Published intervals use exp(log estimate ± 2 log-SE), rather than the exact normal 0.975 quantile.

## Inventory

SAM 0.12.0, RemoteSha 1cc464b80f6f. Model ages 0–8+. One total-catch fleet (486 observations, 1972–2025). Two fitted delta-GAM survey indices: Q1 ages 1–8+ (352 observations, 1983–2026, sampling fraction 0.125); Q3+Q4 ages 0–8+ (315 observations, 1991–2025, fraction 0.75). Supplied survey relative weights are material inputs and must be retained; catch weights are missing by design. Native inputs, final object and detailed chapter are cached in source_cache/ices_haddock_north_sea_2026/.

Final configuration has age-correlated F random-walk increments (corFlag=2), one F innovation variance, separate N process variances for recruitment, intermediate ages and plus group, and independent observation errors (obsCorStruct=ID). Parameter-sharing keys and survey weight semantics require implementation review before full assumption extraction.

## Forecast distinction

The native fitted 2026 SSB is 699,809.6 tonnes. The advice forecast uses different weights and maturity, giving 667,758 tonnes. Recruitment is resampled from 2000–2025 for advice, rather than using the fitted 2026 recruitment. Report section 8.6 explains these differences; preserve fitted and forecast quantities separately.

Next: verify age-specific surfaces; inspect the 2022 benchmark (published 2023), stock annex and 2025 survey/model review; establish original survey units and relative-weight transformation; import all material biological, catch-component and observation-weight inputs.

## Age surfaces and observation-weight checks

All 486 historical F-at-age and 495 N-at-age values match tables 8.3.4–8.3.5 within report rounding. The fitted object contains N through 2026 and F states through 2026, but historical reported F ends in 2025.

All 667 native survey weights match 1/log(1+CV^2) using each raw CV-file row's leading value as the first modeled age and omitting the trailing value. The leading value is 1 in the inspected rows. Do not silently discard the leading column as a legacy flag: doing so shifts the weights relative to those used in the accepted object. The exact file/reader convention remains to be reviewed before declaring the observation-weight input complete. Preserve fit$data$weight and its year/fleet/age alignment as authoritative for the accepted run.

## Resolved reader and plus-group construction

The saved SAM revision's read.ices() reads matrix files and selects the first declared number of age columns, with no survey-effort-column removal on this matrix-file path. The published src/datascript.R reads the two CV matrices through that function and assigns weight = 1/log(CV^2+1). This explains the exact native weight alignment; it is reproduced by all 667 row-level checks. Code sources are cached as SAM_reading.R (revision 1cc464b80f6f) and SAM_datascript.R. The apparent leading/trailing-column discrepancy must not be silently corrected in the database.

The stock script constructs age 8+ before fitting: sums catch/landings/discard numbers; weights stock weight, maturity and M by catch-number composition (terminal composition reused for the added model year); weights catch, landings and discard mean weights by the corresponding component numbers. These final input matrices are available directly in fit$data, so there is no need to reconstruct them from report values. land.frac is derived from aggregated landings numbers divided by total catch numbers. The original discard-number file cited by this script is also cached.

## Canonical import

5,249 inputs, 2,353 outputs and 24 assumptions. Catch, landings fractions, stock/catch/landings/discard weights, maturity, M, both surveys and all 667 relative log-SD factors are represented. The SD factor is 1/sqrt(native precision weight), with the source weight retained in each row note; this is a relative input, not the fitted final observation SD. All 648 summary/interval and 981 N/F surface crosschecks passed. Structural validation and scripts/010_validate_north_sea_haddock.R passed. No core model changes.

## Spawning definition resolved

The accepted object's propF and propM matrices each contain 495 cells, spanning 1972–2026 and ages 0–8+, with all values zero. Both source matrices are now retained as biological inputs. These values define beginning-of-year spawning biomass; they are separate from survey sampling times and do not justify substituting a nonzero spawning fraction. The source-specific validator checks every cell against the cached native exports. Coverage is now 6,239 inputs, 2,353 outputs and 25 assumptions. Survey-unit clarification, initial-state semantics, framework follow-up, q and state uncertainty remain unresolved; statuses remain partial.
## Catchability power correction

Inspection of the accepted fit and exact SAM revision revealed an error in the previous assumption text: this model DOES estimate density-dependent survey powers. keyQpow activates separate powers at Q1 age 1 and Q3+Q4 ages 0 and 1. The other survey ages use power 1. The exact prediction code multiplies log abundance after survey-time survival by exp(logQpow), then adds logFpar: I=q*(N at survey time)^power. Thus age-specific q coefficients at affected ages cannot be interpreted as simple fractions caught.

The incorrect no-power assumption has been replaced. Seventeen q coefficients and three powers are exported from native fixed effects and their covariance, with natural-scale delta-method SEs and 95% log-Wald intervals. Source-specific checks confirm parameter-to-age mapping, values and SE scales. The accepted object is unchanged. Coverage is 6,239 inputs, 2,373 outputs and 25 assumptions; all statuses remain partial. This finding supersedes the earlier description of catchability without density dependence.

Code source:
https://github.com/fishfollower/SAM/blob/1cc464b80f6f/stockassessment/inst/include/SAM/predobs.hpp
## Initial-state semantics resolved

Both native initN/initF vectors are empty and initState=0. The accepted SAM revision's n.hpp and f.hpp evaluate ordinary process transitions from the second modeled year onward, with no explicit density for the first logN/logF state. Those first-year states remain latent random effects in the Laplace fit. A broad first-state normal prior is added only when calculating observation residuals; it must not be described as a regular assessment prior. The native configuration checks pass. This resolves the earlier initial-state gap without changing the fitted object or model mathematics.

Exact source files:
https://github.com/fishfollower/SAM/blob/1cc464b80f6f/stockassessment/inst/include/SAM/n.hpp
https://github.com/fishfollower/SAM/blob/1cc464b80f6f/stockassessment/inst/include/SAM/f.hpp
