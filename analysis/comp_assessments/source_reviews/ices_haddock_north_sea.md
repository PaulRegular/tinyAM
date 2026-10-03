# Northern Shelf (North Sea) haddock: 2026 source review

Status: canonical record ices_haddock_north_sea_2026 imported. Inputs, outputs and assumptions remain partial pending benchmark review, survey-unit clarification, and remaining model-detail checks.

Charbonneau–Keith identifier: ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus.
The source stock covers North Sea, West of Scotland and Skagerrak (had.27.46a20).

## Accepted-run evidence

Detailed 2026 chapter: https://ndownloader.figshare.com/files/66104291
Working-group report: https://doi.org/10.17895/ices.pub.32676345
Native final run: https://stockassessment.org/datadisk/stockassessment/userdirs/user3/NShaddock_WGNSSK2026_Run1/run/model.RData

The baserun/model.RData object in this directory is historical (ends 2021), so it is not used. The run/model.RData object spans 1972–2026. All 648 annual recruitment, SSB, Fbar and TSB estimates and interval endpoints for 1972–2025 match the final-model Table 8.3.6 within published rounding. Published intervals use exp(log estimate ± 2 log-SE), rather than the exact normal 0.975 quantile.

## Inventory

SAM 0.12.0, RemoteSha 1cc464b80f6f. Model ages 0–8+. One total-catch fleet (486 observations, 1972–2025). Two fitted delta-GAM survey indices: Q1 ages 1–8+ (352 observations, 1983–2026, sampling fraction 0.125); Q3+Q4 ages 0–8+ (315 observations, 1991–2025, fraction 0.75). Supplied survey relative weights are material inputs and are now retained explicitly; catch weights are missing by design. Native inputs, final object and detailed chapter are cached in source_cache/ices_haddock_north_sea_2026/.

Final configuration has age-correlated F random-walk increments (corFlag=2), one F innovation variance, separate N process variances for recruitment, intermediate ages and plus group, and independent observation errors (obsCorStruct=ID). Parameter-sharing keys remain to be reconciled; native survey-weight semantics and age mapping are verified.

## Forecast distinction

The native fitted 2026 SSB is 699,809.6 tonnes. The advice forecast uses different weights and maturity, giving 667,758 tonnes. Recruitment is resampled from 2000–2025 for advice, rather than using the fitted 2026 recruitment. Report section 8.6 explains these differences; preserve fitted and forecast quantities separately.

Remaining: inspect the 2022 benchmark (published 2023), stock annex and 2025 survey/model review; establish original survey units; reconcile the remaining parameter-sharing assumptions.

## Age surfaces and observation-weight checks

All 486 historical F-at-age and 495 N-at-age values match tables 8.3.4–8.3.5 within report rounding. The fitted object contains N through 2026 and F states through 2026, but historical reported F ends in 2025.

All 667 native survey weights match 1/log(1+CV^2) using the first declared number of age columns from each raw CV matrix. This includes the leading value (1) as the first modeled age and omits the trailing column, exactly as the accepted read.ices() code specifies. The raw weights are now recorded in inputs.csv as relative_precision_weight rows by survey/year/age. The matching log_index_sd rows retain the derived relative SD factors, 1/sqrt(weight); catch-fleet weights are missing.

## Resolved reader and plus-group construction

The saved SAM revision's read.ices() selects the first declared number of age columns from each matrix, without removing a survey-effort column on this file path. The published src/datascript.R assigns weight = 1/log(CV^2+1). The first CV column is therefore used for the first model age, and the trailing extra column is omitted. This native convention matches all 667 accepted-run weights and is retained in the database. The pinned reader and data script are cached as SAM_reading.R and SAM_datascript.R.

The stock script constructs age 8+ before fitting: sums catch/landings/discard numbers; weights stock weight, maturity and M by catch-number composition (terminal composition reused for the added model year); weights catch, landings and discard mean weights by the corresponding component numbers. These final input matrices are available directly in fit$data, so there is no need to reconstruct them from report values. land.frac is derived from aggregated landings numbers divided by total catch numbers. The original discard-number file cited by this script is also cached.

## Canonical import

6,906 inputs, 2,373 outputs and 25 assumptions. Catch, landings fractions, stock/catch/landings/discard weights, maturity, M, both surveys, all 667 native relative precision weights and all 667 derived relative log-SD factors are represented. The SD factor is 1/sqrt(weight); it is a relative input, not the fitted final observation SD. All 648 summary/interval and 981 N/F surface crosschecks passed. Structural validation and scripts/database/010_validate_north_sea_haddock.R passed. No core model changes.

## Spawning definition resolved

The accepted object's propF and propM matrices each contain 495 cells, spanning 1972–2026 and ages 0–8+, with all values zero. Both source matrices are now retained as biological inputs. These values define beginning-of-year spawning biomass; they are separate from survey sampling times and do not justify substituting a nonzero spawning fraction. The source-specific validator checks every cell against the cached native exports. Coverage is now 6,906 inputs, 2,373 outputs and 25 assumptions. Survey-unit clarification, framework follow-up, parameter-sharing and state-uncertainty questions remain unresolved; statuses remain partial.
## Catchability power correction

Inspection of the accepted fit and exact SAM revision revealed an error in the previous assumption text: this model DOES estimate density-dependent survey powers. keyQpow activates separate powers at Q1 age 1 and Q3+Q4 ages 0 and 1. The other survey ages use power 1. The exact prediction code multiplies log abundance after survey-time survival by exp(logQpow), then adds logFpar: I=q*(N at survey time)^power. Thus age-specific q coefficients at affected ages cannot be interpreted as simple fractions caught.

The incorrect no-power assumption has been replaced. Seventeen q coefficients and three powers are exported from native fixed effects and their covariance, with natural-scale delta-method SEs and 95% log-Wald intervals. Source-specific checks confirm parameter-to-age mapping, values and SE scales. The accepted object is unchanged. Coverage is 6,906 inputs, 2,373 outputs and 25 assumptions; all statuses remain partial. This finding supersedes the earlier description of catchability without density dependence.

Code source:
https://github.com/fishfollower/SAM/blob/1cc464b80f6f/stockassessment/inst/include/SAM/predobs.hpp
## Initial-state semantics resolved

Both native initN/initF vectors are empty and initState=0. The accepted SAM revision's n.hpp and f.hpp evaluate ordinary process transitions from the second modeled year onward, with no explicit density for the first logN/logF state. Those first-year states remain latent random effects in the Laplace fit. A broad first-state normal prior is added only when calculating observation residuals; it must not be described as a regular assessment prior. The native configuration checks pass. This resolves the earlier initial-state gap without changing the fitted object or model mathematics.

Exact source files:
https://github.com/fishfollower/SAM/blob/1cc464b80f6f/stockassessment/inst/include/SAM/n.hpp
https://github.com/fishfollower/SAM/blob/1cc464b80f6f/stockassessment/inst/include/SAM/f.hpp
