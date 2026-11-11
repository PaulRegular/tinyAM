# Curated assessment refresh

The review compared recipes from tinyAM `9a7eb2b6` using the same pinned SAIO snapshot, `54474f12d039a757595ac08236c01b7b53905a98`. Thirty of 35 curated recipes translate and prepare. Five retain their documented source blockers. For each fitted trial, the baseline and revision used the same starting-value policy and a 20-minute limit per fit. Comparisons use finite historical rows with matching units and definitions. Mixed results stay in the local review outputs.

## Retained revisions

| Assessment | Change | Result |
|---|---|---|
| Northern Shelf haddock | Use latent initial recruitment, three `q_power()` groups, native catch and survey SD groups, separate catch weights and Fbar ages 2–4. | Median absolute differences fell from 3.47% for SSB, 17.80% for recruitment and 6.41% for F-at-age to below 0.00008% each. The abundance-weighted Fbar comparison at ages 2–4 fell from 2.67% to 0.000049%. Both fits converged. The revised default-start objective differed from native SAM by `7.41e-9`. |
| North Sea herring | Share F at ages 7–8, use the three source F innovation SD groups and add age AR(1) errors for IBTS-Q3. | Median absolute differences fell from 4.49% to 1.09% for F-at-age, 2.71% to 0.91% for Fbar, 1.51% to 0.80% for SSB and 2.49% to 0.61% for recruitment. Terminal SSB difference moved from +4.20% to −0.38%. Both fits converged. |
| Southern Gulf spring herring | Scope the CPUE catchability random walk with `for_survey()` and remove two indicator columns. | Same 1990–2021 process years and first-year anchor. Both fits converged, with identical objectives, catch and survey predictions and official-comparison rows. |
| Gulf of Alaska pollock | Scope the Shelikof environmental effect with `for_survey()`. | Both fits converged, with identical objectives and predictions. The ADF&G random walk keeps its original years and anchor. A separate double-logistic F trial failed the Hessian check. The reference F-at-age series is reconstructed from rounded official numbers-at-age. |
| Northern cod | Scope the Smith Sound availability interaction to the fall RV survey and keep the shared baseline catchability. | Both fits converged. Objectives and predictions matched across 741 catch and 591 survey rows. The no-Smith-Sound and juvenile design checks also passed. |

Haddock comparisons retain 55 SSB and recruitment years, 486 F-at-age cells and 54 Fbar years. North Sea herring retains 79 years and 711 F-at-age cells. The three scoped rewrites preserve comparison coverage. Southern Gulf spring herring's mature biomass comparison uses January 1. Haddock's Fbar comparison weights by abundance, while published SAM Fbar is an arithmetic mean.

## Case review

| Curated recipe | Disposition |
|---|---|
| `afsc_alaska_plaice_bsai_2024` | Retained. Fishery ages are available for three years and the native selectivity is sex-specific. A separable age curve with annual intensity would be poorly supported by the available catches. |
| `afsc_arrowtooth_flounder_bsai_2022` | Blocked. Missing numerical compositions, age–length key and survey timing. |
| `afsc_atka_mackerel_bsai_2024` | Mixed. The age-10+ F block converged. Median F-at-age difference rose from 21.60% to 23.30% while terminal F agreement improved. |
| `afsc_cod_ai_2024` | Blocked. The public native run predates the production fit and lacks verified matching compositions and weights. |
| `afsc_cod_ebs_2024` | Retained. The native length selectivity and composition likelihoods need more than age-based selectivity functions. |
| `afsc_cod_goa_2026` | Retained. The aggregate fit converged. The three-fleet trial with complete reconstructed grids failed with a nonfinite gradient. |
| `afsc_northern_rock_sole_ebs_2024` | Retained. The baseline failed numerical checks. A trial with shared female F from age 17 converged but estimated M near zero and had large discrepancies from official SSB and F. |
| `afsc_pollock_ai_2024` | Retained. Baseline and age-8+ F sharing trial both failed numerical checks. Published F-at-age is unavailable and Fbar differs in definition. |
| `afsc_pollock_ebs_2024` | Retained. The native covariance of the biomass indices and age compositions does not establish an age AR(1) error model for reconstructed ages. |
| `afsc_pollock_goa_2024` | Revised. Shelikof environmental effect scoped. Scoping the ADF&G random walk would shift its historical anchor from 1990 to 2000. |
| `afsc_pop_bsai_2024` | Retained. A tensor smooth would approximate bicubic selectivity and also replace the F process. The earlier 65-year by 38-age fit reached its time limit. |
| `afsc_yellowfin_sole_ebs_2024` | Retained. The age-17+ F block trial failed numerical checks from both default and source starts. The unchanged recipe also failed its source-start retry, so neither fit provides an eligible comparison. |
| `dfo_canary_bc_2022` | Retained. No age composition for non-trawl catch, so available data do not support separate fleet F. |
| `dfo_cod_2j3kl_2025` | Revised. The RV-only Smith Sound interaction preserves shared baseline catchability, with identical fitted catch and survey predictions. The age-0 sensitivity remains a separate model. |
| `dfo_cod_4t4vn_2019` | Retained. Age-2 catches and the native composition likelihood are missing. Fixed M increment variance alone would still omit the source's initial-M priors. |
| `dfo_herring_4tvn_spring_2024` | Revised. CPUE scope preserves the fitted likelihood. Separate gear observations cannot reproduce the native total-catch and composition likelihoods. |
| `dfo_herring_sog_2025` | Blocked. The survey reports spawning biomass without age samples. |
| `dfo_petrale_bc_2024` | Blocked. Native compositions, biology, timing and final model configuration are unresolved. |
| `dfo_pop_bc_2023` | Blocked. Comparable coastwide official population outputs are unavailable. |
| `ices_bluewhiting_northeast_atlantic_2026` | Mixed. Source F sharing and correlated increments improved median F-at-age and SSB agreement but worsened recruitment and terminal SSB. The age AR(1) approximation also worsened several series. |
| `ices_cod_northeast_arctic_2026` | Mixed. Age-3 F SD worsened terminal F and recruitment. Approximate survey AR(1) improved median SSB difference from 7.15% to 3.83% and terminal difference from +30.24% to +2.06%, while SSB trend correlation fell from 0.963 to 0.950. SAM uses several correlation groups and shares some across surveys. |
| `ices_cod_north_sea_2025` | Retained. The native model couples three substocks. Fleet-specific F does not represent substock dynamics. |
| `ices_haddock_iceland_2025` | Retained. Native process sharing and survey timing lack a verified mapping for the same assessment. |
| `ices_haddock_north_sea_2026` | Revised. Native SAM likelihood and reporting compared using the pinned fit and default-start refit. |
| `ices_herring_north_sea_2026` | Revised. F sharing and IBTS-Q3 correlation improved agreement. The larval spawning-component likelihoods remain outside the translation. |
| `ices_herring_norwegian_spring_2025` | Retained. Native F sharing and observation covariance are unresolved for this assessment. |
| `ices_herring_western_baltic_2026` | Retained. Shared F at ages 7–8 is outside the modeled age-5+ range. The segmented stock–recruit and aggregate index likelihoods remain unmatched. |
| `ices_norway_pout_north_sea_2026_benchmark` | Retained. Annual F SD proxy for quarterly source groups failed numerical checks. |
| `ices_plaice_north_sea_2026` | Retained. Correlated F increments increased the median F-at-age difference from 9.28% to 12.56% and SSB difference from 7.84% to 15.12%. |
| `ices_saithe_north_sea_2026` | Retained. Baseline and documented F SD grouping both failed numerical checks. |
| `ices_sprat_baltic_2026` | Mixed. The shared age-7–8 F state and correlated increments lowered median SSB and F-at-age differences, but SSB and biomass trend correlations worsened. |
| `ices_whiting_north_sea_2026` | Mixed. Native catch and survival SD groups and approximate survey AR(1) were tried separately. Catch SD groups worsened historical SSB agreement. AR(1) improved median SSB agreement and worsened terminal SSB, F and recruitment. Native survey plus ages and correlation groups remain unmatched. |
| `joint_hake_coastal_2025` | Retained. Shared selectivity increments do not imply equal absolute F. A supplied SD multiplier cannot reproduce the aggregate survey's additive extra variance. |
| `nefsc_atlantic_mackerel_2018` | Baseline failed. A separable F mean with an age-6+ plateau converged, with median official gaps of 24.05% in SSB and 39.98% in F-at-age. The aggregate egg index and native composition likelihood are still absent, so there is no matched baseline or native-equivalence result. |
| `nefsc_summer_flounder_2018` | Retained. Both the aggregate model and a trial with four disjoint native catch streams failed numerical checks. |

## Validation

The comparative suite passed 135 of 137 scripts, including all curated preparation checks. Both failures reproduced on unchanged tinyAM `9a7eb2b6` with the same SAIO snapshot. The direct reduced-age western Baltic herring fit failed the convergence check. The preliminary Rockall haddock recipe found duplicate weight-at-age rows. Targeted checks passed for all five revisions, and each revised dashboard rendered with the official assessment alongside it.

The local review outputs contain the paired fits, warnings, comparison metrics, and dashboards. The retained recipes run through the standard assessment workflow using the pinned SAIO snapshot.
