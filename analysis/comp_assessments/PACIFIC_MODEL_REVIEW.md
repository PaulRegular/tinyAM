# Review of the twelve new Pacific cases

Baseline: commit 86e44bf. This review checks the canonical inputs, importers, source reviews, translation code and numerical fits. Native prior equations were revisited directly for the M trials. It is not a claim of exhaustive independent verification of every published cell. Canonical observations are not tuned to match fitted outputs.

Each trial changes a documented assumption or precision treatment while retaining observations, years and ages. Normal package starts run first; one existing source-start construction is allowed after numerical failure. All failed and incomplete trials remain visible. Priors use external source assumptions, not fitted official M values.

The rejected northern rock sole q-prior trial expresses the index in units relative to the fixed source prior center 1.5, retaining original values in source_obs. Its log-index residuals are mathematically unchanged. This is an explicit unit conversion, not an adjustment using fitted q; the retained recipe and canonical observations do not change.

Median differences below are median absolute percent differences, with the existing annual/mean-within-age correlation definitions. No agreement threshold is optimized. Retained changes can have trade-offs, which are reported rather than hidden.

## Case-by-case review

| assessment_id | data_review | assumptions_and_fit_review | decision |
| --- | --- | --- | --- |
| afsc_cod_ebs_2024 | VAST is a number index in thousands, not biomass. Active ages and empirical fishery age-length allocation retained; uncovered length mass stays unallocated. Fitted weights/proxy maturity remain explicit. | Test shared annual F, correlated F and catch SD. Retain the baseline: the shared-level alternative improves F/SSB but weakens N/recruitment trends. Its dashboard remains available as an alternative. | Retain baseline; alternatives not adopted |
| afsc_cod_ai_2024 | Final November and preliminary September runs do not match. Original total catch/biomass inputs retained; final likelihood arrays/weights remain unverified. | Keep the source-data guard. Do not mix preliminary fitted states with final observations. | Blocked; no fit invented |
| afsc_pollock_ai_2024 | Original biomass, number proportions and fleet weights retained. Survey weights proxy unavailable population weights; later aggregate-only years are not given artificial ages. | Test common annual F, correlated F and the verified M log-penalty. Fixed-q normalization, BH settings and source population weights remain limitations. | Retain baseline; alternatives not adopted |
| afsc_atka_mackerel_bsai_2024 | Original kg weights, maturity, number compositions and biomass CVs retained. Actual composition weights are not replaced by raw otolith counts. Whole-million N rounding limits precision. | Keep baseline: all tested F/precision alternatives weaken important comparisons. Joint composition likelihood, q prior and source selectivity penalties remain different. | Retain baseline; alternatives not adopted |
| afsc_yellowfin_sole_ebs_2024 | Joint-sex observed compositions are converted using grams before retaining females. Malformed native tokens use verified observed report echoes, not fitted predictions. Covariates retain native centering. | Test correlated/common-level F and age-dependent catch SD. Failed numerical checks remain failures; different biomass/composition reference times and source composition likelihood remain material. | Retain baseline; alternatives not adopted |
| afsc_northern_rock_sole_ebs_2024 | Joint-sex conversion precedes female restriction; kg weights differ from yellowfin grams. Rounded-derived F stays approximate. Framework confirms the female log-M penalty; male-prior activation is not established. | Test correlated/common-level F, catch SD, random initial N and verified female M/q penalties. The joint M/q trial still fails (gradient 3.22e10, non-PD Hessian). Random N passes numerical checks but produces implausible abundance/q and is rejected. Source joint composition likelihood remains different. | Retain baseline; alternatives not adopted |
| afsc_alaska_plaice_bsai_2024 | Joint-sex conversions retain fixed source q calibration and computed growth weights. Only three fishery age years are available; absent compositions remain missing. | Test correlated/common-level F and catch SD; preserve the source fixed-M/q meanings. Sex pooling and omitted lengths constrain interpretation. | Retain baseline; alternatives not adopted |
| afsc_arrowtooth_flounder_bsai_2022 | ATF-only catch, sex-specific fixed biology and active surveys retained. Numerical length/age likelihood inputs, observed age-length key and survey times remain unverified. | Keep the source-data guard. Fitted N/F cannot manufacture catch or survey observations. | Blocked; no fit invented |
| afsc_pop_bsai_2024 | Original grams, ages 40+, two surveys and active biomass likelihoods retained. Fitted maturity is a labelled plug-in; source M density includes its documented bias adjustment. | Test common/correlated F and the source M density. Incomplete/interrupted results are not numerical successes. | Retain baseline; alternatives not adopted |
| joint_hake_coastal_2025 | Native percentage compositions normalize jointly; original weights/fecundity, age-2 index timing and posterior-median definitions retained. No fitted M is used to change observations. | Retain one shared log-M effect with native fixed center/SD. M is estimated and its density counted once. A BH alternative gives essentially the same comparisons and a weakly estimated curve, so constant recruitment is retained. Laplace integration differs from native Bayesian reporting; bounds are omitted. | Adopt source_M_prior |
| dfo_herring_sog_2025 | Native survey is aggregate spawning biomass with no observed survey ages. Five-OM ensemble, model-year timing, both-sex SSB and posterior point conventions remain distinct. | Keep the observation-model guard. A fishery-only fit or fabricated survey composition would not represent the source abundance information. | Blocked; no fit invented |
| dfo_pop_bc_2023 | Original sex/area biomass conversions and zero-evaluation fixed biology retained. Catch pools only complete three-area years. Historical official N/F/R/SSB and native reproductive-output definition remain unresolved. | Retain a moment-matched positive approximation to the native normal M prior. It is not the same density or bound treatment. Only M has a usable historical official comparison. | Adopt source_M_prior |

## Trial diagnostics

| assessment_id | candidate | status | optimizer_code | objective | max_gradient | elapsed_seconds | sdreport_created | positive_definite_hessian | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| afsc_cod_ebs_2024 | age_dependent_catch_sd | converged |  0 |   244.1958 | 9.567040e-03 |   31.98873 | TRUE | TRUE | not retained |
| afsc_cod_ebs_2024 | annual_F_selectivity | converged |  0 |   280.7143 | 1.134746e-03 |   10.71482 | TRUE | TRUE | not retained |
| afsc_cod_ebs_2024 | correlated_F | converged |  0 |   245.0561 | 5.525401e-04 |   14.05842 | TRUE | TRUE | not retained |
| afsc_pollock_ai_2024 | annual_F_selectivity | not_converged |  1 |   775.6687 | 5.583191e+07 |  512.66172 | TRUE | FALSE | not retained |
| afsc_pollock_ai_2024 | correlated_F | not_converged |  1 |   574.3588 | 2.597941e-05 |  250.50232 | TRUE | TRUE | not retained |
| afsc_pollock_ai_2024 | source_M_prior | not_converged |  1 |   628.7145 | 1.185797e-04 |  126.91757 | TRUE | TRUE | not retained |
| afsc_atka_mackerel_bsai_2024 | age_dependent_catch_sd | converged |  0 |   612.9946 | 6.386351e-03 |  441.92145 | TRUE | TRUE | not retained |
| afsc_atka_mackerel_bsai_2024 | annual_F_selectivity | converged |  0 |   564.4640 | 1.056299e-03 |   30.99187 | TRUE | TRUE | not retained |
| afsc_atka_mackerel_bsai_2024 | correlated_F | converged |  0 |   588.4734 | 1.433639e-03 |   32.55083 | TRUE | TRUE | not retained |
| afsc_yellowfin_sole_ebs_2024 | age_dependent_catch_sd | not_converged |  1 |  8298.9689 | 1.311823e+12 | 6464.10445 | TRUE | FALSE | not retained |
| afsc_yellowfin_sole_ebs_2024 | annual_F_selectivity | no completed result | NA | NA | NA | NA | NA | NA | not retained |
| afsc_yellowfin_sole_ebs_2024 | correlated_F | not_converged |  1 |  8951.0075 | 3.698039e+11 |  875.91019 | TRUE | FALSE | not retained |
| afsc_yellowfin_sole_ebs_2024 | regularized_initial_abundance | no completed result | NA | NA | NA | NA | NA | NA | not retained |
| afsc_northern_rock_sole_ebs_2024 | age_dependent_catch_sd | not_converged |  1 | 88483.5041 | 1.074307e+05 | 1978.90062 | TRUE | FALSE | not retained |
| afsc_northern_rock_sole_ebs_2024 | annual_F_selectivity | not_converged |  0 |  3638.0419 | 2.502611e-02 |  528.08116 | TRUE | FALSE | not retained |
| afsc_northern_rock_sole_ebs_2024 | correlated_F | not_converged |  0 |  3929.6209 | 4.367527e-02 |  559.13213 | TRUE | FALSE | not retained |
| afsc_northern_rock_sole_ebs_2024 | regularized_initial_abundance | converged |  0 |  4055.9862 | 1.398970e-03 |  319.49118 | TRUE | TRUE | not retained |
| afsc_northern_rock_sole_ebs_2024 | source_M_prior | not_converged |  1 |  4313.0128 | 2.711434e+03 |  384.24010 | TRUE | FALSE | not retained |
| afsc_northern_rock_sole_ebs_2024 | source_M_q_priors | not_converged |  1 |  5298.8668 | 3.224368e+10 |  463.89969 | TRUE | FALSE | not retained |
| afsc_alaska_plaice_bsai_2024 | age_dependent_catch_sd | not_converged |  1 |  1325.6070 | 1.954213e+04 |  169.81632 | TRUE | FALSE | not retained |
| afsc_alaska_plaice_bsai_2024 | annual_F_selectivity | not_converged |  1 |   626.0983 | 3.085989e+09 |   30.42684 | TRUE | FALSE | not retained |
| afsc_alaska_plaice_bsai_2024 | correlated_F | not_converged |  1 |   350.8768 | 6.850977e+05 |   24.71134 | TRUE | FALSE | not retained |
| afsc_pop_bsai_2024 | annual_F_selectivity | no completed result | NA | NA | NA | NA | NA | NA | not retained |
| afsc_pop_bsai_2024 | correlated_F | no completed result | NA | NA | NA | NA | NA | NA | not retained |
| afsc_pop_bsai_2024 | source_M_prior | no completed result | NA | NA | NA | NA | NA | NA | not retained |
| joint_hake_coastal_2025 | age_dependent_catch_sd | not_converged |  0 |  1117.0053 | 1.675599e-02 | 1275.46467 | TRUE | TRUE | not retained |
| joint_hake_coastal_2025 | annual_F_selectivity | converged |  0 |  1032.2037 | 7.589024e-04 |  129.75908 | TRUE | TRUE | not retained |
| joint_hake_coastal_2025 | correlated_F | not_converged |  1 |  1131.0278 | 1.131393e-04 |  164.29316 | TRUE | FALSE | not retained |
| joint_hake_coastal_2025 | source_M_prior | converged |  0 |  1185.7297 | 8.678038e-04 |  140.30939 | TRUE | TRUE | retained |
| joint_hake_coastal_2025 | source_M_prior_BH | converged |  0 |  1183.7688 | 1.112179e-04 |  261.16326 | TRUE | TRUE | not retained |
| dfo_pop_bc_2023 | annual_F_selectivity | converged |  0 |  2194.9376 | 2.774667e-03 |   32.84153 | TRUE | TRUE | not retained |
| dfo_pop_bc_2023 | correlated_F | converged |  0 |  2178.2325 | 1.194460e-03 |   60.72606 | TRUE | TRUE | not retained |
| dfo_pop_bc_2023 | source_M_prior | converged |  0 |  2206.6083 | 1.640159e-03 |   17.48687 | TRUE | TRUE | retained |

## Retained changes: every comparable core metric

Hake's M treatment improves all available core scale/trend comparisons. BC POP improves M only; other historical core comparisons remain unavailable or non-equivalent. The EBS cod shared-level F alternative improves F/SSB but weakens N/recruitment trends, so it is not adopted automatically.

| Case | Metric | Median absolute % before | After | r before | r after | Status |
| --- | --- | --- | --- | --- | --- | --- |
| joint_hake_coastal_2025 | M |  100 | 21.5 | — | — | matched |
| joint_hake_coastal_2025 | N |   71 | 54.8 | 0.883 | 0.886 | approximate |
| joint_hake_coastal_2025 | abundance | 71.8 | 42.8 | 0.785 | 0.904 | approximate |
| joint_hake_coastal_2025 | recruitment | 77.8 | 55.2 | 0.959 | 0.969 | matched |
| joint_hake_coastal_2025 | ssb | 69.6 | 44.8 | 0.619 | 0.751 | approximate |
| dfo_pop_bc_2023 | M |  100 | 22.5 | — | — | matched |

M prior implementation: the existing fixed-SD Gaussian mean effect represents one shared log-M level, with one density contribution across all years/ages. Hake matches the source quadratic log-M penalty; BC POP matches the mean and variance of its natural-scale normal prior with a lognormal approximation. Neither fixes M to an official estimate. Native hard bounds and integration/point-estimate conventions remain different. Mean-only M uncertainty is available for the shared log-M effect; the current full M surface is reported without an automatically propagated interval.

The full baseline/candidate comparison table is results/pacific_model_review.csv. Fits, exact specifications, initialization histories and rejected trials remain in gitignored results/cache/model_review_20261016/. Selected dashboards are refreshed in the usual per-stock cache. Missing statistics are not passes; summaries of failed fits are exploratory.

## Validation and boundaries

Canonical numerical inputs and outputs are unchanged. The northern rock sole source note now distinguishes the verified female M penalty from unresolved male-prior activation. Core package mathematics, APIs and likelihood code are unchanged.

The normal-start/conditional-retry policy applies separately to each documented model specification. A failed initialization does not justify changing observations. Ten-minute wall limits were used for expensive later trials; missing completed results remain incomplete, not convergence failures inferred from their last gradient.

The full package suite and all assessment-analysis files were run. Current validation results and existing analysis-fixture failures are recorded in PACIFIC_EXPANSION.md. Targeted tests cover prior scales, one shared density contribution, unchanged source observations, conditional starts, common-definition SSB and projection-safe alternative F design.

## Related existing-case corrections

North Sea saithe now estimates free initial older-age abundance, with normal package starts. Every comparable core median difference and trend improves; observation conversion and plus-group comparisons are unchanged.

| Metric | Median absolute % before | After | r before | r after |
| --- | --- | --- | --- | --- |
| F | 15.7 | 14.1 | 0.807 | 0.841 |
| N | 13.4 | 12.9 | 0.938 | 0.951 |
| abundance | 7.65 | 6.52 | 0.911 | 0.935 |
| recruitment | 11.3 | 9.79 | 0.926 | 0.94 |
| ssb | 10.8 | 9.37 | 0.97 | 0.981 |

Baltic sprat now applies the source 40% pre-spawning fraction for F and M. The derived official SSB uses official N/F/M and shared biology, with no reconstructed uncertainty. This corrects the comparison definition; its statistics are not presented as a like-for-like improvement over the previous January-1 comparison. Native reported SSB remains separate.
