# Pacific expansion review

The database has 36 current detailed stock records. The original twelve-case Pacific expansion has four numerical passes, three non-converged fits, one interrupted run and four source blockers. Nine of those cases historically received optimization attempts; BC POP is now source-output blocked despite a numerical pass. Seven active expansion dashboards and the historical BC trial are cached locally. This falls short of twelve runnable translations; blocked guards are not counted as fitted recipes.

The subsequent BC Canary Rockfish case is labelled Additional below and does not replace or complete an original expansion case. Its female translation passes after one source-informed initialization retry, but SSB has mean absolute difference 27.9%, terminal difference -59.1% and correlation 0.295; estimated M is 36.5% lower. Its source review and dashboard document the approximations and unavailable N/F/recruitment comparisons. No observations were adjusted to improve agreement.

The initial expansion retained existing stock specifications. The subsequent source/model review updates hake and BC POP mortality treatment, North Sea saithe initial abundance and Baltic sprat spawning time; other existing comparisons remain retained results. New diagnostics and comparisons come from the documented default-first attempts. Source additions were committed before their translation/readiness records and pushed separately. No core model mathematics changed.

Review thresholds are mean absolute percent difference >20, absolute terminal mean percent difference >20, or trend correlation <0.90. Annual aggregate correlation and mean within-age temporal correlation are retained, not squared. Missing statistics are not passes. These are review prompts, not tuning targets or scientific-equivalence tests.

## Source and model review

The twelve new cases were reviewed against their cached source inputs, native prior equations and translation code. Per-stock source reviews retain the unit, timing, sex/area, likelihood and completeness findings. Canonical numerical inputs and official outputs were not changed.

Hake's restored source log-M penalty improves median absolute SSB difference from 69.6% to 44.8% (r 0.619 to 0.751), total abundance from 71.8% to 42.8% (r 0.785 to 0.904), and age-2 recruitment from 77.8% to 55.2% (r 0.959 to 0.969). Hake M differs by 21.5%, with finite shared log-M effect uncertainty. The historical BC POP M-prior trial also had finite uncertainty, but is now retired from comparison because complete matching coastwide population series have not been recovered. Native bounds and integration conventions remain different.

EBS cod, Atka and Alaska plaice retain their baselines because alternatives fail checks or worsen other core comparisons. Northern rock sole still fails with verified M/q penalties restored; a random-N initializer passes numerical checks but yields implausible states and is rejected. AI pollock and yellowfin sole remain non-converged; later expensive yellowfin/POP trials are incomplete. The original three source-blocked cases remain blocked, and BC POP is now a fourth.

The related saithe free-initial-abundance change improves all five comparable core metrics; sprat now applies the source 40% pre-spawning fraction. Sprat's revised SSB definition is not a like-for-like agreement improvement. The full one-off trial log and table are retained in [the review checkpoint](https://github.com/PaulRegular/tinyAM/blob/e7b2ac6/analysis/comp_assessments/PACIFIC_MODEL_REVIEW.md); trial fits and initialization histories remain in the local gitignored cache. The experimental runners and duplicate table are retired from the routine workflow.

## All flagged core comparisons

| Case | Metric | Mean absolute % | Terminal % | r | Status | Limitation |
| --- | --- | --- | --- | --- | --- | --- |
| Existing: Northern cod | Total N | 27.6 | 9.88 | 0.943 | matched | Reduced fleets and M/observation dependence |
| Existing: Northern cod | F-at-age | 40.9 | 2.82 | 0.808 | matched | Reduced fleets and M/observation dependence |
| Existing: Northern cod | M | 87.3 | -46.3 | 0.853 | matched; estimated | Reduced fleets and M/observation dependence |
| Existing: Northern cod | N-at-age | 41.2 | 5.09 | 0.923 | matched | Reduced fleets and M/observation dependence |
| Existing: Northern cod | Recruitment | 35.7 | 3.63 | 0.921 | matched | Reduced fleets and M/observation dependence |
| Existing: Northern cod | SSB | 36.8 |  1.8 | 0.956 | matched | Reduced fleets and M/observation dependence |
| Existing: Southern Gulf of St. Lawrence cod | Total N | 22.2 | 3.01 | 0.948 | matched | Mortality priors and observation likelihood differ |
| Existing: Southern Gulf of St. Lawrence cod | F-at-age | 43.8 | 34.1 | 0.915 | matched | Mortality priors and observation likelihood differ |
| Existing: Southern Gulf of St. Lawrence cod | M | 39.4 | -39.4 | — | matched; estimated | Mortality priors and observation likelihood differ |
| Existing: Southern Gulf of St. Lawrence cod | N-at-age | 23.2 | -17.9 | 0.937 | matched | Mortality priors and observation likelihood differ |
| Existing: Southern Gulf of St. Lawrence cod | Recruitment | 36.3 | 19.7 | 0.903 | matched | Mortality priors and observation likelihood differ |
| Existing: Southern Gulf of St. Lawrence cod | SSB | 18.3 | -24.9 | 0.958 | matched | Mortality priors and observation likelihood differ |
| Existing: Northern Shelf cod (formerly North Sea cod) | F-at-age |   19 | -39.5 | 0.799 | matched | Substocks pooled; source-fitted biology simplified |
| Existing: Northern Shelf cod (formerly North Sea cod) | N-at-age | 21.5 |   23 | 0.958 | matched | Substocks pooled; source-fitted biology simplified |
| Existing: North Sea autumn-spawning herring | Recruitment | 3.44 | 21.3 | 0.997 | matched | Source observation dependence and sharing simplified |
| Existing: Eastern Bering Sea walleye pollock | Total N | 22.9 | 20.8 | 0.649 | matched | Composition likelihood and selectivity constraints omitted |
| Existing: Eastern Bering Sea walleye pollock | F-at-age | 50.9 | -8.23 | 0.481 | matched | Composition likelihood and selectivity constraints omitted |
| Existing: Eastern Bering Sea walleye pollock | N-at-age | 41.5 | 23.2 | 0.867 | matched | Composition likelihood and selectivity constraints omitted |
| Existing: Eastern Bering Sea walleye pollock | Recruitment | 33.7 | 21.5 | 0.638 | matched | Composition likelihood and selectivity constraints omitted |
| Existing: Eastern Bering Sea walleye pollock | SSB | 39.6 |   25 | 0.896 | approximate | Composition likelihood and selectivity constraints omitted |
| Existing: Gulf of Alaska Pacific cod | SSB | 32.1 | 7.73 | 0.945 | approximate | Catch-length conversion; source-estimated M fixed |
| Existing: Northern Shelf haddock | Recruitment | 24.8 | 34.3 | 0.968 | matched | Survey precision and observation dependence differ |
| Existing: Northeast Arctic cod | F-at-age | 9.63 | -13.1 | 0.891 | matched | Observation covariance and iterative cannibalism M |
| Existing: Northeast Arctic cod | N-at-age | 10.7 | 25.9 | 0.987 | matched | Observation covariance and iterative cannibalism M |
| Existing: Northeast Arctic cod | SSB | 12.3 | 30.2 | 0.963 | matched | Observation covariance and iterative cannibalism M |
| Existing: Gulf of Alaska walleye pollock | Total N | 29.2 | 83.9 | 0.937 | matched | Fleet aggregation and source q/selectivity priors |
| Existing: Gulf of Alaska walleye pollock | N-at-age | 46.6 |  104 | 0.95 | matched | Fleet aggregation and source q/selectivity priors |
| Existing: Gulf of Alaska walleye pollock | Recruitment | 50.4 | 79.2 | 0.949 | matched | Fleet aggregation and source q/selectivity priors |
| Existing: Gulf of Alaska walleye pollock | SSB | 29.7 |   11 | 0.924 | matched | Fleet aggregation and source q/selectivity priors |
| Existing: Southern Gulf of St. Lawrence spring-spawning Atlantic herring | Total N | 44.9 | -72.3 | 0.792 | matched | CPUE/acoustic conversion and mortality priors |
| Existing: Southern Gulf of St. Lawrence spring-spawning Atlantic herring | F-at-age | 75.6 | 68.9 | 0.876 | matched | CPUE/acoustic conversion and mortality priors |
| Existing: Southern Gulf of St. Lawrence spring-spawning Atlantic herring | N-at-age | 52.3 | -55.4 | 0.797 | matched | CPUE/acoustic conversion and mortality priors |
| Existing: Southern Gulf of St. Lawrence spring-spawning Atlantic herring | Recruitment | 81.1 |  -90 | 0.653 | matched | CPUE/acoustic conversion and mortality priors |
| Existing: Southern Gulf of St. Lawrence spring-spawning Atlantic herring | SSB | 28.6 | -50.4 | 0.864 | matched | CPUE/acoustic conversion and mortality priors |
| Existing: Northeast Atlantic blue whiting | F-at-age | 23.8 | 22.8 | 0.747 | matched | Age pooling and native covariance omitted |
| Existing: Norwegian spring-spawning herring | Total N |   20 | -8.93 | 0.791 | matched | Native process/sharing configuration unresolved |
| Existing: Norwegian spring-spawning herring | F-at-age |   29 | -11.3 | 0.817 | matched | Native process/sharing configuration unresolved |
| Existing: Norwegian spring-spawning herring | N-at-age | 21.5 | 6.71 | 0.943 | matched | Native process/sharing configuration unresolved |
| Existing: Norwegian spring-spawning herring | Recruitment | 45.7 | -23.4 | 0.735 | matched | Native process/sharing configuration unresolved |
| Existing: North Sea plaice | F-at-age | 22.9 | -0.209 | 0.929 | matched | Fleet aggregation and F sharing/correlation differ |
| Existing: North Sea saithe | Total N | 11.1 | 46.1 | 0.935 | matched | Source sharing and observation covariance differ |
| Existing: North Sea saithe | F-at-age | 30.4 |  8.2 | 0.841 | matched | Source sharing and observation covariance differ |
| Existing: North Sea saithe | Recruitment | 15.6 | 48.8 | 0.94 | matched | Source sharing and observation covariance differ |
| Existing: North Sea whiting | F-at-age | 24.3 | 6.57 | 0.929 | matched | Source GMRF M replaced by supplied inputs |
| Existing: North Sea whiting | M | 5.44 | -2.4 | 0.889 | matched; fixed approximation | Source GMRF M replaced by supplied inputs |
| Existing: North Sea whiting | Recruitment | 15.9 | -17.8 | 0.868 | matched | Source GMRF M replaced by supplied inputs |
| Existing: Norway pout | F-at-age | 87.8 | 56.4 | 0.937 | matched | Quarterly source model represented annually |
| Existing: Norway pout | N-at-age | 37.3 | -50.3 | 0.862 | matched | Quarterly source model represented annually |
| Existing: Norway pout | Recruitment |   28 | -47.3 | 0.987 | matched | Quarterly source model represented annually |
| Existing: Norway pout | SSB | 26.7 | -39.3 | 0.989 | matched | Quarterly source model represented annually |
| Existing: Atlantic mackerel | Total N |   26 | -38.8 | 0.877 | matched | Native observation likelihood and selectivity differ |
| Existing: Atlantic mackerel | F-at-age | 77.4 | 35.2 | 0.83 | matched | Native observation likelihood and selectivity differ |
| Existing: Atlantic mackerel | N-at-age | 41.6 | -28.2 | 0.886 | matched | Native observation likelihood and selectivity differ |
| Existing: Atlantic mackerel | Recruitment | 36.9 | -42.9 | 0.614 | matched | Native observation likelihood and selectivity differ |
| Existing: Atlantic mackerel | SSB | 26.4 | -29.2 | 0.946 | matched | Native observation likelihood and selectivity differ |
| Existing: Icelandic haddock | F-at-age | 29.1 | 17.5 | 0.717 | matched | Native sharing settings not fully recovered |
| New: Bering Sea/Aleutian Islands Atka mackerel | F-at-age |  266 | 33.8 | 0.647 | matched | Source selectivity/precision constraints omitted |
| New: Bering Sea/Aleutian Islands Atka mackerel | N-at-age | 30.5 | -22.2 | 0.874 | matched | Source selectivity/precision constraints omitted |
| New: Bering Sea/Aleutian Islands Atka mackerel | SSB | 25.7 | -23.6 | 0.966 | approximate | Source selectivity/precision constraints omitted |
| New: Aleutian Islands walleye pollock | M |   38 |   38 | — | matched; estimated; exploratory not_converged | Non-convergence and extreme N/q scale trade-off |
| New: Aleutian Islands walleye pollock | N-at-age | 5.31e+09 | 5.95e+09 | 0.697 | matched; exploratory not_converged | Non-convergence and extreme N/q scale trade-off |
| New: Aleutian Islands walleye pollock | Total N | 5.1e+09 | 6.5e+09 | 0.535 | matched; exploratory not_converged | Non-convergence and extreme N/q scale trade-off |
| New: Aleutian Islands walleye pollock | Recruitment | 8.33e+09 | 9.31e+09 | 0.759 | matched; exploratory not_converged | Non-convergence and extreme N/q scale trade-off |
| New: Aleutian Islands walleye pollock | SSB | 3.85e+09 | 4.05e+09 | 0.104 | approximate; exploratory not_converged | Non-convergence and extreme N/q scale trade-off |
| New: Bering Sea/Aleutian Islands Alaska plaice | F-at-age |  587 | 1.2e+03 | -0.177 | matched | Sex pooling and selectivity constraints omitted |
| New: Bering Sea/Aleutian Islands Alaska plaice | N-at-age | 29.7 | 10.2 | 0.816 | matched | Sex pooling and selectivity constraints omitted |
| New: Bering Sea/Aleutian Islands Alaska plaice | Total N | 45.2 | 77.2 | 0.277 | matched | Sex pooling and selectivity constraints omitted |
| New: Bering Sea/Aleutian Islands Alaska plaice | Recruitment |  111 | 18.2 | 0.698 | matched | Sex pooling and selectivity constraints omitted |
| New: Eastern Bering Sea yellowfin sole | F-at-age | 3.09e+05 | 1.02e+05 | 0.13 | matched; exploratory not_converged | Non-convergence; female subset and source timing differ |
| New: Eastern Bering Sea yellowfin sole | N-at-age | 3.63e+03 | 2.35e+03 | 0.677 | matched; exploratory not_converged | Non-convergence; female subset and source timing differ |
| New: Eastern Bering Sea yellowfin sole | Total N | 6.8e+03 | 6.76e+03 | 0.139 | matched; exploratory not_converged | Non-convergence; female subset and source timing differ |
| New: Eastern Bering Sea yellowfin sole | Recruitment | 6.41e+04 | 3.67e+04 | -0.0786 | matched; exploratory not_converged | Non-convergence; female subset and source timing differ |
| New: Eastern Bering Sea yellowfin sole | SSB | 47.9 | -63.9 | 0.67 | approximate; exploratory not_converged | Non-convergence; female subset and source timing differ |
| New: Eastern Bering Sea northern rock sole | F-at-age | 1.65e+04 | 4.32e+03 | 0.0623 | approximate; exploratory not_converged | Non-convergence; female subset and joint composition likelihood differ |
| New: Eastern Bering Sea northern rock sole | M |  100 | -100 | — | matched; estimated; exploratory not_converged | Non-convergence; female subset and joint composition likelihood differ |
| New: Eastern Bering Sea northern rock sole | N-at-age | 2.02e+04 | 2.79e+04 | 0.784 | matched; exploratory not_converged | Non-convergence; female subset and joint composition likelihood differ |
| New: Eastern Bering Sea northern rock sole | Total N | 6.12e+04 | 9.13e+04 | 0.308 | matched; exploratory not_converged | Non-convergence; female subset and joint composition likelihood differ |
| New: Eastern Bering Sea northern rock sole | Recruitment | 1.66e+05 | 1.68e+05 | 0.505 | matched; exploratory not_converged | Non-convergence; female subset and joint composition likelihood differ |
| New: Eastern Bering Sea northern rock sole | SSB |  720 | 1.28e+03 | -0.274 | approximate; exploratory not_converged | Non-convergence; female subset and joint composition likelihood differ |
| New: Joint U.S.-Canada coastal Pacific hake | M | 21.5 | -21.5 | — | matched; estimated | Source M penalty restored; recruitment/selectivity and likelihood differ |
| New: Joint U.S.-Canada coastal Pacific hake | N-at-age |  134 |  128 | 0.886 | approximate | Source M penalty restored; recruitment/selectivity and likelihood differ |
| New: Joint U.S.-Canada coastal Pacific hake | Total N | 41.7 | -30.8 | 0.904 | approximate | Source M penalty restored; recruitment/selectivity and likelihood differ |
| New: Joint U.S.-Canada coastal Pacific hake | Recruitment |  137 |  346 | 0.969 | matched | Source M penalty restored; recruitment/selectivity and likelihood differ |
| New: Joint U.S.-Canada coastal Pacific hake | SSB | 42.2 | -30.1 | 0.751 | approximate | Source M penalty restored; recruitment/selectivity and likelihood differ |
| New: Eastern Bering Sea Pacific cod | F-at-age | 70.1 |  107 | -0.2 | approximate | Pooled empirical catch-age key and selectivity simplification |
| New: Eastern Bering Sea Pacific cod | N-at-age | 42.1 | -43.7 | 0.681 | matched | Pooled empirical catch-age key and selectivity simplification |
| New: Eastern Bering Sea Pacific cod | Total N |   19 | -28.8 | 0.861 | matched | Pooled empirical catch-age key and selectivity simplification |
| New: Eastern Bering Sea Pacific cod | Recruitment | 35.6 | 20.9 | 0.697 | matched | Pooled empirical catch-age key and selectivity simplification |
| New: Eastern Bering Sea Pacific cod | SSB | 21.7 | -56.7 | 0.581 | approximate | Pooled empirical catch-age key and selectivity simplification |
| Additional: British Columbia Canary Rockfish | M | 36.5 | -36.5 | — | matched; estimated | Female approximation; joint likelihood, aggregate-only inputs and source recruitment differ |
| Additional: British Columbia Canary Rockfish | SSB | 27.9 | -59.1 | 0.295 | approximate | Female approximation; joint likelihood, aggregate-only inputs and source recruitment differ |

Constant M series have no defined temporal correlation. A fixed approximation to source-estimated M is labelled explicitly; equality is not evidence that tinyAM estimated M successfully.

## New-case numerical outcomes

SD-report creation is separate from usable uncertainty. AI pollock, yellowfin sole and northern rock sole have failed curvature/covariance checks. Hake passes numerical checks with its restored source mortality information. BC POP's historical trial passed too, but its current diagnostic is a source-output blocker, not an optimizer failure. A numerical pass does not establish scientific agreement.

| Case | Status | Code | Objective | Max gradient | SD report created | PD Hessian |
| --- | --- | --- | --- | --- | --- | --- |
| New: Eastern Bering Sea Pacific cod | converged |  0 |  372 | 0.00136 | TRUE | TRUE |
| New: Aleutian Islands Pacific cod | translation_failed | NA | — | — | NA | NA |
| New: Aleutian Islands walleye pollock | not_converged |  1 |  626 | 2.18e-05 | TRUE | FALSE |
| New: Bering Sea/Aleutian Islands Atka mackerel | converged |  0 |  512 | 0.000267 | TRUE | TRUE |
| New: Eastern Bering Sea yellowfin sole | not_converged |  1 | 1.05e+04 | 1.49e+11 | TRUE | FALSE |
| New: Eastern Bering Sea northern rock sole | not_converged |  0 | 3.96e+03 | 0.0788 | TRUE | FALSE |
| New: Bering Sea/Aleutian Islands Alaska plaice | converged |  0 |  296 | 0.000904 | TRUE | TRUE |
| New: Bering Sea and Aleutian Islands arrowtooth flounder | translation_failed | NA | — | — | NA | NA |
| New: Bering Sea/Aleutian Islands Pacific ocean perch | execution_failed | NA | — | — | NA | NA |
| New: Joint U.S.-Canada coastal Pacific hake | converged |  0 | 1.19e+03 | 0.000868 | TRUE | TRUE |
| New: Strait of Georgia Pacific herring | translation_failed | NA | — | — | NA | NA |
| New: British Columbia Pacific ocean perch | translation_failed | NA | — | — | NA | NA |

Initial expansion fits drove M near zero (hake log-M SE 298.1; BC POP 1924.9). Restoring a shared source log-M penalty gives hake M 0.18335 with effect SE 0.09591. The retired BC POP trial used a positive moment approximation to the source normal prior, giving M 0.04053 with effect SE 0.22890. Native hard bounds and integration conventions remain different. Each model specification uses normal starts first and at most one retry after numerical failure; poor agreement alone never triggers an initialization retry. Atka, Alaska plaice and hake used one conditional source-start optimization; EBS cod and the historical BC POP trial passed from defaults. Both allowed attempts failed numerical checks for AI pollock, yellowfin sole and northern rock sole. The interrupted BSAI POP run lacks a completed attempt history; its stage is not guessed.

## Existing numerical or source blockers

| assessment_id | status | max_abs_gradient | positive_definite_hessian |
| --- | --- | --- | --- |
| nefsc_haddock_georges_bank_2026 | blocked_no_stock_spec | NA | NA |
| dfo_cod_3pn4rs_2025 | blocked_no_stock_spec | NA | NA |
| nefsc_summer_flounder_2018 | not_converged | 0.004174769 | FALSE |
| dfo_petrale_bc_2024 | translation_failed | NA | NA |

## Unavailable or non-equivalent comparisons

Each entry below lists missing or incompatible core comparisons. Constant-series correlation is also undefined even when scale comparisons are available. New source-blocked cases and the interrupted POP case have no fitted comparison statistics.

| Case | Unavailable |
| --- | --- |
| Existing: Georges Bank haddock | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| Existing: Gulf of Alaska Pacific cod | Total N: non_equivalent; F-at-age: unavailable; N-at-age: unavailable; Recruitment: non_equivalent |
| Existing: Gulf of Alaska walleye pollock | F-at-age: unavailable |
| Existing: Northern Gulf of St. Lawrence cod | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| Existing: Southern Gulf of St. Lawrence spring-spawning Atlantic herring | M: unavailable |
| Existing: Norway pout | Total N: non_equivalent |
| Existing: Summer flounder | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| New: Bering Sea/Aleutian Islands Pacific ocean perch | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| New: Aleutian Islands walleye pollock | F-at-age: unavailable |
| New: Joint U.S.-Canada coastal Pacific hake | F-at-age: unavailable |
| New: Strait of Georgia Pacific herring | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| New: British Columbia Pacific ocean perch | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| New: Aleutian Islands Pacific cod | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| New: Bering Sea and Aleutian Islands arrowtooth flounder | SSB: no fitted statistics; Recruitment: no fitted statistics; Total N: no fitted statistics; N-at-age: no fitted statistics; F-at-age: no fitted statistics; M: no fitted statistics |
| Additional: British Columbia Canary Rockfish | F-at-age: unavailable; N-at-age: unavailable; Total N: non_equivalent; Recruitment: non_equivalent |

## Fixed-M input checks

This checks aligned values only; it does not credit fixed M as successful estimation. Differences can reflect age pooling or biological approximations and require source-definition review. Missing checks remain unavailable.

| Case | Mean absolute difference / year | Check |
| --- | --- | --- |
| Existing: North Sea autumn-spawning herring | 7.89e-18 | No detected difference in aligned values |
| Existing: Eastern Bering Sea walleye pollock |    0 | No detected difference in aligned values |
| Existing: Northern Shelf haddock | 6.73e-18 | No detected difference in aligned values |
| Existing: Northeast Arctic cod | 3.63e-18 | No detected difference in aligned values |
| Existing: Gulf of Alaska walleye pollock | 5.55e-18 | No detected difference in aligned values |
| Existing: Northeast Atlantic blue whiting |    0 | No detected difference in aligned values |
| Existing: Norwegian spring-spawning herring |    0 | No detected difference in aligned values |
| Existing: Baltic sprat | 3.07e-18 | No detected difference in aligned values |
| Existing: North Sea plaice | 2.78e-18 | No detected difference in aligned values |
| Existing: North Sea saithe | 0.00243 | Review input matching/aggregation |
| Existing: Western Baltic spring-spawning herring | 0.0036 | Review input matching/aggregation |
| Existing: Norway pout |    0 | No detected difference in aligned values |
| Existing: Atlantic mackerel |    0 | No detected difference in aligned values |
| Existing: Icelandic haddock |    0 | No detected difference in aligned values |
| New: Bering Sea/Aleutian Islands Atka mackerel |    0 | No detected difference in aligned values |
| New: Bering Sea/Aleutian Islands Alaska plaice |    0 | No detected difference in aligned values |
| New: Eastern Bering Sea yellowfin sole | 1.39e-17 | No detected difference in aligned values |
| New: Eastern Bering Sea Pacific cod |    0 | No detected difference in aligned values |

## Limits and source blockers

- SoG herring: aggregate spawning-biomass survey with no observed survey ages; density-dependent M and discrete fleets also differ. No artificial age allocation or catch-only substitute.
- AI cod: public September native fit differs from final November production run; final composition/weight inputs and SSB scaling unresolved.
- Arrowtooth: no verified numerical likelihood compositions/age–length key or survey timing; recovered fitted N/F are not observations.
- BC POP: a follow-up recovered regional biomass draws for 3CD and 5DE, but the public 5ABC object ends in 2017. They cannot form the matching coastwide comparison. N/F/recruitment and the native reproductive-output/mature-biomass definition remain unresolved; the case stays blocked. The source review also records the U.S. West Coast petrale alternative.
- BC petrale: an additional partial 2024 summary/output record preserves 1938–2024 SSB posterior summaries, separate published terminal/reference values, three M scalars and 33 assumptions. Original observations, biology and configuration remain missing. It is an applied but non-current summary record, excluded from the detailed-case batch and blocked before optimization. No replacement was fitted and no agreement statistics are available. Mixture intent is documented; component membership of the 2000-draw subset remains unresolved.
- Other new fits omit source composition covariance, ageing error, normalized/sex/fleet selectivity constraints, and in some cases informative M/q priors. Zero composition categories retain source zero values but tinyAM lognormal observations treat them as missing.
- Every new source file is cached with provenance. Completeness remains partial where likelihood arrays, activation, point conventions, biological inputs or joint posterior outputs are unresolved.

See the per-stock source reviews, compact result CSVs and gitignored results/cache/<assessment_id>/dashboard.html for diagnostics and visual comparisons. Source-blocked records have no current comparison dashboard or invented optimizer; the prior BC trial remains archived as historical work. CAPABILITY_GAPS.md distinguishes supported features, approximations and source-data blockers.
