# Gulf of Alaska pollock: accepted-assessment source review

Status: canonical record afsc_pollock_goa_2024 imported; inputs, outputs and assumptions remain partial.
Charbonneau identifier: AFSC_GOA_Gadus_chalcogrammus.

## Accepted assessment

The November 2024 detailed SAFE chapter is cached as source_cache/afsc_pollock_goa_2024/2024_GOA_pollock_SAFE.pdf (118 pages).
Source: https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=GOA+Pollock+2024.pdf&p=bd9a9aa2-3a33-4841-9133-d9ccf593e051.pdf

Final December 2024 SSC report, page 40, explicitly accepts Model 23d for the Western/Central/West Yakutat stock. It endorses revised survey CVs/sample sizes, a Shelikof catchability covariate, removal of Shelikof age-1/age-2 indices, and Dirichlet-multinomial composition likelihoods. The separate Southeast Outside Tier 5 assessment must not be combined with this age-structured stock.
SSC source: https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=SSC+Report+Dec+2024_FINAL.pdf&p=2587e1ac-2a37-4a80-b68b-7a8c7e7a02ca.pdf

The 2025 catch-only notice is cached separately; its rollover was verified against the accepted 2024 run as described below.
Source: https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=1.+Pollock.pdf&p=b8e2f063-e4d4-4749-a5c3-ac04a7bbc612.pdf

## Native files

Official repository: https://github.com/afsc-assessments/GOApollock
Inspected commit: aefe1692520510d55fd25db60121b292a53840a3.
Cached metadata, recursive tree, README, data/2024/pk24_12.txt, data/2024/goa_pk.cpp, data/2024/run_assessment.R, R/tmb_fns.R, R/utils.R and data/READMD.md.

The README explains that annual releases reproduce assessments, while some SAFE inputs/outputs are not publicly shared. No fitted 2024 object appears in the inspected tree. Earlier ADMB/WHAM runs are context, not replacements for Model 23d.

The 2024 script explicitly labels the fit '23d: 2024 final' and requests D-M composition likelihoods. Its first comment incorrectly says 2023; the actual path, data filename, version label and terminal year identify 2024. It fixes log_q4/log_q5 and doubles multN_srv1, multN_srv3 and multN_srv6 after reading. Those post-read transformations must be retained in model-ready composition weights.

The inspected official read_dat function was evaluated independently, without fitting or executing the assessment workflow. It validated the final -999 sentinel, start 1970, terminal 2024, recruitment age 1 and terminal modeled age 10. Parsed values are cached as native_inputs_raw.rds and a lossless key/row/column/value CSV. These are raw parsed inputs, not yet fully transformed inputs consumed by the model.

## Inventory and pending checks

Native total catch: 55 years, tonnes, 1970–2024; terminal catch is assumed 131,000 t. Native fishery age compositions: 49 years, 1975–2023, 10 fitted columns, with lower/upper accumulation ages. Final SAFE Table 1.6 also publishes catch numbers at ages 1–15 for those years; those detailed values are now stored separately in the canonical inputs with their report provenance. Annual fishery weights: 55 by 10. Fishery length-composition placeholders have zero sample sizes and must not be mistaken for fitted observations.

Survey blocks 1, 2, 3 and 6 contain respectively 31, 16, 37 and 6 total-index observations, and 31, 16, 20 and 6 age-composition rows. Survey identities, units, fitted accumulation ages, positive-weight length compositions and annual timing must be matched to the report/source before import. Blocks 4 and 5 remain in the file, but their log-SD values must be inspected: implementation only contributes their likelihood when log SD is positive. Merely finding these arrays does not establish that they are fitted.

The source uses annual survey timing fractions directly in N * exp(-timing * Z), with survey-specific weight matrices. Preserve the timing for each observation year. The environmental catchability component has 40 observed years plus a latent AR1 surface; it is a material input, not optional descriptive context.

Natural mortality is coded as the age vector 1.39, 0.69, 0.48, 0.37, 0.34, 0.30, 0.30, 0.29, 0.28, 0.29, multiplied by natMscalar. The default parameter map fixes natMscalar at 1. The final report confirms externally supplied age-specific M rescaled to an older-age average of 0.3; the scalar is not a substitute for the vector.

The default reader/preparation code fixes sigmaR at 1.3. Final report page 16 confirms this value; the September proposal of 1.0 was not retained. The accepted 2024 repository revision is pinned above.

The 2025 rollover notice, fitted-stream and biology inventory, native/report cross-checks, historical outputs and uncertainty have since been reviewed below. Remaining source gaps are listed at the end.

## Finalized report and staging checks

The finalized 120-page chapter is now cached as 2024_GOA_pollock_SAFE_final.pdf, from https://files.npfmc.org/SAFE/2024/GOApollock.pdf. Prefer it to the earlier 118-page Plan Team draft for extraction.

The 2025 catch notice explicitly confirms no new assessment in 2025. Its GOA-wide totals include Southeast Outside: 2025 ABC 181,022 + 9,749 = 190,771 t; 2026 ABC 133,075 + 9,749 = 142,824 t; 2025 OFL 210,111 + 12,998 = 223,109 t; 2026 OFL 153,971 + 12,998 = 166,969 t. Thus the notice supports production use of the 2024 result; its combined totals must not be stored as outputs of the W/C/WYK age-structured model.

Final report page 16 confirms recruitment sigmaR=1.3, agreeing with the native default fixed parameter. September's proposed 1.0 was not the final convention. Pages 19–20 explain fixed external age-specific M, rescaled to average 0.3 for older fish, and a constant maturity vector based on 1983–2024 female observations. The native maturity vector matches the printed all-years row to rounding. The annual maturity estimates in Table 1.16 are supporting observations, not the maturity surface consumed by this model.

Native arrays confirm all Shelikof age-1/age-2 index log SDs are zero, so their likelihood contributions are disabled. All fishery and survey length-composition sample sizes are zero. Do not import these disabled streams as fitted input coverage.

Active indices are Shelikof winter acoustic (block 1, 31 observations, 1992–2024), NMFS bottom trawl (block 2, 16, 1990–2023), ADF&G crab/groundfish (block 3, 37, 1988–2024), and summer acoustic (block 6, 6, 2013–2023). Source timing fractions are 0.209, annual 0.543–0.584, 0.60989 and 0.519 respectively. Some report prose/tabulated likelihood summaries have stale end years; native arrays control the actual fitted input dates.

Cached stage_report_outputs.py extracted 550 historical numbers-at-age values from Table 1.22 and 110 recruitment/SSB estimates with published 95% intervals and CVs from Table 1.24. Checks require 55 years, 10 ages, valid interval ordering, and exact agreement between printed age-1 N and recruitment. They passed. Units are million fish and thousand tonnes. No SE has been inferred from the rounded CVs. These source-staged outputs have been imported into the canonical records.

## Canonical import and validation

The canonical record contains 5,909 input values, 708 historical output values and 49 assumptions after the catch-table update. scripts/database/017_validate_goa_pollock.R passed, together with full database structural validation. Coverage includes one combined fishery, four active biomass/composition surveys, annual survey/fishery/population/spawning weights, and constant supplied maturity and M vectors. Model years are 1970–2024, ages 1–10+, recruitment age 1, combined-sex abundance with female SSB fraction 0.5. Spawning survival timing is 0.21, distinct from the winter acoustic observation timing 0.209.

The accepted fishery likelihood accumulates ages 1–2 into its first fitted bin and ages 10+ into the terminal bin. The canonical database retains these native grouped proportions and separately stores the report's finer catch numbers at ages 1–15. tinyAM now uses the finer published catch table, summing reported ages 10–15 into its model age-10 plus group. Shelikof compositions accumulate ages 1–3; these model-ready aggregations are labeled reconstructed_source_input with explicit formulas. Remaining compositions preserve native numerical values, without silent renormalization. Survey total biomass units are million tonnes on the native scale. SD rows match the same observation years and sampling times, with blank fish age for aggregate indices. Supplied spawning weights use measure spawning_weight_at_age, distinct from population weight_at_age and survey-specific weights; biological purposes are not treated as spatial regions.

Remaining gaps are the physical environmental-covariate units and age-specific F outputs, which are not available in the cached report or model files. Active selectivity priors and the Dirichlet-multinomial parameter prior are now recorded from the pinned source, alongside the temporal penalty scales, environmental observations, observation SD, age-error matrix and composition sample sizes. Completeness statuses remain partial. No model fit or numerical substitutions were used.

## Native process and initial-state review

The accepted 2024 script uses `prepare_pk_input()` and its default parameter map. Extra initial-age deviations are fixed at zero. Initial ages 2–10 are constructed from first-year recruitment and M, with the exact zero-based mortality indexing recorded in the assumptions table; subsequent survival is deterministic conditional on F and M. Fishery selectivity is double-logistic, normalized at age 7, with penalized annual ascending-limb changes and fixed descending-limb deviations.

Catchability differs by survey: Shelikof uses a baseline plus an integrated latent environmental effect, bottom trawl uses constant estimated q with a log-scale prior, ADF&G uses penalized annual log-q changes, and summer acoustic uses constant estimated q. Thus the implementation does not support treating all survey catchabilities as interchangeable random walks. Evidence is the cached, pinned `R/tmb_fns.R`, `data/2024/run_assessment.R`, and `data/2024/goa_pk.cpp`; the environmental observations and age-error matrix are now represented as described below.

The 40 native environmental observations and their exact observation years are now canonical covariate inputs. Values are preserved without restandardization; physical units remain unresolved. The accepted map fixes their observation SD at 0.02 on the native scale. The fitted latent annual environmental process is distinct from these supplied observations. Validation compares all 40 values and years directly with the cached native export.

The supplied composition sample sizes are now explicit numeric fields, including the accepted script’s doubling for surveys 1, 3 and 6. They remain distinct from D-M effective sample sizes. The constant 10-by-10 `age_trans` matrix is preserved row by row in the biological assumptions record. Its orientation is verified from prediction code: true-age row compositions multiply the matrix to obtain observed-age compositions before bin accumulation. All four active surveys and the fishery use this matrix.

Native age-transition row sums differ from one by up to 0.0001 because of supplied numerical precision. Values are retained exactly, without renormalization; validation allows this source precision.

Survey selectivity shapes and fixed/estimated limbs are now recorded separately for all four active surveys. All 54 fishery ascending log-slope increment SDs equal 0.05; ascending inflection SDs are 0.2. The 54 ADF&G log-q increment SDs are preserved in transition order and checked directly against native inputs. Supplied q1/q2 penalty arrays are not treated as active annual processes because their corresponding deviations are fixed at zero.

Final Table 1.23 adds 48 accepted-run beginning-year age-3+ biomass estimates for 1977–2024, in thousand tonnes. The previous-assessment columns are excluded. SSB and recruitment in the same rows agree exactly with independently extracted Table 1.24. Harvest rate is catch divided by age-3+ biomass, not instantaneous F; it is not used to manufacture F-at-age.

## Earlier tinyAM translation review

The following records earlier exploratory models. The scale audit and latest
model-review settings are recorded in the subsequent sections.

The pinned selectivity equations do not pool older survey ages explicitly. They
use a descending logistic for Shelikof, an ascending logistic for ADF&G, and
double logistics for bottom trawl and summer acoustic. Bottom-trawl descending
slope and inflection are fixed at exp(1) and 20, so its descending multiplier is
indistinguishable from one over ages 1-10. Summer acoustic has its ascending
slope and inflection fixed at exp(4.9) and 0.5, making its ascending multiplier
indistinguishable from one over these ages. Thus increasing trawl and decreasing
acoustic curves are defensible approximations to the implemented source shapes.
The revised tinyAM formula uses independent monotone steps by survey, with
reversed age order for the acoustics. Plateaus are possible, rather than an
arbitrary fixed threshold for pooling old ages. These are not fitted source
logistic curves. The bottom-trawl absolute-q prior remains unrepresented, so
shape constraints alone do not reproduce the source abundance-scale constraint.

ADF&G annual effects use a single year factor in place of manually created
dummy columns, retaining the first composition year as the baseline. The
observed environmental covariate still applies only to Shelikof. Source latent
environmental uncertainty, annual-q penalties and age-reading error are not
implemented in tinyAM.

A quadratic log-SD was tested to allow higher variability at both ends of the
fitted age range, but the combined model did not converge from either the
original or staged starts. The retained model uses a common catch log-SD.
The detailed Table 1.6 catch numbers-at-age are now used by tinyAM for ages
1-15, with ages 10-15 summed into the model's age-10 plus group. This restores
separate young-age observations and avoids reconstructing their counts from the
pooled composition. Published values are rounded to 0.01 million fish; a printed
zero is below that reporting precision and is treated as missing by tinyAM's
lognormal likelihood. The accepted assessment itself continues to fit its
ages-1-2 and 10+ pooled composition with age-reading error.

The yield panel uses the original annual total catch biomass (tonnes converted
to kg) and predicted catches at all model ages multiplied by corresponding
catch weights. This remains a reporting comparison, not an additional
total-catch likelihood. Predictions are conditional medians, without a
lognormal mean correction. Plotly's missing-observation trace indexing is
corrected separately; the underlying catch predictions were present at every
age.

Accepted SSB is `sum(N * exp(-0.21 * Z) * wt_spawn * 0.5 * mat)`;
tinyAM SSB is `sum(N * wt_pop * 0.5 * mat)` at the start of the year.
The source preparation sets `wt_spawn = wt_srv1` and `wt_pop = wt_srv2`.
Differences in weights and timing, as well as q's omitted absolute-scale prior,
must be considered before attributing all SSB disagreement to age-specific q.

### Numerical and reporting checks

The survey-shape fits below were exploratory runs made before the report's
age-specific catch table was connected to tinyAM. Their fit statistics describe
that earlier, incomplete catch translation and are retained only as context.
The final translation uses the same reviewed survey curves, F/N/M settings and
common catch-SD formula, now fitted to the reported catch numbers at ages 1-15.
Its likelihood contains 490 age-year rows, including 49 observations each for
ages 1 and 2. Printed zeroes in the report are below its rounding precision and
are omitted by the lognormal likelihood.

| Fit | Objective | Maximum projected gradient | Positive-definite Hessian | Outcome |
|---|---:|---:|---|---|
| Original | 1018.232 | 0.00121 | Yes | Converged |
| Survey shape, common catch SD | 1056.099 | 0.000724 | Yes | Converged; retained |
| Polynomial catch SD, free age q | 1007.277 | 0.0170 | Yes | Above gradient tolerance |
| Survey shape + polynomial catch SD | 2303.757 | 139.0 | No | Evaluation limit; rejected |
| Same combined model, started from the converged survey-shape fit | 1043.167 | 3.57 | Yes | Evaluation limit; rejected |
| Detailed catch-at-age data, retained settings | 1208.747 | 0.0008 | Yes | Converged; earlier fit |

In the earlier age-incomplete fits, the common catch log-SD was 0.417 and the
survey-shape fit had 18 exactly-zero monotone increments. The polynomial-only
fit estimated larger SD at ages 3 and 10 (0.689 and 0.584), with its lowest
value near age 7 (0.303), but this pattern did not yield a stable combined
model. These exploratory results predate the catch-table correction.

The earlier finding that ages 1-2 predictions were unchecked no longer applies:
the later fits use the published young-age catch numbers. The median
predicted/observed ratio for age 1 is 1.48 and for age 2 is 0.92. Overall
trajectory correlations and scale differences are in the refreshed comparison
table; SSB remains definition-mismatched because spawning time and weight
definitions differ as noted above.
Source total-catch and grouped-age likelihoods, source fishery selectivity,
and the bottom-trawl absolute-q prior remain important structural differences.
The revision does not establish that any one of these explains the remaining
abundance/SSB discrepancy.

A quadratic mean-log-F sensitivity was also attempted with the survey-shape
and polynomial-SD changes, but was stopped after 20 minutes without a completed
fit. Its convergence and Hessian are unknown. It was not retained at that time.

## Scale audit and pre-refinement translation (7 October 2026)

The recipe at the start of this review used exponential initial abundance, IID N
process deviations, age-specific temporal random walks in F, a common catch
log-SD, and survey-specific paired-age q blocks (1-2, 3-4, 5-6, 7-8, 9-10).
The Shelikof observed environmental effect and ADF&G annual effects remain.
Survey log-SDs combine supplied aggregate-index log-SDs with an estimated level
per survey. These settings supersede the earlier exploratory fits above.

The pinned C++ source constructs catch biomass as
`1e6 * sum(C * wt_fsh)` in tonnes, with C and population N in billions of fish
and weights in kg per fish. Its survey biomass equation omits that factor:
`sum(q * N * exp(-timing * Z) * selectivity * wt_srv)` is therefore in million
tonnes. The canonical survey totals correctly retain that native unit. The
translation converts million tonnes to kg with a factor of 1e9, then obtains
individual fish as `biomass_kg * p_age / sum(p_age * weight_kg)`.
The recent conversion fix in commit 4dd97fd corrected the former factor of
1000 for these totals, which understated survey numbers by a factor of one
million. No further unit error was found.

Catch numbers from final SAFE Table 1.6 remain in million fish in the database
and become individual fish with a factor of 1e6; total catch tonnes become kg
with a factor of 1000. Annual stock, catch and survey weights remain kg per
fish. Accepted N and recruitment remain million fish, and SSB/biomass remain
thousand tonnes; comparison factors of 1e-6 convert tinyAM fish/kg to those
native reporting units without altering source values.

Before omitting the pooled Shelikof bin, reconstructed numbers multiplied by
matching survey weights reproduce every corresponding source biomass total
for all four surveys (maximum relative error 2.3e-16). Omission of that bin
leaves the retained ages unchanged. As an independent scale check, the 1992
Shelikof ages 4-9 reproduce final SAFE Table 1.11 to its printed precision
(0.1 million fish). Native-input reconstructions for bottom-trawl and summer
surveys differ slightly from the published Table 1.9 counts; those tables are
not interchangeable with the native composition/weight inputs. No correction
was applied to force agreement. Age-reading error and grouped young ages are
still translation approximations, not unit conversions.

That unchanged baseline model converged: optimizer code 0, objective 1261.842,
maximum absolute gradient 0.000919, and positive-definite Hessian. Median
observation-level q was 0.409 (ADF&G), 1.906 (bottom trawl), 1.251 (Shelikof),
and 1.041 (summer acoustic). These are the full observation q values including
annual/environmental effects, rather than exponentiated covariate coefficients.
q is an estimated index-to-population multiplier and is not constrained to
0-1. The source's separately normalized selectivity curves and bottom-trawl q
prior are not reproduced by this recipe. The remaining abundance/q scale
tradeoff has not been resolved by this unit audit; no additional scaling or
model adjustment was applied merely to lower q.

## Optional logit-q sensitivity (7 October 2026)

The initial link sensitivity fitted the same observations and settings with
`index_settings$q_link` set to `"log"` or `"logit"`. The review script is now
`scripts/translation/review_goa_pollock.R` and separates the link change from
the older-age N-process change. Fits, summaries and dashboards are cached
locally, without replacing aggregate batch outputs.

Both links converged with positive-definite Hessians. The logit fit had objective
1266.537 versus 1261.842 for the log link, and maximum gradient 0.00106 versus
0.00092. All logit q predictions were below one, with survey maxima 0.398
(ADF&G), 0.703 (bottom trawl), 0.677 (Shelikof), and 0.420 (summer acoustic).
Catch/index standardized residual SDs were 0.679/0.904, versus 0.690/0.906.

Terminal 2024 SSB increased from 134,050 to 467,441 tonnes. Its log-scale SE
increased from 0.262 to 0.715; the logit fit's approximate 95% interval was
115,019-1,899,690 tonnes. A second fit starting other parameters from the
previous log-link fit and projecting capped baseline q onto the logit design
reproduced the same objective and SSB, with maximum gradient 0.000192.
The cap (0.99) was used only to construct finite starting values, not as a
model bound or adjustment to observations.

The logit option changes environmental/year effects from constant multipliers
of q to additive effects on logit-q. It is therefore more than a simple cap
on the original log-link predictor. Its restriction does not anchor q near
one or reproduce the source's bottom-trawl baseline-q prior. The trial is
numerically valid but does not establish a better abundance scale, and its
wider uncertainty warrants review before replacing the existing translation.

## Population-process and catchability review (7 October 2026)

All trials retain the original observations, fixed age-specific M, ages 1-10+,
1970-2024 period, F temporal random walks, survey timing, weights and maturity.
No package feature was added during this model review. The retained changes
use existing settings: deterministic older-age survival (`N` process off)
and a logit q link. Exponential initial abundance, paired-age survey q blocks,
environmental/year effects and common catch SD are retained.

The four controlled fits separate the effects of the link and N process:

| Model | Objective | Maximum gradient | N mean absolute % difference | Recruitment trend correlation | Shared-biology SSB mean absolute % difference |
|---|---:|---:|---:|---:|---:|
| Previous: IID N, log q | 1261.842 | 0.000919 | 57.0 | 0.629 | 62.7 |
| IID N, logit q | 1266.537 | 0.001059 | 87.8 | 0.611 | 37.7 |
| N off, log q | 1316.473 | 0.000854 | 58.5 | 0.938 | 65.1 |
| N off, logit q | 1358.922 | 0.000227 | 46.6 | 0.949 | 29.7 |

All four have optimizer code 0 and positive-definite Hessians. These are
descriptive comparisons, not a likelihood ranking or scientific validation.
SSB differences above use accepted N with the same translated stock weights
and maturity as tinyAM. The dashboard preserves native accepted SSB, whose
spawning weights and survival timing differ from tinyAM's start-year definition.

With N off and logit q, 2024 tinyAM SSB is 351,022 tonnes (approximate 95%
interval 244,988-502,950), versus native accepted SSB of 302,000 tonnes
(published interval 241,000-379,000). On the shared-biology definition the
terminal SSB difference is +11.0%, compared with -57.6% previously. Total
abundance and recruitment remain high in 2024 (+83.9% and +79.2%); their
historical trends improve, but this is not a uniform improvement in scale.
Survey standardized-residual SDs are 0.965-0.993, with means -0.036 to 0.050;
catch residual means by age are -0.128 to approximately zero. The common catch
log-SD is 0.527. Catch predictions remain conditional medians, so aggregate
yield need not equal a sum of arithmetic means.

Uncertainty checks identify weakly estimated bottom-trawl q blocks: ages 5-6
and 7-8 approach one, with logit-coefficient SEs about 71 and 612. A positive
Hessian and small gradient do not make their symmetric Wald intervals reliable.
Pooling bottom-trawl ages 5+ did not resolve this: its q approached one and
the coefficient SE increased to about 2568. That pooling was not retained.
The logit restriction supplies no source-style q prior; it cannot establish
an absolute abundance scale by itself.

Additional sensitivities were not retained. With IID N, an AR1 F process with
age-specific means produced false convergence and a non-positive Hessian;
AR1 N converged but strongly inflated abundance. Monotone q trials with either
link and either IID or deterministic N had non-positive Hessians. Quadratic
catch SD converged with the previous log link but worsened agreement; its
logit trials had curvature/convergence failures. Survey-specific linear and
quadratic q formulas under IID N improved SSB but left recruitment trends
weak. Under deterministic N, linear q converged but its acoustic curves rose
with age, unlike the source shapes; quadratic q on raw age retained a gradient
of 0.030. The combined deterministic-N, quadratic-q/catch-SD trial was stopped
after several minutes without completion; no convergence result is claimed.

A quadratic q formula using `survey:poly(age, 2)` was also checked. Its default
start converged to a different solution (objective 1340.513) with substantially
worse abundance agreement. Transforming the raw-quadratic fit's full predictor
exactly into this basis and reusing its other estimates gave objective 1331.791,
gradient 0.000802 and a positive Hessian. This was a change of basis and starting
values only; predicted q at the start was verified equal to within 1e-9, without
clipping or altering observations. Neither solution was retained: the extra
curve assumptions and sensitivity to starts do not offer a clear advantage
over the simpler paired-age representation. The retained fit remains a
provisional illustrative translation, particularly because its near-one
bottom-trawl q is weakly estimated.
