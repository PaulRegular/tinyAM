#' The tinyAM population and observation model
#'
#' @description
#' tinyAM follows fish as they recruit, survive, and move into older age groups.
#' Catch and survey data help estimate abundance and mortality. Process error
#' describes variation in population dynamics; observation error describes
#' variation in the measurements. This reference specifies the model fitted by
#' [fit_tam()] and simulated by [sim_tam()].
#'
#' Choose process assumptions to express biological hypotheses. `"iid"` allows
#' independent departures, `"ar1"` allows departures to be similar in adjacent
#' years and ages, and `"rw"` allows departures to accumulate over time without
#' returning to a mean. These choices make different assumptions, even when
#' they produce similar fitted curves. N and M can also have their process
#' switched `"off"`. Recruitment always follows a temporal random walk.
#'
#' @details
#' ## Indices, units, and states
#'
#' Let \eqn{t=1,\ldots,T} index consecutive modeled years, with historical years
#' ending at \eqn{H \le T}. Let \eqn{a=1,\ldots,A} index consecutive modeled
#' ages (these indices need not equal the biological age labels). The final
#' age is a plus group. At least two years and two ages are required.
#' \eqn{N_{t,a}} is abundance at the start of the year, \eqn{F_{t,a}} and
#' \eqn{M_{t,a}} are annual instantaneous mortality rates, and
#' \eqn{Z_{t,a}=F_{t,a}+M_{t,a}}. Supplied weight \eqn{W_{t,a}} and mature
#' proportion \eqn{P_{t,a}} are treated as known without uncertainty.
#'
#' `log_r`, `log_n0`, `log_n`, `log_f`, and `log_m` are absolute latent log
#' states. Full surfaces `log_N`, `log_F`, and `log_M` are constructed from
#' them by recursion, age-block expansion, and projection. An `eta_*` quantity
#' is a deviation from an expected log state, never an additional absolute state.
#' All matrices have years in rows. Vectorization stacks columns, as in R.
#'
#' Numbers and weights must use consistent units. There is no automatic unit
#' conversion. For [cod_obs], numbers in thousands and weights in kg give
#' biomass in tonnes. Survey catchability absorbs differences in survey units.
#'
#' ## Recruitment and initial abundance
#'
#' Recruitment is \eqn{R_t=N_{t,1}}. The first log recruitment `log_r0` is an
#' estimated fixed state with no process penalty. Subsequent `log_r` states obey
#' \deqn{\log R_t=\log R_{t-1}+\epsilon^R_t,\qquad
#' \epsilon^R_t\stackrel{ind}{\sim}N(0,\sigma_R^2),\quad t=2,\ldots,T.}
#' There is no stock-recruitment relationship or drift term.
#'
#' Initial older-age abundance is chosen independently of the subsequent N
#' process. Write \eqn{b_a=\log N_{1,a-1}-Z_{1,a-1}}, for \eqn{a=2,\ldots,A}.
#' With `init = "exp"`, \eqn{\log N_{1,a}=b_a}. With `"free"`, these log
#' abundances (`log_n0`) are separate unpenalized fixed parameters. With
#' `"random"`, \eqn{\log N_{1,a}=b_a+\eta^0_a}, where the
#' \eqn{\eta^0_a} are independent \eqn{N(0,\sigma_0^2)} residuals and
#' \eqn{\sigma_0} is estimated separately from recruitment and N process SDs.
#' The baseline always uses the realized preceding initial-age state.
#' No equilibrium plus-group correction is applied in this first year.
#' This is an initial-state assumption, not a reconstruction of past cohorts.
#'
#' ## Cohort survival
#'
#' For \eqn{t\ge2}, define expected log abundance for non-plus older ages by
#' \deqn{g_{t,a}=\log N_{t-1,a-1}-Z_{t-1,a-1},\quad 2\le a<A.}
#' For the plus group,
#' \deqn{g_{t,A}=\log\{N_{t-1,A-1}e^{-Z_{t-1,A-1}}+
#' N_{t-1,A}e^{-Z_{t-1,A}}\}.}
#' With the N process off, \eqn{\log N_{t,a}=g_{t,a}}. Otherwise `log_n`
#' contains the realized states and \eqn{\eta^N_{t,a}=\log N_{t,a}-g_{t,a}}
#' follows the selected process below. Its matrix is \eqn{(T-1)\times(A-1)}.
#' An N random walk acts on these cohort residuals, not on absolute abundance.
#'
#' ## Mortality means and latent states
#'
#' Formula design matrices define
#' \deqn{\log\mu^F=X_F\beta_F,\qquad
#' \log\mu^M=\log M^{supplied}+X_M\beta_M.}
#' An absent F formula contributes zero. Either M contribution may be absent
#' (then contributes zero on the log scale), but at least one is required.
#' The M formula intercept is removed when supplied M is also present.
#' Design rows follow the catch table for F and the weight table for M.
#'
#' Historical F has one `log_f` state per year and age:
#' \deqn{\log F_{t,a}=\log\mu^F_{t,a}+\eta^F_{t,a},\quad t\le H.}
#' `log_f` is the left-hand side; `eta_log_f` is the deviation on the right.
#' For an active M process, `log_m` has one state per selected year and age
#' block. Its deviation is `eta_log_m = log_m - log_mu_M`, evaluated at the
#' block's first age. The full M surface is the exponential of that state
#' repeated across the block; the mean is not added again. Each mean component
#' must be constant across ages within a fitted block. Blocks share actual
#' states, not just correlated deviations.
#'
#' The M process starts at `M_settings$first_dev_year`, by default the second
#' historical year, and extends through projection years. Before that start,
#' outside selected age blocks, or with M off, \eqn{M=\mu^M}. Default blocks use
#' `cut_ages(ages[-1], unique(range(ages[-1])))`: normally one shared block
#' excluding the youngest age (two older ages form two singletons under the
#' [cut_ages()] endpoint convention). F has no
#' corresponding delayed start or coupled boundary state.
#'
#' `mu_F` and `mu_M` are exponentiated mean *log* surfaces. For zero-mean
#' Gaussian deviations they are medians, not arithmetic means on the natural
#' scale. No lognormal mean correction is applied to states or observations.
#'
#' ## Process distributions
#'
#' The following distributions apply separately to the matrices
#' \eqn{\eta^N}, \eqn{\eta^F}, and \eqn{\eta^M}; each has its own SD and,
#' for AR1, two correlations. Processes are independent of one another and
#' of recruitment and initial-age innovations before population recursion.
#'
#' * **IID:** all entries are independent \eqn{N(0,\sigma^2)}.
#' * **RW:** successive-row differences are independent
#'   \eqn{N(0,\sigma^2)}, independently across ages or age blocks. The first
#'   row is unpenalized. Conditional on starting residual \eqn{e_{1,a}},
#'   \eqn{E(e_{i,a})=e_{1,a}} and
#'   \eqn{\operatorname{Cov}(e_{i,a},e_{j,b})=
#'   \mathbb{1}(a=b)\sigma^2\min(i-1,j-1)}.
#'   This defines a density of increments, not a proper prior over all levels.
#' * **AR1:** the whole matrix has a stationary zero-mean Gaussian density with
#'   \deqn{\operatorname{Cov}\{\operatorname{vec}(\eta)\}=
#'   \frac{\sigma^2}{(1-\phi_a^2)(1-\phi_t^2)}\;C_a\otimes C_t,}
#'   where \eqn{(C_a)_{ij}=\phi_a^{|i-j|}} and
#'   \eqn{(C_t)_{ij}=\phi_t^{|i-j|}}. Thus \eqn{\sigma} is an innovation
#'   scale; marginal SD is larger when correlations are positive.
#'   Correlation between age blocks uses their index distance, not block widths.
#'
#' Scalar SDs use `exp(log_sd_*)`. Fitted AR1 correlations use
#' `plogis(logit_phi_*)` and are restricted to positive correlations less than
#' one. A singleton matrix axis has its correlation fixed at zero, since it
#' cannot be distinguished from the variance scale. The standalone AR1 helpers
#' also accept negative correlations with absolute value below one.
#' RW has no correlation parameters. If a RW matrix has only one row, its SD
#' is held fixed because there are no increments. Mean coefficients that cancel
#' from all RW increments and do not control any other mortality cells are
#' held at their starting values rather than estimated from a flat likelihood.
#'
#' ## Catchability and observations
#'
#' For survey observation row \eqn{i},
#' \deqn{\log q_i=X_{q,i}\beta_q+B_i d,\qquad d_j\ge0.}
#' Ordinary terms in `q_form` build \eqn{X_q}. An additive [mono()] term
#' builds cumulative step indicators \eqn{B}; its optimized increments are
#' `dq`, on the log-q scale. For level \eqn{k}, its contribution is
#' \eqn{\sum_{j<k}d_j}. A zero step gives an exact plateau. Separate `by`
#' groups have independent steps; ordinary terms supply their baselines.
#' Monotonicity holds with other covariates held constant.
#'
#' Catch and index predictions are
#' \deqn{\widetilde C_{t,a}=N_{t,a}\frac{F_{t,a}}{Z_{t,a}}
#' (1-e^{-Z_{t,a}}),\qquad
#' \widetilde I_i=q_i N_{t_i,a_i}e^{-s_i Z_{t_i,a_i}},}
#' where \eqn{s_i\in[0,1]} is survey sampling time within the year.
#' For each positive observation \eqn{Y_i}, independently conditional on states,
#' \deqn{\log Y_i\sim N(\log\widetilde Y_i,\tau_i^2),\qquad
#' \tau_i=\tau_i^{supplied}\exp(X_{sd,i}\beta_{sd}).}
#' Catch and survey SDs have separate designs and coefficients. Without a supplied
#' SD the baseline is one; with a supplied SD the formula intercept is removed.
#' All supplied SDs and M offsets must be positive and finite; formula covariates
#' must produce one finite design row per observation.
#' Predictions \eqn{\widetilde Y} are conditional medians. Arithmetic means
#' are \eqn{\widetilde Y\exp(\tau^2/2)}.
#'
#' Zero and NA catch/index values are treated as missing, not as censored values
#' or counts. With `fill_missing = TRUE`, their log observations are Gaussian
#' random effects with the same observation equation; integrating each out
#' contributes one, hence no data information. Otherwise they contribute no
#' density. Complete catch, weight, and maturity year-age grids are required;
#' surveys may be sparse. Reducing the maximum age sums catch and index within
#' year (and within survey for indices), and takes unweighted arithmetic means
#' of older-age weight and maturity. Other covariates come from the plus age
#' or the first available older age. Check that this aggregation suits the data.
#'
#' ## Likelihood and inference
#'
#' The joint negative log-likelihood is minus the sum of recruitment innovation,
#' active process, random initial-age, and observation log densities. Gaussian
#' normalizing constants are retained; there are no added hyperparameter priors.
#' Observation densities are evaluated on log data. The response Jacobian
#' \eqn{-\sum_i\log Y_i} is omitted because it is constant in the parameters
#' for the same observed data. Latent parameters are defined in log coordinates;
#' their residual transformations have unit triangular Jacobians.
#'
#' [fit_tam()] integrates `log_r`, active `log_n`, `log_f`, active `log_m`,
#' filled missing values, and random `log_n0` using RTMB's Laplace approximation.
#' All remaining estimated parameters, including free `log_n0`, are fixed
#' effects optimized by `nlminb`; `dq` has a zero lower bound. Uncertainty is a
#' local Gaussian approximation from curvature, with delta-method uncertainty
#' for reported quantities. Boundary steps can have unreliable symmetric Wald
#' intervals. Convergence checks assess gradients (accounting for active step
#' bounds) and positive-definite curvature, not biological identifiability.
#' Inspect optimizer status as well. Intrinsic RW starting levels require care
#' when comparing likelihoods with models having proper initial distributions.
#'
#' ## Projection and simulation
#'
#' For projection year \eqn{H+h}, \eqn{F_{H+h,a}=k_hF_{H,a}}. Multipliers are
#' relative to terminal historical F, not compounded across years. There are
#' no projected F process states. Recruitment, N, and active M processes extend
#' forward under the same equations. A zero F multiplier is replaced by
#' `1e-12` with a warning because the catch equation is evaluated on the log scale.
#' Weight and maturity use age-specific arithmetic
#' means over the requested recent historical years; other covariates are copied
#' from the last available table year. Future catch and index observations are
#' missing. New formula levels cannot be estimated solely from unobserved future
#' data. These projections do not forecast changes in supplied covariates.
#'
#' [nll_fun()] simulation draws process states before constructing N/F/M/Z,
#' then predictions and observation draws using each row's SD. Initial RW
#' states are retained: F/M walks start from the supplied first latent row;
#' the N walk retains its first `log_n` row and computes that row's residual
#' against the newly generated cohort prediction. Random N0 is redrawn; free
#' N0 and first recruitment remain supplied fixed states.
#'
#' [sim_tam()] can hold parameters fixed, draw fixed effects with their estimated
#' covariance while retaining fitted random states, or draw all parameters with
#' the approximate joint precision. These are local Gaussian uncertainty draws,
#' not an exact posterior sampler. Negative drawn `dq` values are clipped to
#' zero. `redraw_random = TRUE` then regenerates process histories across the
#' whole modeled period; it is not a future-only forecast conditional on all
#' historical states. With `FALSE`, observations alone are redrawn conditional
#' on the fitted or sampled states. Filled missing log observations can be
#' returned, but originally missing cells remain excluded from reported catch
#' and yield totals.
#'
#' ## Derived quantities and uncertainty tables
#'
#' Abundance is \eqn{\sum_a N_{t,a}}, biomass is \eqn{\sum_a W_{t,a}N_{t,a}},
#' and SSB is \eqn{\sum_a W_{t,a}P_{t,a}N_{t,a}} at the start of the year.
#' No spawning-time survival or sex-ratio adjustment is added. Fully selected F
#' is \eqn{\max_a F_{t,a}} and selectivity is F divided by that maximum.
#' `F_bar` and `M_bar` are abundance-weighted means over their selected ages.
#' Total catch sums numbers; yield uses the supplied stock weights, not a
#' separate fishery weight series. Observed totals omit missing cells; predicted
#' totals sum the conditional median predictions over all catch cells.
#'
#' Tidied estimates and confidence limits are on reported scales. `se` stays on
#' the estimation scale, identified by `se_scale`. For log quantities, limits are
#' \eqn{\exp(\widehat\ell\pm z\,SE_\ell)}. A small log-scale SE approximates
#' relative uncertainty (CV), but it is not universally a CV. Use the asymmetric
#' confidence interval when uncertainty is large. Zero derived totals cannot
#' have finite log-scale intervals. See [tidy_par()] and [tidy_sdrep()].
#'
#' @seealso [make_dat()], [make_par()], [nll_fun()], [fit_tam()], [sim_tam()],
#'   [dprocess_ar1()], [dprocess_rw()], [mono()]
#' @name tinyAM-model
NULL
