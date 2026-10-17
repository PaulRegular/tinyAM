source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)

## Observations ----
years <- 2000:2023
ages <- 5:25
inputs <- source$inputs
outputs <- source$outputs
weights <- outputs[outputs$measure == "weight_at_age", ]
catch_weights <- outputs[outputs$measure == "catch_weight_at_age", ]
key <- function(x) paste(x$sex, x$age)
compositions <- inputs[inputs$measure == "proportion_at_age", ]
observations <- lapply(split(compositions, compositions$observation_id), function(p) {
  total <- inputs[inputs$type == p$type[1] & inputs$measure == "total_biomass" &
    inputs$year == p$year[1], ]
  stopifnot(nrow(total) == 1)
  w <- catch_weights$value[match(key(p), key(catch_weights))] / 1000
  p$value <- total$value * 1000 * p$value / sum(p$value * w)
  p <- p[order(p$sex, p$age), ]
  pooled <- p[p$sex == "female", ]
  pooled$value <- pooled$value + p$value[p$sex == "male"]
  pooled$sex <- "combined"
  pooled$measure <- "numbers_at_age"
  pooled$basis <- "numbers"
  pooled$unit <- "fish"
  pooled$notes <- "Translation only: total biomass times joint-sex number proportions divided by their sex-specific mid-year weighted mean. Growth-derived weights held at the official estimate."
  pooled
})
observations <- do.call(rbind, observations)
# Fixed source q is a known calibration, not a fitted adjustment to observations.
observations$value[observations$type == "index"] <-
  observations$value[observations$type == "index"] / exp(.182322)
biology <- inputs[inputs$type == "M" & inputs$sex == "female", ]
biology$sex <- "combined"
wf <- weights$value[weights$sex == "female"] / 1000
wm <- weights$value[weights$sex == "male"] / 1000
mat <- inputs$value[inputs$type == "maturity"]
stock_weight <- maturity <- biology
stock_weight$type <- "weight"
stock_weight$measure <- "weight_at_age"
stock_weight$basis <- "kg_per_fish"
stock_weight$unit <- "kg/fish"
stock_weight$value <- (wf + wm) / 2
stock_weight$notes <- "Translation: equal-sex average of computed beginning-year weights."
maturity$type <- "maturity"
maturity$measure <- "maturity_at_age"
maturity$basis <- maturity$unit <- "proportion"
maturity$value <- mat * wf / (wf + wm)
maturity$notes <- "Translation: effective female maturity for equal abundance of the sexes; weight times maturity equals half female weight times female maturity."
obs <- database_to_tam_obs(source$assessment$assessment_id,
  rbind(observations, biology, stock_weight, maturity), years = years, ages = ages,
  assumptions = source$assumptions)
obs$maturity$spawn_time <- 0

## Background and comparisons ----
age_plus_group <- 25
n <- outputs[outputs$measure == "numbers_at_age", ]
f <- outputs[outputs$measure == "fishing_mortality_at_age", ]
f$value <- f$value * n$value[match(paste(f$sex,f$year,f$age),paste(n$sex,n$year,n$age))]
pooled_n <- aggregate(value ~ year + age, n, sum)
pooled_f <- aggregate(value ~ year + age, f, sum)
pooled_f$value <- pooled_f$value / pooled_n$value
pool <- function(x, template) {
  z <- template[template$sex == "female", ]
  z$value <- x$value[match(paste(z$year,z$age),paste(x$year,x$age))]
  z$sex <- "combined"
  z$se <- z$lwr <- z$upr <- NA_real_
  z$source_type <- "derived_source_output"
  z$notes <- "Translation comparison: sum female/male N; mortality is weighted by native N. No covariance or pooled uncertainty available."
  z
}
ssb_age <- n[n$sex == "female", ]
ssb_age$value <- ssb_age$value * weights$value[match(key(ssb_age),key(weights))] *
  mat[match(ssb_age$age,biology$age)]
ssb_age$sex <- "combined"
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- "January 1 female mature biomass from sex-specific official N, computed weight and original maturity; no uncertainty reconstructed."
comparison_outputs <- rbind(outputs[!outputs$measure %in% c("numbers_at_age","fishing_mortality_at_age"), ],
  pool(pooled_n,n), pool(pooled_f,f), ssb_age)
comparison_definitions <- list(ssb = list(status="approximate",
  definition="January 1 female mature biomass at modeled ages",
  reason="The official sex-specific calculation is compared with a pooled-sex model using equal-sex biological weights."))
background <- c(
  "### Alaska plaice: official 2024 Model 24.1b", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Years and ages | 1975–2024, ages 0–25+, separate sexes. | 2000–2023, ages 5–25+, pooled sexes. | Start with the first fishery age sample; end with the last survey age sample. Young ages without useful catch data are omitted. Totals are converted before age restriction. |",
  "| N and recruitment | Deterministic survival, early recruitment deviations; constant recruitment mean (fixed steepness 1), fixed log SD 0.4243, terminal bias rules. | Deterministic survival, free initial ages, IID age-5 recruitment with estimated SD. | Initial distribution, juvenile survival and terminal fixing are approximated. |",
  "| F | Annual intensity and time-invariant sex-specific selectivity. | Independent F random walks at age. | Sex and common-intensity constraints are omitted; the random-walk penalty is an approximation. |",
  "| M | Fixed 0.13 for both sexes. | Same fixed M. | Input consistency, not estimated agreement. |",
  "| Catch and survey | Biomass plus separate joint-sex length/age compositions; survey time 0.5. | Age numbers reconstructed from biomass, joint compositions and computed mid-year weights, then pooled by sex; separate estimated observation SDs. | Only the three fishery age-sample years and survey years with ages enter. Length likelihood, ageing error and source composition weights are omitted. |",
  "| q | Fixed scalar exp(0.182322), with ascending sex-specific age selectivity. | Divide survey numbers by this known fixed calibration; fit one rising logistic selectivity with asymptote 1. | The source fixed q is retained; pooled logistic shape approximates two sex-specific curves. No estimated q is used to rescale observations. |",
  "| Biology | Growth estimated; fixed weight-length parameters and female maturity. January 1 female SSB. | Hold computed weights at official estimates; equal-sex population weights and effective female maturity, spawning time 0. | Growth uncertainty and changing adult sex ratios are omitted. Official sex-specific mature biomass supplies the approximate modeled-age comparison. |",
  "", "Native age-0 recruitment and all-age SSB remain distinct from age-5 entry and modeled-age SSB. Default initialization is tried first; one conditional retry uses aligned pooled N and N-weighted F without changing assumptions."
)

## Model ----
fit_model <- function(...) tinyAM::fit_tam(
  data = obs, years = years, ages = ages,
  N_settings = list(process = "off", init = "free"),
  rec_settings = list(process = "iid", mu_form = ~ 1),
  F_settings = list(process = "rw"),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(q_form = ~ 0 + logistic(age), sd_form = ~ 1),
  silent = silent, ...
)
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- fit_default_first(fit_model, retry_start = function() {
    nn <- .translation_start_surface(comparison_outputs,"population","numbers_at_age",years,ages,1000)
    ff <- .translation_start_surface(comparison_outputs,"mortality","fishing_mortality_at_age",years,ages)
    list(log_r0=log(nn[1,1]),log_r=log(nn[-1,1]),log_n0=log(nn[1,-1]),log_f=log(ff))
  }, attempt_file=file.path("analysis/comp_assessments/results/cache",source$assessment$assessment_id,"initialization_attempts.rds"))
}
