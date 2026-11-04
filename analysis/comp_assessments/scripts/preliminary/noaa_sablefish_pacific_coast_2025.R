## Source observations and biology ----
years <- 2006:2023
ages <- 0:30
survey <- "WCGBTS"
inputs <- source$inputs
age_samples <- inputs[inputs$measure == "frequency_at_age" & inputs$year %in% years &
  !grepl("FltSvy code -", inputs$notes), ]
age_samples$sex_code <- sub(".*sex code ([0-9]+).*", "\\1", age_samples$notes)
if (any(!is.na(age_samples$length_bin_lower) | !is.na(age_samples$length_bin_upper))) {
  cli::cli_abort("Sablefish reconstruction requires marginal age samples, not conditional age-at-length data.")
}
weights <- inputs[inputs$measure %in% c("weight_at_age", "start_weight_at_age",
  "catch_weight_at_age", "mature_weight_at_age") & inputs$year %in% years, ]
weight_at <- function(year, age, sex, fleet, measure) {
  sex <- rep_len(sex, length(age))
  x <- weights[weights$year == year & weights$fleet == fleet &
    weights$measure == measure, ]
  if (anyDuplicated(paste(x$age, x$sex)) || any(!grepl("^kg", x$unit))) {
    cli::cli_abort("Sablefish weights need unique age/sex rows in kg per fish.")
  }
  lookup <- function(sex) x$value[match(paste(pmin(age, 30), sex), paste(x$age, x$sex))]
  value <- ifelse(sex == "combined", (lookup("female") + lookup("male")) / 2, lookup(sex))
  if (any(!is.finite(value) | value < 0)) cli::cli_abort("Sablefish source weights do not cover the requested cells.")
  value
}
age_share <- function(type, stream, year) {
  x <- age_samples[age_samples$type == type & age_samples$year == year &
    (if (type == "catch") age_samples$fleet %in% stream else age_samples$survey %in% stream), ]
  # Prefer the one joint female/male sample; do not also add unsexed samples.
  if (any(x$sex_code == "3")) x <- x[x$sex_code == "3", ]
    else x <- x[x$sex_code == "0", ]
  if (!nrow(x) || sum(x$value) <= 0) return(NULL)
  sample <- sub("-(female|male|combined)$", "", x$observation_id)
  if (length(unique(sample)) != 1L || anyDuplicated(paste(x$age, x$sex)) ||
    any(!is.finite(x$value) | x$value < 0)) {
    cli::cli_abort("Sablefish needs one valid joint age sample per stream and year.")
  }
  x$share <- x$value / sum(x$value)
  x
}

## Reconstruct only complete annual removals ----
catch <- expand.grid(year = years, age = ages)
catch$obs <- NA_real_
catch$unit <- "thousand fish"
catch_coverage <- data.frame(year = years, complete = FALSE)
catch_reconstruction <- list()
for (year in years) {
  totals <- inputs[inputs$type == "catch" & inputs$measure == "total_biomass" &
    inputs$year == year & inputs$value > 0, ]
  if (!nrow(totals) || anyDuplicated(totals$fleet)) cli::cli_abort("Sablefish catch totals need unique fleets.")
  numbers <- list()
  for (i in seq_len(nrow(totals))) {
    total <- totals[i, ]
    if (!grepl("^metric ton", total$unit)) cli::cli_abort("Sablefish catch totals must be metric tonnes.")
    share <- age_share("catch", total$fleet, year)
    if (is.null(share)) next
    w <- weight_at(year, share$age, share$sex, total$fleet, "catch_weight_at_age")
    value <- total$value * share$share / sum(share$share * w)
    numbers[[total$fleet]] <- data.frame(age = pmin(share$age, max(ages)), obs = value)
    catch_reconstruction[[length(catch_reconstruction) + 1L]] <- data.frame(
      year = year, fleet = total$fleet, source_tonnes = total$value,
      reconstructed_tonnes = sum(value * w))
  }
  complete <- length(numbers) == nrow(totals)
  catch_coverage$complete[catch_coverage$year == year] <- complete
  if (complete) {
    annual <- aggregate(obs ~ age, do.call(rbind, numbers), sum)
    rows <- catch$year == year
    catch$obs[rows] <- annual$obs[match(catch$age[rows], annual$age)]
  }
}
catch_reconstruction <- do.call(rbind, catch_reconstruction)

## Relative survey indices and supplied biological schedules ----
totals <- inputs[inputs$type == "index" & inputs$measure == "total_biomass" &
  inputs$survey %in% survey & inputs$year %in% years, ]
if (anyDuplicated(totals$year)) cli::cli_abort("Sablefish survey totals need unique years.")
index <- do.call(rbind, lapply(seq_len(nrow(totals)), function(i) {
  total <- totals[i, ]
  share <- age_share("index", survey, total$year)
  if (is.null(share)) return(NULL)
  w <- weight_at(total$year, share$age, share$sex, survey, "weight_at_age")
  rows <- data.frame(age = pmin(share$age, max(ages)),
    obs = total$value * share$share / sum(share$share * w))
  rows <- aggregate(obs ~ age, rows, sum)
  transform(rows, year = total$year, survey = survey,
    samp_time = total$sampling_time, unit = "native survey index")
}))
# Evaluate each original annual biological schedule separately.
biology <- do.call(rbind, lapply(years, function(year) {
  w <- weight_at(year, ages, "female", "beginning-season population", "start_weight_at_age")
  mid <- weight_at(year, ages, "female", "mid-season population", "weight_at_age")
  mature <- weight_at(year, ages, "female", "mature weight schedule", "mature_weight_at_age")
  data.frame(year = year, age = ages, weight = w, maturity = 0.5 * mature / mid)
}))
if (any(!is.finite(biology$maturity) | biology$maturity < 0 | biology$maturity > 1)) {
  cli::cli_abort("Sablefish mature-weight conversion does not yield valid proportions.")
}
obs <- tinyAM::as_tam_data(list(catch = catch, index = index,
  weight = data.frame(year = biology$year, age = biology$age, obs = biology$weight, unit = "kg"),
  maturity = data.frame(year = biology$year, age = biology$age, obs = biology$maturity)))
attr(obs, "translation") <- list(
  notes = "Joint number frequencies normalized once over both sexes; inactive samples excluded. Combined-sex samples use mean female/male weights. Ages above 30 use the empirical age-30 terminal weight. Years with any unsampled positive catch stream remain wholly missing.",
  catch_coverage = catch_coverage)

N_settings <- list(process = "off", init = "exp")
rec_settings <- list(process = "iid", mu_form = ~1)
F_settings <- list(process = "ar1", mu_form = ~logistic(age), mean_ages = 3:15)
M_settings <- list(process = "off", mu_form = ~1, mu_supplied = NULL)
catch_settings <- list(sd_form = ~1)
index_settings <- list(q_form = ~logistic(age), sd_form = ~1)
age_plus_group <- 30L
source("analysis/comp_assessments/R/pool_source_states.R", local = TRUE)
comparison_outputs <- pool_source_states(source$outputs, plus_age = 30L)
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(
  N = list(status = "approximate", reason = "Sexes and terminal ages pooled; ageing error is omitted."),
  F = list(status = "approximate", reason = "Fleet F summed, then N-weighted across sexes and terminal ages; pooled dynamics remain approximate."),
  ssb = list(status = "approximate", reason = "Shared female weight proxy and a fixed 1:1 sex ratio replace the native two-sex biology."))

## Background ----
background <- c(
  "### U.S. West Coast sablefish: preliminary comparison with the Official 2025 assessment", "",
  print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM representation | Reason for difference |",
  "|---|------|------|------|",
  "| Coverage | Two sexes, ages 0–70, empirical weights end at age 30 and observed age bins at 50. | Fit 2006–2023, ages 0–30+, pooling sexes. | This is the recent overlap with empirical biology and mostly complete removal age samples. Older-bin weights use the terminal age-30 schedule, explicitly approximating the tail. |",
  "| Removals | Six retained/discard fleets with biomass totals and separate compositions. | Reconstruct each stream separately with its own sex/age weights; sum only when every positive stream is represented. | 2014 has missing pot ages and remains wholly missing. Unsampled removals are never silently omitted from a partial total. |",
  "| Composition | Raw marginal frequencies, including joint sexes and inactive rows. | Exclude inactive samples; normalize each joint sample once across sexes. Use a joint sample in preference to a separate unsexed sample. | Age error, sample weights and extra unsexed likelihoods are omitted. Combined-sex samples use the mean of female/male weights, an explicit sex-mixture approximation. |",
  "| Index | Four biomass indices and a recruitment-deviation index. | Use WCGBTS total and its own age sample and weights at the source timing. | Earlier surveys fall outside this fit. Recruitment deviations are not abundance data. Unknown physical index multipliers are not invented; q stays relative. |",
  "| N and recruitment | Two-sex SS3 with Beverton–Holt recruitment and recruitment bias adjustment. | Deterministic pooled survival, IID recruitment around a mean and exponential initial abundance. | Spawner–recruitment dependence, sex dynamics and source bias adjustments are not reproduced. |",
  "| F and M | Fleet-specific selectivity and one estimated constant M with a prior. | Aggregate F has AR1 age/year residuals about a logistic mean; estimate constant M without the prior. | F dynamics and regularization differ. M is estimated, not fixed at its official estimate. |",
  "| Biology | Original sex-specific beginning/mid-season weights and maturity-times-weight schedules. | Derive female maturity by dividing mature weight by mid-season female weight; use beginning-year female weight as a shared proxy and multiply maturity by 0.5. | A fixed 1:1 sex ratio approximates source sex structure; biology is not obtained from fitted N/F. Catch/index conversions still use stream-specific sex weights. |",
  "| Comparison | Native two-sex N and fleet F; native SSB and recruitment. | Sum source N; add fleet F before N-weighting rates across sexes and terminal ages. | Pooled rates are descriptive common-definition quantities, not equivalent sex-specific dynamics. No uncertainty is invented for these sums. Native aggregate outputs remain separate. |",
  "", "Default starts are attempted first. One numerical retry may use aligned, explicitly pooled official states; poor agreement is not a reason to retry.",
  "", "The first validated run passed numerical checks but estimated M close to zero with extreme uncertainty. It is flagged for review: omitting the M prior can leave mortality weakly identified. Numerical convergence is not evidence of a successful replication."
)

## Model ----
fit <- NULL
if (do_fit) {
  source("analysis/comp_assessments/R/fit_default_first.R", local = TRUE)
  source("analysis/comp_assessments/R/source_start_parameters.R", local = TRUE)
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit_model <- function(start_par = NULL) tinyAM::fit_tam(
    data = obs, years = years, ages = ages, N_settings = N_settings,
    rec_settings = rec_settings, F_settings = F_settings, M_settings = M_settings,
    catch_settings = catch_settings, index_settings = index_settings,
    start_par = start_par, silent = silent)
  fit <- fit_default_first(fit_model, retry_start = function() {
    dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages,
      N_settings = N_settings, rec_settings = rec_settings, F_settings = F_settings,
      M_settings = M_settings, catch_settings = catch_settings, index_settings = index_settings)
    states <- source; states$outputs <- comparison_outputs
    source_start_parameters(states, dat)
  }, attempt_file = file.path(.assessment_cache(source$assessment$assessment_id,
    translation_tier, database = source), "initialization_attempts.rds"))
}
