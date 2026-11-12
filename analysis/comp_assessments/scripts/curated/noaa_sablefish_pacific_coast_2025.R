## Source observations and biology ----
years <- 1983:2024
ages <- 0:30
fleets <- c("Trawl", "Hook-and-line", "Pot",
  "Trawl discard", "Hook-and-line discard", "Pot discard")
surveys <- c("Early Triennial survey", "Late Triennial survey", "NWFSC Slope Survey", "WCGBTS")
inputs <- source$inputs
age_samples <- subset(inputs, measure == "frequency_at_age" & year %in% years &
  !grepl("FltSvy code -", notes))
age_samples$sex_code <- sub(".*sex code ([0-9]+).*", "\\1", age_samples$notes)
if (any(!is.na(age_samples$length_bin_lower) | !is.na(age_samples$length_bin_upper))) {
  cli::cli_abort("Sablefish reconstruction requires marginal age samples.")
}

# SS3 uses the -1890 input schedule before annual empirical weights begin in 2003.
weights <- subset(inputs, measure %in% c("weight_at_age", "start_weight_at_age",
  "catch_weight_at_age", "mature_weight_at_age") & year %in% c(-1890L, years))
stopifnot(all(grepl("^kg", weights$unit)), all(is.finite(weights$value) & weights$value >= 0))
weight_key <- with(weights, paste(year, fleet, measure, age, sex))
stopifnot(!anyDuplicated(weight_key))
weight_at <- function(year, age, sex, fleet, measure) {
  weight_year <- ifelse(year < 2003L, -1890L, year)
  sex <- rep_len(sex, length(age))
  lookup <- function(sex) weights$value[match(
    paste(weight_year, fleet, measure, pmin(age, max(ages)), sex), weight_key)]
  value <- ifelse(sex == "combined", (lookup("female") + lookup("male")) / 2, lookup(sex))
  if (any(!is.finite(value) | value < 0)) {
    cli::cli_abort("Sablefish source weights do not cover the requested cells.")
  }
  value
}

# Normalize a joint female/male sample once, before collapsing sexes or ages.
age_share <- function(type, stream, year) {
  stream_name <- if (type == "catch") age_samples$fleet else age_samples$survey
  x <- age_samples[age_samples$type == type & age_samples$year == year &
    stream_name %in% stream, ]
  x <- x[x$sex_code == if (any(x$sex_code == "3")) "3" else "0", ]
  if (!nrow(x) || sum(x$value) <= 0) return(NULL)
  sample <- sub("-(female|male|combined)$", "", x$observation_id)
  stopifnot(length(unique(sample)) == 1L, !anyDuplicated(paste(x$age, x$sex)),
    all(is.finite(x$value) & x$value >= 0))
  x$share <- x$value / sum(x$value)
  x
}

## Separate retained and dead-discard catches ----
catch <- expand.grid(year = years, age = ages, fleet = fleets, KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE)
catch$obs <- NA_real_
catch$unit <- "thousand fish"
catch_weight <- transform(catch[c("year", "age", "fleet")], obs = weight_at(
  catch$year, catch$age, "combined", catch$fleet, "catch_weight_at_age"), unit = "kg", type = "catch")
catch_totals <- subset(inputs, type == "catch" & measure == "total_biomass" &
  year %in% years & fleet %in% fleets)
stopifnot(all(grepl("^metric ton", catch_totals$unit)),
  all(is.finite(catch_totals$value) & catch_totals$value > 0),
  nrow(catch_totals) == length(years) * length(fleets),
  !anyDuplicated(paste(catch_totals$year, catch_totals$fleet)))
catch_reconstruction <- transform(catch_totals[c("year", "fleet")],
  source_tonnes = catch_totals$value, reconstructed_tonnes = NA_real_)

for (i in seq_len(nrow(catch_totals))) {
  total <- catch_totals[i, ]
  share <- age_share("catch", total$fleet, total$year)
  if (is.null(share)) next
  w <- weight_at(total$year, share$age, share$sex, total$fleet, "catch_weight_at_age")
  numbers <- total$value * share$share / sum(share$share * w)
  annual <- aggregate(cbind(numbers, tonnes = numbers * w) ~ age,
    data.frame(age = pmin(share$age, max(ages)), numbers = numbers, w = w), sum)
  rows <- catch$year == total$year & catch$fleet == total$fleet
  annual <- annual[match(catch$age[rows], annual$age), ]
  stopifnot(!anyNA(annual$age))
  catch$obs[rows] <- annual$numbers
  positive <- annual$numbers > 0
  catch_weight$obs[which(rows)[positive]] <- with(annual, tonnes[positive] / numbers[positive])
  catch_reconstruction$reconstructed_tonnes[i] <- sum(annual$tonnes)
}

# Group the oldest ages where selectivity approaches a tail plateau. Keep the
# source changes in time as age-shape predictors, with no fitted SS3 values.
F_breaks <- setNames(list(c(0:15, 20, 31), c(0:10, 13, 16, 20, 31),
  c(0:10, 13, 16, 20, 31), c(0:4, 8, 15, 31),
  c(0:6, 31), c(0:3, 31)), fleets)
catch$F_age <- catch$age
for (fleet in fleets) {
  rows <- catch$fleet == fleet
  breaks <- F_breaks[[fleet]]
  catch$F_age[rows] <- breaks[findInterval(catch$age[rows], breaks)]
}
catch$retained_period <- factor(ifelse(catch$year < 2002, "early", "late"))
shared_gear <- catch$fleet %in% c("Hook-and-line", "Pot")
# Mirror the retained shape across gears, allowing a steep rise, a broader
# decline and an older-age plateau. The 2002 shift changes the younger limb.
catch$fixed_gear_young <- shared_gear * pmin(catch$F_age, 2) / 2
catch$fixed_gear_rise <- shared_gear * pmin(pmax(catch$F_age - 2, 0), 4) / 4
catch$fixed_gear_decline <- shared_gear * pmin(pmax(catch$F_age - 6, 0), 14) / 14
catch$fixed_gear_shift <- shared_gear * (catch$year >= 2002) * pmax(0, 1 - catch$F_age / 6)
catch$trawl_tail <- as.integer(catch$F_age >= 15)
catch$trawl_discard_age <- pmin(catch$F_age, 4)
catch$discard_young_2011 <- (catch$year >= 2011) * pmax(0, 1 - catch$F_age / 7)
catch$discard_young_2019 <- (catch$year >= 2019) * pmax(0, 1 - catch$F_age / 7)

# Approximate the source discard curves with age groups that retain their
# changes in the youngest ages. Shift the hook-and-line groups by one age in
# 2019 and pool the nearly flat older ages in pot discards.
catch$discard_group <- "other"
hkl_discard <- catch$fleet == "Hook-and-line discard"
hkl_period <- cut(catch$year[hkl_discard], c(-Inf, 2010, 2018, Inf), labels = FALSE)
hkl_groups <- rbind(
  c("0", "early 1", "early 2–3", "early 2–3", "tail", "tail", "tail"),
  c("0", "rising 1", "rising 2", "rising 3", "peak", "descending", "tail"),
  c("0", "late 1", "rising 1", "rising 2", "rising 3", "peak", "tail"))
catch$discard_group[hkl_discard] <- hkl_groups[cbind(hkl_period, catch$F_age[hkl_discard] + 1L)]
pot_discard <- catch$fleet == "Pot discard"
pot_period <- 1L + (catch$year[pot_discard] >= 2011)
pot_groups <- rbind(
  c("0", "early 1", "tail", "tail"),
  c("0", "late 1", "late 2", "tail"))
catch$discard_group[pot_discard] <- pot_groups[cbind(pot_period, catch$F_age[pot_discard] + 1L)]

## Survey biomass expanded with observed age compositions ----
totals <- subset(inputs, type == "index" & measure == "total_biomass" &
  survey %in% surveys & year %in% years)
stopifnot(!anyDuplicated(paste(totals$survey, totals$year)),
  all(grepl("^native biomass index units", totals$unit)), all(totals$value > 0))
index_reconstruction <- transform(totals[c("survey", "year")], source_biomass_index = totals$value,
  reconstructed_biomass_index = NA_real_)
index_rows <- vector("list", nrow(totals))
for (i in seq_len(nrow(totals))) {
  total <- totals[i, ]
  share <- age_share("index", total$survey, total$year)
  if (is.null(share)) next
  w <- weight_at(total$year, share$age, share$sex, total$survey, "weight_at_age")
  numbers <- total$value * share$share / sum(share$share * w)
  rows <- aggregate(obs ~ age,
    data.frame(age = pmin(share$age, max(ages)), obs = numbers), sum)
  index_rows[[i]] <- transform(rows, year = total$year, survey = total$survey,
    samp_time = total$sampling_time, unit = "native survey index")
  index_reconstruction$reconstructed_biomass_index[i] <- sum(numbers * w)
}
index <- do.call(rbind, index_rows)
index$q_age <- with(index, ifelse(survey == "Early Triennial survey",
  ifelse(age %in% 1:3, 1, pmin(age, 5)), ifelse(survey == "Late Triennial survey",
  pmin(age, 8), ifelse(survey == "NWFSC Slope Survey", pmin(age, 6), pmin(age, 2)))))
index$slope_age <- pmin(index$age, 4)
index$slope_tail <- pmin(pmax(index$age - 4, 0), 2)
index$wcgbts_order <- factor(index$age, levels = c(0L, 2:30, 1L),
  labels = c("0", rep("2+", 29), "1"))

## Pooled stock biology and official reference ----
biology <- expand.grid(year = years, age = ages, KEEP.OUT.ATTRS = FALSE)
biology$weight <- with(biology, weight_at(year, age, "combined",
  "beginning-season population", "start_weight_at_age"))
female_mature <- with(biology, weight_at(year, age, "female",
  "mature weight schedule", "mature_weight_at_age"))
# January spawning. A 1:1 sex ratio gives 0.5 * female mature mass per pooled fish.
biology$maturity <- 0.5 * female_mature / biology$weight
stopifnot(all(is.finite(biology$maturity) & biology$maturity >= 0 & biology$maturity <= 1))
obs <- tinyAM::as_tam_data(list(catch = catch, index = index,
  weight = rbind(data.frame(year = biology$year, age = biology$age, fleet = NA_character_,
    obs = biology$weight, unit = "kg", type = "stock"), catch_weight),
  maturity = data.frame(year = biology$year, age = biology$age, obs = biology$maturity,
    spawn_time = 0),
  pop_covariates = transform(biology[c("year", "age")], M_prior_group = factor("stock"))))

age_plus_group <- max(ages)
source("analysis/comp_assessments/R/pool_source_states.R", local = TRUE)
source_history <- subset(source$outputs, year %in% years)
comparison_outputs <- pool_source_states(source_history, plus_age = age_plus_group)
source_n <- subset(source_history, type == "population" & measure == "numbers_at_age")
source_f <- subset(source_history, type == "mortality" & measure == "fishing_mortality_at_age")
fleet_f <- lapply(fleets, function(fleet) {
  states <- pool_source_states(rbind(source_n, source_f[source_f$fleet %in% fleet, ]),
    plus_age = age_plus_group)
  rates <- subset(states, measure == "fishing_mortality_at_age")
  rates$fleet <- fleet
  rates$notes <- "Source fleet F weighted by numbers across sexes and terminal ages."
  rates
})
comparison_outputs <- rbind(comparison_outputs, do.call(rbind, fleet_f))
comparison_aggregates <- c("abundance", "ssb")
comparison_definitions <- list(
  N = list(status = "approximate", reason = "Sexes and terminal ages pooled. Ageing error is omitted."),
  F = list(status = "approximate", reason = "Fleet F summed, then N-weighted across sexes and terminal ages."),
  ssb = list(status = "approximate", reason = "Female mature mass weighted by an assumed 1:1 sex ratio."))

## Background ----
background <- c(
  "### Assessment and translation", "",
  print_sources(source$assessment), "",
  "| Component | Official assessment (Stock Synthesis) | tinyAM translation |",
  "|---|---|---|",
  "| History and population | One stock area, separate females and males, 1890–2024, ages 0–70. | Pool sexes and fit 1983–2024, ages 0–30+. 1983 is the first year with active catch and survey age samples. Early biomass-only removals do not enter the catch-at-age likelihood. |",
  "| Fisheries | Trawl, hook-and-line and pot. Each gear has a retained and a dead-discard stream, giving six fleets with biomass totals and separate marginal age compositions. | Fit six fleet surfaces. Reconstruct catch numbers from each stream’s own age composition, weights and total tonnes. Gear-years without active age samples remain missing. Discard biomass already accounts for discard mortality. |",
  "| Selectivity and F | Joined double-normal age curves with time blocks. Retained hook-and-line has sex-specific selectivity, mirrored by pot. | Retained hook-and-line and pot share a piecewise log-age curve with an older-age plateau and a change in 2002. Hook-and-line discards shift from grouped young ages before 2011 to a dome shifted one age older in 2019. Pot discards pool ages 2+ before 2011 and 3+ thereafter. Fleet-specific AR1 year effects change F intensity across ages. |",
  "| Natural mortality | A common constant M for females and males. Its prior has log mean −2.631 and log SD 0.31. Bounds are 0.04–0.13 per year. | Estimate a single constant M using a one-level Gaussian IID effect. Use the bounded prior's log centre (−2.63062) and a tighter log SD (0.19036), with 0.1% prior probability below 0.04. The effect is integrated over positive M. |",
  "| Recruitment and cohort survival | Beverton–Holt recruitment with steepness fixed at 0.75, log recruitment SD fixed at 1.4, annual deviations from 1975 and a bias adjustment. Older cohorts follow survival implied by F and M. An environmental index informs 2020–2024 recruitment deviations. | IID log recruitment and independent log-abundance deviations around expected cohort survivors, each with an estimated SD. Use exponential initialization. The environmental recruitment series stays outside the index likelihood. |",
  "| Biology and reference | Empirical weights through age 30, baseline schedule before 2003, annual schedules thereafter and January spawning. SSB is female mature biomass. | Use both sexes’ average stock weight and sex-specific sample weights for catch/index reconstruction. Use the official female mature weight and assume half the population is female at each age. Reference N sums the sexes, and reference F is N-weighted after summing fleets. |",
  "", "### Survey catchability", "",
  "| Series | Official biomass-index years | tinyAM years with active ages | Official relative catchability at age | tinyAM age form |",
  "|---|---|---|---|---|",
  "| Early Triennial | 1980, 1983, 1986, 1989, 1992 | 1983, 1989, 1992 | High at ages 1–3, with a steep drop by age 5 | Fit levels at 0, 1–3, 4 and 5+ |",
  "| Late Triennial | 1995, 1998, 2001, 2004 | 1995, 2004 | Peak at age 2, declining to a low older-age plateau | Double logistic with an 8+ plateau |",
  "| NWFSC Slope | 1998–2002 | 1998–2002 | Rises to age 4, then levels off near 0.78 | Logistic rise through age 4 with a fitted shift through age 6+ |",
  "| WCGBTS | 2003–2024, excluding 2020 | Same as official | Peak at age 1, age 0 about 0.077, ages 2+ about 0.80 | Fit 0, 1 and 2+ levels with q(0) ≤ q(2+) ≤ q(1), allowing ties |",
  "",
  "Stock Synthesis derives one constant q per survey against vulnerable biomass. Its joined double-normal shapes are the same across years and sexes. The age-0 settings imply low positive selectivity. The model also allows positive terminal plateaus. The AFSC Slope Survey was excluded from the 2025 base fit.",
  "",
  "The tinyAM indices combine the source biomass and age samples at midyear. Catchability scales and age shapes are estimated from these data. The official joined curve and its fixed widths are approximated.",
  "",
  "SS3 fits the biomass indices and age compositions separately. The tinyAM translation uses survey-specific lognormal SDs for reconstructed age rows and omits ageing error, composition weights and within-year dependence. Biomass observations without active age samples remain outside the likelihood.",
  "", "### Coverage and comparison", "",
  "| Removals | Fleet-years with ages | First and last sampled years |",
  "|---|---:|---|",
  vapply(fleets, function(fleet) {
    covered <- catch_reconstruction$fleet == fleet &
      is.finite(catch_reconstruction$reconstructed_tonnes)
    span <- paste(range(catch_reconstruction$year[covered]), collapse = "–")
    paste0("| ", fleet, " | ", sum(covered), " | ", span, " |")
  }, character(1)),
  "",
  "Known annual biomass scales the age samples and appears on the aggregate yield panel. There is no separate likelihood for total catch biomass. Fleet yield panels show the biomass represented by observed age samples.",
  "",
  "The IID cohort deviations allow apparent survival to vary beyond the fitted fishing and natural mortality. They can absorb differences introduced by pooling sexes and reconstructing observations at age.",
  "",
  "Agreement is stronger for SSB and average F in recent years. The translation starts with lower SSB, and total abundance remains lower despite a similar trend. The source review reports the full-period and terminal-year differences."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    rec_settings = list(process = "iid", mu_form = ~1),
    F_settings = list(process = "off", mean_ages = 3:15, age_breaks = F_breaks,
      mu_form = ~0 + ar1(year, by = fleet) +
        fixed_gear_young + fixed_gear_rise + fixed_gear_decline + fixed_gear_shift +
        for_fleet("Trawl", ~ logistic(F_age, by = retained_period) + trawl_tail) +
        for_fleet("Hook-and-line", ~1) +
        for_fleet("Pot", ~1) +
        for_fleet("Trawl discard", ~ normal(trawl_discard_age) + discard_young_2011 + discard_young_2019) +
        for_fleet("Hook-and-line discard", ~ factor(discard_group)) +
        for_fleet("Pot discard", ~ factor(discard_group))),
    M_settings = list(process = "off", mu_supplied = ~ I(exp(-2.63061660)),
      mu_form = ~0 + iid(M_prior_group, sd = 0.19036084)),
    catch_settings = list(sd_form = ~0 + fleet),
    index_settings = list(sd_form = ~0 + survey, q_form = ~0 +
      for_survey("Early Triennial survey", ~0 + factor(q_age)) +
      for_survey("Late Triennial survey", ~ double_logistic(q_age)) +
      for_survey("NWFSC Slope Survey", ~logistic(slope_age) + slope_tail) +
      for_survey("WCGBTS", ~mono(wcgbts_order))),
    silent = silent
  )
  # Show the native annual catch total even where catch-at-age is absent.
  annual_totals <- aggregate(value ~ year, catch_totals, sum)
  fit$pop$total_yield$est <- annual_totals$value[match(fit$pop$total_yield$year, annual_totals$year)]
  fit$rep$total_yield <- fit$pop$total_yield$est
}
