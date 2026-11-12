## Source observations and biology ----
years <- 1979:2018
ages <- 0:45
fishery <- "TRAWL_FISHERY_BC"
surveys <- c("NMFS_TRIENNIAL", "QCS_SYNOPTIC", "WCVI_SYNOPTIC")
inputs <- source$inputs

# The official model supplies fixed sex-specific weight schedules. Convert
# biomass to numbers over the joint sample before summing the sexes by age.
age_numbers <- function(composition, totals, weights, stream) {
  stream_year <- function(x) paste(x[[stream]], x$year)
  sex_age <- function(x) paste(x[[stream]], x$sex, x$age)
  total_key <- stream_year(totals)
  weight_key <- sex_age(weights)
  total_row <- match(stream_year(composition), total_key)
  w <- weights$value[match(sex_age(composition), weight_key)]
  stopifnot(nrow(composition) > 0L, !anyDuplicated(total_key), !anyDuplicated(weight_key),
    !anyDuplicated(paste(stream_year(composition), composition$sex, composition$age)),
    all(composition$basis == "proportion_numbers"), all(composition$unit == "proportion"),
    all(totals$unit == "t"), all(is.finite(totals$value) & totals$value > 0),
    all(grepl("^kg", weights$unit)), all(is.finite(w) & w > 0),
    all(is.finite(composition$value) & composition$value >= 0), !anyNA(total_row))
  groups <- split(seq_len(nrow(composition)), total_row)
  stopifnot(all(vapply(groups, function(i) {
    length(i) == 2L * (length(ages) - 1L) &&
      setequal(composition$sex[i], c("female", "male")) &&
      setequal(composition$age[i], ages[-1L])
  }, logical(1))))
  mean_mass <- ave(composition$value * w, total_row, FUN = sum)
  stopifnot(all(mean_mass > 0))
  composition$numbers <- totals$value[total_row] * composition$value / mean_mass
  composition$tonnes <- composition$numbers * w
  composition$samp_time <- totals$sampling_time[total_row]
  composition
}

catch_total <- subset(inputs, type == "catch" & measure == "total_biomass" &
  fleet == fishery & year %in% years)
catch_composition <- subset(inputs, type == "catch" & measure == "proportion_at_age" &
  fleet == fishery & year %in% years)
trawl_weight <- subset(inputs, type == "catch_weight" & measure == "weight_at_age" &
  fleet == fishery & season == "mid_year")
catch_numbers <- age_numbers(catch_composition, catch_total, trawl_weight, "fleet")
catch_age <- aggregate(cbind(numbers, tonnes) ~ year + age, catch_numbers, sum)
catch <- expand.grid(year = years, age = ages, KEEP.OUT.ATTRS = FALSE)
year_age <- function(x) paste(x$year, x$age)
i <- match(year_age(catch), year_age(catch_age))
catch$obs <- catch_age$numbers[i]
catch$unit <- "thousand fish"

index_total <- subset(inputs, type == "index" & measure == "total_biomass" &
  survey %in% surveys & year %in% years)
index_composition <- subset(inputs, type == "index" & measure == "proportion_at_age" &
  survey %in% surveys & year %in% years)
survey_weight <- subset(inputs, type == "weight" & measure == "weight_at_age" &
  survey %in% surveys & season == "mid_year")
index_numbers <- age_numbers(index_composition, index_total, survey_weight, "survey")
index <- aggregate(numbers ~ year + age + survey + samp_time, index_numbers, sum)
names(index)[names(index) == "numbers"] <- "obs"
index$unit <- "thousand fish"
index_sd <- subset(inputs, type == "index" & measure == "log_index_sd" &
  survey %in% surveys & year %in% years)
survey_year <- function(x) paste(x$survey, x$year)
stopifnot(!anyDuplicated(survey_year(index_sd)))
# Preserve the native aggregate precision as metadata. The age-specific SDs
# estimated below belong to the reconstructed indices.
index$relative_sd <- index_sd$value[match(survey_year(index), survey_year(index_sd))]
survey_time <- c(NMFS_TRIENNIAL = 0.5, QCS_SYNOPTIC = 0.5, WCVI_SYNOPTIC = 1 / 3)
stopifnot(all(is.finite(index$relative_sd) & index$relative_sd > 0),
  all(abs(index$samp_time - survey_time[index$survey]) < 1e-12),
  all(abs(index_composition$sampling_time - survey_time[index_composition$survey]) < 1e-12))

stock_weight <- subset(inputs, type == "weight" & measure == "weight_at_age" &
  season == "begin_year" & (is.na(survey) | !nzchar(survey)))
female_reproductive_weight <- subset(inputs, type == "biology" &
  measure == "maturity_fecundity_at_age" & sex == "female" & season == "spawning")
average_weight <- function(x) {
  stopifnot(nrow(x) == 2L * length(ages), !anyDuplicated(paste(x$age, x$sex)),
    setequal(x$age, ages), setequal(x$sex, c("female", "male")),
    all(grepl("^kg", x$unit)), all(is.finite(x$value) & x$value > 0))
  by_age <- tapply(x$value, x$age, mean)
  unname(by_age[as.character(catch$age)])
}
stock_kg <- average_weight(stock_weight)
catch_kg <- average_weight(trawl_weight)
positive <- which(is.finite(catch$obs) & catch$obs > 0)
catch_kg[positive] <- catch_age$tonnes[i[positive]] / catch$obs[positive]
stopifnot(nrow(female_reproductive_weight) == length(ages),
  setequal(female_reproductive_weight$age, ages),
  all(grepl("^kg", female_reproductive_weight$unit)))
female_kg <- female_reproductive_weight$value[match(catch$age, female_reproductive_weight$age)]
# Pooled abundance assumes half the fish at each age are female.
obs <- tinyAM::as_tam_data(list(
  catch = catch,
  index = index,
  weight = rbind(
    transform(catch[c("year", "age")], obs = stock_kg, unit = "kg", type = "stock"),
    transform(catch[c("year", "age")], obs = catch_kg, unit = "kg", type = "catch")),
  maturity = transform(catch[c("year", "age")], obs = 0.5 * female_kg / stock_kg,
    spawn_time = 0),
  pop_covariates = transform(catch[c("year", "age")], M_prior_group = factor("stock"))
))
catch_reporting <- list(weights = data.frame(year = catch$year, age = catch$age, weight = catch_kg),
  totals = data.frame(year = catch_total$year, yield = 1000 * catch_total$value))

## Official reference and dashboard background ----
source("analysis/comp_assessments/R/pool_source_states.R", local = TRUE)
age_plus_group <- max(ages)
source_history <- subset(source$outputs, year %in% years &
  measure %in% c("SSB", "numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age"))
# Use the source female numbers here so that the reference retains its native
# sex ratio. Otherwise the reference builder applies the pooled fit's biology.
source_ssb_at_age <- subset(source_history, measure == "numbers_at_age" & sex == "female")
source_ssb_at_age$value <- source_ssb_at_age$value * female_reproductive_weight$value[
  match(source_ssb_at_age$age, female_reproductive_weight$age)]
source_ssb_at_age$type <- "biomass"
source_ssb_at_age$measure <- "mature_biomass_at_age"
source_ssb_at_age$unit <- "t"
source_ssb_at_age$source_type <- "derived_source_output"
source_ssb_at_age$source_reference <- paste(source_ssb_at_age$source_reference,
  "and native effective reproductive weights at age")
source_ssb_at_age$notes <- "Native January female numbers multiplied by native effective reproductive weight at age."
native_ssb <- subset(source_history, measure == "SSB")
ssb_sum <- tapply(source_ssb_at_age$value, source_ssb_at_age$year, sum)
stopifnot(length(ssb_sum) == length(years), nrow(native_ssb) == length(years),
  all(abs(ssb_sum[as.character(native_ssb$year)] / native_ssb$value - 1) < 1e-5))
comparison_outputs <- rbind(pool_source_states(source_history, plus_age = age_plus_group),
  source_ssb_at_age)
comparison_aggregates <- "abundance"
comparison_definitions <- list(
  N = list(status = "approximate", reason = "Official N is summed over sexes."),
  F = list(status = "approximate", reason = "Official continuous F is weighted by numbers over sexes."),
  F_bar = list(status = "approximate", reason = "Both estimates average continuous F across ages 9–20, weighted by abundance."),
  ssb = list(status = "approximate", reason = paste(
    "Official January female SSB uses the native point fit.",
    "The tinyAM approximation assumes equal sex composition at each age."))
)

background <- c(
  "### Assessment and translation", "",
  print_sources(source$assessment), "",
  "The authors provide the [research document source and appendices](https://github.com/pbs-software/pbs-synth/tree/16f1d4368ad7d3e854219ca1c2d9d6353252c861/models/sweave/2024/YTR/RD). The [DFO Science Advisory Report 2025/026](https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs70-6/Fs70-6-2025-026-eng.pdf#page=7) cites the research document by Starr and Haigh as in preparation.", "",
  "| Component | Official assessment | tinyAM translation |",
  "|---|---|---|",
  "| Coverage | One coastwide stock in areas 3CD5ABCDE, two sexes, ages 0–45+, and 1935–2024. | One pooled population at ages 0–45+. Fit 1979–2018, the first to last year with trawl age compositions. |",
  "| Fishery | One pooled trawl fishery, with separate biomass totals and marginal age samples. | Expand the joint-sex compositions using the trawl biomass and midyear weights, then sum numbers by age. The two unsampled catch years retain missing ages. The dashboard shows annual reported totals alongside predicted yield. |",
  "| F | Static asymptotic age curves with male age offsets. The hybrid catch method solves annual F from biomass. | Fit a shared rising logistic curve with annual AR1 effects on log F. Baranov catches use reconstructed numbers and one estimated catch SD. Compare abundance-weighted F at ages 9–20, starting at the native summary-biomass age. |",
  "| Population | Sex-specific deterministic survival, Beverton–Holt recruitment and fixed log recruitment SD of 0.9. | Deterministic pooled survival, exponential initial age structure and IID recruitment with fixed SD of 0.9 about a constant mean. |",
  "| M | Constant M estimated separately by sex with the same natural-scale Normal prior (mean 0.08 and SD 0.024 per year) and bounds of 0.02–0.20. | Estimate one positive M with a Gaussian effect on log M. The lognormal prior matches the stated mean and SD on the natural scale and has no upper bound. |",
  "| Biology and SSB | Fixed growth and weight relationships, January female spawning biomass. | Use source-derived sex-specific weight schedules to expand age samples. Stock weight averages both sexes. SSB uses the female reproductive weights and a 1:1 sex ratio at age. |",
  "",
  "### Surveys and comparison", "",
  "Retain the six NMFS Triennial, eight Queen Charlotte Sound and six West Coast Vancouver Island index years with matched age samples in the fitted window. Estimate a separate logistic age pattern, scale and age-level lognormal SD for each survey. The West Coast Vancouver Island midpoint lies well beyond the observed ages and is poorly determined. The original aggregate log-SDs remain in the input table.",
  "",
  "The reconstructed age observations use the same annual biomass within each stream. The native model fits biomass and compositions separately and accounts for ageing error and composition weighting. The tinyAM likelihood treats the reconstructed ages as independent.",
  "",
  "The official comparison uses the native point fit for January female SSB, recruitment, abundance and F. Official F and M are weighted across the sexes by native numbers-at-age. The public SSB posterior is a separate reference."
)

## Model ----
# The lognormal prior approximation matches the native natural-scale mean and SD.
M_log_sd <- sqrt(log1p((0.024 / 0.08)^2))
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs, years = years, ages = ages,
    N_settings = list(process = "off", init = "exp"),
    rec_settings = list(process = "iid", mu_form = ~ 1, sd = 0.9),
    F_settings = list(process = "off", mu_form = ~ logistic(age) + ar1(year), mean_ages = 9:20),
    M_settings = list(process = "off", mu_supplied = ~ I(0.08 / sqrt(1 + (0.024 / 0.08)^2)),
      mu_form = ~ 0 + iid(M_prior_group, sd = M_log_sd)),
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(q_form = ~ 0 + survey + logistic(age, by = survey),
      sd_form = ~ 0 + survey),
    silent = silent
  )
}
