root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "database_to_tam_obs.R"))
source(file.path("R", "obs.R"))
source(file.path(root, "R", "audit_assumptions.R"))
if (!requireNamespace("cli", quietly = TRUE)) {
  .translation_abort <- function(message) stop(message, call. = FALSE)
}

expect_equal <- function(actual, expected, tolerance = 1e-10) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop("Values differ for ", deparse(substitute(actual)), ": ",
         paste(capture.output(str(actual)), collapse = " "),
         " expected ", paste(capture.output(str(expected)), collapse = " "),
         call. = FALSE)
  }
}

row <- function(type, measure, basis, survey = NA_character_, fleet = NA_character_,
                year, age, value, unit, sampling_time = NA_real_,
                year_basis = "calendar_year") {
  data.frame(
    assessment_id = "translation_fixture", type, measure, basis, fleet, survey,
    sex = NA_character_, region = NA_character_, season = NA_character_,
    year, year_basis, age, value, unit, sampling_time,
    source_type = "official_table", source_reference = "fixture table",
    transformation = NA_character_, notes = NA_character_,
    stringsAsFactors = FALSE
  )
}

inputs <- do.call(rbind, list(
  row("weight", "weight_at_age", "kg_per_fish", survey = "RV",
      year = 2000:2001, age = 1, value = c(1, 1), unit = "kg/fish"),
  row("weight", "weight_at_age", "kg_per_fish", survey = "RV",
      year = 2000:2001, age = 2, value = c(2, 2), unit = "kg/fish"),
  row("maturity", "maturity_at_age", "proportion", year = 2000:2001,
      age = 1, value = c(0.1, 0.1), unit = "proportion"),
  row("maturity", "maturity_at_age", "proportion", year = 2000:2001,
      age = 2, value = c(0.8, 0.8), unit = "proportion"),
  row("catch", "numbers_at_age", "numbers", fleet = "fishery",
      year = 2000, age = 1:2, value = c(10, 20), unit = "fish"),
  row("catch", "numbers_at_age", "numbers", fleet = "fishery",
      year = 2001, age = 1, value = 30, unit = "fish"),
  row("index", "proportion_at_age", "proportion_numbers", survey = "RV",
      year = 2000, age = 1:2, value = c(0.25, 0.75), unit = "proportion",
      sampling_time = 0.75),
  row("index", "total_biomass", "biomass", survey = "RV",
      year = 2000, age = NA_real_, value = 1, unit = "t per tow",
      sampling_time = 0.75),
  row("index", "proportion_at_age", "proportion_numbers", survey = "Longline",
      year = 2000, age = 1:2, value = c(0.25, 0.75), unit = "proportion"),
  row("index", "total_biomass", "biomass", survey = "Longline",
      year = 2000, age = NA_real_, value = 1, unit = "t per tow"),
  row("index", "numbers_at_age", "numbers", survey = "RV",
      year = 2001, age = 1:2, value = c(4, 0), unit = "fish per 1,000 hooks",
      sampling_time = 0.75),
  row("index", "numbers_at_age", "numbers", survey = "Acoustic",
      year = 2000, age = 1:2, value = c(40, 50), unit = "native survey index",
      sampling_time = 0.25),
  row("index", "numbers_at_age", "numbers", survey = "Relative",
      year = 2000, age = 1:2, value = c(2, 3), unit = "survey_index",
      sampling_time = 0.4)
))

sampling_times <- c(Longline = 0.6)
obs <- database_to_tam_obs(
  "translation_fixture", inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expected_age_numbers <- 1000 * c(0.25, 0.75) / sum(c(0.25, 0.75) * c(1, 2))
expect_equal(obs$index$obs[obs$index$year == 2000 & obs$index$survey == "RV"],
             expected_age_numbers)
expect_equal(obs$index$obs[obs$index$year == 2000 & obs$index$survey == "Longline"],
             expected_age_numbers)
expect_equal(obs$index$samp_time[obs$index$survey == "RV"], rep(0.75, 4))
expect_equal(obs$index$samp_time[obs$index$survey == "Longline"], rep(0.6, 2))
expect_equal(as.character(obs$index$q_block), as.character(obs$index$age))
expect_equal(length(unique(obs$index$q_key)), 8L)
expect_equal(obs$index$obs[obs$index$survey == "Acoustic"], c(40, 50))
expect_equal(obs$index$obs[obs$index$survey == "Relative"], c(2, 3))
provenance <- attr(obs, "translation")$source_provenance
expect_equal(setequal(provenance$component, c("catch", "index", "weight", "maturity")), TRUE)
expect_equal(all(grepl("fixture table", provenance$source_reference)), TRUE)
expect_equal(grepl("retained on their source scale",
                   provenance$method[provenance$component == "index" &
                                       provenance$survey == "Acoustic"]), TRUE)
expect_equal(grepl("retained on their source scale",
                   provenance$method[provenance$component == "index" &
                                       provenance$survey == "Relative"]), TRUE)
expect_equal(obs$catch$obs[obs$catch$year == 2001 & obs$catch$age == 2], NA_real_)
expect_equal(obs$index$obs[obs$index$year == 2001], c(4, 0))
constant_biology <- database_to_tam_obs(
  "translation_fixture", inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times,
  weight_reference_year = 2000, maturity_reference_year = 2000,
  maturity_multiplier = 0.5
)
expect_equal(constant_biology$weight$obs[constant_biology$weight$year == 2001], c(1, 2))
expect_equal(constant_biology$maturity$obs[constant_biology$maturity$year == 2001], c(0.05, 0.4))
expect_equal(attr(constant_biology, "translation")$maturity_multiplier, 0.5)

static_inputs <- inputs[!(inputs$type == "weight" |
                            inputs$type == "maturity"), , drop = FALSE]
static_inputs <- rbind(
  static_inputs,
  row("weight", "weight_at_age", "kg_per_fish", survey = "RV",
      year = NA_real_, year_basis = NA_character_, age = 1:2,
      value = c(1, 2), unit = "kg/fish"),
  row("maturity", "maturity_at_age", "proportion",
      year = NA_real_, year_basis = NA_character_, age = 1:2,
      value = c(0.1, 0.8), unit = "proportion"),
  row("M", "natural_mortality_at_age", "per_year",
      year = NA_real_, year_basis = NA_character_, age = 1:2,
      value = c(0.2, 0.3), unit = "per year")
)
static_obs <- database_to_tam_obs(
  "translation_fixture", static_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(static_obs$weight$obs[static_obs$weight$year == 2001], c(1, 2))
expect_equal(static_obs$maturity$obs[static_obs$maturity$year == 2001], c(0.1, 0.8))
expect_equal(static_obs$weight$M_assumption[static_obs$weight$age == 1], rep(0.2, 2))
expect_equal(static_obs$weight$M_assumption[static_obs$weight$age == 2], rep(0.3, 2))
expect_equal(grepl("Time-invariant source maturity vector expanded",
                   attr(static_obs, "translation")$source_provenance$method[
                     attr(static_obs, "translation")$source_provenance$component == "maturity"]), TRUE)

stock_weight_inputs <- rbind(
  inputs,
  row("weight", "weight_at_age", "kg_per_fish", year = 2000:2001,
      age = 1, value = c(10, 10), unit = "kg/fish"),
  row("weight", "weight_at_age", "kg_per_fish", year = 2000:2001,
      age = 2, value = c(20, 20), unit = "kg/fish")
)
stock_weight_obs <- database_to_tam_obs(
  "translation_fixture", stock_weight_inputs, years = 2000:2001, ages = 1:2,
  sampling_times = sampling_times
)
expect_equal(stock_weight_obs$weight$obs[stock_weight_obs$weight$year == 2000],
             c(10, 20))
expect_equal(stock_weight_obs$index$obs[
  stock_weight_obs$index$survey == "RV" & stock_weight_obs$index$year == 2000
], expected_age_numbers)
expected_longline <- 1000 * c(0.25, 0.75) / sum(c(0.25, 0.75) * c(10, 20))
expect_equal(stock_weight_obs$index$obs[
  stock_weight_obs$index$survey == "Longline" & stock_weight_obs$index$year == 2000
], expected_longline)
stock_index_weight_obs <- database_to_tam_obs(
  "translation_fixture", stock_weight_inputs, years = 2000:2001, ages = 1:2,
  sampling_times = sampling_times, index_weight_source = "stock"
)
expect_equal(stock_index_weight_obs$index$obs[
  stock_index_weight_obs$index$survey == "RV" & stock_index_weight_obs$index$year == 2000
], expected_longline)

index_sd_inputs <- rbind(
  inputs,
  row("index", "index_sd", "index_scale", survey = "RV", year = 2000,
      age = NA_real_, value = 0.1, unit = "t per tow")
)
index_sd_obs <- database_to_tam_obs(
  "translation_fixture", index_sd_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(index_sd_obs$index$obs, obs$index$obs)

relative_error_inputs <- rbind(
  inputs,
  row("index", "relative_standard_error", "relative_scale", survey = "RV",
      year = 2000, age = 1:2, value = c(0.1, 0.2), unit = "relative standard error")
)
relative_error_obs <- database_to_tam_obs(
  "translation_fixture", relative_error_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(relative_error_obs$index$obs, obs$index$obs)

log_index_sd_inputs <- rbind(
  inputs,
  row("index", "log_index_sd", "log_scale", survey = "RV", year = 2000,
      age = NA_real_, value = 0.1, unit = "log scale")
)
log_index_sd_obs <- database_to_tam_obs(
  "translation_fixture", log_index_sd_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(log_index_sd_obs$index$relative_sd[
  log_index_sd_obs$index$survey == "RV" & log_index_sd_obs$index$year == 2000
], rep(0.1, 2))

plus_obs <- database_to_tam_obs(
  "translation_fixture", inputs, years = 2000:2001, ages = 1,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(plus_obs$catch$obs[plus_obs$catch$year == 2000], 30)
expect_equal(plus_obs$index$obs[plus_obs$index$year == 2000 &
                                 plus_obs$index$survey == "RV"],
             sum(expected_age_numbers))
expect_equal(plus_obs$index$obs[plus_obs$index$survey == "Acoustic"], 90)
expect_equal(plus_obs$index$obs[plus_obs$index$survey == "Relative"], 5)
expect_equal(plus_obs$index$obs[plus_obs$index$year == 2001], 4)

catch_proportion_inputs <- inputs[inputs$type != "catch", , drop = FALSE]
catch_proportion_inputs <- rbind(
  catch_proportion_inputs,
  row("catch", "proportion_at_age", "proportion_numbers", fleet = "fishery",
      year = 2000, age = 1:2, value = c(0.2, 0.7), unit = "proportion"),
  row("catch", "total_numbers", "numbers", fleet = "fishery",
      year = 2000, age = NA_real_, value = 400, unit = "thousand fish")
)
catch_proportion_obs <- database_to_tam_obs(
  "translation_fixture", catch_proportion_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(catch_proportion_obs$catch$obs[catch_proportion_obs$catch$year == 2000],
             c(80000, 280000))
expect_equal(all(is.na(catch_proportion_obs$catch$obs[
  catch_proportion_obs$catch$year == 2001])), TRUE)
expect_equal(grepl("without renormalizing",
                   attr(catch_proportion_obs, "translation")$catch_method), TRUE)
missing_catch_total <- catch_proportion_inputs[
  catch_proportion_inputs$measure != "total_numbers", , drop = FALSE]
catch_error <- tryCatch(database_to_tam_obs(
  "translation_fixture", missing_catch_total, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
), error = identity)
expect_equal(grepl("matching annual total_numbers", conditionMessage(catch_error)), TRUE)

missing_index_total <- inputs[!(inputs$type == "index" &
                                  inputs$measure == "total_biomass"), , drop = FALSE]
index_error <- tryCatch(database_to_tam_obs(
  "translation_fixture", missing_index_total, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
), error = identity)
expect_equal(grepl("no matching total index values", conditionMessage(index_error)), TRUE)

biomass_catch <- catch_proportion_inputs
biomass_catch$measure[biomass_catch$measure == "total_numbers"] <- "total_biomass"
biomass_catch$unit[biomass_catch$measure == "total_biomass"] <- "thousand t"
biomass_catch$value[biomass_catch$measure == "total_biomass"] <- 0.001
biomass_catch <- rbind(biomass_catch,
  row("catch_weight", "weight_at_age", "kg_per_fish", fleet = "fishery",
      year = 2000, age = 1:2, value = c(1000, 2000), unit = "g/fish"))
bc <- .translation_catch_at_age(biomass_catch, 2000:2001, 1:2)$catch
expect_equal(bc$obs[bc$year == 2000], 1000 * c(0.2, 0.7) / 1.6)
expect_equal(sum(bc$obs[bc$year == 2000] * c(1, 2)), 1000)
restricted_bc <- .translation_catch_at_age(biomass_catch, 2000:2001, 2)$catch
expect_equal(restricted_bc$obs[restricted_bc$year == 2000], 1000 * 0.7 / 1.6)
biomass_catch$basis[biomass_catch$type == "catch" &
                      biomass_catch$measure == "proportion_at_age"] <- "proportion_biomass"
bc <- .translation_catch_at_age(biomass_catch, 2000:2001, 1:2)$catch
expect_equal(bc$obs[bc$year == 2000], c(200, 350))
biomass_catch$measure[biomass_catch$measure == "total_biomass"] <- "total_numbers"
biomass_catch$unit[biomass_catch$measure == "total_numbers"] <- "fish"
biomass_catch$value[biomass_catch$measure == "total_numbers"] <- 550
bc <- .translation_catch_at_age(biomass_catch, 2000:2001, 1:2)$catch
expect_equal(bc$obs[bc$year == 2000], c(200, 350))
missing_catch_weights <- biomass_catch[biomass_catch$type != "catch_weight", ]
weight_error <- tryCatch(.translation_catch_at_age(missing_catch_weights, 2000:2001, 1:2),
                         error = identity)
expect_equal(grepl("matching catch weights", conditionMessage(weight_error)), TRUE)
stock_weights_catch <- biomass_catch
stock_weights_catch$type[stock_weights_catch$type == "catch_weight"] <- "weight"
stock_weight_error <- tryCatch(
  .translation_catch_at_age(stock_weights_catch, 2000:2001, 1:2),
  error = identity
)
expect_equal(grepl("stock weights are not substituted",
                   conditionMessage(stock_weight_error)), TRUE)

thousand_inputs <- inputs
thousand_inputs$unit[thousand_inputs$type == "catch"] <- "thousand fish"
thousand_obs <- database_to_tam_obs(
  "translation_fixture", thousand_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(thousand_obs$catch$obs[thousand_obs$catch$year == 2000],
             c(10000, 20000))
expect_equal(attr(thousand_obs, "translation")$catch_units$multiplier_to_fish,
             1000)
expect_equal(unname(vapply(c("native survey index", "survey_index"),
                           .translation_index_multiplier, numeric(1))), c(1, 1))

m_inputs <- inputs[0, ]
m_inputs <- rbind(
  m_inputs,
  row("M", "natural_mortality_at_age", "per_year", year = NA_real_,
      age = 1:2, value = c(0.2, 0.3), unit = "per_year", year_basis = NA_character_)
)
m <- database_to_tam_M("translation_fixture", m_inputs, years = 2000:2001, ages = 1:2)
expect_equal(m$status, "fixed_numerical_input")
expect_equal(m$surface$M_assumption, c(0.2, 0.2, 0.3, 0.3))
expect_equal(m$M_settings$process, "off")

fixed_m_inputs <- rbind(
  inputs,
  row("M", "natural_mortality_at_age", "per_year", year = NA_real_,
      age = 1:2, value = c(0.2, 0.3), unit = "per year",
      year_basis = NA_character_)
)
fixed_m_inputs_before <- fixed_m_inputs
obs_with_m <- database_to_tam_obs(
  "translation_fixture", fixed_m_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(obs_with_m$weight$M_assumption, c(0.2, 0.2, 0.3, 0.3))
expect_equal(attr(obs_with_m, "translation")$M$status,
             "fixed_numerical_input")
expect_equal(sum(attr(obs_with_m, "translation")$source_provenance$component == "M"),
             1L)
expect_equal(fixed_m_inputs, fixed_m_inputs_before)

read_table <- function(name) {
  read.csv(file.path(root, "database", name), stringsAsFactors = FALSE,
           na.strings = c("", "NA"), check.names = FALSE)
}
assumptions_db <- read_table("assumptions.csv")
inputs_db <- read_table("inputs.csv")
audit <- audit_assumptions("dfo_cod_4t4vn_2019", assumptions_db)
expect_equal(audit$tinyam_support[audit$component == "M" &
                                    audit$setting == "process"],
             "partially_supported")
expect_equal(audit$tinyam_support[audit$component == "M" &
                                    audit$setting == "initial_priors"],
             "unsupported")
expect_equal(audit$tinyam_support[audit$component == "index" &
                                    audit$setting == "sampling_time" &
                                    grepl("unknown", audit$value, ignore.case = TRUE)],
             "unsupported")

m_estimated <- database_to_tam_M("dfo_cod_4t4vn_2019", inputs_db,
                                  assumptions_db, years = 1971:2018, ages = 2:11)
expect_equal(m_estimated$status, "estimated_in_source")
expect_equal(is.null(m_estimated$surface), TRUE)

# Raw landings require the documented stock-recipe catch proxy.
sg_error <- tryCatch(database_to_tam_obs(
  "dfo_cod_4t4vn_2019", inputs_db, years = 1986:2018, ages = 2:11,
  weight_survey = "DFO September RV survey",
  sampling_times = c("DFO September RV survey" = 0.75),
  surveys = "DFO September RV survey"
), error = identity)
expect_equal(inherits(sg_error, "error"), TRUE)
expect_equal(grepl("No source catch-at-age", conditionMessage(sg_error)), TRUE)

nea_id <- "ices_cod_northeast_arctic_2026"
nea_obs <- database_to_tam_obs(nea_id, inputs_db,
                                years = 1946:2026, ages = 3:15)
nea_m <- database_to_tam_M(nea_id, inputs_db, assumptions_db,
                            years = 1946:2026, ages = 3:15)
expect_equal(nrow(nea_obs$catch), 1053L)
expect_equal(sum(!is.na(nea_obs$catch$obs)), 1028L)
expect_equal(nrow(nea_obs$index), 1349L)
expect_equal(sort(as.integer(table(nea_obs$index$survey))),
             sort(c(290L, 130L, 403L, 336L, 190L)))
expect_equal(range(nea_obs$index$age), c(3, 12))
expect_equal(sum(nea_obs$index$age %in% 13:15), 0L)
expect_equal(range(nea_obs$catch$age), c(3, 15))
nea_catch <- inputs_db[inputs_db$assessment_id == nea_id &
                         inputs_db$type == "catch" &
                         inputs_db$measure == "numbers_at_age", , drop = FALSE]
nea_catch_match <- match(
  paste(nea_catch$year, nea_catch$age),
  paste(nea_obs$catch$year, nea_obs$catch$age)
)
expect_equal(nea_obs$catch$obs[nea_catch_match],
             as.numeric(nea_catch$value) * 1000)
nea_direct <- inputs_db[inputs_db$assessment_id == nea_id &
                           inputs_db$type == "index" &
                           inputs_db$measure == "numbers_at_age", , drop = FALSE]
nea_index_match <- match(
  paste(nea_direct$year, nea_direct$age, nea_direct$survey),
  paste(nea_obs$index$year, nea_obs$index$age, nea_obs$index$survey)
)
expect_equal(nea_obs$index$obs[nea_index_match], as.numeric(nea_direct$value))
nea_provenance <- attr(nea_obs, "translation")$source_provenance
expect_equal(all(grepl("retained on their source scale",
                       nea_provenance$method[nea_provenance$component == "index"])), TRUE)
expect_equal(nea_m$status, "fixed_numerical_input")
expect_equal(nea_m$M_settings$process, "off")
expect_equal(nrow(nea_m$surface), 1053L)
nea_m_source <- inputs_db[inputs_db$assessment_id == nea_id &
                            inputs_db$type == "M" &
                            inputs_db$measure == "natural_mortality_at_age", , drop = FALSE]
nea_m_match <- match(paste(nea_m_source$year, nea_m_source$age),
                     paste(nea_m$surface$year, nea_m$surface$age))
expect_equal(nea_m$surface$M_assumption[nea_m_match],
             as.numeric(nea_m_source$value))
for (survey in unique(nea_obs$index$survey)) {
  source_time <- unique(inputs_db$sampling_time[
    inputs_db$assessment_id == nea_id & inputs_db$type == "index" &
      inputs_db$survey == survey & !is.na(inputs_db$sampling_time)])
  expect_equal(unique(nea_obs$index$samp_time[nea_obs$index$survey == survey]),
               as.numeric(source_time))
}

herring_inputs <- inputs_db[inputs_db$assessment_id == "ices_herring_north_sea_2026", , drop = FALSE]
herring_error <- tryCatch(.translation_index(
  herring_inputs, data.frame(), data.frame(), years = 1947:2026, ages = 0:8,
  sampling_times = NULL
), error = identity)
expect_equal(grepl("Unsupported index measure", conditionMessage(herring_error)), TRUE)

herring_surveys <- c("HERAS", "IBTS0", "IBTS-Q1", "IBTS-Q3")
herring_missing_time <- inputs_db
herring_missing_time$sampling_time[
  herring_missing_time$assessment_id == "ices_herring_north_sea_2026" &
    herring_missing_time$survey %in% herring_surveys
] <- NA_real_
herring_timing_error <- tryCatch(database_to_tam_obs(
  "ices_herring_north_sea_2026", herring_missing_time,
  years = 1947:2025, ages = 0:8, surveys = herring_surveys
), error = identity)
expect_equal(grepl("Survey timing is unknown", conditionMessage(herring_timing_error)), TRUE)

herring_before <- herring_inputs
herring_obs <- database_to_tam_obs(
  "ices_herring_north_sea_2026", inputs_db, years = 1947:2025, ages = 0:8,
  surveys = herring_surveys
)
expect_equal(setequal(unique(herring_obs$index$survey), herring_surveys), TRUE)
expect_equal(sort(unique(herring_obs$index$samp_time)), c(0.125, 0.5, 0.625))
expect_equal(setequal(attr(herring_obs, "translation")$excluded_surveys,
                      c("LAI-SNS", "LAI-CNS", "LAI-BUN", "LAI-ORSH")), TRUE)
herring_direct <- herring_inputs[herring_inputs$type == "index" &
                                   herring_inputs$survey %in% herring_surveys &
                                   herring_inputs$year <= 2025, , drop = FALSE]
herring_match <- match(paste(herring_direct$year, herring_direct$age, herring_direct$survey),
                       paste(herring_obs$index$year, herring_obs$index$age,
                             herring_obs$index$survey))
expect_equal(herring_obs$index$obs[herring_match], as.numeric(herring_direct$value))
expect_equal(herring_inputs, herring_before)
expect_equal(sum(is.na(herring_obs$catch$obs)), 18L)

source(file.path(root, "R", "read_database.R"))
working_db <- read_database()
db <- read_committed_database()
expect_equal(working_db$database_source, "working_tree")
expect_equal(db$database_source, "committed")
expect_equal(grepl("^[[:xdigit:]]{40}$", working_db$commit), TRUE)

fixture_dir <- tempfile("assessment_database_")
dir.create(fixture_dir)
database_files <- c("stocks.csv", "assessments.csv", "assumptions.csv",
                    "inputs.csv", "outputs.csv")
for (name in database_files) {
  write.csv(data.frame(source_file = name), file.path(fixture_dir, name),
            row.names = FALSE)
}
fixture_db <- read_database(fixture_dir)
expect_equal(unname(vapply(fixture_db[c("stocks", "assessments", "assumptions",
                                       "inputs", "outputs")],
                           function(table) table$source_file[[1]], character(1))),
             database_files)

working_ebs <- read_assessment("afsc_pollock_ebs_2024", working_db)
ebs <- read_committed_assessment("afsc_pollock_ebs_2024", database = db)
expect_equal(nrow(ebs$assessment), 1L)
expect_equal(all(ebs$inputs$assessment_id == "afsc_pollock_ebs_2024"), TRUE)
expect_equal(grepl("^[[:xdigit:]]{40}$", ebs$commit), TRUE)
expect_equal(working_ebs$inputs, ebs$inputs)
expect_equal(read_assessment("afsc_pollock_ebs_2024", db)$inputs, ebs$inputs)
ebs_surveys <- c("NMFS bottom-trawl VAST", "NMFS acoustic-trawl",
                 "NMFS acoustic-trawl age-1 index")
ebs_obs <- database_to_tam_obs(
  "afsc_pollock_ebs_2024", ebs$inputs, years = 1964:2024, ages = 1:15,
  surveys = ebs_surveys, maturity_reference_year = 1964,
  maturity_multiplier = 0.5,
  index_weight_source = "stock",
  assumptions = ebs$assumptions
)
expect_equal(dim(ebs_obs$catch), c(915L, 3L))
expect_equal(sum(is.na(ebs_obs$catch$obs)), 15L)
expect_equal(dim(ebs_obs$maturity), c(915L, 3L))
expect_equal(ebs_obs$weight$M_assumption[ebs_obs$weight$age == 1],
             rep(0.9, 61))
expect_equal(ebs_obs$weight$M_assumption[ebs_obs$weight$age == 2],
             rep(0.45, 61))
expect_equal(ebs_obs$weight$M_assumption[ebs_obs$weight$age >= 3],
             rep(0.3, 61 * 13))
expect_equal(attr(ebs_obs, "translation")$M$status,
             "fixed_numerical_input")
ebs_maturity <- ebs$inputs[ebs$inputs$type == "maturity", , drop = FALSE]
expect_equal(ebs_obs$maturity$obs[ebs_obs$maturity$year == 1964],
             ebs_maturity$value[match(1:15, ebs_maturity$age)] * 0.5)
expect_equal(setequal(unique(ebs_obs$index$survey), ebs_surveys), TRUE)
expect_equal(sum(ebs_obs$index$survey == "NMFS bottom-trawl VAST"), 630L)
expect_equal(sum(ebs_obs$index$survey == "NMFS acoustic-trawl"), 266L)
expect_equal(sum(ebs_obs$index$survey == "NMFS acoustic-trawl age-1 index"), 18L)
expect_equal(any(ebs_obs$index$survey == "NMFS acoustic-trawl age-1 index" &
                   ebs_obs$index$year == 2024), FALSE)
ebs_index_method <- attr(ebs_obs, "translation")$index
expect_equal(all(grepl("stock weights as an approximation",
                       ebs_index_method$method[ebs_index_method$survey %in%
                                                 c("NMFS bottom-trawl VAST",
                                                   "NMFS acoustic-trawl")])), TRUE)
ebs_catch_weight <- ebs$inputs[ebs$inputs$type == "catch_weight" &
                                 ebs$inputs$year < 2024, , drop = FALSE]
ebs_catch_weight$key <- paste(ebs_catch_weight$year, ebs_catch_weight$age)
ebs_catch <- ebs_obs$catch
ebs_catch$key <- paste(ebs_obs$catch$year, ebs_obs$catch$age)
weighted_catch <- ebs_obs$catch$obs[match(ebs_catch_weight$key, ebs_catch$key)] *
  ebs_catch_weight$value
catch_biomass <- tapply(weighted_catch, ebs_catch_weight$year, sum)
total_biomass <- ebs$inputs[ebs$inputs$type == "catch" &
                              ebs$inputs$measure == "total_biomass" &
                              ebs$inputs$year < 2024, , drop = FALSE]
expect_equal(as.numeric(catch_biomass),
             as.numeric(total_biomass$value[
               match(names(catch_biomass), total_biomass$year)] * 1e6))

source(file.path(root, "R", "database_to_tam_ref.R"))
ebs_reference <- database_to_tam_ref(
  "afsc_pollock_ebs_2024", ebs$outputs, obs = ebs_obs,
  years = 1964:2024, ages = 1:15, terminal_year = 2024
)
expect_equal(inherits(ebs_reference, "tam_ref"), TRUE)
expect_equal("source_fit" %in% names(ebs_reference), FALSE)
expect_equal(nrow(ebs_reference$pop$ssb), 61L)
expect_equal(nrow(ebs_reference$pop$M), 915L)
expect_equal(all(is.na(ebs_reference$pop$M$se)), TRUE)
expect_equal(all(is.na(ebs_reference$obs_pred$catch$pred)), TRUE)
expect_equal(all(is.na(ebs_reference$pop$ssb$se_scale)), TRUE)

ns_haddock_id <- "ices_haddock_north_sea_2026"
ns_haddock <- read_committed_assessment(ns_haddock_id, database = db)
ns_haddock_obs <- database_to_tam_obs(
  ns_haddock_id, ns_haddock$inputs,
  years = 1972:2026, ages = 0:8, assumptions = ns_haddock$assumptions
)
expect_equal(check_obs(ns_haddock_obs), TRUE)
expect_equal(nrow(ns_haddock_obs$index), 667L)
expect_equal(all(is.finite(ns_haddock_obs$index$relative_sd)), TRUE)
expect_equal(all(abs(ns_haddock_obs$index$relative_sd^2 *
                       ns_haddock_obs$index$relative_precision_weight - 1) < 1e-10), TRUE)

nea_reference <- database_to_tam_ref(
  "ices_cod_northeast_arctic_2026", db$outputs,
  obs = nea_obs, years = 1946:2026, ages = 3:12,
  terminal_year = 2026, age_plus_group = 12
)
source_n_plus <- db$outputs[db$outputs$assessment_id == "ices_cod_northeast_arctic_2026" &
                              db$outputs$type == "population" &
                              db$outputs$measure == "numbers_at_age" &
                              db$outputs$year == 2026 & db$outputs$age >= 12, , drop = FALSE]
expect_equal(nea_reference$pop$N$est[
  nea_reference$pop$N$year == 2026 & nea_reference$pop$N$age == 12
], sum(as.numeric(source_n_plus$value)))

template_years <- 2000:2002
template_ages <- 1:2
template_grid <- expand.grid(year = template_years, age = template_ages)
template_pop <- list(
  ssb = data.frame(year = template_years, est = NA_real_, se = NA_real_,
                   se_scale = NA_character_, lwr = NA_real_, upr = NA_real_,
                   is_proj = FALSE),
  N = data.frame(template_grid, est = NA_real_, is_proj = FALSE),
  F = data.frame(template_grid, est = NA_real_, is_proj = FALSE),
  M = data.frame(template_grid, est = NA_real_, is_proj = FALSE),
  recruitment = data.frame(year = template_years, est = NA_real_,
                           is_proj = FALSE)
)
template_obs <- list(
  catch = data.frame(year = rep(template_years, each = 2),
                     age = rep(template_ages, length(template_years)),
                     fleet = "fishery", obs = 10:15),
  index = data.frame(year = rep(template_years, each = 2),
                     age = rep(template_ages, length(template_years)),
                     survey = "RV", obs = 20:25),
  weight = data.frame(year = rep(template_years, each = 2), age = rep(template_ages, 3),
                      M_assumption = 0.2),
  maturity = data.frame(year = rep(template_years, each = 2), age = rep(template_ages, 3),
                        obs = 0.5)
)
attr(template_obs, "translation") <- list(M = list(status = "fixed_numerical_input"))

template_fit <- list(
  call = quote(fit_tam()), refit_args = list(),
  dat = list(obs = template_obs, years = template_years,
             ages = template_ages, is_proj = rep(FALSE, 3)),
  obj = list(fn = function(x) x), opt = list(par = 1), rep = list(
    ssb = setNames(rep(NA_real_, 3), template_years),
    recruitment = setNames(rep(NA_real_, 3), template_years),
    N = matrix(NA_real_, 3, 2,
               dimnames = list(year = template_years, age = template_ages)),
    M = matrix(NA_real_, 3, 2,
               dimnames = list(year = template_years, age = template_ages))
  ),
  sdrep = NULL,
  fixed_par = data.frame(par = "log_sd_catch", coef = "fishery", est = 1,
                         se = 0.1, se_scale = "sd", lwr = 0.8, upr = 1.2),
  random_par = list(log_f = data.frame(year = template_years[1:2], age = 1L,
                                       est = 1:2)),
  obs_pred = list(
    catch = data.frame(year = template_obs$catch$year, age = template_obs$catch$age,
                       fleet = template_obs$catch$fleet, obs = template_obs$catch$obs,
                       pred = 1:6, sd = 1, std_res = 0, osa_res = 0),
    index = data.frame(year = template_obs$index$year, age = template_obs$index$age,
                       survey = template_obs$index$survey, obs = template_obs$index$obs,
                       pred = 1:6, sd = 1, std_res = 0, osa_res = 0, q = 1)
  ),
  pop = template_pop, is_converged = TRUE, grad_tol = 0.01
)
source_outputs <- data.frame(
  assessment_id = "translation_fixture",
  type = c("biomass", "population", "mortality", "mortality"),
  measure = c("SSB", "numbers_at_age", "fishing_mortality_at_age",
              "fishing_mortality_at_age"),
  year = c(2000, 2000, 2000, 2000), age = c(NA, 1, 1, 2), age_group = NA_character_,
  fleet = c(NA_character_, NA_character_, "Residual catch", "Residual catch"),
  value = c(900, 100, 0.1, 0.2), se = NA_real_, lwr = NA_real_, upr = NA_real_,
  unit = c("t", "fish", "per year", "per year"), source_type = "official_table",
  source_reference = "fixture", notes = NA_character_
)
template_reference <- database_to_tam_ref(
  "translation_fixture", source_outputs, obs = template_obs,
  years = template_years, ages = template_ages, terminal_year = 2002,
  comparison_scales = c(ssb = 1e-3), template = template_fit
)
expect_equal(names(template_reference), c(names(template_fit), "comparison_scales"))
expect_equal(names(template_reference$pop), names(template_fit$pop))
expect_equal(names(template_reference$rep), names(template_fit$rep))
expect_equal(names(template_reference$obs_pred), names(template_fit$obs_pred))
expect_equal(names(template_reference$random_par), names(template_fit$random_par))
expect_equal(nrow(template_reference$pop$N), nrow(template_fit$pop$N))
expect_equal(nrow(template_reference$random_par$log_f),
             nrow(template_fit$random_par$log_f))
expect_equal(template_reference$random_par$log_f$year,
             template_fit$random_par$log_f$year)
expect_equal(all(is.na(template_reference$random_par$log_f$est)), TRUE)
expect_equal(template_reference$pop$N$est[template_reference$pop$N$year == 2000 &
                                            template_reference$pop$N$age == 1], 100)
expect_equal(template_reference$pop$F$est[template_reference$pop$F$year == 2000 &
                                            template_reference$pop$F$age == 1], 0.1)
expect_equal(template_reference$pop$F$est[template_reference$pop$F$year == 2000 &
                                            template_reference$pop$F$age == 2], 0.2)
expect_equal(all(is.na(template_reference$pop$F$est[
  template_reference$pop$F$year != 2000
])), TRUE)
expect_equal(all(is.na(template_reference$pop$N$est[
  !(template_reference$pop$N$year == 2000 & template_reference$pop$N$age == 1)
])), TRUE)
expect_equal(template_reference$pop$M$est, rep(0.2, 6))
expect_equal(template_reference$rep$ssb, c(`2000` = 900, `2001` = NA, `2002` = NA))
expect_equal(template_reference$rep$N["2000", "1"], 100)
expect_equal(template_reference$rep$N["2001", "1"], NA_real_)
expect_equal(template_reference$rep$M[, "2"],
             setNames(rep(0.2, 3), as.character(template_years)))
expect_equal(template_reference$dat$obs, template_obs)
expect_equal(template_reference$obs_pred$catch$pred, rep(NA_real_, 6))
expect_equal(template_reference$fixed_par$est, NA_real_)
expect_equal(template_fit$fixed_par$est, 1)
expect_equal(template_reference$comparison_scales, c(ssb = 1e-3))
cat("Translation helper tests passed.\n")
