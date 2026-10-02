root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "database_to_tiny_obs.R"))
source(file.path(root, "R", "database_to_tiny_M.R"))
source(file.path(root, "R", "audit_assumptions.R"))
if (!requireNamespace("cli", quietly = TRUE)) {
  .translation_abort <- function(message) stop(message, call. = FALSE)
}

expect_equal <- function(actual, expected, tolerance = 1e-10) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop("Values differ: ", paste(capture.output(str(actual)), collapse = " "),
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
      sampling_time = 0.75)
))

sampling_times <- c(Longline = 0.6)
obs <- database_to_tiny_obs(
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
expect_equal(length(unique(obs$index$q_key)), 4L)
provenance <- attr(obs, "translation")$source_provenance
expect_equal(setequal(provenance$component, c("catch", "index", "weight", "maturity")), TRUE)
expect_equal(all(grepl("fixture table", provenance$source_reference)), TRUE)
expect_equal(obs$catch$obs[obs$catch$year == 2001 & obs$catch$age == 2], NA_real_)
expect_equal(obs$index$obs[obs$index$year == 2001], c(4, 0))

plus_obs <- database_to_tiny_obs(
  "translation_fixture", inputs, years = 2000:2001, ages = 1,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(plus_obs$catch$obs[plus_obs$catch$year == 2000], 30)
expect_equal(plus_obs$index$obs[plus_obs$index$year == 2000 &
                                 plus_obs$index$survey == "RV"],
             sum(expected_age_numbers))
expect_equal(plus_obs$index$obs[plus_obs$index$year == 2001], 4)

thousand_inputs <- inputs
thousand_inputs$unit[thousand_inputs$type == "catch"] <- "thousand fish"
thousand_obs <- database_to_tiny_obs(
  "translation_fixture", thousand_inputs, years = 2000:2001, ages = 1:2,
  weight_survey = "RV", sampling_times = sampling_times
)
expect_equal(thousand_obs$catch$obs[thousand_obs$catch$year == 2000],
             c(10000, 20000))
expect_equal(attr(thousand_obs, "translation")$catch_units$multiplier_to_fish,
             1000)

m_inputs <- inputs[0, ]
m_inputs <- rbind(
  m_inputs,
  row("M", "natural_mortality_at_age", "per_year", year = NA_real_,
      age = 1:2, value = c(0.2, 0.3), unit = "per_year", year_basis = NA_character_)
)
m <- database_to_tiny_M("translation_fixture", m_inputs, years = 2000:2001, ages = 1:2)
expect_equal(m$status, "fixed_numerical_input")
expect_equal(m$surface$M_assumption, c(0.2, 0.2, 0.3, 0.3))
expect_equal(m$M_settings$process, "off")

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

m_estimated <- database_to_tiny_M("dfo_cod_4t4vn_2019", inputs_db,
                                  assumptions_db, years = 1971:2018, ages = 2:11)
expect_equal(m_estimated$status, "estimated_in_source")
expect_equal(is.null(m_estimated$surface), TRUE)

sg_obs <- database_to_tiny_obs(
  "dfo_cod_4t4vn_2019", inputs_db, years = 1986:2018, ages = 2:11,
  weight_survey = "DFO September RV survey",
  sampling_times = c("DFO September RV survey" = 0.75,
                     "Mobile Sentinel August survey" = 0.625,
                     "Longline Sentinel survey" = 0.67),
  surveys = c("DFO September RV survey", "Mobile Sentinel August survey",
              "Longline Sentinel survey"),
  exclude_index_years = list("DFO September RV survey" = 2003)
)
expect_equal(nrow(sg_obs$weight), 330L)
expect_equal(sum(is.na(sg_obs$catch$obs)), 33L)
expect_equal(any(!is.finite(sg_obs$index$obs)), FALSE)
expect_equal(as.integer(table(sg_obs$index$survey)), c(320L, 161L, 160L))
expect_equal(any(sg_obs$index$survey == "DFO September RV survey" &
                   sg_obs$index$year == 2003), FALSE)
expect_equal(all(sg_obs$index$samp_time %in% c(0.75, 0.625, 0.67)), TRUE)

cat("Translation helper tests passed.\n")
