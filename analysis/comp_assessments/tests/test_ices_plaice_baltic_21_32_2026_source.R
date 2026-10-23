source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))

assessment_id <- "ices_plaice_baltic_21_32_2026"
database <- read_database()
source_record <- read_assessment(assessment_id, database)
inputs <- source_record$inputs
outputs <- source_record$outputs
fit_env <- new.env(parent = baseenv())
load(file.path(root, "source_cache", "ices_plaice_baltic_21_32_2026", "run__model.RData"), fit_env)
fit <- fit_env$fit
years <- as.integer(fit$data$years)
ages <- seq_len(ncol(fit$data$stockMeanWeight))
aux <- as.data.frame(fit$data$aux)
names(aux) <- c("year", "fleet", "age")

observation_key <- function(year, age, fleet) paste(year, age, fleet, sep = ":")
same_matrix <- function(rows, matrix, years, ages) {
  rows <- rows[order(as.integer(rows$year), as.integer(rows$age)), , drop = FALSE]
  expected <- matrix[cbind(match(as.integer(rows$year), years), match(as.integer(rows$age), ages))]
  !anyDuplicated(observation_key(rows$year, rows$age, rows$fleet)) &&
    isTRUE(all.equal(as.numeric(rows$value), as.numeric(expected), tolerance = 1e-10))
}

stock <- source_record$stock
assessment <- source_record$assessment
catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age", ]
q1 <- inputs[inputs$type == "index" & inputs$measure == "numbers_at_age" & inputs$survey == "Q1IBTS+BITS1-7+", ]
q4 <- inputs[inputs$type == "index" & inputs$measure == "numbers_at_age" & inputs$survey == "Q4IBTS+BITS1-7+", ]
weights <- inputs[inputs$measure == "relative_precision_weight", ]
n_age <- outputs[outputs$type == "population" & outputs$measure == "numbers_at_age", ]
f_age <- outputs[outputs$type == "mortality" & outputs$measure == "fishing_mortality_at_age", ]
total_n <- outputs[outputs$measure == "total_numbers", ]
summary <- outputs[outputs$source_type == "official_table" &
  outputs$measure %in% c("recruitment", "SSB", "Fbar"), ]
predictions <- outputs[outputs$measure %in% c("predicted_catch", "predicted_index"), ]
q <- outputs[outputs$type == "catchability" & outputs$measure == "q", ]
survey_names <- c("Q1IBTS+BITS1-7+", "Q4IBTS+BITS1-7+")
expected_q_keys <- rbind(
  c(0L, 1L, 2L, 2L, 2L, 3L, 4L),
  c(5L, 6L, 7L, 8L, 8L, 8L, 8L)
)
q_group_ages <- function(i) {
  group <- q$age_group[i]
  if (!is.na(group) && grepl("^[0-9]+-[0-9]+$", group)) {
    bounds <- as.integer(strsplit(group, "-", fixed = TRUE)[[1L]])
    return(seq.int(bounds[1L], bounds[2L]))
  }
  if (!is.na(group) && grepl("^[0-9]+\\+$", group)) {
    return(seq.int(as.integer(sub("\\+$", "", group)), max(ages)))
  }
  as.integer(q$age[i])
}
recruitment_age <- source_record$assumptions[
  source_record$assumptions$component == "N" &
    source_record$assumptions$setting == "recruitment_age", ]
plus_group <- source_record$assumptions[
  source_record$assumptions$component == "population" &
    source_record$assumptions$setting == "plus_group", ]

stopifnot(
  nrow(stock) == 1L,
  stock$authority_stock_id == "ple.27.21-32",
  stock$ocean == "Northeast Atlantic",
  stock$region == "Baltic Sea",
  nrow(assessment) == 1L,
  assessment$is_current,
  assessment$assessment_year == 2026L,
  assessment$terminal_year == 2025L,
  assessment$inputs_status == "complete",
  assessment$assumptions_status == "partial",
  assessment$outputs_status == "partial",
  as.integer(fit$conf$minAge) == 1L,
  as.integer(fit$conf$maxAge) == 7L,
  nrow(recruitment_age) == 1L,
  recruitment_age$value == "1",
  nrow(plus_group) == 1L,
  plus_group$value == "7+",
  length(years) == 24L,
  identical(years, 2002:2025),
  identical(ages, 1:7),
  nrow(catch) == 161L,
  nrow(q1) == 168L,
  nrow(q4) == 161L,
  all(catch$unit == "thousand fish"),
  all(q1$basis == "index_scale"),
  all(q4$basis == "index_scale"),
  all(q1$unit == "native survey index (unit unresolved)"),
  all(q4$unit == "native survey index (unit unresolved)"),
  all(q1$sampling_time == 0.1507111),
  all(q4$sampling_time == 0.8173179),
  length(weights$value) == sum(is.finite(fit$data$weight)),
  all(weights$source_type == "native_model"),
  all(weights$basis == "relative_precision"),
  nrow(n_age) == 24L * 7L,
  nrow(f_age) == 24L * 7L,
  all(n_age$age_group[n_age$age == 7] == "7+"),
  all(f_age$age_group[f_age$age == 7] == "7+"),
  nrow(total_n) == 24L,
  max(abs(as.numeric(total_n$value) -
    vapply(2002:2025, function(y) sum(as.numeric(n_age$value[n_age$year == y])), numeric(1)))) < 1e-8,
  nrow(summary) == 24L * 3L,
  as.numeric(summary$value[summary$measure == "recruitment" & summary$year == 2025]) == 8796107,
  as.numeric(summary$value[summary$measure == "SSB" & summary$year == 2025]) == 96773,
  as.numeric(summary$value[summary$measure == "Fbar" & summary$year == 2025]) == 0.02,
  as.numeric(summary$lwr[summary$measure == "SSB" & summary$year == 2025]) == 79727,
  nrow(predictions) == 490L
)

stopifnot(
  nrow(q) == 9L,
  identical(
    unname(expected_q_keys),
    rbind(
      as.integer(fit$conf$keyLogFpar[2L, ]),
      as.integer(fit$conf$keyLogFpar[3L, ])
    )
  ),
  all(q$source_type == "native_model"),
  all(is.na(q$year)),
  all(q$unit == "native survey index units per thousand fish (physical index unit unresolved)"),
  setequal(
    paste(q$survey[!is.na(q$age_group)], q$age_group[!is.na(q$age_group)], sep = ":"),
    c("Q1IBTS+BITS1-7+:3-5", "Q1IBTS+BITS1-7+:7+", "Q4IBTS+BITS1-7+:4-7")
  )
)
for (survey in survey_names) {
  survey_rows <- which(q$survey == survey)
  covered_ages <- sort(unlist(lapply(survey_rows, q_group_ages)))
  stopifnot(identical(as.integer(covered_ages), 1:7L))
}
for (i in seq_len(nrow(q))) {
  survey_row <- match(q$survey[i], survey_names) + 1L
  mapped_ages <- q_group_ages(i)
  parameter_key <- unique(as.integer(fit$conf$keyLogFpar[survey_row, mapped_ages]))
  stopifnot(length(parameter_key) == 1L, parameter_key >= 0L)
  parameter <- parameter_key + 1L
  log_estimate <- fit$pl$logFpar[parameter]
  log_se <- fit$plsd$logFpar[parameter]
  estimate <- exp(log_estimate)
  stopifnot(
    abs(as.numeric(q$value[i]) - estimate) < 1e-12,
    abs(as.numeric(q$se[i]) - estimate * log_se) < 1e-12,
    abs(as.numeric(q$lwr[i]) - exp(log_estimate - qnorm(0.975) * log_se)) < 1e-12,
    abs(as.numeric(q$upr[i]) - exp(log_estimate + qnorm(0.975) * log_se)) < 1e-12
  )
}

observation_rows <- rbind(catch, q1, q4)
expected <- data.frame(
  year = aux$year,
  age = aux$age,
  fleet = aux$fleet,
  value = exp(fit$data$logobs)
)
observed_fleet <- c(
  "Residual catch", "Q1IBTS+BITS1-7+", "Q4IBTS+BITS1-7+"
)
expected$survey <- ifelse(expected$fleet == 1L, "", observed_fleet[expected$fleet])
observation_rows$survey[is.na(observation_rows$survey)] <- ""
expected_key <- observation_key(expected$year, expected$age, expected$survey)
input_key <- observation_key(observation_rows$year, observation_rows$age, observation_rows$survey)
expected <- expected[match(input_key, expected_key), ]
stopifnot(
  !anyNA(expected$value),
  max(abs(as.numeric(observation_rows$value) - expected$value)) < 1e-10
)

nat_m <- inputs[inputs$measure == "natural_mortality_at_age", ]
maturity <- inputs[inputs$measure == "maturity_at_age", ]
stopifnot(
  nrow(nat_m) == length(years) * length(ages),
  nrow(maturity) == length(years) * length(ages),
  all(nat_m$year_basis == "calendar_year"),
  all(nat_m$unit == "per year"),
  all(maturity$unit == "proportion"),
  all(is.finite(as.numeric(nat_m$value))),
  all(as.numeric(maturity$value) >= 0 & as.numeric(maturity$value) <= 1)
)

review <- paste(readLines(file.path(root, "source_reviews", "ices_plaice_baltic_21_32_2026.md"), warn = FALSE), collapse = " ")
stopifnot(
  grepl("fbarRange=3:6", review, fixed = TRUE),
  grepl("14,910,084", review, fixed = TRUE),
  grepl("95% confidence", review, fixed = TRUE)
)

cat("ICES Baltic plaice source records passed.\n")
