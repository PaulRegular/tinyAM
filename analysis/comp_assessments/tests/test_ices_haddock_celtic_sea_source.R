source("analysis/comp_assessments/tests/helper_stock.R")
root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_haddock_celtic_sea_7b_k_2026"
database <- read_database()
source_record <- read_assessment(assessment_id, database)
inputs <- source_record$inputs
outputs <- source_record$outputs
years <- 1993:2025
ages <- 0:8

fit_env <- new.env(parent = baseenv())
load(file.path(root, "source_cache", "ices_haddock_7b_k_2026", "run", "model.RData"),
  fit_env)
fit <- fit_env$fit

catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age", , drop = FALSE]
index <- inputs[inputs$type == "index" & inputs$measure == "numbers_at_age", , drop = FALSE]
n_at_age <- outputs[outputs$measure == "numbers_at_age", , drop = FALSE]
f_at_age <- outputs[outputs$measure == "fishing_mortality_at_age", , drop = FALSE]
total_numbers <- outputs[outputs$measure == "total_numbers", , drop = FALSE]
predicted_catch <- outputs[outputs$measure == "predicted_catch", , drop = FALSE]
predicted_index <- outputs[outputs$measure == "predicted_index", , drop = FALSE]
q <- outputs[outputs$type == "catchability" & outputs$measure == "q", , drop = FALSE]
output_key <- do.call(paste, c(lapply(
  outputs[c("type", "measure", "fleet", "survey", "sex", "region", "season", "year", "age", "age_group")],
  function(x) ifelse(is.na(x), "", as.character(x))
), sep = "\r"))

fit_aux <- as.data.frame(fit$data$aux, stringsAsFactors = FALSE)
names(fit_aux) <- c("year", "fleet", "age")
fit_aux$year <- as.integer(fit_aux$year)
fit_aux$fleet <- as.integer(fit_aux$fleet)
fit_aux$age <- as.integer(fit_aux$age)
catch_rows <- which(fit_aux$fleet == 1L)
index_rows <- which(fit_aux$fleet == 2L)
expected_catch <- exp(fit$rep$predObs[catch_rows])
expected_index <- exp(fit$rep$predObs[index_rows])
catch_match <- match(
  paste(predicted_catch$year, predicted_catch$age),
  paste(fit_aux$year[catch_rows], fit_aux$age[catch_rows])
)
index_match <- match(
  paste(predicted_index$year, predicted_index$age),
  paste(fit_aux$year[index_rows], fit_aux$age[index_rows])
)
q <- q[order(as.integer(q$age)), , drop = FALSE]
q_expected <- exp(fit$pl$logFpar)
q_log_se <- fit$plsd$logFpar

stopifnot(
  nrow(source_record$stock) == 1L,
  source_record$stock$authority_stock_id == "had.27.7b-k",
  source_record$stock$common_name == "Celtic Sea haddock",
  source_record$stock$ocean == "Northeast Atlantic",
  nrow(source_record$assessment) == 1L,
  source_record$assessment$is_current,
  source_record$assessment$terminal_year == 2025L,
  source_record$assessment$inputs_status == "complete",
  source_record$assessment$outputs_status == "complete",
  identical(as.integer(fit$data$years), years),
  nrow(catch) == length(years) * length(ages),
  nrow(index) == 23L * 8L,
  sum(is.na(index$value)) == 1L,
  any(index$year == 2005L & index$age == 7L & is.na(index$value)),
  nrow(n_at_age) == length(years) * length(ages),
  nrow(f_at_age) == length(years) * length(ages),
  nrow(total_numbers) == length(years),
  nrow(predicted_catch) == length(catch_rows),
  nrow(predicted_index) == length(index_rows),
  nrow(q) == 8L,
  all(predicted_catch$source_type == "native_model"),
  all(predicted_index$source_type == "native_model"),
  all(q$source_type == "native_model"),
  all(predicted_catch$fleet == "Residual catch"),
  all(predicted_index$survey == "FRA_IRL_IBTS"),
  identical(as.integer(q$age), 0:7),
  all(q$survey == "FRA_IRL_IBTS"),
  all(predicted_catch$age_group[predicted_catch$age == 8L] == "8+"),
  !any(predicted_catch$age == 0L),
  sum(is.na(fit$data$logobs)) == 3L,
  all(is.finite(fit$rep$predObs)),
  all(is.finite(catch_match)),
  all(is.finite(index_match)),
  isTRUE(all.equal(
    as.numeric(predicted_catch$value),
    as.numeric(expected_catch[catch_match]),
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    as.numeric(predicted_index$value),
    as.numeric(expected_index[index_match]),
    tolerance = 1e-10
  )),
  isTRUE(all.equal(as.numeric(q$value), as.numeric(q_expected), tolerance = 1e-10)),
  isTRUE(all.equal(
    as.numeric(q$se), as.numeric(q_expected * q_log_se), tolerance = 1e-10
  )),
  isTRUE(all.equal(
    as.numeric(q$lwr),
    as.numeric(exp(fit$pl$logFpar - stats::qnorm(0.975) * q_log_se)),
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    as.numeric(q$upr),
    as.numeric(exp(fit$pl$logFpar + stats::qnorm(0.975) * q_log_se)),
    tolerance = 1e-10
  )),
  all(outputs$type %in% c(
    "population", "mortality", "biomass", "catch", "index", "catchability"
  )),
  all(outputs$source_type %in% c("native_model", "derived_source_output")),
  !anyDuplicated(output_key),
  isTRUE(all.equal(as.numeric(total_numbers$value),
    as.numeric(colSums(exp(fit$pl$logN))), tolerance = 1e-10)),
  isTRUE(all.equal(as.numeric(n_at_age$value),
    as.numeric(exp(fit$pl$logN)), tolerance = 1e-10))
)

obs <- database_to_tam_obs(
  assessment_id,
  inputs,
  years = years,
  ages = ages,
  surveys = "FRA_IRL_IBTS",
  assumptions = source_record$assumptions
)

stopifnot(
  nrow(obs$catch) == length(years) * length(ages),
  nrow(obs$index) == 23L * 8L - 1L,
  !any(obs$index$year == 2005L & obs$index$age == 7L),
  all(is.finite(obs$index$obs))
)

cat("ICES Celtic Sea haddock source records passed.\n")
