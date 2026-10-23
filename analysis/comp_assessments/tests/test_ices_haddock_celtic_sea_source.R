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

stopifnot(
  nrow(source_record$stock) == 1L,
  source_record$stock$authority_stock_id == "had.27.7b-k",
  source_record$stock$common_name == "Celtic Sea haddock",
  source_record$stock$ocean == "Northeast Atlantic",
  nrow(source_record$assessment) == 1L,
  source_record$assessment$is_current,
  source_record$assessment$terminal_year == 2025L,
  source_record$assessment$inputs_status == "complete",
  source_record$assessment$outputs_status == "partial",
  identical(as.integer(fit$data$years), years),
  nrow(catch) == length(years) * length(ages),
  nrow(index) == 23L * 8L,
  sum(is.na(index$value)) == 1L,
  any(index$year == 2005L & index$age == 7L & is.na(index$value)),
  nrow(n_at_age) == length(years) * length(ages),
  nrow(f_at_age) == length(years) * length(ages),
  nrow(total_numbers) == length(years),
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
