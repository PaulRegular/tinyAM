root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "read_database.R"))
assessment_id <- "jrn_afwg_haddock_northeast_arctic_2026"
source_data <- read_assessment(assessment_id, read_database())
inputs <- source_data$inputs
outputs <- source_data$outputs

count <- function(type, measure) {
  sum(inputs$type == type & inputs$measure == measure)
}
one_output <- function(measure, year, age = NA_integer_) {
  x <- outputs[outputs$measure == measure & outputs$year == year, , drop = FALSE]
  if (!is.na(age)) x <- x[x$age == age, , drop = FALSE]
  stopifnot(nrow(x) == 1L)
  x$value[[1L]]
}

stopifnot(
  nrow(source_data$stocks) == 1L,
  source_data$stocks$stock_id == "ck_ices_afwg_aeglefinus_nea_1_2",
  source_data$stocks$authority == "JRN-AFWG",
  source_data$stocks$authority_stock_id == "had.27.1-2",
  source_data$assessment$assessment_year == 2026L,
  source_data$assessment$terminal_year == 2026L,
  source_data$assessment$estimate_terminal_year == 2026L,
  is.na(source_data$assessment$framework_year),
  source_data$stocks$ocean == "Northeast Atlantic",
  count("catch", "numbers_at_age") == 836L,
  count("catch", "numbers_at_age_for_M_reconstruction") == 836L,
  count("index", "numbers_at_age") == 678L,
  count("index", "index_by_age_group") == 112L,
  count("weight", "weight_at_age") == 847L,
  count("maturity", "maturity_at_age") == 847L,
  count("M", "natural_mortality_at_age") == 847L,
  count("biology", "fraction_F_before_spawning") == 847L,
  count("biology", "fraction_M_before_spawning") == 847L,
  length(unique(inputs$survey[inputs$type == "index"])) == 4L,
  one_output("numbers_at_age", 2025L, 3L) == 494847,
  one_output("recruitment", 2025L) == 494847,
  one_output("SSB", 2025L) == 157942,
  one_output("Fbar", 2025L) == 0.381,
  one_output("total_biomass", 2025L) == 494692,
  one_output("natural_mortality_at_age", 2025L, 7L) == 0.2,
  one_output("natural_mortality_at_age", 2025L, 13L) == 0.2,
  all(is.na(outputs$se)),
  all(is.na(outputs$lwr[outputs$measure %in% c("numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age")]))
)

cat("JRN-AFWG Northeast Arctic haddock source records passed.\n")
