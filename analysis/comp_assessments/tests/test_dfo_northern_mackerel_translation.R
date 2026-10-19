source("analysis/comp_assessments/tests/helper_stock.R")

root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "dfo_northern_mackerel_2025"
source_data <- read_assessment(assessment_id, read_database())
original <- source_data
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat

stopifnot(
  identical(translated$years, 1969:2025),
  identical(translated$ages, 1:10),
  translated$age_plus_group == 10,
  tinyAM::check_obs(obs),
  nrow(obs$catch) == 570L,
  nrow(obs$index) == 0L,
  nrow(obs$weight) == 570L,
  nrow(obs$maturity) == 570L,
  identical(attr(obs, "translation")$selected_surveys, character()),
  "Total egg production (TEP)" %in% attr(obs, "translation")$excluded_surveys,
  all(obs$pop_covariates$M_assumption == 0.28),
  dat$F_settings$process == "cor_rw",
  identical(dat$F_settings$mean_ages, 5:10),
  translated$settings$rec_settings$process == "iid",
  grepl("bh(ssb)", paste(deparse(translated$settings$rec_settings$mu_form), collapse = ""), fixed = TRUE),
  dat$M_settings$process == "off",
  ncol(dat$q_modmat) == 0L,
  ncol(dat$sd_index_modmat) == 0L,
  identical(translated$settings$N_settings$init, "exp"),
  identical(paste(deparse(translated$settings$index_settings$q_form), collapse = ""), "~0"),
  identical(paste(deparse(translated$settings$index_settings$sd_form), collapse = ""), "~0"),
  setequal(unique(translated$comparison_outputs$year), 1969:2025),
  all(c("numbers_at_age", "fishing_mortality_at_age", "SSB", "recruitment") %in%
    translated$comparison_outputs$measure),
  identical(translated$comparison_aggregates, "ssb"),
  identical(translated$comparison_definitions$M$status, "matched"),
  identical(translated$comparison_definitions$ssb$status, "approximate"),
  any(grepl("catch-only", translated$background, fixed = TRUE)),
  is.null(translated$fit),
  is.null(translated$start_par),
  identical(source_data, original)
)

cat("Northern mackerel translation preparation and comparison settings passed.\n")
