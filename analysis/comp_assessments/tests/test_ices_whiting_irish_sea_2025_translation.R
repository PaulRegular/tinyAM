source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/read_database.R")
source("analysis/comp_assessments/R/database_to_tam_obs.R")

assessment_id <- "ices_whiting_irish_sea_2025_benchmark"
source_record <- read_assessment(assessment_id)
translated <- .test_stock(source_record)
obs <- translated$obs

stopifnot(
  identical(translated$years, 1980:2023),
  identical(translated$ages, 0:6),
  translated$age_plus_group == 6L,
  nrow(obs$catch) == 308L,
  nrow(obs$index) == 444L,
  setequal(unique(obs$index$survey), c("NIGFS-Q1", "NIGFS-Q4", "NIMIK")),
  sum(obs$index$survey == "NIGFS-Q1") == 192L,
  sum(obs$index$survey == "NIGFS-Q4") == 224L,
  sum(obs$index$survey == "NIMIK") == 28L,
  setequal(unique(obs$index$q_key[obs$index$survey == "NIGFS-Q1"]),
    interaction("NIGFS-Q1", c("age1", "ages2_5", "age6"), drop = TRUE, lex.order = TRUE)),
  length(unique(obs$index$q_key)) == 9L,
  sum(is.finite(obs$index$relative_precision_weight)) == 416L,
  all(is.finite(obs$index$relative_sd)),
  all(obs$index$relative_sd[!is.finite(obs$index$relative_precision_weight)] == 1),
  !2020L %in% obs$index$year[obs$index$survey == "NIMIK"],
  !2022L %in% obs$index$year[obs$index$survey == "NIMIK"],
  identical(translated$settings$N_settings$process, "iid"),
  identical(translated$settings$N_settings$init, "exp"),
  identical(translated$settings$rec_settings$process, "rw"),
  identical(translated$settings$F_settings$process, "cor_rw"),
  identical(translated$settings$M_settings$process, "off"),
  identical(deparse(translated$settings$catch_settings$sd_form), "~1"),
  length(translated$warm_call) > 0L,
  !"start_par" %in% names(translated$warm_call),
  "start_par" %in% names(translated$fit_call),
  any(grepl("AR-correlated across ages", translated$background, fixed = TRUE)),
  any(grepl("independent index residuals", translated$background, fixed = TRUE))
)

cat("Irish Sea whiting tinyAM preparation and source-to-model mappings passed validation.\n")
