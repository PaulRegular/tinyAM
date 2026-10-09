root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_herring_norwegian_spring_2025"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
recipe <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     paste0(assessment_id, ".R")), envir = recipe)
translated <- recipe$translate_stock(source_data)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

inputs <- source_data$inputs
maturity <- inputs[inputs$type == "maturity" &
                     inputs$measure == "maturity_at_age", ]
cohort_value <- maturity$value[
  as.integer(maturity$year) == 1986 & as.integer(maturity$age) == 2
]
translated_value <- translated$obs$maturity$obs[
  translated$obs$maturity$year == 1988 & translated$obs$maturity$age == 2
]
catch <- inputs[inputs$assessment_id == assessment_id & inputs$type == "catch" &
                  inputs$measure == "numbers_at_age" & inputs$year == "1988" &
                  as.integer(inputs$age) >= 12, ]
plus_catch <- translated$obs$catch$obs[
  translated$obs$catch$year == 1988 & translated$obs$catch$age == 12
]
n_source <- source_data$outputs[source_data$outputs$type == "population" &
                                  source_data$outputs$measure == "numbers_at_age" &
                                  source_data$outputs$year == "1988" &
                                  source_data$outputs$age == "2", ]

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1988:2024),
  identical(translated$ages, 2:12),
  nrow(translated$obs$catch) == 37L * 11L,
  nrow(translated$obs$weight) == 37L * 11L,
  nrow(translated$obs$maturity) == 37L * 11L,
  all(translated$obs$weight$M_assumption[
    translated$obs$weight$age == 2
  ] == 0.9),
  all(translated$obs$weight$M_assumption[
    translated$obs$weight$age > 2
  ] == 0.15),
  identical(as.numeric(translated_value), as.numeric(cohort_value)),
  identical(as.numeric(plus_catch), sum(as.numeric(catch$value) * 1000)),
  !"RFID" %in% unique(translated$obs$index$survey),
  !any(translated$obs$index$survey == "IESNS_Barents" &
         translated$obs$index$year == 2008),
  all(is.finite(par$log_f)),
  all(is.finite(translated$start_par$log_q)),
  identical(dim(translated$start_par$log_f), c(37L, 11L)),
  isTRUE(all.equal(translated$start_par$log_r0,
                   log(as.numeric(n_source$value) * 1e6))),
  identical(dat$N_settings$process, "iid"),
  length(unique(translated$obs$index$survey)) == 5L
)

cat("Norwegian spring-spawning herring translation structure passed.\n")
