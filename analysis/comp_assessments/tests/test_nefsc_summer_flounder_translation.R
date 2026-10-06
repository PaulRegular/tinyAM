root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

assessment_id <- "nefsc_summer_flounder_2018"
database <- read_database()
source_data <- read_assessment(assessment_id, database)
recipe <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     paste0(assessment_id, ".R")), envir = recipe)
translated <- recipe$translate_stock(source_data)
dat <- do.call(tinyAM::make_dat, c(
  list(obs = translated$obs, years = translated$years, ages = translated$ages),
  translated$settings
))
par <- tinyAM::make_par(dat)

published_catch <- source_data$inputs[
  source_data$inputs$type == "catch" &
    source_data$inputs$measure == "numbers_at_age" &
    source_data$inputs$fleet == "Published total" &
    source_data$inputs$year %in% translated$years &
    source_data$inputs$age %in% translated$ages,
  c("year", "age", "value")
]
expected_catch <- published_catch$value[match(
  paste(translated$obs$catch$year, translated$obs$catch$age),
  paste(published_catch$year, published_catch$age)
)] * 1000

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1982:2016),
  identical(translated$ages, 0:7),
  nrow(translated$obs$catch) == 35L * 8L,
  nrow(translated$obs$weight) == 35L * 8L,
  nrow(translated$obs$maturity) == 35L * 8L,
  nrow(translated$obs$index) == 16L * 8L,
  isTRUE(all.equal(translated$obs$catch$obs, expected_catch)),
  setequal(unique(translated$obs$index$survey),
           c("NEFSC Bigelow spring trawl", "NEFSC Bigelow fall trawl")),
  setequal(unique(translated$obs$index$samp_time), c(.25, .75)),
  identical(dat$N_settings$process, "iid"),
  identical(dat$F_settings$process, "ar1"),
  identical(dat$F_settings$mean_ages, 4),
  identical(dat$M_settings$process, "off"),
  identical(dim(par$log_n), c(34L, 7L)),
  all(is.finite(translated$obs$weight$obs)),
  all(is.finite(translated$obs$weight$M_assumption)),
  grepl("published total numbers-at-age once", paste(translated$background, collapse = " "))
)

cat("Summer flounder translation structure passed.\n")
