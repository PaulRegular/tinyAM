root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_whiting_north_sea_2026"
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

catch_1978 <- source_data$inputs[
  source_data$inputs$type == "catch" &
    source_data$inputs$measure == "numbers_at_age" &
    source_data$inputs$year == 1978,
]
m_surface <- translated$obs$weight
m_2022 <- m_surface[m_surface$year == 2022, c("age", "M_assumption")]

stopifnot(
  tinyAM::check_obs(translated$obs),
  identical(translated$years, 1978:2026),
  identical(translated$ages, 0:8),
  nrow(translated$obs$catch) == 49L * 9L,
  nrow(translated$obs$index) == 509L,
  all(translated$obs$index$samp_time[
    translated$obs$index$survey == "IBTS-Q1"
  ] == .125),
  all(translated$obs$index$samp_time[
    translated$obs$index$survey == "IBTS-Q3"
  ] == .625),
  length(unique(translated$obs$index$q_key)) == 13L,
  translated$obs$catch$obs[
    translated$obs$catch$year == 1978 & translated$obs$catch$age == 0
  ] == catch_1978$value[catch_1978$age == 0] * 1000,
  all(is.na(translated$obs$catch$obs[translated$obs$catch$year == 2026])),
  all(vapply(2023:2026, function(year) {
    m_year <- m_surface[m_surface$year == year, c("age", "M_assumption")]
    identical(m_year$M_assumption[match(m_2022$age, m_year$age)],
              m_2022$M_assumption)
  }, logical(1))),
  identical(dat$N_settings$process, "rw"),
  identical(dat$N_settings$init, "exp"),
  identical(dat$F_settings$process, "rw"),
  identical(dat$F_settings$mean_ages, 2:5),
  identical(dat$M_settings$process, "off"),
  identical(dim(translated$start_par$log_f), c(49L, 9L)),
  identical(dim(translated$start_par$log_f), dim(par$log_f)),
  identical(dim(translated$start_par$log_n), dim(par$log_n)),
  !any(c("log_n0", "log_m") %in% names(translated$start_par)),
  identical(names(translated$start_par$log_q), names(par$log_q)),
  all(is.finite(translated$start_par$log_f)),
  all(is.finite(translated$start_par$log_n)),
  all(is.finite(translated$start_par$log_q)),
  grepl("carried forward", paste(translated$background, collapse = " "))
)

cat("North Sea whiting translation structure passed.\n")
