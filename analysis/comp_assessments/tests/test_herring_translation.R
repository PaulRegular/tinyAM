root <- "analysis/comp_assessments"
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "ices_herring_north_sea_2026.R"), envir = stock)
database <- read_committed_database()
source_data <- read_committed_assessment("ices_herring_north_sea_2026", database)
translated <- stock$translate_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))

stopifnot(
  identical(translated$years, 1947:2025),
  identical(translated$ages, 0:8),
  nrow(obs$catch) == 711L,
  nrow(obs$index) == 532L,
  nrow(obs$weight) == 711L,
  nrow(obs$maturity) == 711L,
  setequal(unique(obs$index$survey), c("HERAS", "IBTS-Q1", "IBTS0", "IBTS-Q3")),
  all(obs$index$age %in% 0:8),
  !any(obs$index$survey %in% c("LAI-SNS", "LAI-CNS", "LAI-BUN", "LAI-ORSH")),
  all(obs$weight$M_assumption > 0),
  all(obs$index$samp_time == c(HERAS = 0.5, `IBTS-Q1` = 0.125,
                              IBTS0 = 0.125, `IBTS-Q3` = 0.625)[obs$index$survey]),
  setequal(levels(obs$index$q_key), c("HERAS_1_2", "HERAS_3_8", "IBTS-Q1",
                                      "IBTS-Q3_0", "IBTS-Q3_1", "IBTS-Q3_2",
                                      "IBTS-Q3_3", "IBTS-Q3_4", "IBTS-Q3_5",
                                      "IBTS0")),
  dat$N_settings$init == "free",
  dat$F_settings$process == "rw",
  dat$M_settings$process == "off",
  max(dat$years) == 2025,
  translated$comparison_scales[["F_bar"]] == 1,
  all(is.finite(translated$start_par$log_f)),
  all(is.finite(translated$start_par$log_r)),
  all(is.finite(translated$start_par$log_n)),
  all(is.finite(translated$start_par$log_n0))
)

m_source <- source_data$inputs[
  source_data$inputs$type == "M" &
    source_data$inputs$measure == "natural_mortality_at_age", , drop = FALSE
]
m_key <- paste(obs$weight$year, obs$weight$age)
source_key <- paste(m_source$year, m_source$age)
expected_m <- as.numeric(m_source$value[match(m_key, source_key)]) + 0.02
stopifnot(isTRUE(all.equal(obs$weight$M_assumption, expected_m)))

cat("North Sea herring translation tests passed.\n")
