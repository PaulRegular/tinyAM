source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))
database <- read_database()
for (id in c(
  "ices_haddock_north_sea_2026", "ices_haddock_iceland_2025",
  "ices_norway_pout_north_sea_2026_benchmark",
  "nefsc_summer_flounder_2018"
)) {
  source_data <- read_assessment(id, database)
  stock <- new.env(parent = globalenv())
  translated <- .test_stock(source_data)
  dat <- do.call(tinyAM::prepare_tam, c(
    list(data = translated$obs, years = translated$years, ages = translated$ages),
    translated$settings
  ))
  par <- tinyAM::make_par(dat)
  stationary <- id == "nefsc_summer_flounder_2018"
  stopifnot(
    dat$F_settings$process == if (stationary) "ar1" else "cor_rw",
    if (stationary) length(par$log_mu_f) == length(translated$ages) else is.null(par$log_mu_f),
    if (stationary) length(par$logit_phi_f) == 2L else length(par$atanh_rho_f) == 1L,
    is.finite(tinyAM::nll_fun(par, dat))
  )
}
translated <- .test_stock(read_assessment("ices_saithe_north_sea_2026", database))
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = translated$obs, years = translated$years, ages = translated$ages), translated$settings
))
par <- tinyAM::make_par(dat)
stopifnot(
  dat$F_settings$process == "cor_rw", is.null(par$log_mu_f),
  is.finite(tinyAM::nll_fun(par, dat))
)
cat("Retained F-process parameterization tests passed.\n")
