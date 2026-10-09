root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))
database <- read_database()
for (id in c("ices_haddock_north_sea_2026", "ices_haddock_iceland_2025",
             "ices_norway_pout_north_sea_2026_benchmark",
             "ices_saithe_north_sea_2026", "nefsc_summer_flounder_2018")) {
  source_data <- read_assessment(id, database)
  stock <- new.env(parent = globalenv())
  sys.source(file.path(root, "scripts", "translation", "stocks",
                       paste0(id, ".R")), stock)
  translated <- stock$translate_stock(source_data)
  dat <- do.call(tinyAM::prepare_tam, c(
    list(data = translated$obs, years = translated$years, ages = translated$ages),
    translated$settings
  ))
  par <- tinyAM::make_par(dat)
  stopifnot(dat$F_settings$process == "ar1",
            length(par$log_mu_f) == length(translated$ages),
            is.finite(tinyAM::nll_fun(par, dat)))
}
cat("Stationary F translation mean tests passed.\n")
