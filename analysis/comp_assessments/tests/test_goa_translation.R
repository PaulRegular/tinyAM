root <- "analysis/comp_assessments"
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_committed_assessment.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "afsc_pollock_goa_2024.R"), envir = stock)
database <- read_committed_database()
source_data <- read_committed_assessment("afsc_pollock_goa_2024", database)
translated <- stock$translate_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
stopifnot(all(is.na(obs$catch$obs[obs$catch$age == 2])))
shelikof <- obs$index$survey == "Shelikof winter acoustic"
stopifnot(!any(shelikof & obs$index$age == 3),
          any(is.finite(obs$index$obs[shelikof & obs$index$age >= 4])))

dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))
stopifnot(identical(dat$years, translated$years),
          identical(dat$ages, translated$ages))
cat("GOA pollock translation tests passed.\n")
