root <- "analysis/comp_assessments"
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "afsc_pollock_ebs_2024.R"), envir = stock)
database <- read_committed_database()
source_data <- read_committed_assessment("afsc_pollock_ebs_2024", database)
translated <- stock$translate_stock(source_data)
obs <- translated$obs

tinyAM::check_obs(obs)
dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  translated$settings
))
stopifnot(identical(dat$years, 1964:2024),
          identical(dat$ages, 1:15),
          dat$N_settings$init == "exp",
          all(obs$catch$age %in% dat$ages),
          all(obs$index$age %in% dat$ages))
cat("EBS pollock translation tests passed.\n")
