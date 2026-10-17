source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/run_assessment.R")
pkgload::load_all(".", quiet = TRUE)

id <- "ices_herring_western_baltic_2026"
source_data <- read_assessment(id, read_database())
before <- source_data
translated <- .test_stock(source_data)
obs <- translated$obs
dat <- translated$dat
input <- source_data$inputs

catch <- subset(input, type == "catch" & measure == "numbers_at_age")
catch$age <- pmin(catch$age, 5)
expected <- aggregate(value ~ year + age, catch, sum)
key <- function(d) paste(d$year, d$age)
stopifnot(
  tinyAM::check_tam_data(obs),
  identical(source_data, before),
  identical(translated$years, 1991:2025),
  identical(translated$ages, 0:5),
  translated$age_plus_group == 5,
  nrow(obs$catch) == 35L * 6L,
  nrow(obs$weight) == 35L * 6L,
  isTRUE(all.equal(obs$catch$obs, expected$value[match(key(obs$catch), key(expected))])),
  setequal(obs$index$survey, c("HERAS", "GERAS", "N20", "IBTS/BITSQ1")),
  max(obs$index$age) == 5,
  length(unique(obs$index$sd_block)) == 4L,
  length(unique(obs$catch$sd_block)) == 3L,
  !any(obs$index$survey == "HERAS" & obs$index$year == 1999),
  all(obs$index$relative_sd[obs$index$survey != "IBTS/BITSQ1"] == 1),
  identical(dat$N_settings$process, "iid"),
  identical(dat$N_settings$init, "free"),
  identical(dat$F_settings$process, "cor_rw"),
  identical(dat$F_settings$mean_ages, 2:5),
  identical(dat$M_settings$process, "off"),
  is.finite(tinyAM::nll_fun(tinyAM::make_par(dat), dat))
)

timing <- c(HERAS = .625, GERAS = .8, N20 = .4, `IBTS/BITSQ1` = .136365)
stopifnot(all(obs$index$samp_time == unname(timing[obs$index$survey])))
precision <- subset(input, measure == "relative_precision_weight")
ibts <- subset(obs$index, survey == "IBTS/BITSQ1")
stopifnot(isTRUE(all.equal(
  ibts$relative_sd, 1 / sqrt(precision$value[match(key(ibts), key(precision))])
)))

heras <- subset(input, type == "index" & measure == "numbers_at_age" & survey == "HERAS" & age >= 5)
heras$value <- heras$value * vapply(heras$unit, .translation_index_multiplier, numeric(1))
expected <- aggregate(value ~ year, heras, sum)
plus <- subset(obs$index, survey == "HERAS" & age == 5)
stopifnot(isTRUE(all.equal(plus$obs, expected$value[match(plus$year, expected$year)])))
cat("Western Baltic herring translation structure passed.\n")
