pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
record <- read_assessment("ices_cod_west_scotland_2022", read_database())
before <- record
stock <- .assessment_stock(record, do_fit = FALSE, silent = TRUE)
stopifnot(is.null(stock$error), is.null(stock$env$fit), identical(record, before))
env <- stock$env
obs <- env$obs
dat <- tinyAM::prepare_tam(
  data = obs, years = env$years, ages = env$ages,
  N_settings = list(process = "iid", init = "free", sd_form = ~1),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "cor_rw", mu_form = NULL,
    sd_form = ~I(age == 1), mean_ages = 2:5),
  M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
  catch_settings = list(sd_form = ~0 + sd_block, sd_supplied = ~relative_sd),
  index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey,
    sd_supplied = ~relative_sd)
)
par <- tinyAM::make_par(dat)
stopifnot(
  identical(env$years, 1981:2021), identical(env$ages, 1:7),
  identical(dim(par$log_f), c(41L, 7L)),
  identical(dim(par$log_n), c(40L, 6L)), length(par$log_n0) == 6,
  ncol(dat$q_modmat) == 21, ncol(dat$sd_catch_modmat) == 3,
  ncol(dat$sd_index_modmat) == 5,
  all(obs$catch$relative_sd[obs$catch$year < 2006] == 1),
  all(obs$catch$relative_sd[obs$catch$year >= 2006] == 2),
  all(obs$maturity$spawn_time == 0),
  all(is.na(obs$catch$obs[obs$catch$year == 2021 & obs$catch$age %in% 1:2])),
  abs(obs$index$obs[obs$index$survey == "SCO.Q4" & obs$index$year == 2021 & obs$index$age == 2] - 0.008912) < 1e-12,
  env$comparison_definitions$F$status == "approximate",
  env$comparison_definitions$ssb$status == "matched"
)
for (survey in unique(obs$index$survey)) {
  x <- obs$index[obs$index$survey == survey, ]
  source_q <- record$outputs[record$outputs$measure == "q" & record$outputs$survey == survey, ]
  q <- source_q$value[match(x$age, source_q$age)]
  stopifnot(all(is.finite(q)), all(vapply(split(q, x$q_key, drop = TRUE),
    function(values) diff(range(values)) < 1e-12, logical(1))))
}
native_index <- record$inputs[record$inputs$type == "index" & record$inputs$measure == "numbers_at_age", ]
key <- function(x) paste(x$year, x$age, x$survey, sep = "/")
at <- match(key(obs$index), key(native_index))
stopifnot(!anyNA(at), isTRUE(all.equal(obs$index$obs, native_index$value[at] * 0.001)))
native_catch <- record$inputs[record$inputs$type == "catch" & record$inputs$measure == "numbers_at_age", ]
at <- match(paste(obs$catch$year, obs$catch$age), paste(native_catch$year, native_catch$age))
stopifnot(isTRUE(all.equal(obs$catch$obs, native_catch$value[at])))
q <- record$outputs[record$outputs$measure == "q", ]
scaled_q <- env$comparison_outputs[env$comparison_outputs$measure == "q", ]
for (field in c("value", "se", "lwr", "upr")) {
  stopifnot(isTRUE(all.equal(scaled_q[[field]], q[[field]] * 0.001)))
}
cat("West of Scotland cod translation preparation, q sharing and precision passed.\n")
