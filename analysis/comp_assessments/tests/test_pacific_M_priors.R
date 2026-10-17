source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/tests/helper_stock.R")
database <- read_database()
for (id in "joint_hake_coastal_2025") {
  record <- read_assessment(id, database)
  original <- record
  spec <- .test_stock(record)
  baseline_script <- tempfile(fileext = ".R")
  status <- system2("git", c("show", shQuote(paste0(
    "86e44bf:analysis/comp_assessments/scripts/translation/stocks/", id, ".R"
  ))), stdout = baseline_script)
  stopifnot(status == 0)
  baseline <- .test_stock(record, script = baseline_script)
  unlink(baseline_script)
  for (table in c("catch", "index", "weight", "maturity")) {
    stopifnot(identical(spec$obs[[table]], baseline$obs[[table]]))
  }
  stopifnot(attr(spec$obs, "translation")$M$status != "fixed_numerical_input")
  dat <- spec$dat
  par <- tinyAM::make_par(dat)
  term <- dat$M_terms[[1L]]
  stopifnot(identical(record, original), is.null(spec$fit),
    dat$M_settings$process == "off", ncol(dat$M_modmat) == 0L,
    length(dat$M_terms) == 1L, length(par[[term$parameter]]) == 1L,
    is.null(term$sd_parameter), length(unique(dat$log_mu_supplied_m)) == 1L)
  for (effect in c(-.07, 0, .03)) {
    par[[term$parameter]][] <- effect
    m <- tinyAM:::.mean_effects(par, dat, "M")
    stopifnot(all(m$contribution == effect),
      abs(m$nll + stats::dnorm(effect, 0, term$sd, log = TRUE)) < 1e-12)
  }
  center <- exp(unique(dat$log_mu_supplied_m))
  stopifnot(abs(log(center) + 1.60944) < 1e-12, term$sd == .1,
    !"mu_m" %in% names(spec$retry_start()))
}
cat("One shared M level receives one fixed-SD Gaussian penalty; source prior scales and input preservation passed.\n")

# A fixed change of index units centers the native q penalty without changing
# any log-index likelihood residual or using fitted values to adjust data.
record <- read_assessment("afsc_northern_rock_sole_ebs_2024", database)
spec <- .test_stock(record)
obs <- spec$obs
original_index <- obs$index$obs
obs$index$obs <- original_index / 1.5
obs$index$q_age <- pmin(obs$index$age, 17)
obs$index$q_prior_group <- factor("stock")
settings <- spec$settings
settings$index_settings$q_form <- ~0+logistic(q_age)+iid(q_prior_group,sd=.2)
dat <- do.call(tinyAM::prepare_tam,c(list(data=obs,years=spec$years,ages=spec$ages),settings))
par <- tinyAM::make_par(dat)
term <- Filter(function(x)x$type=="iid",dat$q_terms)[[1L]]
stopifnot(length(par[[term$parameter]])==1L,length(par$log_q)==0L)
par[[term$parameter]][] <- .03
effect <- tinyAM:::.q_effects(par,dat)
stopifnot(abs(effect$nll+dnorm(.03,0,.2,log=TRUE))<1e-12,
  identical(spec$obs$index$obs,original_index))
keep <- is.finite(original_index) & original_index>0
log_q <- effect$contribution+effect$log_selectivity
native <- dnorm(log(original_index[keep]),log(1e6)+log(1.5)+log_q[keep],.4,log=TRUE)
standardized <- dnorm(log(obs$index$obs[keep]),log(1e6)+log_q[keep],.4,log=TRUE)
stopifnot(max(abs(native-standardized))<1e-10)
cat("Source q penalty is counted once; fixed index-unit conversion preserves every log-likelihood residual.\n")
