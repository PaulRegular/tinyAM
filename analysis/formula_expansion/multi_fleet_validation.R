# Run from the repository root: Rscript analysis/formula_expansion/multi_fleet_validation.R core
# The artificial catch split below exercises workflows; it is not an assessment
# trial or a recovery experiment. Scientific trials use separately recorded truth.
pkgload::load_all(".", quiet = TRUE)

fleet_workflow_inputs <- function() {
  d <- tinyAM::cod_obs
  d$catch <- rbind(transform(d$catch, fleet = "Trawl", obs = .6 * obs),
    transform(d$catch, fleet = "Gillnet", obs = .4 * obs))
  d
}

run_fleet_core_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/core") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  warnings <- character(); stage <- "fit"
  results <- withCallingHandlers(tryCatch({
    fit <- tinyAM::fit_tam(fleet_workflow_inputs(), years = 1983:2000, ages = 2:8,
      N_settings = list(process = "off"),
      rec_settings = list(process = "off", mu_form = ~ factor(year)),
      F_settings = list(process = "off", mu_form = ~ fleet + factor(year) + factor(age)),
      M_settings = list(process = "off", mu_supplied = ~ I(.3)), silent = TRUE)
    diagnostics <- tinyAM::check_tam(fit)
    stage <- "projection"
    projected <- update(fit, start_par = fit,
      proj_settings = list(n_proj = 2, n_mean = 2, catch = 1000), silent = TRUE)
    stopifnot(max(abs(projected$rep$total_yield_pred[projected$dat$is_proj] - 1000)) < 1e-6)
    stage <- "simulation"
    set.seed(6515)
    simulations <- tinyAM::sim_tam(fit, n = 2, par_uncertainty = "none")
    stage <- "hindcast"
    hindcast <- tinyAM::fit_hindcast(fit, folds = 1, progress = FALSE)
    stage <- "dashboard"
    tinyAM::vis_tam(fit, open_file = FALSE, output_file = file.path(output, "dashboard.html"),
      render_args = list(quiet = TRUE))
    list(fit = fit, diagnostics = diagnostics, projected = projected,
      simulations = simulations, hindcast = hindcast, error = NULL)
  }, error = function(e) list(error = conditionMessage(e), stage = stage)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  results$warnings <- warnings
  results$session <- utils::sessionInfo()
  saveRDS(results, file.path(output, "workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Fleet workflow failed during {results$stage}: {results$error}")
  print(results$diagnostics)
  print(results$hindcast)
  invisible(results)
}

run_fleet_active_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/active") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  warnings <- character(); results <- list()
  withCallingHandlers(tryCatch({
    results$fit <- tinyAM::fit_tam(fleet_workflow_inputs(), years = 1983:2000, ages = 2:8,
      F_settings = list(process = "cor_rw", sd_form = ~fleet, cor_form = ~fleet,
        age_breaks = list(Gillnet = c(2,4,6,9), Trawl = NULL)),
      M_settings = list(process = "off", mu_supplied = ~I(.3)), silent = TRUE)
    results$diagnostics <- tinyAM::check_tam(results$fit)
    tinyAM::vis_tam(results$fit, open_file = FALSE,
      output_file = file.path(output,"dashboard.html"),render_args = list(quiet = TRUE))
  }, error = function(e) { results$error <<- conditionMessage(e) }),
    warning = function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })
  results$warnings <- warnings; results$session <- utils::sessionInfo()
  saveRDS(results,file.path(output,"workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Active-F rendering check failed: {results$error}")
  print(results$diagnostics)
  invisible(results)
}

run_fleet_scoped_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/scoped") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  d <- fleet_workflow_inputs()
  d$index <- rbind(transform(d$index, survey = "RV"),
    transform(d$index, survey = "Sentinel", obs = .8 * obs))
  settings <- list(data = d, years = 1983:2000, ages = 2:8,
    N_settings = list(process = "off"),
    rec_settings = list(process = "off", mu_form = ~ factor(year)),
    M_settings = list(process = "off", mu_supplied = ~ I(.3)), silent = TRUE)
  warnings <- character(); results <- list(); stage <- "fit"
  withCallingHandlers(tryCatch({
    results$base <- do.call(tinyAM::fit_tam, c(settings, list(
      F_settings = list(process = "off", mu_form = ~ fleet + factor(year) + factor(age)),
      index_settings = list(q_form = ~ 0 + survey:q_block))))
    results$scoped <- do.call(tinyAM::fit_tam, c(settings, list(
      F_settings = list(process = "off",
        mu_form = ~ 0 + factor(age) + factor(year) + for_fleet("Trawl", ~ 1)),
      index_settings = list(q_form = ~ 0 + for_survey("RV", ~ q_block) +
        for_survey("Sentinel", ~ q_block)))))
    results$diagnostics <- tinyAM::check_tam(results$base, results$scoped)
    results$AIC <- stats::AIC(results$base, results$scoped)
    results$cAIC <- tinyAM::cAIC(results$base, results$scoped, details = TRUE)
    stage <- "hindcast"
    results$hindcast <- tinyAM::fit_hindcast(results$scoped, folds = 1, progress = FALSE)
    stage <- "osa"
    before <- results$scoped$parameter_values
    results$osa <- stats::residuals(results$scoped, type = "osa", subset = 1:3, trace = FALSE)
    stopifnot(identical(before, results$scoped$parameter_values))
    stage <- "dashboard"
    tinyAM::vis_tam(Base = results$base, Scoped = results$scoped, open_file = FALSE,
      output_file = file.path(output, "dashboard.html"), render_args = list(quiet = TRUE))
  }, error = function(e) { results$error <<- conditionMessage(e); results$stage <<- stage }),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  results$warnings <- warnings; results$session <- utils::sessionInfo()
  saveRDS(results, file.path(output, "workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Scoped workflow failed during {results$stage}: {results$error}")
  print(results$diagnostics); print(results$cAIC); print(results$hindcast)
  invisible(results)
}

# Fully simulated SCAA example: known F shapes, supplied M, stable survey
# selectivity, and finite lognormal observation errors. No assessment is changed.
fleet_curve_example <- function(seed = 47281, family = "curves") {
  family <- match.arg(family,c("curves","smooths"))
  grid <- expand.grid(year = 2000:2019, age = 1:12)
  stock <- transform(grid, obs = .2*age^1.5, type = "stock", fleet = NA_character_, unit = "kg")
  d <- list(catch = rbind(transform(grid, obs = 1, fleet = "Trawl", unit = "thousand fish"),
      transform(grid, obs = 1, fleet = "Gillnet", unit = "thousand fish")),
    index = transform(grid, obs = 1, survey = "RV", samp_time = .5, unit = "fish per tow"),
    weight = rbind(stock, transform(stock, obs = .8*obs, type = "catch", fleet = "Trawl"),
      transform(stock, obs = 1.2*obs, type = "catch", fleet = "Gillnet")),
    maturity = transform(grid, obs = plogis(age-4)))
  settings <- list(years = 2000:2019, ages = 1:12,
    N_settings = list(process = "off"), rec_settings = list(process = "iid", mu_form = ~1),
    F_settings = list(process = "off", mu_form = ~ 0 + fleet:factor(year) +
      for_fleet("Trawl", ~ 0 + logistic(age)) + for_fleet("Gillnet", ~ 0 + warped_normal(age))),
    M_settings = list(process = "off", mu_supplied = ~ I(.2)), index_settings = list(q_form = ~1))
  if (family == "smooths") {
    settings$F_settings$mu_form <- ~0+fleet:factor(year)+
      for_fleet("Trawl",~0+s(age,k=5))+for_fleet("Gillnet",~0+s(age,k=5))
    settings$index_settings$q_form <- ~0+for_survey("RV",~s(age,k=4))
  }
  dat <- do.call(tinyAM::prepare_tam, c(list(data = d), settings))
  p <- tinyAM::make_par(dat)
  p$log_r0[] <- log(1000); p$rec_beta[] <- log(1000); p$log_r[] <- log(1000)
  p$log_sd_r[] <- log(.25); p$log_sd_catch[] <- p$log_sd_index[] <- log(.12)
  p$log_q[] <- 0
  p$log_q[grepl("Intercept",names(p$log_q))] <- log(.03)
  mean <- log(ifelse(dat$F_data$fleet == "Trawl", .25, .12)) + .2*sin((dat$F_data$year-2000)/3)
  p$log_mu_f[] <- qr.solve(dat$F_modmat, mean)
  for (term in dat$F_terms) {
    if (term$type == "smooth") {
      p[[term$sd_parameter]][] <- log(.2)
      p$log_mu_f[colnames(term$smooth_fixed)] <- if (term$scope$name == "Trawl") .3 else -.15
    } else if (term$type == "logistic") {
      p[[term$a50_parameter]][] <- 3
      p[[term$slope_parameter]][] <- log(1.2)
    } else {
      p[[term$a50_parameter]][] <- 3.5
      p[[term$scale_parameter]][] <- log(1.5)
      p[[term$shape_parameter]][] <- .35
    }
  }
  for (term in dat$q_terms) if (term$type == "smooth") {
    p[[term$sd_parameter]][] <- log(.15)
    p$log_q[colnames(term$smooth_fixed)] <- .1
  }
  set.seed(seed)
  draw <- tinyAM::nll_fun(p, dat, simulate = TRUE)
  p$log_r <- draw$log_r
  names <- intersect(tinyAM:::.formula_random_parameters(dat),names(draw))
  p[names] <- draw[names]
  objective <- RTMB::MakeADFun(function(p) tinyAM::nll_fun(p,dat),p,silent = TRUE)
  truth <- objective$report()
  for (type in c("catch","index")) {
    d[[type]] <- dat$obs[[type]]
    d[[type]]$obs <- exp(draw$log_obs[dat$obs_map$type == type])
  }
  list(data = d, settings = settings, parameters = p, truth = truth, seed = seed, family = family)
}

run_fleet_curve_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/curves") {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  warnings <- character(); results <- list(); stage <- "generation"
  withCallingHandlers(tryCatch({
    results$example <- fleet_curve_example()
    stage <- "fit"
    results$fit <- do.call(tinyAM::fit_tam,c(list(data=results$example$data,silent=TRUE),results$example$settings))
    results$diagnostics <- tinyAM::check_tam(results$fit)
    results$AIC <- stats::AIC(results$fit)
    results$cAIC <- tinyAM::cAIC(results$fit, details = TRUE)
    stage <- "projection"
    results$projected <- update(results$fit, start_par = results$fit,
      proj_settings = list(n_proj = 2, n_mean = 3, catch = 100), silent = TRUE)
    stage <- "simulation"
    set.seed(5831)
    results$simulations <- tinyAM::sim_tam(results$fit,n=2,par_uncertainty="none",redraw_random=TRUE)
    stage <- "hindcast"
    results$hindcast <- tinyAM::fit_hindcast(results$fit,folds=1,progress=FALSE)
    stage <- "dashboard"
    tinyAM::vis_tam(results$fit,open_file=FALSE,output_file=file.path(output,"dashboard.html"),
      render_args=list(quiet=TRUE))
  },error=function(e){results$error <<- conditionMessage(e); results$stage <<- stage}),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")})
  results$warnings <- warnings; results$session <- utils::sessionInfo()
  saveRDS(results,file.path(output,"workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Curve workflow failed during {results$stage}: {results$error}")
  print(results$diagnostics); print(results$cAIC); print(results$hindcast)
  invisible(results)
}

run_fleet_smooth_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/smooths") {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  warnings <- character(); results <- list(); stage <- "generation"
  withCallingHandlers(tryCatch({
    results$example <- fleet_curve_example(family="smooths")
    stage <- "fit"
    results$fit <- do.call(tinyAM::fit_tam,c(list(data=results$example$data,silent=TRUE),results$example$settings))
    results$diagnostics <- tinyAM::check_tam(results$fit)
    results$cAIC <- tinyAM::cAIC(results$fit,details=TRUE)
    stage <- "projection"
    results$projected <- update(results$fit,start_par=results$fit,
      proj_settings=list(n_proj=2,n_mean=3,catch=100),silent=TRUE)
    stage <- "simulation"
    set.seed(5831)
    results$simulations <- tinyAM::sim_tam(results$fit,n=2,par_uncertainty="none",redraw_random=TRUE)
    stage <- "hindcast"
    results$hindcast <- tinyAM::fit_hindcast(results$fit,folds=1,progress=FALSE)
    stage <- "residuals"
    before <- results$fit$parameter_values
    results$osa <- stats::residuals(results$fit,type="osa",subset=1:3,trace=FALSE)
    stopifnot(identical(before,results$fit$parameter_values))
    stage <- "dashboard"
    tinyAM::vis_tam(results$fit,open_file=FALSE,output_file=file.path(output,"dashboard.html"),
      render_args=list(quiet=TRUE))
  },error=function(e){results$error <<- conditionMessage(e); results$stage <<- stage}),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")})
  results$warnings <- warnings; results$session <- utils::sessionInfo()
  saveRDS(results,file.path(output,"workflow.rds"))
  if (!is.null(results$error)) cli::cli_abort("Smooth workflow failed during {results$stage}: {results$error}")
  print(results$diagnostics); print(results$cAIC); print(results$hindcast)
  invisible(results)
}

if (sys.nframe() == 0L) {
  mode <- commandArgs(TRUE)
  if (!length(mode)) mode <- "core"
  run <- list(core = run_fleet_core_workflow, active = run_fleet_active_workflow,
    scoped = run_fleet_scoped_workflow, curves = run_fleet_curve_workflow,
    smooths = run_fleet_smooth_workflow)
  if (!mode[1L] %in% names(run)) cli::cli_abort("Available validation modes: {paste(names(run), collapse = ', ')}.")
  run[[mode[1L]]]()
}
