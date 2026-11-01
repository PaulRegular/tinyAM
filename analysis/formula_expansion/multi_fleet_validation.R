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

# Fully simulated smooth example with supplied M, stable survey selectivity
# and finite lognormal observation errors. No assessment is changed.
fleet_smooth_example <- function(seed = 47281) {
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
    F_settings = list(process = "off", mu_form = ~0+fleet:factor(year)+
      for_fleet("Trawl",~0+s(age,k=5))+for_fleet("Gillnet",~0+s(age,k=5))),
    M_settings = list(process = "off", mu_supplied = ~ I(.2)),
    index_settings = list(q_form = ~0+for_survey("RV",~s(age,k=4))))
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
  list(data = d, settings = settings, parameters = p, truth = truth, seed = seed, family = "smooths")
}

run_fleet_smooth_workflow <- function(output = "analysis/formula_expansion/results/multi_fleet/smooths") {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  warnings <- character(); results <- list(); stage <- "generation"
  withCallingHandlers(tryCatch({
    results$example <- fleet_smooth_example()
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

run_fleet_herring_trial <- function(database_dir,
  output="analysis/formula_expansion/results/multi_fleet/herring") {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  source("analysis/comp_assessments/R/read_database.R",local=TRUE)
  source("analysis/comp_assessments/R/database_to_tam_obs.R",local=TRUE)
  source("analysis/comp_assessments/R/database_to_tam_ref.R",local=TRUE)
  id <- "ices_herring_west_scotland_ireland_2020"
  paths <- file.path(database_dir,.assessment_database_files)
  if (!all(file.exists(paths))) cli::cli_abort("Supply the source database containing the historical 2020 herring assessment.")
  hashes <- tools::md5sum(paths)
  read_source_table <- function(path) {
    lines <- if (basename(path) == "assumptions.csv") readLines(path,warn=FALSE) else character()
    if (length(lines) == 1L && grepl('"\\r\\n"',lines,fixed=TRUE)) {
      # A local unfinished import encoded record separators literally. Decode
      # in memory for this trial only; the user's source file stays untouched.
      return(utils::read.csv(text=gsub("\\r\\n","\n",lines,fixed=TRUE),
        stringsAsFactors=FALSE,na.strings=c("","NA"),check.names=FALSE))
    }
    utils::read.csv(path,stringsAsFactors=FALSE,na.strings=c("","NA"),check.names=FALSE)
  }
  source_cache <- file.path(output,"source.rds")
  if (file.exists(source_cache) && identical(readRDS(source_cache)$hashes,hashes)) {
    source <- readRDS(source_cache)$source
  } else {
    database <- .assessment_database(function(name) {
      cli::cli_inform("Reading source {name}.")
      x <- read_source_table(file.path(database_dir,name))
      if ("assessment_id" %in% names(x)) x <- x[!is.na(x$assessment_id) & x$assessment_id == id,,drop=FALSE]
      x
    },source=normalizePath(database_dir,winslash="/"))
    source <- read_assessment(id,database)
    saveRDS(list(source=source,hashes=hashes),source_cache)
  }
  # Also repair the metadata of an earlier source-only cache, never its inputs.
  assumptions <- read_source_table(file.path(database_dir,"assumptions.csv"))
  source$assumptions <- assumptions[assumptions$assessment_id == id,,drop=FALSE]
  result <- list(source=source,source_hashes=hashes); warnings <- character(); stage <- "preparation"
  withCallingHandlers(tryCatch({
    raw <- source$inputs
    # Missing source index cells are unavailable observations, not zeros. The
    # canonical converter expects observed rows and expands missing grid cells.
    missing_index <- raw$type == "index" & is.na(raw$value)
    result$missing_source_index <- raw[missing_index,,drop=FALSE]
    raw <- raw[!missing_index,,drop=FALSE]
    raw$unit[raw$type == "index"] <- "native survey index"
    fleets <- unique(raw$fleet[raw$type == "catch" & raw$measure == "numbers_at_age"])
    catch_row <- raw$type == "catch" & raw$measure == "numbers_at_age"
    # The canonical converter selects one stream at a time. Preserve both streams
    # here; do not alter the settled aggregate recipe or canonical database.
    converted <- lapply(fleets,function(f)database_to_tam_obs(id,
      raw[!catch_row | (!is.na(raw$fleet) & raw$fleet == f),,drop=FALSE],
      years=1957:2019,ages=1:9,surveys=c("MS_HERAS","IBTS_Q1","IBTS_Q4"),
      assumptions=source$assumptions))
    obs <- converted[[1]]
    missing <- result$missing_source_index
    for (j in seq_len(nrow(missing))) {
      key <- obs$index$year == as.numeric(missing$year[j]) & obs$index$age == as.numeric(missing$age[j]) & obs$index$survey == missing$survey[j]
      if (any(key)) {obs$index$obs[key] <- NA_real_; next}
      row <- obs$index[which(obs$index$survey == missing$survey[j])[1L],,drop=FALSE]
      row$year <- as.numeric(missing$year[j]); row$age <- as.numeric(missing$age[j]); row$obs <- NA_real_
      obs$index <- rbind(obs$index,row)
    }
    obs$catch <- do.call(rbind,lapply(seq_along(fleets),function(i)
      transform(converted[[i]]$catch,fleet=fleets[i])))
    stock <- transform(obs$weight,type="stock",fleet=NA_character_)
    catch_weights <- lapply(fleets,function(f) {
      x <- raw[raw$type == "catch_weight" & !is.na(raw$fleet) & raw$fleet == f,,drop=FALSE]
      data.frame(year=as.numeric(x$year),age=as.numeric(x$age),
        obs=mapply(.translation_weight_to_kg,as.numeric(x$value),x$unit),
        unit="kg",type="catch",fleet=f)
    })
    obs$weight <- do.call(rbind,c(list(stock),catch_weights))
    obs$maturity$spawn_time <- .67
    i <- obs$index
    i$q_key <- interaction(i$survey,ifelse(i$age >= 8,"8plus",paste0("age",i$age)),drop=TRUE)
    sd_age <- ifelse(i$survey == "MS_HERAS",ifelse(i$age == 1,"age1",ifelse(i$age == 2,"age2",ifelse(i$age <= 6,"ages3_6","ages7_9"))),
      ifelse(i$survey == "IBTS_Q1",ifelse(i$age == 2,"age2",ifelse(i$age <= 5,"ages3_5",ifelse(i$age <= 7,"ages6_7","ages8_9"))),
        ifelse(i$age == 2,"age2",ifelse(i$age <= 6,"ages3_6","ages7_9"))))
    i$sd_block <- interaction(i$survey,sd_age,drop=TRUE)
    # Unknown survey units remain unspecified, rather than pretending they are counts.
    i$unit <- NULL; obs$index <- i
    result$data <- obs
    tinyAM::check_tam_data(obs)
    stage <- "fit"
    result$fit <- tinyAM::fit_tam(obs,years=1957:2019,ages=1:9,
      N_settings=list(process="iid",init="exp"),rec_settings=list(process="iid",mu_form=~1),
      F_settings=list(process="cor_rw",sd_form=~fleet,cor_form=~fleet,mean_ages=3:6),
      M_settings=list(process="off",mu_supplied=~M_assumption),
      catch_settings=list(sd_form=~fleet),index_settings=list(q_form=~0+q_key,sd_form=~0+sd_block),silent=TRUE)
    result$diagnostics <- tinyAM::check_tam(result$fit)
    result$population <- tinyAM::tidy_pop(result$fit)
    stage <- "dashboard"
    result$reference <- render_fleet_herring(result,output)
  },error=function(e){result$error <<- conditionMessage(e); result$stage <<- stage}),
    warning=function(w){warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")})
  result$warnings <- warnings; result$session <- utils::sessionInfo()
  previous <- file.path(output,"trial.rds")
  if (file.exists(previous)) {
    prior <- readRDS(previous)
    result$earlier_attempts <- c(prior$earlier_attempts,list(prior[c("stage","error","warnings")]))
  }
  saveRDS(result,file.path(output,"trial.rds"))
  if (!is.null(result$error)) cli::cli_abort("Herring trial stopped during {result$stage}: {result$error}")
  print(result$diagnostics)
  invisible(result)
}

render_fleet_herring <- function(result,output="analysis/formula_expansion/results/multi_fleet/herring") {
  source("analysis/comp_assessments/R/database_to_tam_obs.R",local=TRUE)
  source("analysis/comp_assessments/R/compare_assessment.R",local=TRUE)
  source("analysis/comp_assessments/R/database_to_tam_ref.R",local=TRUE)
  outputs <- result$source$outputs
  area_to_fleet <- c("ICES Division 6.a North"="Catch 6aN",
    "ICES divisions 6.a South, 7.b-c"="Catch 6aS/7bc")
  f <- outputs$measure == "fishing_mortality_at_age"
  if (any(!outputs$region[f] %in% names(area_to_fleet))) stop("Unknown source F area; do not infer its fleet.")
  outputs$fleet[f] <- unname(area_to_fleet[outputs$region[f]])
  reference <- database_to_tam_ref("ices_herring_west_scotland_ireland_2020",outputs,
    obs=result$data,years=1957:2019,ages=1:9,age_plus_group=9,template=result$fit,
    assumptions=result$source$assumptions)
  tinyAM::vis_tam(tinyAM=result$fit,Official=reference,open_file=FALSE,
    output_file=file.path(output,"dashboard.html"),render_args=list(quiet=TRUE),
    background=c("### Historical west-of-Scotland/Ireland herring trial",
      "Two source catch fleets and their own catch weights are retained. Separate correlated F random walks share one population.",
      "This is a parsimonious fleet trial, not an exact official specification: fleet SD/correlation sharing is explicit, source F bindings are unresolved, survey observation covariance is unavailable, and no executed source fit was recovered.",
      "Catch units are inferred from the source record; survey units remain unspecified. No rescue tuning or canonical recipe changes were applied."))
  invisible(reference)
}

if (sys.nframe() == 0L) {
  mode <- commandArgs(TRUE)
  if (!length(mode)) mode <- "core"
  if (mode[1L] == "herring-dashboard") {
    output <- if(length(mode)>1L)mode[2L] else "analysis/formula_expansion/results/multi_fleet/herring"
    result <- readRDS(file.path(output,"trial.rds"))
    result$reference <- render_fleet_herring(result,output)
    result$earlier_attempts <- c(result$earlier_attempts,list(result[c("stage","error","warnings")]))
    result$error <- result$stage <- NULL
    saveRDS(result,file.path(output,"trial.rds"))
    quit(status=0L)
  }
  run <- list(core = run_fleet_core_workflow, active = run_fleet_active_workflow,
    scoped = run_fleet_scoped_workflow,
    smooths = run_fleet_smooth_workflow,herring=run_fleet_herring_trial)
  if (!mode[1L] %in% names(run)) cli::cli_abort("Available validation modes: {paste(names(run), collapse = ', ')}.")
  if (mode[1L] == "herring") {
    if (length(mode) < 2L) cli::cli_abort("Supply the source assessment database directory after herring.")
    run_fleet_herring_trial(mode[2L])
  } else run[[mode[1L]]]()
}
