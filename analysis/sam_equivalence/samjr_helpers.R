# Analysis-only conversion and checks. Never used by package likelihood code.
samjr_fixture_hash <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  while (length(lines) && !nzchar(tail(lines, 1L))) lines <- head(lines, -1L)
  path <- tempfile(); on.exit(unlink(path))
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), path)
  unname(tools::md5sum(path))
}

samjr_fixture <- function(case, folder) {
  pin <- "0384774a3f2f8ec3dc3330375998bd73ef6bbe7c"
  files <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat", "mo.dat", "nm.dat", "pf.dat", "pm.dat", "sw.dat", "survey.dat", "script.R")
  if (case %in% c("nscodXtraSd", "codIN3Ben")) files <- c(files, "model.cfg")
  if (case == "nscodXtraSd") files <- c(files, "model2.cfg")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  urls <- paste0("https://raw.githubusercontent.com/fishfollower/samjr/", pin, "/testmore/", case, "/", files)
  for (i in seq_along(files)) if (!file.exists(file.path(folder, files[i]))) {
    utils::download.file(urls[i], file.path(folder, files[i]), mode = "wb", quiet = TRUE)
  }
  manifest <- read.csv("analysis/sam_equivalence/fixtures.csv", stringsAsFactors = FALSE)
  expected <- manifest$canonical_md5[match(paste(case, files, sep = "/"), manifest$path)]
  actual <- vapply(file.path(folder, files), samjr_fixture_hash, character(1))
  stopifnot(!anyNA(expected), identical(unname(actual), expected))
  write.csv(data.frame(file = files, url = urls,
    md5 = unname(tools::md5sum(file.path(folder, files)))), file.path(folder, "provenance.csv"), row.names = FALSE)
}

samjr_inputs <- function(native) {
  src <- native$data; cf <- native$conf
  stopifnot(cf$stockRecruitmentModelCode == 0, cf$noScaledYears == 0,
    cf$initState == 0, all(cf$logNMeanAssumption == 0),
    all(as.character(cf$obsCorStruct) == "ID"), all(as.vector(src$propF) == as.vector(src$propM)),
    cf$stockWeightModel == 0, cf$catchWeightModel == 0, cf$matureModel == 0,
    cf$mortalityModel == 0, all(cf$fracMixN == 0), cf$fracMixF == 0,
    all(cf$fracMixObs == 0), all(is.finite(src$logobs) | is.na(src$logobs)))
  years <- as.integer(src$years); ages <- cf$minAge:cf$maxAge
  grid <- expand.grid(year = years, age = ages)
  aux <- as.data.frame(src$aux); names(aux)[1:3] <- c("year", "fleet", "age")
  key <- function(v) factor(ifelse(v < 0, NA, paste0("K", v)))
  catches <- aux$fleet == 1L
  j <- match(paste(grid$year, grid$age), paste(aux$year[catches], aux$age[catches]))
  catch <- transform(grid, obs = exp(src$logobs[catches])[j], unit = "thousand fish")
  catch$sd_key <- key(cf$keyVarObs[1L, match(catch$age, ages)])
  catch$extra <- 0
  if (nrow(cf$keyXtraSd)) {
    extra <- cf$keyXtraSd
    stopifnot(all(extra[, 1L] == 1), length(unique(extra[, 4L])) == 1L)
    catch$extra <- as.integer(paste(catch$year, catch$age) %in% paste(extra[, 2L], extra[, 3L]))
  }
  i <- which(!catches); ia <- match(aux$age[i], ages)
  index <- data.frame(year = aux$year[i], age = aux$age[i],
    survey = paste0("Survey ", aux$fleet[i] - 1L), obs = exp(src$logobs[i]),
    samp_time = src$sampleTimes[aux$fleet[i]])
  cells <- cbind(aux$fleet[i], ia)
  index$q_key <- key(cf$keyLogFpar[cells]); index$sd_key <- key(cf$keyVarObs[cells])
  stopifnot(all(cf$keyQpow[cells] < 0), all(cf$predVarObsLink[!is.na(cf$predVarObsLink)] < 0))
  index$relative_sd <- ifelse(is.finite(src$weight[i]), 1 / sqrt(src$weight[i]), 1)
  weights <- matrix(src$catchMeanWeight, nrow = dim(src$catchMeanWeight)[1L])
  if (nrow(weights) == length(years) - 1L) {
    stopifnot(!any(aux$year[catches] == max(years)))
    # Survey-only final year: padding is unused by fitting; never compare its yield.
    weights <- rbind(weights, weights[nrow(weights), ])
  }
  stopifnot(identical(dim(weights), dim(src$stockMeanWeight)))
  weight <- rbind(transform(grid, obs = as.vector(src$stockMeanWeight), unit = "kg", type = "stock"),
    transform(grid, obs = as.vector(weights), unit = "kg", type = "catch"))
  maturity <- transform(grid, obs = as.vector(src$propMat), spawn_time = as.vector(src$propF))
  n_keys <- cf$keyVarLogN; n_keys[1L] <- n_keys[2L]
  f_keys <- cf$keyVarF[1L, ]; f_states <- cf$keyLogFsta[1L, ]
  active <- which(f_states >= 0)
  stopifnot(identical(as.integer(active), seq.int(min(active), length(ages))))
  blocks <- rle(f_states[active])
  stopifnot(!anyDuplicated(blocks$values))
  starts <- ages[active[c(1L, head(cumsum(blocks$lengths), -1L) + 1L)]]
  breaks <- if (identical(as.numeric(f_states), seq_along(ages) - 1 + 0)) NULL else c(starts, max(ages) + 1L)
  pop <- transform(grid, M = as.vector(src$natMor), N_sd_key = key(n_keys[match(age, ages)]),
    F_sd_key = key(f_keys[match(age, ages)]))
  obs <- list(catch = catch, index = index, weight = weight, maturity = maturity, pop_covariates = pop)
  catch_form <- if (length(unique(stats::na.omit(catch$sd_key))) == 1L) ~1 else ~0 + sd_key
  if (any(catch$extra != 0)) catch_form <- update(catch_form, ~. + extra)
  settings <- list(N_settings = list(process = "iid", init = "latent",
      sd_form = if (length(unique(n_keys[-1L])) == 1L) ~1 else ~0 + N_sd_key),
    rec_settings = list(process = "rw", init = "latent"),
    F_settings = c(list(process = if (cf$corFlag == 0) "rw" else "cor_rw",
      mu_form = NULL, sd_form = if (length(unique(f_keys[active])) == 1L) ~1 else ~0 + F_sd_key,
      mean_ages = seq(cf$fbarRange[1], cf$fbarRange[2])),
      if (!is.null(breaks)) list(age_breaks = breaks)),
    M_settings = list(process = "off", mu_supplied = ~M),
    catch_settings = list(sd_form = catch_form),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_key, sd_supplied = ~relative_sd))
  list(obs = obs, settings = settings, years = years, ages = ages, aux = aux)
}

samjr_compare <- function(native, label) {
  input <- samjr_inputs(native)
  dat <- do.call(prepare_tam, c(list(data = input$obs, years = input$years, ages = input$ages), input$settings))
  par <- make_par(dat); src <- native$data; cf <- native$conf
  key <- function(x) as.integer(sub(".*K", "", x)) + 1L
  par$log_r0 <- native$pl$logN[1L, 1L]; par$log_n0[] <- native$pl$logN[-1L, 1L]
  par$log_r[] <- native$pl$logN[1L, -1L]; par$log_n[] <- t(native$pl$logN[-1L, -1L])
  par$log_f[] <- t(native$pl$logF)
  par$log_sd_r <- native$pl$logSdLogN[1L]
  for (component in c("N", "F")) {
    name <- dat$process_sd[[component]]$parameter
    source <- if (component == "N") "logSdLogN" else "logSdLogFsta"
    indices <- if (dat$process_sd[[component]]$default) {
      if (component == "N") cf$keyVarLogN[2L] + 1L else cf$keyVarF[1L, which(cf$keyLogFsta[1L, ] >= 0)[1L]] + 1L
    } else key(names(par[[name]]))
    par[[name]][] <- native$pl[[source]][indices]
  }
  if (cf$corFlag == 2) par$atanh_rho_f <- native$pl$itrans_rho[1L]
  par$log_q[] <- native$pl$logFpar[key(names(par$log_q))]
  for (name in c("log_sd_catch", "log_sd_index")) {
    grouped <- grepl("K", names(par[[name]]))
    par[[name]][grouped] <- native$pl$logSdLogObs[key(names(par[[name]])[grouped])]
    intercept <- names(par[[name]]) == "(Intercept)"
    if (any(intercept)) par[[name]][intercept] <- native$pl$logSdLogObs[cf$keyVarObs[1L, 1L] + 1L]
    extra <- names(par[[name]]) == "extra"
    if (any(extra)) par[[name]][extra] <- native$pl$logXtraSd
  }
  tiny <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  sampar <- stockassessment::defpar(src, cf)
  for (name in names(sampar)) sampar[[name]][] <- native$pl[[name]]
  sampar$missing <- as.numeric(native$pl$missing)
  samdata <- c(src, list(forecast = list(), referencepoints = list()), cf,
    list(simFlag = c(1L, 1L), reportingLevel = 0L, resFlag = 0L))
  maps <- list(logFScaleMSY = factor(NA), implicitFunctionDelta = factor(NA), splinePenalty = factor(NA))
  sam <- TMB::MakeADFun(samdata, sampar, DLL = "stockassessment", map = maps, silent = TRUE)
  # Independent missing Normal observations integrate to one. For joint-density
  # checks only, remove their auxiliary densities and compare common states.
  # The original native fit and its full Laplace objective remain unchanged.
  missing_term <- NULL
  if (length(sampar$missing)) {
    missing_dat <- dat
    aux <- input$aux[is.na(src$logobs), ]
    fleet <- ifelse(dat$obs_map$type == "catch", 1L,
      as.integer(sub("Survey ", "", dat$obs_map$survey)) + 1L)
    rows <- match(paste(aux$year, aux$fleet, aux$age),
      paste(dat$obs_map$year, fleet, dat$obs_map$age))
    stopifnot(!anyNA(rows), !anyDuplicated(rows))
    missing_dat$is_observed[] <- FALSE; missing_dat$is_observed[rows] <- TRUE
    missing_dat$observed <- sampar$missing[order(rows)]
    missing_dat$.likelihood_component <- "observation"
    missing_term <- RTMB::MakeADFun(function(p) nll_fun(p, missing_dat), par, silent = TRUE)
  }
  pack <- function(p, template) {
    for (name in unique(names(template))) template[names(template) == name] <- as.vector(p[[name]])
    template
  }
  to_sam <- function(p) {
    out <- sampar
    out$logN[1L, 1L] <- p$log_r0; out$logN[-1L, 1L] <- p$log_n0
    out$logN[1L, -1L] <- p$log_r; out$logN[-1L, -1L] <- t(p$log_n)
    out$logF[] <- t(p$log_f); out$logSdLogN[1L] <- p$log_sd_r
    for (component in c("N", "F")) {
      name <- dat$process_sd[[component]]$parameter
      source <- if (component == "N") "logSdLogN" else "logSdLogFsta"
      indices <- if (dat$process_sd[[component]]$default) {
        if (component == "N") cf$keyVarLogN[2L] + 1L else cf$keyVarF[1L, which(cf$keyLogFsta[1L, ] >= 0)[1L]] + 1L
      } else key(names(p[[name]]))
      out[[source]][indices] <- p[[name]]
    }
    if (cf$corFlag == 2) out$itrans_rho[] <- p$atanh_rho_f
    out$logFpar[key(names(p$log_q))] <- p$log_q
    for (name in c("log_sd_catch", "log_sd_index")) {
      grouped <- grepl("K", names(p[[name]]))
      out$logSdLogObs[key(names(p[[name]])[grouped])] <- p[[name]][grouped]
      intercept <- names(p[[name]]) == "(Intercept)"
      if (any(intercept)) out$logSdLogObs[cf$keyVarObs[1L, 1L] + 1L] <- p[[name]][intercept]
      extra <- names(p[[name]]) == "extra"
      if (any(extra)) out$logXtraSd[] <- p[[name]][extra]
    }
    out
  }
  marker <- par; n <- 0L
  for (name in names(marker)) { marker[[name]][] <- n + seq_along(marker[[name]]); n <- n + length(marker[[name]]) }
  common <- which(names(sam$par) != "missing")
  permutation <- as.integer(pack(to_sam(marker), sam$par)[common])
  stopifnot(identical(sort(permutation), seq_along(tiny$par)))
  order <- common[order(permutation)]
  checks <- list()
  add <- function(test, value, tolerance) {
    checks[[length(checks) + 1L]] <<- data.frame(case = label, test = test,
      difference = value, tolerance = tolerance, passed = is.finite(value) && value <= tolerance)
  }
  for (attempt in 1:2) {
    v <- tiny$par
    if (attempt == 2L) { set.seed(761); v <- v + rnorm(length(v), sd = .001) }
    y <- pack(to_sam(tiny$env$parList(v)), sam$par)
    correction <- if (is.null(missing_term)) 0 else missing_term$fn(v)
    add(paste0("joint_", attempt), abs(tiny$fn(v) - (sam$fn(y) - correction)), 1e-7)
    g <- as.numeric(sam$gr(y))[order]
    if (!is.null(missing_term)) g <- g - as.numeric(missing_term$gr(v))
    add(paste0("gradient_", attempt), max(abs(tiny$gr(v) - g) / (1 + abs(g))), 1e-6)
    h <- sam$he(y)[order, order]
    if (!is.null(missing_term)) h <- h - missing_term$he(v)
    add(paste0("hessian_", attempt), max(abs(tiny$he(v) - h) / (1 + abs(h))), 1e-6)
  }
  fit <- do.call(fit_tam, c(list(data = input$obs, years = input$years, ages = input$ages,
    start_par = par, silent = TRUE), input$settings))
  marginal <- TMB::MakeADFun(samdata, sampar, DLL = "stockassessment", map = maps,
    random = c("logF", "logN", if (length(sampar$missing)) "missing"), silent = TRUE)
  tiny_marginal <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par,
    random = fit$obj$env$.random, map = fit$parameter_map, silent = TRUE)
  add("laplace", abs(tiny_marginal$fn(tiny_marginal$par) - marginal$fn(marginal$par)), 1e-5)
  add("n_fixed", abs(length(fit$opt$par) - length(marginal$par)), 0)
  add("n_random_common", abs(sum(fit$obj$env$lrandom()) - (sum(marginal$env$lrandom()) - length(sampar$missing))), 0)
  add("n_observed", abs(sum(dat$is_observed) - sum(is.finite(src$logobs))), 0)
  add("optimized_objective", abs(fit$opt$objective - native$opt$objective), 1e-4)
  add("numerical_success", as.numeric(!isTRUE(fit$is_converged) || native$opt$convergence != 0), 0)
  outer <- pack(to_sam(fit$obj$env$parList(par = fit$parameter_values)), marginal$par)
  add("fixed_parameters", max(abs(outer - marginal$par) / (1 + abs(marginal$par))), 1e-4)
  add("marginal_AIC", abs((-2 * as.numeric(suppressWarnings(logLik(fit))) + 2 * length(fit$opt$par)) -
    (2 * native$opt$objective + 2 * length(marginal$par))), 1e-4)
  add("relative_N", max(abs(fit$rep$N / exp(t(native$pl$logN)) - 1)), 1e-4)
  active <- cf$keyLogFsta[1L, ] >= 0
  nativeF <- matrix(0, length(input$years), length(input$ages))
  nativeF[, active] <- exp(t(native$pl$logF[cf$keyLogFsta[1L, active] + 1L, , drop = FALSE]))
  add("F_surface", max(abs(fit$rep$F - nativeF) / (1 + abs(nativeF))), 1e-4)
  fleet <- ifelse(dat$obs_map$type == "catch", 1L,
    as.integer(sub("Survey ", "", dat$obs_map$survey)) + 1L)
  observation_order <- match(paste(input$aux$year, input$aux$fleet, input$aux$age),
    paste(dat$obs_map$year, fleet, dat$obs_map$age))
  add("log_predictions", max(abs(fit$rep$log_pred[observation_order] - native$rep$predObs), na.rm = TRUE), 1e-4)
  sd <- tinyAM:::.tam_parameter_summary(fit, "Std. Error")
  sd <- pack(to_sam(sd), sam$par)[order]
  native_sd <- pack(native$plsd, sam$par)[order]
  add("parameter_SE", max(abs(sd - native_sd) / (1 + abs(native_sd))), 1e-4)
  reports <- summary(native$sdrep, "report")
  report_names <- c(ssb = "logssb", biomass = "logtsb", recruitment = "logR", total_yield_pred = "logCatch")
  for (metric in names(report_names)) {
    name <- report_names[[metric]]
    values <- reports[rownames(reports) == name, , drop = FALSE]
    table <- fit$pop[[metric]][seq_len(nrow(values)), ]
    expected <- cbind(exp(values[, 1L]), exp(values[, 1L] - qnorm(.975) * values[, 2L]),
      exp(values[, 1L] + qnorm(.975) * values[, 2L]))
    if (all(c("lwr", "upr") %in% names(table))) {
      actual <- as.matrix(table[c("est", "lwr", "upr")])
    } else { actual <- as.matrix(table["est"]); expected <- expected[, 1L, drop = FALSE] }
    add(paste0("table_", metric), max(abs(actual - expected) / (1 + abs(expected))), 1e-4)
  }
  selected <- match(dat$F_settings$mean_ages, input$ages)
  nativeN <- exp(t(native$pl$logN))
  fbar <- rowSums(nativeF[, selected, drop = FALSE] * nativeN[, selected, drop = FALSE]) /
    rowSums(nativeN[, selected, drop = FALSE])
  add("aligned_Fbar", max(abs(fit$pop$F_bar$est - fbar) / (1 + abs(fbar))), 1e-4)
  q_table <- fit$obs_pred$index
  q_fleet <- as.integer(sub("Survey ", "", q_table$survey)) + 1L
  q_key <- cf$keyLogFpar[cbind(q_fleet, match(q_table$age, input$ages))]
  available <- !is.na(q_key) & q_key >= 0
  q_key <- q_key[available] + 1L
  q <- exp(native$pl$logFpar[q_key]); q_se <- native$plsd$logFpar[q_key]
  expected <- cbind(q, q * exp(-qnorm(.975) * q_se), q * exp(qnorm(.975) * q_se))
  actual <- as.matrix(q_table[available, c("q", "q_lwr", "q_upr")])
  add("table_catchability", max(abs(actual - expected) / (1 + abs(expected))), 1e-4)
  stopifnot(identical(sort(observation_order[is.finite(src$logobs)]), sort(which(dat$is_observed))))
  # Expose a matched native-parameter mapping for subsequent controlled trials.
  perturbation <- tiny$par
  set.seed(762); perturbation <- perturbation + rnorm(length(perturbation), sd = .01)
  start <- tiny$env$parList(perturbation)
  source_start <- stockassessment::defpar(src, cf)
  converted <- to_sam(start)
  for (name in names(source_start)) source_start[[name]][] <- converted[[name]]
  if (length(source_start$missing)) source_start$missing <- converted$missing
  fitted_native <- stockassessment::sam.fit(src, cf, source_start)
  fitted_tiny <- do.call(fit_tam, c(list(data = input$obs, years = input$years, ages = input$ages,
    start_par = start, silent = TRUE), input$settings))
  add("perturbed_optimized_objective", abs(fitted_tiny$opt$objective - fitted_native$opt$objective), 1e-4)
  add("perturbed_N", max(abs(fitted_tiny$rep$N / exp(t(fitted_native$pl$logN)) - 1)), 1e-4)
  add("perturbed_numerical_success", as.numeric(!isTRUE(fitted_tiny$is_converged) || fitted_native$opt$convergence != 0), 0)
  # A common calendar cutoff, rather than peeling each fleet's own last year.
  native_peels <- stockassessment::retro(native, year = max(input$years) - 0:1, ncores = 1)
  retro <- fit_retro(fit, folds = 2, start_from_fit = TRUE, progress = FALSE)
  trials <- list()
  for (peel in native_peels) {
    terminal <- max(peel$data$years)
    if (terminal == max(input$years)) next
    attempt <- retro$fits[[1L]][[as.character(terminal)]]
    if (is.null(attempt)) {
      add(paste0("retro_", terminal, "_available"), Inf, 0)
    } else {
      add(paste0("retro_", terminal, "_objective"), abs(attempt$opt$objective - peel$opt$objective), 1e-4)
      add(paste0("retro_", terminal, "_N"), max(abs(attempt$rep$N / exp(t(peel$pl$logN)) - 1)), 1e-4)
      add(paste0("retro_", terminal, "_numerical_success"), as.numeric(!isTRUE(attempt$is_converged) || peel$opt$convergence != 0), 0)
    }
    trials[[as.character(terminal)]] <- peel
  }
  list(checks = do.call(rbind, checks), fit = fit, native = native, input = input,
    perturbed = list(tinyAM = fitted_tiny, SAM = fitted_native), retro = list(tinyAM = retro, SAM = trials))
}
