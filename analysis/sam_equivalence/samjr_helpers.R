# Analysis-only conversion and checks. Never used by package likelihood code.
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
  write.csv(data.frame(file = files, url = urls,
    md5 = unname(tools::md5sum(file.path(folder, files)))), file.path(folder, "provenance.csv"), row.names = FALSE)
}

samjr_inputs <- function(native) {
  src <- native$data; cf <- native$conf
  stopifnot(cf$stockRecruitmentModelCode == 0, cf$noScaledYears == 0,
    all(as.character(cf$obsCorStruct) == "ID"), all(as.vector(src$propF) == as.vector(src$propM)),
    cf$stockWeightModel == 0, cf$catchWeightModel == 0, cf$matureModel == 0,
    cf$mortalityModel == 0, all(cf$fracMixN == 0), cf$fracMixF == 0,
    all(cf$fracMixObs == 0), all(is.finite(src$logobs)))
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
  stopifnot(identical(active, seq.int(min(active), length(ages))))
  blocks <- rle(f_states[active])
  stopifnot(!anyDuplicated(blocks$values))
  starts <- ages[active[c(1L, head(cumsum(blocks$lengths), -1L) + 1L)]]
  breaks <- if (identical(as.numeric(f_states), seq_along(ages) - 1 + 0)) NULL else c(starts, max(ages) + 1L)
  pop <- transform(grid, M = as.vector(src$natMor), N_sd_key = key(n_keys[match(age, ages)]),
    F_sd_key = key(f_keys[match(age, ages)]))
  obs <- list(catch = catch, index = index, weight = weight, maturity = maturity, pop_covariates = pop)
  settings <- list(N_settings = list(process = "iid", init = "latent",
      sd_form = if (length(unique(n_keys[-1L])) == 1L) ~1 else ~0 + N_sd_key),
    rec_settings = list(process = "rw", init = "latent"),
    F_settings = c(list(process = if (cf$corFlag == 0) "rw" else "cor_rw",
      mu_form = NULL, sd_form = if (length(unique(f_keys[active])) == 1L) ~1 else ~0 + F_sd_key,
      mean_ages = seq(cf$fbarRange[1], cf$fbarRange[2])),
      if (!is.null(breaks)) list(age_breaks = breaks)),
    M_settings = list(process = "off", mu_supplied = ~M),
    catch_settings = list(sd_form = if (any(catch$extra != 0)) ~0 + sd_key + extra else ~0 + sd_key),
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
    if (any(!grouped)) par[[name]][!grouped] <- native$pl$logXtraSd
  }
  tiny <- RTMB::MakeADFun(function(p) nll_fun(p, dat), par, silent = TRUE)
  sampar <- stockassessment::defpar(src, cf)
  for (name in names(sampar)) sampar[[name]][] <- native$pl[[name]]
  sampar$missing <- numeric(0)
  samdata <- c(src, list(forecast = list(), referencepoints = list()), cf,
    list(simFlag = c(1L, 1L), reportingLevel = 0L, resFlag = 0L))
  maps <- list(logFScaleMSY = factor(NA), implicitFunctionDelta = factor(NA), splinePenalty = factor(NA))
  sam <- TMB::MakeADFun(samdata, sampar, DLL = "stockassessment", map = maps, silent = TRUE)
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
      if (any(!grouped)) out$logXtraSd[] <- p[[name]][!grouped]
    }
    out
  }
  marker <- par; n <- 0L
  for (name in names(marker)) { marker[[name]][] <- n + seq_along(marker[[name]]); n <- n + length(marker[[name]]) }
  permutation <- as.integer(pack(to_sam(marker), sam$par))
  stopifnot(identical(sort(permutation), seq_along(tiny$par)))
  order <- order(permutation)
  checks <- list()
  add <- function(test, value, tolerance) {
    checks[[length(checks) + 1L]] <<- data.frame(case = label, test = test,
      difference = value, tolerance = tolerance, passed = is.finite(value) && value <= tolerance)
  }
  for (attempt in 1:2) {
    v <- tiny$par
    if (attempt == 2L) { set.seed(761); v <- v + rnorm(length(v), sd = .001) }
    y <- pack(to_sam(tiny$env$parList(v)), sam$par)
    add(paste0("joint_", attempt), abs(tiny$fn(v) - sam$fn(y)), 1e-7)
    g <- as.numeric(sam$gr(y))[order]
    add(paste0("gradient_", attempt), max(abs(tiny$gr(v) - g) / (1 + abs(g))), 1e-6)
    h <- sam$he(y)[order, order]
    add(paste0("hessian_", attempt), max(abs(tiny$he(v) - h) / (1 + abs(h))), 1e-6)
  }
  fit <- do.call(fit_tam, c(list(data = input$obs, years = input$years, ages = input$ages,
    start_par = par, silent = TRUE), input$settings))
  add("optimized_objective", abs(fit$opt$objective - native$opt$objective), 1e-4)
  add("relative_N", max(abs(fit$rep$N / exp(t(native$pl$logN)) - 1)), 1e-4)
  active <- cf$keyLogFsta[1L, ] >= 0
  nativeF <- matrix(0, length(input$years), length(input$ages))
  nativeF[, active] <- exp(t(native$pl$logF[cf$keyLogFsta[1L, active] + 1L, , drop = FALSE]))
  add("F_surface", max(abs(fit$rep$F - nativeF) / (1 + abs(nativeF))), 1e-4)
  list(checks = do.call(rbind, checks), fit = fit, native = native, input = input)
}
