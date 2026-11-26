# Keep the production objective and fitter. Change one summand in private
# copies of the functions so the ablation uses identical parameter maps,
# optimizer, starts, positive density and conditional lognormal SDs.
delta_occurrence_ablated_fit <- function() {
  fit <- tinyAM::fit_tam
  nll <- tinyAM::nll_fun
  occurrence <- list(quote(delta$log_p_zero[dat$delta$zero]), quote(delta$log_p_positive[positive]))
  replaced <- integer(length(occurrence))
  drop_occurrence <- function(code) {
    index <- which(vapply(occurrence, identical, logical(1), y = code))
    if (length(index)) {
      replaced[index] <<- replaced[index] + 1L
      return(0)
    }
    if (!is.call(code)) return(code)
    entries <- as.list(code)
    for (i in seq_along(entries)[-1L]) {
      if (is.call(entries[[i]])) {
        entries[i] <- list(drop_occurrence(entries[[i]]))
      }
    }
    as.call(entries)
  }
  body(nll) <- drop_occurrence(body(nll))
  if (any(replaced != 1L)) stop("The observation likelihood has changed. Audit the ablation.")
  private <- new.env(parent = environment(fit))
  private$nll_fun <- nll
  environment(fit) <- private
  fit
}

delta_spatial_fit <- function(simulation, plus_age, arm) {
  family <- delta_study_family(arm)
  fit <- tinyAM::fit_tam(simulation$data, ages = seq_len(plus_age), silent = TRUE,
    N_settings = list(process = "off", init = "exp"),
    R_settings = list(process = "rw", sd = 0.50),
    F_settings = list(process = "off", mu_form = ~ factor(year) + normal(age), mean_ages = 5:10),
    M_settings = list(process = "off", mu_supplied = ~ I(0.30)),
    catch_settings = list(sd_form = ~ 1, family = family),
    index_settings = list(q_form = ~ logistic(age), sd_form = ~ 1, family = family),
    start_par = list(log_mu_f = setNames(log(0.12), "(Intercept)")))
  checks <- tinyAM::check_tam(fit)
  reference <- simulation$reference
  numbers <- fit$pop$N[c("year", "age", "est")]
  fishing <- tinyAM::get_f(fit)[c("year", "age", "est")]
  ssb <- tinyAM::get_ssb(fit)[c("year", "est", "lwr", "upr")]
  recruitment <- tinyAM::get_recruitment(fit)[c("year", "est", "lwr", "upr")]
  input <- lapply(fit$obs_pred[c("catch", "index")], function(d) d[c("year", "age", "obs")])
  data_metrics <- function(estimate, truth, label) {
    use <- truth > 0 & is.finite(truth) & is.finite(estimate) & estimate > 0
    error <- log(estimate[use] / truth[use])
    data.frame(quantity = label, n = sum(use), available = length(truth),
      absolute = mean(abs(error)), bias = mean(error))
  }
  old <- stats::aggregate(est ~ year, subset(numbers, age >= 14L & year %in% delta_older_window), sum)
  target <- stats::aggregate(N ~ year, subset(reference, age >= 14L), sum)
  young <- reference$age == min(reference$age)
  old_est <- old$est
  old_true <- target$N[match(old$year, target$year)]
  mean20 <- subset(numbers, age == 20L & year %in% delta_older_window)
  true20 <- subset(reference, age == 20L)
  f <- merge(subset(fishing, age %in% 5:10), reference, by = c("year", "age"))
  fishing_annual <- do.call(rbind, lapply(split(f, f$year), function(d) data.frame(
    year = d$year[1L], F_est = stats::weighted.mean(d$est, d$N),
    F_realized = stats::weighted.mean(d$F_realized, d$N),
    F_reference = stats::weighted.mean(d$F_reference, d$N))))
  ssb$truth <- as.numeric(tapply(reference$ssb, reference$year, sum))[match(ssb$year, delta_years)]
  recruitment$truth <- reference$N[young][match(recruitment$year, reference$year[young])]
  quantities <- rbind(
    data_metrics(old_est, old_true, "Age 14+ abundance"),
    if (nrow(mean20)) data_metrics(mean20$est, true20$N[match(mean20$year, true20$year)], "Age 20+ abundance"),
    data_metrics(ssb$est, ssb$truth, "Spawning biomass"),
    data_metrics(recruitment$est, recruitment$truth, "Recruitment"),
    data_metrics(fishing_annual$F_est, fishing_annual$F_realized, "Fishing mortality"))
  quantities$coverage <- NA_real_
  quantities$coverage[quantities$quantity == "Spawning biomass"] <-
    mean(with(ssb, truth >= lwr & truth <= upr))
  quantities$coverage[quantities$quantity == "Recruitment"] <-
    mean(with(recruitment, truth >= lwr & truth <= upr))
  list(converged = isTRUE(checks$is_converged), numerical = checks$numerical,
    gradient = checks$max_gradient, optimizer = fit$opt$convergence, objective = fit$opt$objective,
    merged_observations = input, metrics = quantities, F = fishing_annual, ssb = ssb,
    recruitment = recruitment, N = numbers, parameters = delta_study_parameters(fit))
}
delta_study_family <- function(arm, y50 = NULL) {
  switch(arm,
    lognormal = tinyAM::lognormal(),
    full = tinyAM::delta_lognormal(type = "poisson-link", y50 = y50),
    positive = tinyAM::delta_lognormal(type = "poisson-link", y50 = y50),
    gengamma = tinyAM::delta_gengamma(type = "poisson-link", y50 = y50),
    stop("Unknown study arm: ", arm))
}

delta_study_parameters <- function(fit) {
  fit$fixed_par[fit$fixed_par$par %in% c("y50_catch", "y50_index", "gengamma_Q"), ,
    drop = FALSE]
}
