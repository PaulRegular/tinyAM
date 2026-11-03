# Render a cached comparison without fitting or calculating OSA residuals:
# Rscript analysis/sam_equivalence/render_samjr.R results/samjr/CASE/FIT.rds
args <- commandArgs(TRUE)
if (!length(args)) stop("Supply a comparison RDS path.")
pkgload::load_all(".", quiet = TRUE)
x <- readRDS(args[1L]); f <- x$fit; s <- x$native
years <- as.integer(s$data$years); ages <- s$conf$minAge:s$conf$maxAge
grid <- expand.grid(year = years, age = ages)
N <- t(exp(s$pl$logN)); F <- Fse <- matrix(0, length(years), length(ages))
key <- s$conf$keyLogFsta[1L, ]; active <- key >= 0
F[, active] <- exp(t(s$pl$logF[key[active] + 1L, , drop = FALSE]))
Fse[, active] <- t(s$plsd$logF[key[active] + 1L, , drop = FALSE])
M <- s$data$natMor; Z <- F + M; W <- s$data$stockMeanWeight
at_age <- function(value, unit, se = NULL) {
  d <- transform(grid, est = as.vector(value), unit = unit, is_proj = FALSE)
  if (!is.null(se)) {
    d$se <- as.vector(se); d$se_scale <- "log"
    d$lwr <- d$est * exp(-qnorm(.975) * d$se); d$upr <- d$est * exp(qnorm(.975) * d$se)
    d$se_scale[d$est == 0] <- "reported"
  }
  d
}
annual <- function(value, unit) data.frame(year = years, est = value, unit = unit, is_proj = FALSE)
report <- summary(s$sdrep, "report")
pop <- list(N = at_age(N, "thousand fish", t(s$plsd$logN)),
  F = at_age(F, "per year", Fse), M = at_age(M, "per year"), Z = at_age(Z, "per year"),
  biomass_mat = at_age(N * W, "tonnes"),
  ssb_mat = at_age(N * W * s$data$propMat * exp(-as.vector(s$data$propF) * Z), "tonnes"),
  abundance = annual(rowSums(N), "thousand fish"))
for (metric in c("ssb", "biomass", "recruitment", "total_yield_pred")) {
  name <- c(ssb = "logssb", biomass = "logtsb", recruitment = "logR", total_yield_pred = "logCatch")[[metric]]
  a <- report[rownames(report) == name, , drop = FALSE]
  d <- annual(rep(NA_real_, length(years)), if (metric == "recruitment") "thousand fish" else "tonnes")
  d$se <- d$lwr <- d$upr <- NA_real_; d$se_scale <- "log"
  j <- seq_len(nrow(a)); d$est[j] <- exp(a[, 1L]); d$se[j] <- a[, 2L]
  d$lwr[j] <- exp(a[, 1L] - qnorm(.975) * a[, 2L]); d$upr[j] <- exp(a[, 1L] + qnorm(.975) * a[, 2L])
  pop[[metric]] <- d
}
for (component in c("F", "M")) {
  j <- match(f$dat[[paste0(component, "_settings")]]$mean_ages, ages)
  rate <- if (component == "F") F else M
  pop[[paste0(component, "_bar")]] <- annual(rowSums(rate[, j, drop = FALSE] * N[, j, drop = FALSE]) /
    rowSums(N[, j, drop = FALSE]), "per year")
}
obs <- f$dat$obs[c("catch", "index")]
aux <- as.data.frame(s$data$aux); names(aux)[1:3] <- c("year", "fleet", "age")
for (type in names(obs)) {
  d <- obs[[type]]
  fleet <- if (type == "catch") rep(1L, nrow(d)) else as.integer(sub("Survey ", "", d$survey)) + 1L
  i <- match(paste(d$year, fleet, d$age), paste(aux$year, aux$fleet, aux$age))
  d$pred <- exp(s$rep$predObs[i])
  if (type == "catch") d$pred[d$age < ages[which(active)[1L]]] <- 0
  d$sd <- exp(s$pl$logSdLogObs[s$conf$keyVarObs[cbind(fleet, match(d$age, ages))] + 1L])
  if (type == "catch" && nrow(s$conf$keyXtraSd)) {
    keys <- s$conf$keyXtraSd
    selected <- paste(d$year, d$age) %in% paste(keys[, 2L], keys[, 3L])
    d$sd[selected] <- d$sd[selected] * exp(s$pl$logXtraSd[1L])
  }
  if (type == "index") {
    qkey <- s$conf$keyLogFpar[cbind(fleet, match(d$age, ages))] + 1L
    d$q <- exp(s$pl$logFpar[qkey]); d$q_se <- s$plsd$logFpar[qkey]; d$q_se_scale <- "log"
    d$q_lwr <- d$q * exp(-qnorm(.975) * d$q_se); d$q_upr <- d$q * exp(qnorm(.975) * d$q_se)
    d$q_unit <- NA_character_
    weight <- s$data$weight[i]
    d$sd <- d$sd * ifelse(is.finite(weight), 1 / sqrt(weight), 1)
  }
  d$std_res <- ifelse(is.finite(d$obs) & d$obs > 0, (log(d$obs) - log(d$pred)) / d$sd, NA_real_)
  d$osa_res <- NA_real_; obs[[type]] <- d
}
ref <- structure(list(dat = list(obs = f$dat$obs, years = years, ages = ages), pop = pop,
  obs_pred = obs, fixed_par = f$fixed_par[FALSE, ], random_par = list()), class = c("tam_ref", "list"))
output <- sub("[.]rds$", "-comparison.html", args[1L])
vis_tam(model_list = list(tinyAM = f, SAM = ref), output_file = normalizePath(output, mustWork = FALSE),
  open_file = FALSE, render_args = list(quiet = TRUE), background = c(
    "### Pinned SAMjr fixture comparison",
    "The SAM reference contains native estimates, predictions and available 95% log-Wald intervals. F and M averages use the same abundance weighting and age ranges as tinyAM; intervals for these redefined averages and total abundance are unavailable. Biological inputs and F blocks match the source specification.",
    "Numbers are thousand fish, biomass/yield tonnes and average weights kg. Final survey-only years without native catch weights have missing SAM yield. Diagnostics and parameter pages describe the tinyAM fit; SAM reference state/observation tables are precomputed. No OSA residuals or new fits are calculated.",
    "See analysis/sam_equivalence/README.md for provenance, tolerances and limitations."))
saveRDS(ref, sub("[.]rds$", "-reference.rds", args[1L]))
cat(normalizePath(output), "\n")
