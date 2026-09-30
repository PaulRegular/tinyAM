# Optional integration validation with the installed SAM package.
# Run from the repository root after 001_download.R and 002_nscod.R.
library(tinyAM)
if (!requireNamespace("stockassessment", quietly = TRUE)) {
  cli::cli_abort("Install stockassessment to run this optional SAM validation.")
}
here <- file.path("analysis", "comp_sam")
out <- file.path(here, "results")
bridge <- readRDS(file.path(out, "nscod_bridge.rds"))
read_input <- function(name) stockassessment::read.ices(file.path(here, "source", paste0(name, ".dat")))
equal <- function(a, b) stopifnot(isTRUE(all.equal(a, b, check.attributes = FALSE)))
for (nm in c("sw", "cw", "dw", "lw", "lf", "mo", "nm", "pf", "pm")) {
  equal(bridge$source$data[[nm]], read_input(nm))
}
equal(bridge$source$data$catch[[1]], read_input("cn"))
surveys <- read_input("survey")
for (i in seq_along(surveys)) equal(bridge$source$data$surveys[[i]], surveys[[i]])

# Explicit inputs mirror the published script, rather than relying on defaults
# or the convenience directory reader's filename conventions.
dat <- stockassessment::setup.sam.data(surveys = surveys, residual.fleets = list(read_input("cn")),
  prop.mature = read_input("mo"), stock.mean.weight = read_input("sw"),
  catch.mean.weight = read_input("cw"), dis.mean.weight = read_input("dw"),
  land.mean.weight = read_input("lw"), natural.mortality = read_input("nm"),
  prop.f = read_input("pf"), prop.m = read_input("pm"), land.frac = read_input("lf"))
conf <- stockassessment::loadConf(dat, file.path(here, "nscod.cfg"), patch = FALSE)
for (nm in names(conf)) equal(conf[[nm]], bridge$source$conf[[nm]])
equal(dat$sampleTimes, bridge$source$fleets$samp_time)

# Check ordinary saved values against SAM's own table methods. qtable returns
# log catchability; the bridge intentionally reports q on the natural scale.
check_reference <- function(fit, reference) {
  equal(reference$tables$N$obs, as.vector(stockassessment::ntable(fit)))
  equal(reference$tables$F$obs, as.vector(stockassessment::faytable(fit)))
  for (nm in c("SSB", "recruitment", "Fbar")) {
    method <- switch(nm, SSB = stockassessment::ssbtable,
                     recruitment = stockassessment::rectable, Fbar = stockassessment::fbartable)
    equal(reference$tables[[nm]]$est, method(fit)[, "Estimate"])
  }
  qt <- stockassessment::qtable(fit)
  q <- reference$tables$q
  equal(q$est, exp(qt[cbind(match(q$fleet_id, which(fit$data$fleetTypes != 0)),
                          q$age - fit$conf$minAge + 1L)]))
}
saved <- new.env(parent = baseenv())
load(file.path(here, "source", "fit.expected.Rdata"), saved)
check_reference(saved$fit.exp, bridge$reference)

# A fresh baseline uses the same inputs and configuration, with SAM's default
# optimizer and Newton steps. Keep it separate from the historical saved fit.
elapsed <- system.time(fit <- stockassessment::sam.fit(dat, conf,
  stockassessment::defpar(dat, conf), silent = TRUE))["elapsed"]
desc <- utils::packageDescription("stockassessment")
provenance <- paste("stockassessment", desc$Version, "source", desc$RemoteSha)
reference <- sam_reference(fit, provenance)
check_reference(fit, reference)
saveRDS(fit, file.path(out, "SAM_current_fit.rds"))
for (nm in names(reference$tables)) {
  utils::write.csv(reference$tables[[nm]], file.path(out, paste0("SAM_current_", nm, ".csv")), row.names = FALSE)
}
diagnostics <- data.frame(package_version = desc$Version, source_commit = desc$RemoteSha,
  convergence = fit$opt$convergence, message = fit$opt$message,
  objective = fit$opt$objective, max_abs_gradient = max(abs(fit$obj$gr(fit$opt$par))),
  sdreport_success = all(is.finite(fit$sdrep$sd)), positive_definite_hessian = isTRUE(fit$sdrep$pdHess),
  elapsed_seconds = unname(elapsed), historical_objective = saved$fit.exp$opt$objective,
  max_abs_SSB_difference = max(abs(reference$tables$SSB$est - bridge$reference$tables$SSB$est)))
utils::write.csv(diagnostics, file.path(out, "SAM_current_diagnostics.csv"), row.names = FALSE)
print(diagnostics)
cat("Installed SAM reader, configuration and reference-table checks passed.\n")
