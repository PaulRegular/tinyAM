cache <- commandArgs(TRUE)[1]
env <- new.env()
load(file.path(cache, "baserun_model.RData"), env)
fit <- env$fit
dat <- fit$data
fleet <- attr(dat, "fleetNames")
obs <- data.frame(dat$aux, value = exp(dat$logobs))
obs$fleet_name <- fleet[obs$fleet]
obs$sampling_time <- dat$sampleTimes[obs$fleet]
obs$prediction <- exp(fit[["rep"]]$predObs)
write.csv(obs, file.path(cache, "native_observations.csv"), row.names = FALSE)
for (name in c("logN", "logF")) {
  x <- fit$pl[[name]]
  ages <- seq.int(fit$conf$minAge, fit$conf$maxAge)
  if (name == "logF") x <- x[fit$conf$keyLogFsta[1, ] + 1, , drop = FALSE]
  grid <- expand.grid(age = ages, year = dat$years)
  grid$value <- c(exp(x))
  write.csv(grid, file.path(cache, paste0("native_", name, ".csv")), row.names = FALSE)
}
for (name in c("propMat", "stockMeanWeight", "catchMeanWeight", "natMor", "landFrac", "disMeanWeight", "landMeanWeight", "propF", "propM")) {
  x <- dat[[name]]
  grid <- expand.grid(dimnames(x), stringsAsFactors = FALSE)
  names(grid) <- c("year", "age", "fleet")[seq_len(ncol(grid))]
  grid$value <- c(x)
  write.csv(grid, file.path(cache, paste0("native_", name, ".csv")), row.names = FALSE)
}
write.csv(data.frame(year = dat$years, recruitment = exp(fit$pl$logN[1, ])), file.path(cache, "native_recruitment.csv"), row.names = FALSE)
v <- fit$sdrep$value
s <- fit$sdrep$sd
summary <- do.call(rbind, lapply(c("logssb", "logfbar", "logtsb", "logR"), function(name) {
  take <- which(names(v) == name)
  if (length(take) != length(dat$years)) return(NULL)
  data.frame(metric = name, year = dat$years, value = exp(v[take]), log_se = s[take], lwr = exp(v[take] - qnorm(.975) * s[take]), upr = exp(v[take] + qnorm(.975) * s[take]))
}))
write.csv(summary, file.path(cache, "native_summary.csv"), row.names = FALSE)
capture.output({
  cat("SAM version", attr(fit, "Version"), "revision", attr(fit, "RemoteSha"), "\n")
  print(aggregate(year ~ fleet_name + sampling_time, obs, range))
  print(aggregate(value ~ fleet_name, obs, length))
  print(tail(summary, 6))
  print(fit$opt[c("objective", "convergence", "message")])
  print(fit$sdrep$pdHess)
  print(fit$conf)
}, file = file.path(cache, "native_inventory.txt"))
print(summary[summary$year >= 2025, ])
