# Sourced by 001_compare.R after fitting. No new assessment mathematics.
results <- file.path(dir, "results")
state_matrix <- function(d) {
  m <- matrix(NA_real_, length(years), length(dat$ages), dimnames = list(year = years, age = dat$ages))
  m[cbind(match(d$year, years), match(d$age, dat$ages))] <- d$est
  m
}
input_matrix <- function(d, field) state_matrix(transform(d, est = d[[field]]))
W <- input_matrix(tam_obs$weight, "obs")
P <- input_matrix(tam_obs$maturity, "obs")
CW <- input_matrix(tam_obs$weight, "catch_weight")
pF <- input_matrix(tam_obs$weight, "propF")
pM <- input_matrix(tam_obs$weight, "propM")
common <- lapply(models, function(fit) {
  tables <- tidy_tam(fit)
  N <- state_matrix(tables$pop$N)
  F <- state_matrix(tables$pop$F)
  M <- state_matrix(tables$pop$M)
  C <- tables$obs_pred$catch
  C$est <- C$pred
  C <- state_matrix(C)
  data.frame(year = years,
    SSB_original_biology = rowSums(N * W * P * exp(-F * pF - M * pM)),
    biomass_original_weight = rowSums(N * W), abundance = rowSums(N),
    recruitment_age1 = N[, 1], Fbar_arithmetic = rowMeans(F[, as.character(settings$F_settings$mean_ages), drop = FALSE]),
    predicted_catch_original_weight = rowSums(C * CW))
})
common_long <- do.call(rbind, lapply(names(common), function(nm) {
  d <- common[[nm]]
  do.call(rbind, lapply(setdiff(names(d), "year"), function(metric)
    data.frame(model = nm, year = d$year, metric = metric, est = d[[metric]])))
}))
write.csv(common_long, file.path(results, "common_definition_trends.csv"), row.names = FALSE)
# Common transformations do not have calculated SEs: do not transfer native SEs.
pair <- function(d, keys, value = "est") {
  a <- d[d$model == "SAM", c(keys, value, intersect(c("lwr", "upr", "se", "se_scale"), names(d))), drop = FALSE]
  b <- d[d$model == "tinyAM", c(keys, value, intersect(c("lwr", "upr", "se", "se_scale"), names(d))), drop = FALSE]
  z <- merge(a, b, by = keys, suffixes = c("_SAM", "_tinyAM"))
  z$relative_difference <- ifelse(z[[paste0(value, "_SAM")]] != 0,
    z[[paste0(value, "_tinyAM")]] / z[[paste0(value, "_SAM")]] - 1, NA_real_)
  z
}
comparisons <- pair(common_long, c("year", "metric"))
comparisons$definition <- "common"
for (nm in c("ssb", "F_bar", "recruitment")) {
  z <- pair(tabs$pop[[nm]], "year")
  z$metric <- paste0("native_", nm)
  z$definition <- if (nm == "recruitment") "same age" else "different definitions"
  for (field in setdiff(names(z), names(comparisons))) comparisons[[field]] <- NA
  for (field in setdiff(names(comparisons), names(z))) z[[field]] <- NA
  comparisons <- rbind(comparisons, z[names(comparisons)])
}
write.csv(comparisons, file.path(results, "trend_comparisons.csv"), row.names = FALSE)
metrics <- function(d, value = "est") {
  a <- d[[paste0(value, "_SAM")]]; b <- d[[paste0(value, "_tinyAM")]]
  i <- is.finite(a) & is.finite(b)
  a <- a[i]; b <- b[i]; y <- d$year[i]
  r <- ifelse(a != 0, b / a - 1, NA_real_)
  terminal <- if (length(y)) which.max(y) else integer()
  data.frame(n = length(a), mean_relative_difference = if (any(is.finite(r))) mean(r, na.rm = TRUE) else NA_real_,
    mean_absolute_relative_difference = if (any(is.finite(r))) mean(abs(r), na.rm = TRUE) else NA_real_,
    trend_correlation = if (length(a) > 1 && sd(a) > 0 && sd(b) > 0) cor(a, b) else NA_real_,
    terminal_year = if (length(y)) max(y) else NA_integer_,
    terminal_relative_difference = if (length(terminal)) r[terminal] else NA_real_)
}
summary <- do.call(rbind, lapply(split(comparisons, comparisons$metric), function(d)
  cbind(metric = d$metric[1], definition = d$definition[1], metrics(d))))
write.csv(summary, file.path(results, "trend_agreement.csv"), row.names = FALSE)
for (nm in c("N", "F")) {
  z <- pair(tabs$pop[[nm]], c("year", "age"))
  write.csv(z, file.path(results, paste0(nm, "_at_age_comparison.csv")), row.names = FALSE)
  age_summary <- do.call(rbind, lapply(split(z, z$age), function(d) cbind(age = d$age[1], metrics(d))))
  write.csv(age_summary, file.path(results, paste0(nm, "_at_age_agreement.csv")), row.names = FALSE)
}
for (nm in c("catch", "index")) {
  d <- tabs$obs_pred[[nm]]
  keys <- c("year", "age", if (nm == "index") "survey")
  z <- pair(d, keys, "pred")
  write.csv(z, file.path(results, paste0(nm, "_prediction_comparison.csv")), row.names = FALSE)
  if (nm == "index") {
    z <- pair(d, keys, "q")
    write.csv(z, file.path(results, "q_comparison.csv"), row.names = FALSE)
  }
}
# A second interactive view isolates reporting definitions from state differences.
plots <- lapply(unique(common_long$metric), function(nm) {
  d <- common_long[common_long$metric == nm, ]
  plotly::plot_ly(d, x = ~year, y = ~est, color = ~model, type = "scatter", mode = "lines") |>
    plotly::layout(title = nm, xaxis = list(title = "Year"), yaxis = list(title = "Original source units"))
})
notes <- htmltools::tags$p("Common period: 1983-2022. SSB uses original maturity, weight and spawning timing for both models. Fbar is arithmetic over ages 2-4. Catch biomass uses original catch weights and is unavailable for 2022. No uncertainty is invented for these recomputed summaries. Native outputs and their available intervals remain in the main dashboard.")
htmltools::save_html(htmltools::tagList(htmltools::tags$h1("SAM and tinyAM: common definitions"), notes, plots),
  file = file.path(results, "common_definitions.html"), libdir = "common_definitions_files")
print(summary, row.names = FALSE)
