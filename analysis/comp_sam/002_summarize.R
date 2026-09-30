# Annual differences ----
pop <- tidy_tam(model_list = models)$pop
metrics <- pop[c("ssb", "recruitment", "N", "F", "M")]
metrics$Fbar <- aggregate(est ~ model + year,
  data = pop$F[pop$F$age %in% settings$F_settings$mean_ages, ], FUN = mean)

percent_differences <- do.call(rbind, lapply(names(metrics), function(metric) {
  d <- metrics[[metric]]
  if (!"age" %in% names(d)) d$age <- NA_integer_
  sam <- d[d$model == "SAM", c("year", "age", "est")]
  tam <- d[d$model == "tinyAM", c("year", "age", "est")]
  names(sam)[3] <- "SAM"
  names(tam)[3] <- "tinyAM"
  d <- merge(sam, tam, by = c("year", "age"))
  d$metric <- metric
  d$percent_difference <- ifelse(d$SAM == 0, NA_real_, 100 * (d$tinyAM / d$SAM - 1))
  d[c("metric", "year", "age", "SAM", "tinyAM", "percent_difference")]
}))
write.csv(percent_differences, file.path(results, "percent_differences.csv"), row.names = FALSE)

# Summary by metric and age ----
mean_difference <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
groups <- split(percent_differences, paste(percent_differences$metric, percent_differences$age))
summary <- do.call(rbind, lapply(groups, function(d) {
  terminal <- d[which.max(d$year), ]
  data.frame(metric = d$metric[1], age = d$age[1],
    n_years = sum(!is.na(d$percent_difference)),
    mean_percent_difference = mean_difference(d$percent_difference),
    mean_absolute_percent_difference = mean_difference(abs(d$percent_difference)),
    terminal_year = terminal$year, terminal_percent_difference = terminal$percent_difference)
}))
write.csv(summary, file.path(results, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
