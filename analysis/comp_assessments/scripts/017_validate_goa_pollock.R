source("analysis/comp_assessments/scripts/002_validate_database.R")
id <- "afsc_pollock_goa_2024"
x <- inputs[inputs$assessment_id == id, ]
y <- outputs[outputs$assessment_id == id, ]
a <- assessments[assessments$assessment_id == id, ]
stopifnot(nrow(x) == 5134, nrow(y) == 660, a$model_version == "23d",
          a$inputs_status == "partial", a$outputs_status == "partial")
index <- x[x$type == "index", ]
totals <- index[index$measure == "total_biomass", ]
stopifnot(length(unique(totals$survey)) == 4, nrow(totals) == 90,
          !any(grepl("age-1|age-2", index$survey)),
          all(is.na(index$age[index$measure == "log_index_sd"])))
for (label in unique(totals$survey)) {
  z <- totals[totals$survey == label, ]
  s <- index[index$survey == label & index$measure == "log_index_sd", ]
  stopifnot(identical(z$year, s$year), identical(z$sampling_time, s$sampling_time))
}
stopifnot(all(totals$sampling_time[totals$survey == "Shelikof winter acoustic"] == .209),
          all(totals$sampling_time[totals$survey == "Summer acoustic"] == .519),
          all(totals$sampling_time[totals$survey == "ADF&G crab/groundfish trawl"] == .60989))
comp <- x[x$measure == "proportion_at_age", ]
fsh <- comp[comp$type == "catch", ]
shelikof <- comp[!is.na(comp$survey) & comp$survey == "Shelikof winter acoustic", ]
stopifnot(nrow(fsh) == 441, min(fsh$age) == 2, max(fsh$year) == 2023,
          nrow(shelikof) == 248, min(shelikof$age) == 3,
          all(fsh$source_type == "reconstructed_source_input"),
          all(shelikof$source_type == "reconstructed_source_input"))
stopifnot(nrow(x[x$type == "M", ]) == 10,
          nrow(x[x$type == "maturity", ]) == 10,
          all(x$year[x$type %in% c("M", "maturity")] == 1970),
          nrow(x[x$measure == "spawning_weight_at_age", ]) == 550)
n <- y[y$measure == "numbers_at_age", ]
r <- y[y$measure == "recruitment", ]
ssb <- y[y$measure == "SSB", ]
stopifnot(nrow(n) == 550, nrow(r) == 55, nrow(ssb) == 55,
          identical(n$value[n$age == 1], r$value),
          all(n$age_group[n$age == 10] == "10+"),
          all(r$lwr <= r$value & r$value <= r$upr),
          all(ssb$lwr <= ssb$value & ssb$value <= ssb$upr), all(is.na(y$se)))
message("GOA pollock survey coverage, timing, grouped compositions, biology and uncertainty checks passed.")
