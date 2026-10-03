source("analysis/comp_assessments/scripts/002_validate_database.R")
id <- "ices_herring_north_sea_2026"
x <- inputs[inputs$assessment_id == id, ]
y <- outputs[outputs$assessment_id == id, ]
a <- assessments[assessments$assessment_id == id, ]
stopifnot(nrow(x) == 3744L, nrow(y) == 1840L,
          a$inputs_status == "partial", a$outputs_status == "partial",
          a$assumptions_status == "partial", a$terminal_year == 2026)
expected <- c("HERAS", "IBTS0", "IBTS-Q1", "IBTS-Q3",
              "LAI-SNS", "LAI-CNS", "LAI-BUN", "LAI-ORSH")
stopifnot(setequal(unique(na.omit(x$survey)), expected))
lai <- x[x$measure == "larval_abundance_index", ]
stopifnot(all(is.na(lai$age)), all(lai$sampling_time == 0.67),
          all(lai$year >= 1973), all(!is.na(lai$season)),
          all(!is.na(lai$region)))
catch <- x[x$type == "catch", ]
stopifnot(!any(catch$year %in% 1978:1979), all(catch$unit == "thousand fish"),
          all(catch$year <= 2025), any(catch$value == 0))
n <- y[y$measure == "numbers_at_age", ]
r <- y[y$measure == "recruitment", ]
stopifnot(nrow(n) == 720L, nrow(r) == 80L,
          identical(n$value[n$age == 0], r$value))
f <- y[y$measure == "fishing_mortality_at_age", ]
fbar <- y[y$measure == "Fbar", ]
stopifnot(identical(f$value[f$age == 7], f$value[f$age == 8]))
means <- aggregate(value ~ year, f[f$age %in% 2:6, ], mean)
stopifnot(all(abs(means$value - fbar$value) < 0.0006),
          all(fbar$age_group == "2–6 winter rings"))
message("Herring coverage, native dimensions, exclusions, recruitment and Fbar checks passed.")
