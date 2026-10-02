root <- "analysis/comp_assessments/database"
read <- function(name) {
  x <- read.csv(file.path(root, paste0(name, ".csv")), stringsAsFactors = FALSE)
  x[x$assessment_id == "ices_haddock_north_sea_2026", ]
}
inputs <- read("inputs")
outputs <- read("outputs")
stopifnot(nrow(inputs) == 5249, nrow(outputs) == 2353)
obs <- inputs[inputs$type == "index" & inputs$measure == "numbers_at_age", ]
sd <- inputs[inputs$measure == "log_index_sd", ]
key <- function(x) paste(x$survey, x$year, x$age)
stopifnot(nrow(obs) == 667, setequal(key(obs), key(sd)), all(sd$value > 0))
stopifnot(all(obs$sampling_time[obs$survey == "delta-GAMNS-WCQ1"] == .125))
stopifnot(all(obs$sampling_time[obs$survey == "delta-GAMNS-WCQ3+Q4"] == .75))
stopifnot(setequal(obs$age[obs$survey == "delta-GAMNS-WCQ1"], 1:8))
stopifnot(setequal(obs$age[obs$survey == "delta-GAMNS-WCQ3+Q4"], 0:8))
catch <- inputs[inputs$type == "catch" & inputs$measure == "numbers_at_age", ]
stopifnot(nrow(catch) == 486, max(catch$year) == 2025)
f <- outputs[outputs$measure == "fishing_mortality_at_age", ]
fbar <- outputs[outputs$measure == "Fbar", ]
mean_f <- aggregate(value ~ year, f[f$age %in% 2:4, ], mean)
joined <- merge(fbar, mean_f, by = "year")
stopifnot(nrow(joined) == 54, all(abs(joined$value.x - joined$value.y) < 1e-12))
n <- outputs[outputs$measure == "numbers_at_age" & outputs$type == "population", ]
r <- outputs[outputs$measure == "recruitment", ]
joined <- merge(n[n$age == 0, ], r, by = "year")
stopifnot(nrow(joined) == 55, all(abs(joined$value.x - joined$value.y) < 1e-7))
ssb <- outputs[outputs$measure == "SSB" & outputs$year == 2026, ]
stopifnot(abs(ssb$value - 699809.6) < .1)
stopifnot(abs(ssb$value - 667758) > 30000)
interval <- outputs[!is.na(outputs$lwr), ]
stopifnot(all(interval$lwr <= interval$value), all(interval$value <= interval$upr))
cat("Haddock survey coverage, SD pairing, age surfaces and forecast separation passed.\n")
