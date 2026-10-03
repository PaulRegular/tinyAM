root <- "analysis/comp_assessments/database"
read <- function(name) {
  x <- read.csv(file.path(root, paste0(name, ".csv")), stringsAsFactors = FALSE)
  x[x$assessment_id == "ices_haddock_north_sea_2026", ]
}
inputs <- read("inputs")
outputs <- read("outputs")
stopifnot(nrow(inputs) == 6239, nrow(outputs) == 2373)
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

for (measure in c("fraction_F_before_spawning", "fraction_M_before_spawning")) {
  x <- inputs[inputs$measure == measure, ]
  stopifnot(nrow(x) == 495, all(x$value == 0),
            setequal(x$year, 1972:2026), setequal(x$age, 0:8))
  native <- read.csv(file.path("analysis/comp_assessments/source_cache/ices_haddock_north_sea_2026",
    if (measure == "fraction_F_before_spawning") "native_propF.csv" else "native_propM.csv"))
  joined <- merge(x, native, by = c("year", "age"))
  stopifnot(nrow(joined) == 495, all(joined$value.x == joined$value.y))
}
cat("Both native spawning-fraction matrices verified.\n")

q <- outputs[outputs$type == "catchability", ]
stopifnot(nrow(q) == 20, all(is.na(q$year)), all(q$se > 0),
          all(q$lwr < q$value), all(q$upr > q$value))
power <- q[q$measure == "q_power", ]
stopifnot(nrow(power) == 3,
          setequal(paste(power$survey, power$age), c("delta-GAMNS-WCQ1 1", "delta-GAMNS-WCQ3+Q4 0", "delta-GAMNS-WCQ3+Q4 1")))
native_q <- read.csv("analysis/comp_assessments/source_cache/ices_haddock_north_sea_2026/native_catchability.csv")
joined <- merge(q, native_q, by = c("measure", "survey", "age"))
stopifnot(nrow(joined) == 20, all(abs(joined$value.x-exp(joined$log_estimate)) < 1e-10),
          all(abs(joined$se.x-joined$value.x*joined$log_se) < 1e-10))
cat("Catchability power mapping and uncertainty scales verified.\n")

model_env <- new.env()
load("analysis/comp_assessments/source_cache/ices_haddock_north_sea_2026/run_model.RData", model_env)
stopifnot(model_env$fit$conf$initState == 0,
          length(model_env$fit$pl$initN) == 0, length(model_env$fit$pl$initF) == 0,
          all(c("logN", "logF") %in% names(model_env$fit$sdrep$par.random)))
cat("Native initial-state configuration verified.\n")
