root <- "analysis/comp_assessments/database"
read <- function(name) {
  x <- read.csv(file.path(root, paste0(name, ".csv")), stringsAsFactors = FALSE)
  x[x$assessment_id == "ices_cod_northeast_arctic_2026", ]
}
inputs <- read("inputs")
outputs <- read("outputs")
assessment <- read("assessments")
stopifnot(nrow(inputs) == 6576, nrow(outputs) == 4793)
stopifnot(assessment$terminal_year == 2026, assessment$estimate_terminal_year == 2026)
index <- inputs[inputs$type == "index", ]
timing <- c("FLT15_I:NorBarTrSur_I" = .137, "FLT15_II:NorBarTrSur_II" = .137,
            "FLT16:NorBarLofAcSur" = .1725, "FLT18:RusSweptArea" = .95,
            "FLT007:Ecosystem" = .7)
stopifnot(setequal(unique(index$survey), names(timing)))
stopifnot(all(index$sampling_time == timing[index$survey]))
stopifnot(max(index$year[index$survey == "FLT15_I:NorBarTrSur_I"]) == 2013)
stopifnot(max(index$year[index$survey == "FLT18:RusSweptArea"]) == 2017)
catch <- inputs[inputs$type == "catch", ]
stopifnot(nrow(catch) == 1028, max(catch$year) == 2025)
stopifnot(!any(catch$year == 2011 & catch$age == 15))
f <- outputs[outputs$measure == "fishing_mortality_at_age", ]
stopifnot(max(f$year) == 2025)
f14 <- f[f$age == 14, ]; f15 <- f[f$age == 15, ]
stopifnot(identical(f14$value, f15$value))
fbar <- outputs[outputs$measure == "Fbar", ]
mean_f <- aggregate(value ~ year, f[f$age %in% 5:10, ], mean)
joined <- merge(fbar, mean_f, by = "year")
stopifnot(all(abs(joined$value.x - joined$value.y) < 1e-12))
n <- outputs[outputs$measure == "numbers_at_age" & outputs$type == "population", ]
r <- outputs[outputs$measure == "recruitment", ]
joined <- merge(n[n$age == 3, ], r, by = "year")
stopifnot(all(abs(joined$value.x - joined$value.y) < 1e-7))
stopifnot(r$value[r$year == 2026] < 140000)
interval <- outputs[!is.na(outputs$lwr), ]
stopifnot(all(interval$lwr <= interval$value), all(interval$value <= interval$upr))
cat("Northeast Arctic cod inventory, timing, exclusions, state mapping and reporting checks passed.\n")

