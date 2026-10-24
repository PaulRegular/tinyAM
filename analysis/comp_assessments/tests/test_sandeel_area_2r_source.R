source("analysis/comp_assessments/R/run_assessment.R")
db <- read_database()
s <- read_assessment("ices_sandeel_area_2r_2026", db)
i <- s$inputs; o <- s$outputs
stopifnot(s$stocks$authority_stock_id == "san.sa.2r", is.na(s$stocks$charbonneau_id),
  s$stocks$scientific_name == "Ammodytes spp.", s$stocks$ocean == "Northeast Atlantic",
  s$assessment$is_current, s$assessment$is_applied, s$assessment$terminal_year == 2025,
  all(unlist(s$assessment[c("assumptions_status", "inputs_status", "outputs_status")]) == "partial"),
  nrow(i) == 2192L, nrow(o) == 1032L, nrow(s$assumptions) == 17L,
  all(i$region == "Sandeel Area 2r"), all(o$region == "Sandeel Area 2r"),
  !any(i$year > 2025), !any(o$year > 2025))
c <- i[i$type == "catch" & i$measure == "numbers_at_age", ]
stopifnot(nrow(c) == 430L, all(c$unit == "million fish"),
  sum(c$value[c$age >= 1] > 0) == 321L,
  c$value[c$year == 1983 & c$age == 1 & c$season == "1"] == 4162,
  c$value[c$year == 2025 & c$age == 1 & c$season == "2"] == 43)
m <- i[i$type == "M", ]
stopifnot(nrow(m) == 387L, all(m$basis == "per_season"), all(m$unit == "per half-year"),
  all(lengths(lapply(split(m$value, interaction(m$age, m$season, drop = TRUE)), unique)) == 1L))
effort <- i[i$measure == "effort", ]
stopifnot(nrow(effort) == 127L,
  effort$value[effort$year == 1983 & effort$season == "annual"] == 4696,
  effort$value[effort$year == 1983 & effort$season == "1"] == 3127,
  effort$value[effort$year == 1983 & effort$season == "2"] == 1569,
  !any(effort$year %in% c(2007, 2022) & effort$season == "2"))
index <- i[i$type == "index", ]
stopifnot(nrow(index) == 44L, all(index$survey == "Dredge"),
  setequal(index$year, 2004:2025), setequal(index$age, 0:1),
  all(index$basis == "index_scale"), all(is.na(index$sampling_time)),
  index$value[index$year == 2004 & index$age == 1] == 339263682,
  index$value[index$year == 2025 & index$age == 0] == 420905820,
  any(s$assumptions$setting == "survey_start_discrepancy"))
n <- o[o$measure == "numbers_at_age", ]
stopifnot(nrow(n) == 215L, all(n$season[n$age == 0] == "2"),
  all(n$season[n$age >= 1] == "1"),
  n$value[n$year == 1983 & n$age == 1] == 13950072,
  n$value[n$year == 1983 & n$age == 4] == 43483.09,
  all(is.na(o$se)), !any(o$measure == "natural_mortality_at_age"),
  o$value[o$measure == "SSB" & o$year == 2025] == 25786,
  o$lwr[o$measure == "SSB" & o$year == 2025] == 17011,
  o$value[o$measure == "predicted_catch" & o$year == 2025] == 25354)
r <- o[o$measure == "recruitment", ]
stopifnot(abs(exp(mean(log(r$value))) - 22706561) < 1,
  abs(exp(mean(log(r$value[r$year %in% 2015:2024]))) - 12832432) < 1,
  any(s$assumptions$setting == "recruitment_averaging_windows"))

# The tabled Yield is consistent with seasonal model predictions, not observed C*W.
g <- expand.grid(year = 1983:2025, age = 1:4)
key <- function(z) paste(z$year, z$age)
surface <- function(rows, season) {
  rows <- rows[rows$season == season, ]
  rows$value[match(key(g), key(rows))]
}
f <- o[o$measure == "fishing_mortality_at_age", ]
f1 <- surface(f, "1"); f2 <- surface(f, "2")
m1 <- surface(m, "1"); m2 <- surface(m, "2")
w <- i[i$type == "catch_weight", ]
w1 <- surface(w, "1"); w2 <- surface(w, "2")
N <- surface(n, "1")
fraction <- function(f, m) f / (f + m) * (-expm1(-f - m))
lo1 <- (N - .5) * fraction(pmax(0, f1 - .0005), m1 + .005) * (w1 - .005)
lo2 <- (N - .5) * exp(-f1 - .0005 - m1 - .005) *
  fraction(pmax(0, f2 - .0005), m2 + .005) * (w2 - .005)
hi1 <- (N + .5) * fraction(f1 + .0005, m1 - .005) * (w1 + .005)
hi2 <- (N + .5) * exp(-pmax(0, f1 - .0005) - m1 + .005) *
  fraction(f2 + .0005, m2 - .005) * (w2 + .005)
lower <- tapply((lo1 + lo2) / 1000, g$year, sum)
upper <- tapply((hi1 + hi2) / 1000, g$year, sum)
yield <- o[o$measure == "predicted_catch", ]
stopifnot(all(yield$value + .5 >= lower[as.character(yield$year)]),
  all(yield$value - .5 <= upper[as.character(yield$year)]))
cat("Area 2r identity, source units/coverage, survey-start conflict, seasonal yield bounds and recruitment averaging windows passed.\n")
