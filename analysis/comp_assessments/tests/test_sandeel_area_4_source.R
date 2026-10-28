source("analysis/comp_assessments/R/run_assessment.R")
db <- read_database()
s <- read_assessment("ices_sandeel_area_4_2025", db)
i <- s$inputs; o <- s$outputs
stopifnot(s$stocks$authority_stock_id == "san.sa.4", is.na(s$stocks$charbonneau_id),
  s$stocks$scientific_name == "Ammodytes spp.", s$stocks$ocean == "Northeast Atlantic",
  s$assessment$assessment_year == 2025, s$assessment$terminal_year == 2024,
  s$assessment$is_current, s$assessment$is_applied, s$assessment$framework_year == 2023,
  grepl("28389008.v5", s$assessment$assessment_url, fixed = TRUE),
  all(unlist(s$assessment[c("assumptions_status", "inputs_status", "outputs_status")]) == "partial"),
  nrow(i) == 1642L, nrow(o) == 768L, nrow(s$assumptions) == 20L,
  all(i$region == "Sandeel Area 4"), all(o$region == "Sandeel Area 4"),
  setequal(i$year, 1993:2024), !any(o$year > 2024))
c <- i[i$type == "catch" & i$measure == "numbers_at_age", ]
stopifnot(nrow(c) == 320L, all(c$unit == "million fish"),
  sum(c$value[c$age >= 1] > 0) == 132L,
  all(c$value[c$year == 2024] == 0))
m <- i[i$type == "M", ]
stopifnot(nrow(m) == 288L, all(m$basis == "per_season"),
  all(m$unit == "per half-year"),
  all(lengths(lapply(split(m$value, interaction(m$age, m$season, drop = TRUE)), unique)) == 1L),
  all(m$value[m$age == 0] == .99), all(m$season[m$age == 0] == "2"))
effort <- i[i$measure == "effort", ]
stopifnot(nrow(effort) == 94L,
  effort$value[effort$year == 1996 & effort$season == "annual"] == 1501,
  effort$value[effort$year == 1996 & effort$season == "1"] == 535,
  effort$value[effort$year == 1996 & effort$season == "2"] == 966,
  all(effort$source_type[effort$season == "2"] == "reconstructed_source_input"),
  all(grepl("Area-5r", effort$transformation[effort$season == "2"], fixed = TRUE)),
  !any(effort$year %in% c(2007, 2022) & effort$season == "2"))
e <- reshape(effort[c("year", "season", "value")], idvar = "year", timevar = "season", direction = "wide")
stopifnot(all(abs(e$value.annual - e$value.1 - e$value.2) <= 1, na.rm = TRUE))
index <- i[i$type == "index", ]
stopifnot(nrow(index) == 44L, sum(!is.na(index$value)) == 43L,
  all(index$basis == "index_scale"), all(is.na(index$sampling_time)),
  setequal(index$year[index$survey == "Dredge"], 1999:2003),
  setequal(index$year[index$survey == "Dredge2"], 2008:2024),
  is.na(index$value[index$survey == "Dredge" & index$year == 2003 & index$age == 1]),
  index$value[index$survey == "Dredge2" & index$year == 2024 & index$age == 0] == 58591395,
  any(grepl("boundary estimates, not fixed", s$assumptions$notes, fixed = TRUE)),
  any(grepl("profiles catch log-SDs", s$assumptions$notes, fixed = TRUE)))
n <- o[o$measure == "numbers_at_age", ]
f <- o[o$measure == "fishing_mortality_at_age", ]
stopifnot(nrow(n) == 160L, all(n$season[n$age == 0] == "2"),
  all(n$season[n$age >= 1] == "1"), all(is.na(o$se)),
  all(f$value[f$year %in% c(2009, 2024)] == 0),
  n$value[n$year == 2024 & n$age == 1] == 1341151,
  o$value[o$measure == "SSB" & o$year == 2024] == 102311,
  !any(o$measure == "natural_mortality_at_age"))
r <- o[o$measure == "recruitment", ]
stopifnot(abs(exp(mean(log(r$value))) - 49075026) < 1,
  abs(exp(mean(log(r$value[r$year %in% 2014:2023]))) - 53694268) < 1)

# Native Catchtot is conditional on positive observed catch cells.
g <- expand.grid(year = 1993:2024, age = 1:4)
key <- function(z) paste(z$year, z$age)
surface <- function(rows, season) {
  rows <- rows[rows$season == season, ]
  rows$value[match(key(g), key(rows))]
}
f1 <- surface(f, "1"); f2 <- surface(f, "2")
m1 <- surface(m, "1"); m2 <- surface(m, "2")
w <- i[i$type == "catch_weight", ]
w1 <- surface(w, "1"); w2 <- surface(w, "2")
N <- surface(n, "1"); c1 <- surface(c, "1") > 0; c2 <- surface(c, "2") > 0
fraction <- function(f, m) f / (f + m) * (-expm1(-f - m))
lower <- tapply((N - .5) * (
  fraction(pmax(0, f1 - .0005), m1 + .005) * (w1 - .005) * c1 +
  exp(-f1 - .0005 - m1 - .005) * fraction(pmax(0, f2 - .0005), m2 + .005) * (w2 - .005) * c2) / 1000, g$year, sum)
upper <- tapply((N + .5) * (
  fraction(f1 + .0005, m1 - .005) * (w1 + .005) * c1 +
  exp(-pmax(0, f1 - .0005) - m1 + .005) * fraction(f2 + .0005, m2 - .005) * (w2 + .005) * c2) / 1000, g$year, sum)
yield <- o[o$measure == "predicted_catch", ]
stopifnot(all(yield$value + .5 >= lower[as.character(yield$year)]),
  all(yield$value - .5 <= upper[as.character(yield$year)]),
  yield$value[yield$year == 1996] == 41518,
  all(yield$value[yield$year %in% c(2006:2011, 2013, 2024)] == 0))
later <- read_assessment("ices_sandeel_area_4_2026_summary", db)
stopifnot(!later$assessment$is_current, later$assessment$is_applied,
  later$assessment$assessment_year == 2026, nrow(later$inputs) == 69L,
  nrow(later$outputs) == 0L,
  later$assessment$inputs_status == "partial",
  later$assessment$outputs_status == "not_applicable")
cat("Area 4 identity, historical arrays, missing cell, effort reassignment, conditional yield and separate summary vintage passed.\n")
