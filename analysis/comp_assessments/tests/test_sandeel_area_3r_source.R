source("analysis/comp_assessments/R/run_assessment.R")
db <- read_database()
s <- read_assessment("ices_sandeel_area_3r_2026", db)
i <- s$inputs; o <- s$outputs
stopifnot(s$stocks$authority_stock_id == "san.sa.3r", is.na(s$stocks$charbonneau_id),
  s$stocks$scientific_name == "Ammodytes spp.", s$stocks$ocean == "Northeast Atlantic",
  s$assessment$is_current, s$assessment$is_applied, s$assessment$terminal_year == 2025,
  s$assessment$framework_year == 2023,
  grepl("^Declared smsR", s$assessment$model_version),
  all(unlist(s$assessment[c("assumptions_status", "inputs_status", "outputs_status")]) == "partial"),
  nrow(i) == 2056L, nrow(o) == 936L, nrow(s$assumptions) == 18L,
  all(i$region == "Sandeel Area 3r"), all(o$region == "Sandeel Area 3r"),
  !any(i$year > 2025), !any(o$year > 2025))
c <- i[i$type == "catch" & i$measure == "numbers_at_age", ]
stopifnot(nrow(c) == 390L, all(c$unit == "million fish"),
  sum(c$value[c$age >= 1] > 0) == 240L,
  c$value[c$year == 1987 & c$age == 1 & c$season == "1"] == 33760,
  c$value[c$year == 2025 & c$age == 1 & c$season == "2"] == 0)
m <- i[i$type == "M", ]
stopifnot(nrow(m) == 351L, all(m$basis == "per_season"), all(m$unit == "per half-year"),
  m$value[m$year == 1987 & m$age == 1 & m$season == "1"] == .61,
  m$value[m$year == 2025 & m$age == 1 & m$season == "1"] == .70,
  all(lengths(lapply(split(m$value, interaction(m$age, m$season, drop = TRUE)), unique)) > 1L))

effort <- i[i$measure == "effort", ]
stopifnot(nrow(effort) == 115L, all(effort$value[effort$year == 2025] == 0),
  all(effort$value[effort$year < 2025 & effort$season == "annual"] > 0),
  !any(effort$year %in% c(2007, 2022) & effort$season == "2"))
index <- i[i$type == "index", ]
d <- index[index$survey == "Dredge", ]; a <- index[index$survey == "Acoustic", ]
stopifnot(nrow(index) == 108L, nrow(d) == 40L, nrow(a) == 68L,
  setequal(d$year, 2006:2025), setequal(d$age, 0:1),
  setequal(a$year, 2009:2025), setequal(a$age, 1:4),
  !anyDuplicated(paste(index$survey, index$year, index$age)),
  all(index$basis == "index_scale"), all(is.na(index$sampling_time)),
  a$value[a$year == 2025 & a$age == 3] == 576.3118,
  a$value[a$year == 2025 & a$age == 4] == 1721.7144,
  any(s$assumptions$setting == "known_zero_year"))
n <- o[o$measure == "numbers_at_age", ]
stopifnot(nrow(n) == 195L, all(n$season[n$age == 0] == "2"),
  all(n$season[n$age >= 1] == "1"),
  n$value[n$year == 1987 & n$age == 1] == 88317286,
  n$value[n$year == 1987 & n$age == 4] == 211312,
  all(o$value[o$measure == "fishing_mortality_at_age" & o$year == 2025] == 0),
  all(c$value[c$year == 2025] == 0), all(is.na(o$se)), !any(o$measure == "natural_mortality_at_age"),
  o$value[o$measure == "SSB" & o$year == 2025] == 91332,
  o$value[o$measure == "predicted_catch" & o$year == 2025] == 0)
r <- o[o$measure == "recruitment", ]
stopifnot(abs(exp(mean(log(r$value))) - 75267201) < 1,
  abs(exp(mean(log(r$value[r$year %in% 1987:2024]))) - 79171378) < 1,
  any(s$assumptions$setting == "recruitment_averaging_windows"))
stopifnot(grepl("surveySD/catchSD/estSD", s$assumptions$notes[
  s$assumptions$setting == "remaining_sources"], fixed = TRUE))

# The tabled Yield is consistent with seasonal model predictions, not observed C*W.
g <- expand.grid(year = 1987:2025, age = 1:4)
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
cat("Area 3r identity, seasonal arrays, time-varying supplied M, survey coverage, exact 2025 zeros and yield/forecast checks passed.\n")
