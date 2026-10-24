pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
db <- read_database()
source <- read_assessment("ices_sandeel_area_1r_2026", db)
stopifnot(source$stocks$stock_id == "ices_sandeel_area_1r",
  is.na(source$stocks$charbonneau_id), source$stocks$authority_stock_id == "san.sa.1r",
  source$stocks$scientific_name == "Ammodytes spp.", source$stocks$ocean == "Northeast Atlantic",
  source$assessment$is_current, source$assessment$is_applied,
  source$assessment$terminal_year == 2025, source$assessment$estimate_terminal_year == 2026,
  all(unlist(source$assessment[c("assumptions_status", "inputs_status", "outputs_status")]) == "partial"))
parent <- db$stocks[db$stocks$stock_id == "ck_ices_hawg_dubius_ns", ]
stopifnot(parent$charbonneau_id == "ICES-HAWG_dubius_NS",
  grepl("one-to-one", parent$notes), !any(db$assessments$stock_id == parent$stock_id))
x <- source$inputs; o <- source$outputs
stopifnot(nrow(x) == 2213, nrow(o) == 1032, nrow(source$assumptions) == 17,
  all(x$region == "Sandeel Area 1r"), all(o$region == "Sandeel Area 1r"))
C <- x[x$type == "catch" & x$measure == "numbers_at_age", ]
stopifnot(nrow(C) == 430, all(C$unit == "million fish"),
  all(C$year_basis == "calendar_year"),
  C$value[C$year == 1983 & C$age == 2 & C$season == "1"] == 28971,
  C$value[C$year == 2024 & C$age == 0 & C$season == "1"] == 26025,
  C$value[C$year == 2025 & C$age == 4 & C$season == "2"] == 6,
  sum(C$value > 0 & C$age >= 1) == 320,
  any(grepl("318", source$assumptions$notes)), any(C$value == 0))
M <- x[x$type == "M", ]
stopifnot(nrow(M) == 387, all(M$basis == "per_season"), all(M$unit == "per half-year"),
  !any(M$age == 0 & M$season == "1"),
  M$value[M$year == 1983 & M$age == 1 & M$season == "1"] == .38,
  M$value[M$year == 1983 & M$age == 1 & M$season == "2"] == .42,
  all(x$value[x$type == "weight"] == x$value[x$type == "catch_weight"]),
  all(x$unit[x$type %in% c("weight", "catch_weight")] == "g per fish"))
E <- x[x$measure == "effort", ]
stopifnot(nrow(E) == 127, !anyNA(E$value),
  !any(E$year %in% c(2007L, 2022L) & E$season == "2"),
  sum(E$season == "annual") == 43, sum(E$season == "1") == 43,
  E$value[E$year == 1983 & E$season == "2"] == 2065)
I <- x[x$type == "index", ]
stopifnot(nrow(I) == 65, all(I$basis == "index_scale"),
  all(I$sampling_time[I$survey == "RTM"] == 0), all(is.na(I$sampling_time[I$survey == "Dredge"])),
  !any(I$survey == "RTM" & I$year %in% c(2021L, 2022L, 2024L)),
  I$value[I$year == 2025 & I$survey == "Dredge" & I$age == 0] == 272283759)
N <- o[o$measure == "numbers_at_age", ]
F <- o[o$measure == "fishing_mortality_at_age" & o$season == "annual", ]
stopifnot(nrow(N) == 215, nrow(F) == 215, all(N$season[N$age == 0] == "2"),
  all(N$season[N$age >= 1] == "1"),
  N$value[N$year == 1983 & N$age == 1] == 13333913,
  F$value[F$year == 1983 & F$age == 1] == .4756039,
  all(o$unit[o$measure == "fishing_mortality_at_age" & o$season != "annual"] == "per half-year"),
  all(is.na(o$se)), !any(o$measure == "natural_mortality_at_age"), !any(o$year > 2025))
total <- o[o$measure == "total_abundance", ]
adult <- N[N$age >= 1, ]
summed <- aggregate(adult$value, list(year = adult$year), sum)
stopifnot(all(total$source_type == "derived_source_output"),
  identical(total$value, summed$x[match(total$year, summed$year)]),
  o$value[o$measure == "SSB" & o$year == 2025] == 142799,
  o$lwr[o$measure == "SSB" & o$year == 2025] == 95198,
  o$upr[o$measure == "SSB" & o$year == 2025] == 214202,
  o$value[o$measure == "recruitment" & o$year == 2025] == 66060391,
  o$value[o$measure == "predicted_catch" & o$year == 2025] == 29044,
  x$value[x$measure == "total_biomass" & x$year == 2025] == 67675.175772)
cat("Sandeel Area 1r identity, seasonal hazards, observation roles, mixed state dates, source gaps and uncertainty passed.\n")

validate_seasonal_fixture <- function(tables) {
  env <- new.env(parent = globalenv())
  env$read.csv <- function(file, ...) tables[[sub(".csv", "", basename(file), fixed = TRUE)]]
  tryCatch({
    invisible(capture.output(sys.source(
      "analysis/comp_assessments/scripts/database/002_validate_database.R", envir = env)))
    NULL
  }, error = identity)
}
at <- which(db$inputs$assessment_id == source$assessment$assessment_id & db$inputs$type == "M")[1L]
for (field in c("type", "season")) {
  wrong <- db
  wrong$inputs[[field]][at] <- if (field == "type") "weight" else NA_character_
  error <- validate_seasonal_fixture(wrong)
  stopifnot(inherits(error, "error"), grepl("M-at-age and a source season", conditionMessage(error), fixed = TRUE))
}
cat("Seasonal mortality validation rejects other input types and missing source seasons.\n")
