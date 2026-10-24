pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
db <- read_database()
framework <- read_assessment("nefsc_yellowtail_snema_2024_framework", db)
update <- read_assessment("nefsc_yellowtail_snema_2025", db)
stopifnot(framework$stocks$charbonneau_id == "NEFSC-GARMIII_ferruginea_SNE, MA",
  framework$stocks$ocean == "Northwest Atlantic",
  framework$assessment$is_current, !framework$assessment$is_applied,
  framework$assessment$terminal_year == 2022,
  !update$assessment$is_current, update$assessment$is_applied,
  update$assessment$terminal_year == 2024, update$assessment$estimate_terminal_year == 2028,
  all(unlist(framework$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  all(unlist(update$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"))
x <- framework$inputs
catch <- x[x$type == "catch", ]
stopifnot(nrow(catch) == 300, all(catch$unit == "thousand fish"),
  identical(sort(unique(catch$year)), 1973:2022),
  identical(sort(unique(catch$age)), 1:6),
  sum(is.na(catch$value)) == 18,
  all(is.na(catch$value[catch$year >= 2020])),
  catch$value[catch$year == 1975 & catch$age == 1] == 8037,
  catch$value[catch$year == 1990 & catch$age == 3] == 40271,
  catch$value[catch$year == 2019 & catch$age == 6] == 16,
  all(grepl("6\\+", catch$notes)))
indices <- x[x$type == "index" & x$measure == "numbers_at_age", ]
stopifnot(nrow(indices) == 648 + 580 + 128, all(indices$unit == "fish per tow"),
  all(is.na(indices$sampling_time)), all(indices$year_basis == "calendar_year"),
  !any(grepl("Larval|MARMAP|ECOMON", indices$survey)))
spring <- indices[indices$season == "spring", ]
fall <- indices[indices$season == "fall", ]
winter <- indices[indices$season == "winter", ]
stopifnot(identical(sort(unique(spring$year)), setdiff(1968:2022, 2020L)),
  identical(sort(unique(spring$age)), 0:11),
  spring$value[spring$year == 1968 & spring$age == 3] == 47.65,
  spring$value[spring$year == 2004 & spring$age == 7] == .019,
  identical(sort(unique(fall$year)), setdiff(1963:2022, c(2017L, 2020L))),
  identical(sort(unique(fall$age)), 1:10),
  fall$value[fall$year == 1984 & fall$age == 4] == 1.865,
  fall$value[fall$year == 1997 & fall$age == 4] == 1.459,
  all(grepl("incomplete age-length key", fall$notes[fall$year == 2013 & fall$age %in% c(2:7, 9)])),
  identical(sort(unique(winter$year)), 1992:2007),
  identical(sort(unique(winter$age)), 1:8),
  winter$value[winter$year == 2006 & winter$age == 2] == 18.736,
  winter$value[winter$year == 1996 & winter$age == 8] == .015,
  any(indices$value == 0), !anyNA(indices$value))
totals <- x[x$type == "index" & x$measure == "total_numbers", ]
for (season in c("spring", "fall")) {
  series <- indices[indices$season == season, ]
  annual <- aggregate(series$value, list(year = series$year), sum)
  total <- totals[totals$season == season, ]
  difference <- round(annual$x - total$value[match(annual$year, total$year)], 6)
  problematic <- abs(difference) > if (season == "spring") .0065 else .0055
  stopifnot(identical(annual$year[problematic], if (season == "spring") c(2012L, 2015L, 2017L) else 2013L),
    identical(difference[problematic], if (season == "spring") c(-.016, -.007, -.012) else -.008))
}
cv <- x[x$measure %in% c("abundance_index_cv", "biomass_index_cv"), ]
stopifnot(nrow(cv) == 240, all(is.na(cv$age)), all(cv$basis == "relative_scale"),
  all(cv$unit == "ratio"),
  cv$value[cv$season == "spring" & cv$year == 2022 & cv$measure == "abundance_index_cv"] == .540,
  cv$value[cv$season == "spring" & cv$year == 2022 & cv$measure == "biomass_index_cv"] == .609,
  totals$value[totals$season == "spring" & totals$year == 2020] == 0,
  all(is.na(cv$value[cv$season == "spring" & cv$year == 2020])),
  !any(x$type %in% c("M", "weight", "maturity")),
  any(grepl("Revised Candidate", framework$assumptions$value)),
  any(grepl("age-specific M", framework$assumptions$value, fixed = TRUE)),
  any(grepl("removed", update$assumptions$value)),
  nrow(framework$outputs) == 3,
  !any(framework$outputs$measure %in% c("SSB", "Ffull", "numbers_at_age", "recruitment")),
  framework$outputs$value[framework$outputs$measure == "Bmsy"] == 191,
  framework$outputs$lwr[framework$outputs$measure == "Bmsy"] == 11,
  framework$outputs$upr[framework$outputs$measure == "Bmsy"] == 3345,
  all(grepl("proxy", framework$outputs$notes)),
  all(is.na(framework$outputs$year)), all(is.na(framework$outputs$se)))
o <- update$outputs
stopifnot(nrow(o) == 30, all(o$year %in% 2015:2024), all(is.na(o$se)),
  all(is.na(o$lwr)), all(is.na(o$upr)),
  identical(o$value[o$measure == "SSB"], c(881, 475, 290, 219, 135, 107, 58, 36, 31, 38)),
  identical(o$value[o$measure == "Ffull"], c(.381, .32, .221, .078, .058, .061, .073, .119, .102, .051)),
  identical(o$value[o$measure == "recruitment"], c(1126, 1623, 1907, 2319, 3264, 743, 514, 537, 739, 2022)),
  all(o$age[o$measure == "recruitment"] == 1),
  all(o$unit[o$measure == "recruitment"] == "thousand fish"),
  all(o$age_group[o$measure == "Ffull"] == "5-6+"),
  !any(o$measure %in% c("Fbar", "numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age")),
  nrow(update$inputs) == 10, all(update$inputs$measure == "total_biomass"),
  !anyNA(update$inputs$value), all(is.na(update$inputs$age)))
ref <- database_to_tam_ref(update$assessment$assessment_id, o, years = 2015:2024,
  ages = 1:6, terminal_year = 2024, assumptions = update$assumptions)
stopifnot(ref$pop$ssb$est[ref$pop$ssb$year == 2024] == 38,
  all(is.na(ref$pop$F$est)), all(is.na(ref$pop$F_bar$est)), all(is.na(ref$pop$N$est)))
for (id in c(framework$assessment$assessment_id, update$assessment$assessment_id)) {
  blocked <- run_assessment(id, fit = FALSE)
  stopifnot(blocked$diagnostics$status == "blocked_no_stock_spec", is.null(blocked$fit))
}
validate_in_memory <- function(tables) {
  env <- new.env(parent = globalenv())
  env$read.csv <- function(file, ...) tables[[sub(".csv", "", basename(file), fixed = TRUE)]]
  tryCatch({
    invisible(capture.output(sys.source("analysis/comp_assessments/scripts/database/002_validate_database.R", envir = env)))
    NULL
  }, error = identity)
}
at <- which(db$inputs$assessment_id == framework$assessment$assessment_id &
  db$inputs$measure == "biomass_index_cv" & db$inputs$season == "spring" & db$inputs$year == 2022)
stopifnot(length(at) == 1L)
wrong <- db
wrong$inputs$age[at] <- 2L
error <- validate_in_memory(wrong)
stopifnot(inherits(error, "error"), grepl("without an age", conditionMessage(error), fixed = TRUE))
wrong <- db
wrong$inputs <- wrong$inputs[!(wrong$inputs$assessment_id == framework$assessment$assessment_id &
  wrong$inputs$measure == "total_biomass" & wrong$inputs$season == "spring" & wrong$inputs$year == 2022), ]
error <- validate_in_memory(wrong)
stopifnot(inherits(error, "error"), grepl("corresponding total_biomass", conditionMessage(error), fixed = TRUE))
cat("SNEMA yellowtail source tables, selected-run separation, units, missing values and unresolved translation gaps passed.\n")
