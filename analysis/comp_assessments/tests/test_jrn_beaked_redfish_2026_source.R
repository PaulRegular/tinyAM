pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
record <- read_assessment("jrn_redfish_beaked_2026", read_database())
inputs <- record$inputs
outputs <- record$outputs
stopifnot(
  record$stocks$charbonneau_id == "ICES-AFWG_mentella_DEEP 1,2",
  record$stocks$scientific_name == "Sebastes mentella",
  record$stocks$ocean == "Northeast Atlantic",
  record$assessment$model_family == "SCAA",
  record$assessment$terminal_year == 2025,
  record$assessment$is_current, record$assessment$is_applied,
  all(unlist(record$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  all(inputs$source_type[inputs$type != "M"] == "official_machine_readable"),
  all(inputs$source_type[inputs$type == "M"] == "official_document"),
  all(is.na(outputs$se)), all(is.na(outputs$lwr)), all(is.na(outputs$upr))
)
select_input <- function(measure, survey = NULL, fleet = NULL) {
  x <- inputs[inputs$measure == measure, ]
  if (!is.null(survey)) x <- x[!is.na(x$survey) & x$survey == survey, ]
  if (!is.null(fleet)) x <- x[!is.na(x$fleet) & x$fleet == fleet, ]
  x
}
cell <- function(x, year, age) {
  value <- x$value[!is.na(x$year) & x$year == year & !is.na(x$age) & x$age == age]
  stopifnot(length(value) == 1L)
  value
}
total <- select_input("numbers_at_age", fleet = "All fisheries (source total)")
pelagic <- select_input("numbers_at_age", fleet = "Pelagic outside EEZ (subset)")
stopifnot(
  nrow(total) == 34 * 14, nrow(pelagic) == 20 * 13,
  all(total$unit == "thousand fish"), all(pelagic$unit == "thousand fish"),
  cell(total, 1992, 7) == 2498, cell(total, 2025, 6) == 1971,
  cell(total, 2025, 19) == 42500, cell(pelagic, 2025, 19) == 3366,
  is.na(cell(pelagic, 2025, 7)),
  all(grepl("Do not add", c(total$notes, pelagic$notes), fixed = TRUE))
)
weight <- inputs[inputs$type == "weight", ]
catch_weight <- inputs[inputs$type == "catch_weight", ]
maturity <- inputs[inputs$type == "maturity", ]
M <- inputs[inputs$type == "M", ]
stopifnot(
  identical(weight$age, 6:19), all(is.na(weight$year)),
  identical(weight$value, c(.141, .188, .237, .286, .334, .381, .424, .465, .503, .537, .569, .597, .623, .755)),
  identical(catch_weight$value, weight$value),
  identical(maturity$age, 6:19), all(is.na(maturity$year)),
  identical(maturity$value, c(.02, .04, .08, .14, .24, .37, .49, .58, .66, .73, .77, .81, .85, 1)),
  identical(M$age, 2:19), all(M$value == .05), all(is.na(M$year))
)
russian <- select_input("numbers_at_age", survey = "Russian autumn groundfish")
eco <- select_input("numbers_at_age", survey = "Ecosystem survey")
winter <- select_input("numbers_at_age", survey = "Winter bottom trawl")
grouped <- select_input("index_by_age_group", survey = "Ecosystem survey")
stopifnot(
  identical(sort(unique(russian$year)), 1974:2015),
  all(russian$year_basis == "birth_cohort"), all(russian$sampling_time == .90),
  all(russian$unit == "fish per hour"),
  sum(grepl("source prints '+'", russian$notes, fixed = TRUE)) == 6,
  is.na(cell(russian, 1997, 2)), cell(russian, 1974, 2) == 4.8,
  all(eco$year_basis == "calendar_year"), all(eco$sampling_time == .75),
  all(eco$unit == "million fish"), all(eco$age %in% 2:15),
  all(is.na(eco$value[eco$year == 2010])),
  cell(eco, 2025, 2) == 40, cell(eco, 2025, 3) == 39,
  cell(eco, 2011, 14) == 0,
  all(grouped$age == 16), cell(grouped, 2025, 16) == 130,
  all(winter$sampling_time == .12), all(winter$unit == "million fish"),
  all(is.na(winter$value[winter$year %in% c(2012, 2019)])),
  !any(inputs$measure == "proportion_at_age")
)
select_output <- function(measure) outputs[outputs$measure == measure, ]
N <- select_output("numbers_at_age")
F <- select_output("fishing_mortality_at_age")
stopifnot(
  nrow(N) == 612, nrow(F) == 612,
  cell(N, 1992, 3) == 435071, cell(N, 2025, 6) == 591807,
  all(N$age_group[N$age == 19] == "19+"),
  all(F$source_type == "derived_source_output"),
  !any(outputs$measure == "natural_mortality_at_age")
)
states <- function(x) {
  g <- expand.grid(year = 1992:2025, age = 2:19)
  key <- function(z) paste(z$year, z$age)
  matrix(x$value[match(key(g), key(x))], nrow = 34)
}
n <- states(N); f <- states(F)
coefficients <- select_output("fleet_F_coefficient")
fy_d <- coefficients$value[coefficients$fleet == "Demersal"]
fy_p <- coefficients$value[coefficients$fleet == "Pelagic"]
selectivity <- select_output("fleet_selectivity_at_age")
dem <- states(selectivity[selectivity$fleet == "Demersal", ])
pel <- selectivity$value[selectivity$fleet == "Pelagic"]
stopifnot(isTRUE(all.equal(f, dem * fy_d + rep(pel, each = 34) * fy_p)))
cohort <- log(n[-1, 2:17]) - log(n[-34, 1:16]) + .05 + f[-34, 1:16]
f_rounding <- .0005 * (dem + rep(pel, each = 34) + fy_d + fy_p) + .0000005
n_rounding <- -log1p(-.5 / n)
# M is reported as 0.05; its two-decimal precision also limits this check.
cohort_bound <- .005 + f_rounding[-34, 1:16] + n_rounding[-34, 1:16] + n_rounding[-1, 2:17]
stopifnot(all(abs(cohort) <= cohort_bound))
bound <- .0015 + .0005 * max(fy_d + fy_p) + .0000005
stopifnot(
  max(abs(rowMeans(f[, 11:17]) - select_output("Fbar")$value)) <= bound,
  max(abs(f[, 18] - select_output("F_age19plus")$value)) <= bound,
  identical(rowSums(n), select_output("total_numbers")$value),
  all(select_output("total_numbers")$unit == "thousand fish"),
  all(select_output("recruitment")$age == 2),
  max(abs(n[, 1] / 1000 - select_output("recruitment")$value)) <= .5,
  max(abs(n[, 5] / 1000 - select_output("recruitment_age6_reported")$value)) <= .5,
  all(grepl("inconsistent", select_output("recruitment_age6_reported")$unit)),
  select_output("SSB")$value[select_output("SSB")$year == 2025] == 978795,
  select_output("total_biomass")$value[select_output("total_biomass")$year == 2025] == 1536941
)
prepared <- run_assessment(record$assessment$assessment_id, fit = FALSE)
stopifnot(
  prepared$diagnostics$status == "not_fitted",
  is.null(prepared$fit),
  isTRUE(tinyAM::check_obs(prepared$obs))
)
cat("Beaked redfish source coverage, fixed biology, missing cells, year basis, units and rounded F passed.\n")
