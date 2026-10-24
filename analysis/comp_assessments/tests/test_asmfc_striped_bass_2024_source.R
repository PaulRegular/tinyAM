pkgload::load_all(".", quiet = TRUE)
source("analysis/comp_assessments/R/run_assessment.R")
record <- read_assessment("asmfc_striped_bass_atlantic_2024", read_database())
inputs <- record$inputs; outputs <- record$outputs
stopifnot(
  record$stocks$authority == "ASMFC", record$stocks$ocean == "Northwest Atlantic",
  record$stocks$charbonneau_id == "NEFSC-GARMIII_saxatilis_GOM, CHATT",
  record$stocks$scientific_name == "Morone saxatilis",
  record$assessment$is_current, record$assessment$is_applied,
  record$assessment$terminal_year == 2023, record$assessment$framework_year == 2018,
  all(unlist(record$assessment[c("inputs_status", "outputs_status", "assumptions_status")]) == "partial"),
  nrow(inputs) == 129, nrow(outputs) == 3486,
  all(inputs$source_type == "official_table"), all(outputs$source_type == "official_table"),
  !any(grepl("github|simulation|115", outputs$source_reference)),
  !any(inputs$measure %in% c("numbers_at_age", "proportion_at_age", "weight_at_age")),
  all(is.na(inputs$sampling_time)), all(is.na(outputs$lwr)), all(is.na(outputs$upr))
)
catch <- inputs[inputs$type == "catch", ]
stopifnot(nrow(catch) == 84, all(catch$measure == "total_numbers"),
  all(catch$unit == "fish"), all(is.na(catch$age)),
  identical(sort(unique(catch$year)), 1982:2023),
  catch$value[catch$year == 2023 & catch$fleet == "Chesapeake Bay"] == 1418439,
  catch$value[catch$year == 2023 & catch$fleet == "Ocean"] == 4163671)
bio <- function(measure) inputs[inputs$measure == measure, ]
stopifnot(
  identical(bio("natural_mortality_at_age")$value, c(1.13,.68,.45,.33,.25,.19,rep(.15,9))),
  identical(bio("maturity_at_age")$value, c(0,0,0,.09,.32,.45,.84,.89,rep(1,7))),
  identical(bio("proportion_female_at_age")$value, c(.53,.56,.56,.52,.57,.65,.73,.81,.88,.92,.95,.97,1,1,1)),
  all(bio("maturity_at_age")$sex == "female"),
  all(is.na(inputs$year[inputs$type != "catch"]))
)
select <- function(measure) outputs[outputs$measure == measure, ]
cell <- function(x, year, age) {
  row <- x[!is.na(x$year) & x$year == year & !is.na(x$age) & x$age == age, ]
  stopifnot(nrow(row) == 1L)
  row$value
}
N <- select("numbers_at_age")
F <- select("fishing_mortality_at_age")
total_F <- F[is.na(F$fleet), ]
stopifnot(
  nrow(N) == 630, nrow(F) == 1890, nrow(total_F) == 630,
  cell(N, 1982, 2) == 8340660, cell(N, 2023, 1) == 94898600,
  cell(N, 2023, 15) == 631181,
  all(N$unit == "fish"), all(N$age_group[N$age == 15] == "15+"),
  all(is.na(N$se)), all(is.na(F$se)),
  cell(total_F, 1982, 15) == .1829,
  !any(outputs$measure %in% c("natural_mortality_at_age", "Fbar"))
)
states <- function(x) {
  g <- expand.grid(year = 1982:2023, age = 1:15)
  matrix(x$value[match(paste(g$year,g$age),paste(x$year,x$age))], nrow = 42)
}
n <- states(N); f <- states(total_F)
bay <- states(F[!is.na(F$fleet) & F$fleet == "Chesapeake Bay", ])
ocean <- states(F[!is.na(F$fleet) & F$fleet == "Ocean", ])
rounding <- .5 * 10^(floor(log10(n)) - 5)
stopifnot(
  max(abs(f - bay - ocean)) <= .00015 + 1e-12,
  all(abs(rowSums(n) - select("total_numbers")$value) <= rowSums(rounding) + .5),
  all(abs(rowSums(n[,8:15]) - select("total_numbers_age8plus")$value) <= rowSums(rounding[,8:15]) + .5),
  identical(n[,1], select("recruitment")$value),
  all(select("recruitment")$age == 1)
)
ssb <- select("SSB")
stopifnot(
  ssb$value[ssb$year == 2023] == 86536, ssb$se[ssb$year == 2023] == 9309,
  all(ssb$sex == "female"), all(ssb$unit == "tonnes"),
  all(grepl("Natural-scale SE", ssb$notes, fixed = TRUE)),
  nrow(select("mature_biomass_at_age")) == 630,
  all(select("mature_biomass_at_age")$sex == "female"),
  !any(!is.na(outputs$se) & !outputs$measure %in% c("SSB", "Ffull"))
)
ssb_age <- states(select("mature_biomass_at_age"))
stopifnot(all(abs(rowSums(ssb_age) - select("SSB_age_table_total")$value) <= 15 * .5 + .05),
  max(abs(ssb$value-select("SSB_age_table_total")$value)) <= 1)
full <- select("Ffull")
full_total <- full[is.na(full$fleet), ]
stopifnot(full_total$value[full_total$year == 2023] == .183,
  full_total$se[full_total$year == 2023] == .024,
  all(grepl("Natural-scale standard error",full$notes,fixed=TRUE)))
# A single-stock reference uses the reported total F, rather than adding its overlapping fleets.
reference_outputs <- outputs[!(outputs$measure == "fishing_mortality_at_age" & !is.na(outputs$fleet)), ]
reference <- database_to_tam_ref(record$assessment$assessment_id, reference_outputs,
  years = 1982:2023, ages = 1:15, terminal_year = 2023, assumptions = record$assumptions)
stopifnot(nrow(reference$pop$N) == 630, nrow(reference$pop$F) == 630,
  reference$pop$ssb$se[reference$pop$ssb$year == 2023] == 9309,
  is.null(reference$pop$F_bar))
blocked <- run_assessment(record$assessment$assessment_id, fit = FALSE)
stopifnot(blocked$diagnostics$status == "blocked_no_stock_spec", is.null(blocked$fit))
cat("Striped bass base-run states, female biology, fleet separation, rounding and uncertainty scales passed.\n")
