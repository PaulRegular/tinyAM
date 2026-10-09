root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "database_to_tam_ref.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

years <- 2000:2002
ages <- 1:2
obs <- list(
  catch = data.frame(year = 2000, age = 1, fleet = "fishery", obs = 1),
  index = data.frame(year = 2000, age = 1, survey = "RV", obs = 1),
  weight = expand.grid(year = years, age = 1:2),
  maturity = expand.grid(year = years, age = 1:2)
)
obs$weight$M_assumption <- 0.2
obs$maturity$obs <- 0.5
output_rows <- function(type, measure, year, age, value, se = NA_real_,
                        lwr = NA_real_, upr = NA_real_) {
  data.frame(
    assessment_id = "plus_group_fixture", type = type, measure = measure,
    year = year, age = age, age_group = ifelse(age >= 2, "2+", NA_character_),
    value = value, se = se, lwr = lwr, upr = upr, unit = "per year",
    source_type = "official_table", source_reference = "fixture",
    notes = "Source table row", stringsAsFactors = FALSE
  )
}
outputs <- rbind(
  output_rows("population", "numbers_at_age", 2000, 1:3, c(100, 60, 40)),
  output_rows("mortality", "fishing_mortality_at_age", 2000, 1:3, c(0.1, 0.1, 0.2)),
  output_rows("mortality", "natural_mortality_at_age", 2000, 1:3, c(0.1, 0.3, 0.6)),
  output_rows("mortality", "fishing_mortality_at_age", 2001, 2, 0.3,
              se = 0.03, lwr = 0.2, upr = 0.4),
  output_rows("mortality", "natural_mortality_at_age", 2001, 2, 0.4,
              se = 0.04, lwr = 0.3, upr = 0.5)
)
outputs$unit[outputs$measure == "numbers_at_age"] <- "fish"
reference <- database_to_tam_ref(
  "plus_group_fixture", outputs, obs = obs, years = years, ages = ages,
  terminal_year = 2002, age_plus_group = 2L
)
plus_n <- reference$pop$N[reference$pop$N$year == 2000 & reference$pop$N$age == 2, ]
plus_f <- reference$pop$F[reference$pop$F$year == 2000 & reference$pop$F$age == 2, ]
plus_m <- reference$pop$M[reference$pop$M$year == 2000 & reference$pop$M$age == 2, ]
terminal_f <- reference$pop$F[reference$pop$F$year == 2001 & reference$pop$F$age == 2, ]
terminal_m <- reference$pop$M[reference$pop$M$year == 2001 & reference$pop$M$age == 2, ]
stopifnot(
  nrow(plus_n) == 1L, plus_n$est == 100,
  nrow(plus_f) == 1L, isTRUE(all.equal(plus_f$est, 0.14)),
  nrow(plus_m) == 1L, isTRUE(all.equal(plus_m$est, 0.42)),
  is.na(plus_f$se), is.na(plus_f$lwr), is.na(plus_f$upr),
  is.na(plus_m$se), is.na(plus_m$lwr), is.na(plus_m$upr),
  nrow(terminal_f) == 1L, terminal_f$est == 0.3,
  terminal_f$se == 0.03, terminal_f$lwr == 0.2, terminal_f$upr == 0.4,
  nrow(terminal_m) == 1L, terminal_m$est == 0.4,
  terminal_m$se == 0.04, terminal_m$lwr == 0.3, terminal_m$upr == 0.5,
  grepl("N-weighted", plus_f$notes, fixed = TRUE),
  grepl("retained without aggregation", terminal_m$notes, fixed = TRUE)
)

# Fixed M on the model grid has one terminal-age row per year and needs no N surface.
fixed_obs <- obs
attr(fixed_obs, "translation") <- list(M = list(status = "fixed_numerical_input"))
fixed_outputs <- output_rows("biomass", "SSB", 2000, NA_real_, 1)
fixed_outputs$age_group <- NA_character_
fixed_reference <- database_to_tam_ref(
  "plus_group_fixture", fixed_outputs, obs = fixed_obs, years = years, ages = 1:2,
  terminal_year = 2002, age_plus_group = 2L
)
stopifnot(
  nrow(fixed_reference$pop$M) == length(years) * 2L,
  all(is.finite(fixed_reference$pop$M$est)),
  all(fixed_reference$pop$M$est == 0.2)
)

cat("Plus-group comparison mapping passed.\n")
