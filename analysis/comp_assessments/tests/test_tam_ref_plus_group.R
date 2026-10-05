root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "database_to_tam_ref.R"))

years <- 2000:2002
ages <- 1:2
grid <- expand.grid(year = years, age = ages)
obs <- list(
  catch = data.frame(year = 2000, age = 1, fleet = "fishery", obs = 1),
  index = data.frame(year = 2000, age = 1, survey = "RV", obs = 1),
  weight = data.frame(year = rep(years, each = 2), age = rep(ages, 3),
                      M_assumption = 0.2),
  maturity = data.frame(year = rep(years, each = 2), age = rep(ages, 3),
                        obs = 0.5)
)
template <- list(
  call = quote(fit_tam()), refit_args = list(),
  dat = list(obs = obs, years = years, ages = ages,
             is_proj = rep(FALSE, length(years))),
  pop = list(
    N = data.frame(grid, est = NA_real_, is_proj = FALSE),
    F = data.frame(grid, est = NA_real_, is_proj = FALSE),
    M = data.frame(grid, est = NA_real_, is_proj = FALSE)
  ),
  rep = list(N = matrix(NA_real_, length(years), length(ages),
                        dimnames = list(years, ages))),
  obs_pred = list(catch = data.frame(year = 2000, age = 1),
                  index = data.frame(year = 2000, age = 1)),
  fixed_par = data.frame(par = character(), coef = character(), est = numeric()),
  random_par = list()
)
outputs <- data.frame(
  assessment_id = "plus_group_fixture",
  type = c(rep("population", 3), rep("mortality", 3)),
  measure = c(rep("numbers_at_age", 3), rep("fishing_mortality_at_age", 3)),
  year = rep(2000, 6), age = c(1:3, 1:3),
  age_group = c(NA, "2+", "2+", NA, "2+", "2+"),
  value = c(100, 60, 40, 0.1, 0.1, 0.2),
  se = NA_real_, lwr = NA_real_, upr = NA_real_,
  unit = c(rep("fish", 3), rep("per year", 3)),
  source_type = "official_table", source_reference = "fixture",
  notes = NA_character_, stringsAsFactors = FALSE
)

reference <- database_to_tam_ref(
  "plus_group_fixture", outputs, obs = obs, years = years, ages = ages,
  terminal_year = 2002, age_plus_group = 2L, template = template
)
stopifnot(reference$pop$N$est[
  reference$pop$N$year == 2000 & reference$pop$N$age == 2
] == 100)
stopifnot(isTRUE(all.equal(reference$pop$F$est[
  reference$pop$F$year == 2000 & reference$pop$F$age == 2
], 0.14)))

cat("Plus-group comparison mapping passed.\n")
