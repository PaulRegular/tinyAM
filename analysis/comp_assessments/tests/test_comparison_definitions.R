source("analysis/comp_assessments/R/run_assessment.R")
.compare_fixture <- function(fit, reference, ...) {
  measures <- c(N = "numbers_at_age", F = "fishing_mortality_at_age",
    M = "natural_mortality_at_age", recruitment = "recruitment",
    abundance = "total_numbers", biomass = "total_biomass", ssb = "SSB",
    biomass_at_age = "biomass_by_age_group", F_bar = "Fbar")
  rows <- lapply(names(reference$pop), function(metric) {
    x <- reference$pop[[metric]]
    if (is.null(x)) return(NULL)
    data.frame(assessment_id = "fixture",
      type = if (metric %in% c("N", "abundance")) "population" else
        if (metric %in% c("F", "M", "F_bar")) "mortality" else
        if (metric == "recruitment") "recruitment" else "biomass",
      measure = measures[[metric]], year = x$year,
      age = if (is.null(x$age)) NA else x$age,
      age_group = if (is.null(x$age_group)) "" else x$age_group,
      value = x$est, se = NA_real_, lwr = NA_real_, upr = NA_real_,
      unit = if (!is.null(x$unit)) x$unit else
        if (metric %in% c("N", "abundance", "recruitment")) "fish" else
        if (metric %in% c("F", "M", "F_bar")) "per year" else "kg",
      source_type = "fixture", source_reference = "fixture", notes = "")
  })
  fit$rep <- fit$obs_pred <- fit$random_par <- list()
  fit$fixed_par <- data.frame()
  ref <- database_to_tam_ref("fixture", do.call(rbind, rows),
    template = fit, comparison_scales = reference$comparison_scales, ...)
  .assessment_percent_differences(ref)
}

years <- 2000:2002
n <- expand.grid(year = years, age = 0:3)
n$est <- rep(c(100, 50, 10, 20), each = 3)
f <- n
f$est <- rep(c(0, 0, 0.1, 0.4), each = 3)
tiny_n <- n[n$age %in% 2:3, ]
fit <- list(
  dat = list(years = years, ages = 2:3, W = matrix(2, 3, 2), P = matrix(0.5, 3, 2),
             F_settings = list(mean_ages = 2:3), M_settings = list(mean_ages = 2:3)),
  pop = list(N = tiny_n, F = f[f$age %in% 2:3, ],
             recruitment = data.frame(year = years, est = 10),
             abundance = data.frame(year = years, est = 30),
             biomass = data.frame(year = years, est = 60),
             ssb = data.frame(year = years, est = 30),
             F_bar = data.frame(year = years, est = 0.3))
)
reference <- list(pop = list(
  N = n, F = f, recruitment = data.frame(year = years, age = 0L, est = 100),
  abundance = data.frame(year = years, est = 180),
  ssb = data.frame(year = years, est = 999),
  F_bar = data.frame(year = years, est = 0.25)
), comparison_scales = c(N = 1, F = 1, recruitment = 1, abundance = 1,
                         biomass = 1, ssb = 1, F_bar = 1))
original <- reference
x <- .compare_fixture(fit, reference)
stopifnot(
  identical(reference, original),
  all(x$comparison_status[x$metric == "recruitment"] == "matched"),
  all(x$percent_difference[x$metric == "recruitment"] == 0),
  all(x$source[x$metric == "abundance"] == 30),
  all(x$source[x$metric == "biomass"] == 60),
  all(x$source[x$metric == "ssb"] == 30),
  all(abs(x$source[x$metric == "F_bar"] - 0.3) < 1e-12),
  all(x$percent_difference[x$metric == "N"] == 0)
)
timed_fit <- fit
timed_fit$dat$ssb_settings <- list(spawn_time = .25)
timed_comparison <- .compare_fixture(timed_fit, reference)
stopifnot(all(timed_comparison$comparison_status[timed_comparison$metric == "ssb"] == "non_equivalent"),
          all(is.na(timed_comparison$percent_difference[timed_comparison$metric == "ssb"])))
# Accepted N without accepted mortality cannot supply spawning-time survival.
reference$pop$recruitment$age <- 2L
reference$pop$recruitment$est <- 10
x <- .compare_fixture(fit, reference)
stopifnot(all(x$percent_difference[x$metric == "recruitment"] == 0))
reference$pop$N$unit <- "million fish"
reference$pop$N$est <- reference$pop$N$est / 1e6
x <- .compare_fixture(fit, reference)
stopifnot(all(x$source[x$metric == "abundance"] == 30),
          all(x$percent_difference[x$metric == "N"] == 0))
reference <- original
reference$pop$recruitment$age <- NULL
reference$pop$recruitment$est <- 10
x <- .compare_fixture(fit, reference,
  assumptions = data.frame(setting = "recruitment_age", value = "2"))
stopifnot(all(x$percent_difference[x$metric == "recruitment"] == 0))

# Incomplete age coverage must never produce a partial annual total.
reference$pop$N <- n[!(n$year == 2001 & n$age == 3), ]
x <- .compare_fixture(fit, reference)
stopifnot(!2001 %in% x$year[x$metric == "abundance"])
reference$pop$N <- NULL
x <- .compare_fixture(fit, reference)
stopifnot(x$comparison_status[x$metric == "abundance"] == "non_equivalent",
          x$comparison_status[x$metric == "ssb"] == "non_equivalent",
          x$comparison_status[x$metric == "F_bar"] == "non_equivalent")
x <- .compare_fixture(fit, reference,
                                     comparison_aggregates = "ssb")
stopifnot(x$comparison_status[x$metric == "ssb"] == "matched",
          all(x$source[x$metric == "ssb"] == 999),
          grepl("Reported aggregate female SSB", x$definition[x$metric == "ssb"]),
          x$comparison_status[x$metric == "F_bar"] == "non_equivalent")

# Zero accepted values still contribute to absolute differences, not percentages.
example <- data.frame(metric = "N", year = rep(years, 2), age = rep(2:3, each = 3),
                      source = c(0, 1, 2, 100, 101, 102),
                      tinyAM = c(2, 1, 0, 102, 101, 100),
                      percent_difference = c(NA, 0, -100, 2, 0, -100 * 2/102))
summary <- .assessment_comparison_summary(example, "test")
stopifnot(summary$n == 6L, summary$mean_absolute_difference == 8/6,
          abs(summary$trend_correlation + 1) < 1e-12)
flat <- data.frame(
  metric = "M_bar", year = years, age = NA_integer_,
  source = 0.2 + c(-2, 0, 2) * .Machine$double.eps,
  tinyAM = 0.2 + c(1, -2, 1) * .Machine$double.eps,
  percent_difference = 0, comparison_status = "matched",
  definition = "N-weighted M", reason = "", unit = "per year"
)
stopifnot(is.na(.assessment_comparison_summary(flat, "test")$trend_correlation))
failed <- .assessment_diagnostics("test", list(commit = "test"), "not_converged",
  fit = list(opt = structure("NA/NaN gradient evaluation", class = "try-error"),
             is_converged = FALSE, sdrep = list(pdHess = FALSE),
             obj = list(env = list(random = 1:3))))
stopifnot(failed$status == "fit_failed", is.na(failed$optimizer_code),
          is.na(failed$max_abs_gradient), !failed$positive_definite_hessian,
          grepl("NA/NaN", failed$reason))

# A tinyAM plug-in M baseline must not become an invented accepted M surface.
db <- read_database()
output <- db$outputs[db$outputs$assessment_id == "dfo_cod_2j3kl_2025" &
                       db$outputs$measure == "SSB", ][1, ]
obs <- list(weight = data.frame(year = output$year, age = 2L, M_assumption = 0.3))
for (status in c("estimated_in_source", "no_numerical_input")) {
  attr(obs, "translation") <- list(M = list(status = status))
  ref <- database_to_tam_ref("dfo_cod_2j3kl_2025", output, obs = obs)
  stopifnot(is.null(ref$pop$M))
}
attr(obs, "translation") <- list(M = list(status = "fixed_numerical_input"))
ref <- database_to_tam_ref("dfo_cod_2j3kl_2025", output, obs = obs)
stopifnot(ref$pop$M$est == 0.3)
cat("Matched assessment definition tests passed.\n")
# Explicit display-age mappings keep report groups distinct from the model plus group.
ages <- 1:15
years <- 2000:2001
tiny_n <- expand.grid(year = years, age = ages)
tiny_n$est <- tiny_n$age * 1e9
tiny_bio <- tiny_n
tiny_bio$est <- tiny_bio$age * 1000
accepted_n <- expand.grid(year = years, age = 1:10)
accepted_n$age_group <- NA_character_
accepted_n$est <- NA_real_
for (year in years) {
  for (age in 1:9) {
    row <- accepted_n$year == year & accepted_n$age == age
    accepted_n$est[row] <- tiny_n$est[tiny_n$year == year & tiny_n$age == age] / 1e9
  }
  row <- accepted_n$year == year & accepted_n$age == 10
  accepted_n$age_group[row] <- "10+"
  accepted_n$est[row] <- sum(tiny_n$est[tiny_n$year == year & tiny_n$age %in% 10:15]) / 1e9
}
accepted_n$unit <- "billion fish"
accepted_bio <- data.frame(
  year = years, age_group = "3+", est = vapply(years, function(year) {
    sum(tiny_bio$est[tiny_bio$year == year & tiny_bio$age %in% 3:15]) / 1e6
  }, numeric(1)), unit = "thousand t", stringsAsFactors = FALSE
)
accepted_ssb <- data.frame(year = years, est = 100, unit = "thousand t")
accepted_rec <- data.frame(year = years, age = 1, est = 1, unit = "million fish")
fit_grouped <- list(
  dat = list(years = years, ages = ages),
  pop = list(
    N = tiny_n,
    abundance = aggregate(est ~ year, tiny_n, sum),
    biomass_at_age = tiny_bio,
    ssb = data.frame(year = years, est = 100),
    recruitment = data.frame(year = years, age = 1, est = 1e6)
  )
)
ref_grouped <- list(pop = list(
  N = accepted_n, biomass_at_age = accepted_bio, ssb = accepted_ssb,
  recruitment = accepted_rec
))
attr(ref_grouped, "source_pop") <- ref_grouped$pop
ref_grouped$comparison_scales <- c(
  N = 1e-9, abundance = 1e-9, recruitment = 1e-6,
  biomass_at_age = 1e-6, ssb = 1e-6
)
grouped <- .compare_fixture(
  fit_grouped, ref_grouped, comparison_age_groups = list(
    N = list("10+" = 10:15),
    biomass_at_age = list("3+" = 3:15)
  ), comparison_aggregates = "ssb",
  comparison_definitions = list(ssb = list(
    status = "approximate",
    definition = "Spawning-time source SSB versus begin-year tinyAM SSB",
    reason = "The within-year survival adjustment differs."
  ))
)
n_group <- grouped[grouped$metric == "N" & grouped$age == "10+", ]
n_age <- grouped[grouped$metric == "N" & grouped$age == "10", ]
abundance <- grouped[grouped$metric == "abundance", ]
bio_group <- grouped[grouped$metric == "biomass_at_age" & grouped$age == "3+", ]
ssb_comparison <- grouped[grouped$metric == "ssb", ]
stopifnot(
  nrow(n_group) == length(years),
  all(n_group$comparison_status == "matched"),
  all(n_group$source == n_group$tinyAM),
  grepl("10-15", n_group$definition[1L], fixed = TRUE),
  nrow(n_age) == 0L,
  nrow(abundance) == length(years),
  all(abundance$source == abundance$tinyAM),
  grepl("explicit reported age groups", abundance$definition[1L], fixed = TRUE),
  nrow(bio_group) == length(years),
  all(bio_group$source == bio_group$tinyAM),
  all(ssb_comparison$comparison_status == "approximate"),
  grepl("survival adjustment", ssb_comparison$reason[1L], fixed = TRUE)
)
unmapped <- .compare_fixture(fit_grouped, ref_grouped)
stopifnot(
  unmapped$comparison_status[unmapped$metric == "N"] == "non_equivalent",
  all(is.na(unmapped$source[unmapped$metric == "N"])),
  unmapped$comparison_status[unmapped$metric == "abundance"] == "non_equivalent"
)

grouped_output <- data.frame(
  assessment_id = "grouped",
  type = "biomass",
  measure = "biomass_by_age_group",
  year = 2020L,
  age = NA_integer_,
  age_group = "3+",
  value = 10,
  se = NA_real_,
  lwr = NA_real_,
  upr = NA_real_,
  unit = "thousand t",
  source_type = "official_table",
  source_reference = "Table",
  notes = "Reported 3+ biomass.",
  stringsAsFactors = FALSE
)
template <- list(
  dat = list(years = 2020L, ages = 1:15, is_proj = FALSE),
  pop = list(),
  rep = list(),
  obs_pred = list(),
  fixed_par = data.frame(par = character(), coef = character(), est = numeric()),
  random_par = list()
)
grouped_ref <- database_to_tam_ref(
  "grouped", grouped_output, years = 2020L, ages = 1:15,
  age_plus_group = 15L, template = template
)
stopifnot(nrow(grouped_ref$pop$biomass_at_age) == 15L,
          all(is.na(grouped_ref$pop$biomass_at_age$est)),
          attr(grouped_ref, "native_pop")$biomass_at_age$age_group == "3+")
cat("Reported age-group comparison tests passed.\n")
