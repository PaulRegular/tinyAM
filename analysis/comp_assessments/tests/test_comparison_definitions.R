source("analysis/comp_assessments/R/run_assessment.R")
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
x <- .assessment_percent_differences(fit, reference)
stopifnot(
  identical(reference, original),
  x$comparison_status[x$metric == "recruitment"] == "non_equivalent",
  all(is.na(x$percent_difference[x$metric == "recruitment"])),
  all(x$source[x$metric == "abundance"] == 30),
  all(x$source[x$metric == "biomass"] == 60),
  all(x$source[x$metric == "ssb"] == 30),
  all(abs(x$source[x$metric == "F_bar"] - 0.3) < 1e-12),
  all(x$percent_difference[x$metric == "N"] == 0)
)
reference$pop$recruitment$age <- 2L
reference$pop$recruitment$est <- 10
x <- .assessment_percent_differences(fit, reference)
stopifnot(all(x$percent_difference[x$metric == "recruitment"] == 0))
reference$pop$N$unit <- "million fish"
reference$pop$N$est <- reference$pop$N$est / 1e6
x <- .assessment_percent_differences(fit, reference)
stopifnot(all(x$source[x$metric == "abundance"] == 30),
          all(x$percent_difference[x$metric == "N"] == 0))
reference <- original
reference$pop$recruitment$age <- NULL
reference$pop$recruitment$est <- 10
x <- .assessment_percent_differences(fit, reference,
  assumptions = data.frame(setting = "recruitment_age", value = "2"))
stopifnot(all(x$percent_difference[x$metric == "recruitment"] == 0))

# Incomplete age coverage must never produce a partial annual total.
reference$pop$N <- n[!(n$year == 2001 & n$age == 3), ]
x <- .assessment_percent_differences(fit, reference)
stopifnot(!2001 %in% x$year[x$metric == "abundance"])
reference$pop$N <- NULL
x <- .assessment_percent_differences(fit, reference)
stopifnot(x$comparison_status[x$metric == "abundance"] == "non_equivalent",
          x$comparison_status[x$metric == "ssb"] == "non_equivalent",
          x$comparison_status[x$metric == "F_bar"] == "non_equivalent")
x <- .assessment_percent_differences(fit, reference,
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
