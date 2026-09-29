test_that("tidy uncertainty keeps SE scales explicit without changing estimates or CIs", {
  par <- list(log_q = c(baseline = log(2)), logit_phi_f = c(age = qlogis(.7)),
              mu_m = c(effect = -.2), dq = c(step = .3))
  obj <- RTMB::MakeADFun(function(p) {
    sum((p$log_q - par$log_q)^2 + (p$logit_phi_f - par$logit_phi_f)^2 +
        (p$mu_m - par$mu_m)^2 + (p$dq - par$dq)^2) / (2 * .1^2)
  }, par, silent = TRUE)
  fit <- structure(list(obj = obj, sdrep = RTMB::sdreport(obj), dat = list()),
                   class = "tam_fit")
  tab <- tidy_par(fit)$fixed
  expect_identical(tail(names(tab), 5), c("est", "lwr", "upr", "se", "se_scale"))
  expect_equal(tab$se, rep(.1, 4))
  expect_identical(tab$se_scale, c("log", "logit", "reported", "reported"))
  expect_equal(tab$est, c(2, .7, -.2, .3))
  z <- qnorm(.975)
  expect_equal(tab$lwr, c(exp(log(2) - z * .1), plogis(qlogis(.7) - z * .1),
                          -.2 - z * .1, .3 - z * .1))
  expect_equal(tab$upr, c(exp(log(2) + z * .1), plogis(qlogis(.7) + z * .1),
                          -.2 + z * .1, .3 + z * .1))

  table <- tinyAM:::.coef_table(tab)
  expect_s3_class(table, "data.frame")
  expect_identical(names(table), c("Estimate", "Lower 95%", "Upper 95%", "Std. Error", "SE scale"))
  expect_type(table$Estimate, "double")
  expect_equal(table$`Std. Error`, tab$se)
  display <- tinyAM:::.format_coef_display(table)
  expect_identical(unname(display[, "SE scale"]), tab$se_scale)
  expect_equal(unname(display[, "Std. Error"]), rep("0.100", 4))

  attr(tab, "interval") <- .9
  expect_identical(names(tinyAM:::.coef_table(tab))[2:3], c("Lower 90%", "Upper 90%"))
  expect_equal(nrow(tinyAM:::.coef_table(NULL)), 0)
})

test_that("terminal tables show confidence limits before log-scale SEs", {
  pop <- list(abundance = data.frame(year = c(2024, 2025), est = c(10000, 20000),
    lwr = c(8000, 15000), upr = c(13000, 25000), se = c(.1, .2),
    se_scale = "log", is_proj = c(FALSE, TRUE)))
  attr(pop, "interval") <- .9
  tab <- tinyAM:::.terminal_table(pop, 2024)
  expect_identical(names(tab), c("Estimate", "Lower 90%", "Upper 90%", "Std. Error", "SE scale"))
  expect_equal(tab$Estimate, 10000)
  expect_equal(tab$`Std. Error`, .1)
  display <- tinyAM:::.format_terminal_display(tab)
  expect_equal(unname(display[1, ]), c("10,000", "8,000", "13,000", "0.100", "log"))
  expect_equal(nrow(tinyAM:::.terminal_table(NULL, 2024)), 0)
})

test_that("tidy population and random parameters use the same uncertainty order", {
  pop <- tidy_sdrep(default_fit)
  expect_identical(names(pop$ssb), c("year", "est", "lwr", "upr", "se", "se_scale", "is_proj"))
  expect_true(all(pop$ssb$se_scale == "log"))
  expected <- as.list(default_fit$sdrep, "Std. Error", report = TRUE)$log_ssb
  expect_equal(pop$ssb$se, as.numeric(expected))
  random <- tidy_par(default_fit)$random$log_f
  expect_identical(tail(names(random), 6), c("est", "lwr", "upr", "se", "se_scale", "is_proj"))
  expect_true(all(random$se_scale == "log"))

  printed <- paste(capture.output(print(default_fit)), collapse = "\n")
  expect_match(printed, "Estimate +Lower 95% +Upper 95% +Std. Error +SE scale")
  expect_match(printed, "0.10 is about 10%", fixed = TRUE)
  expect_match(paste(capture.output(print(summary(default_fit))), collapse = "\n"), "SE scale")
  expect_type(summary(default_fit)$terminal$abundance, "double")
})
