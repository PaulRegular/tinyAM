source("analysis/comp_assessments/R/run_assessment.R")

years <- 2000:2002
ages <- 2:3
source_grid <- expand.grid(year = 1999:2003, age = 0:4)
source_grid$value <- source_grid$year * 10 + source_grid$age
source_grid$value[source_grid$year == 2001 & source_grid$age == 4] <- NA
output <- transform(source_grid, assessment_id = "aligned", type = "population",
  measure = "numbers_at_age", age_group = "", se = 2, lwr = value - 4,
  upr = value + 4, unit = "thousand fish", source_type = "official_table",
  source_reference = "N table", notes = "Beginning of year")
recruitment <- output[output$age == 0, ]
recruitment$type <- "recruitment"
recruitment$measure <- "recruitment"
recruitment$value <- 999
recruitment$source_reference <- "Native age-0 recruitment"
outputs <- rbind(output, recruitment)
before <- outputs
model_grid <- expand.grid(year = years, age = ages)
model_grid$est <- 1
annual <- data.frame(year = years, est = 1)
fit <- list(
  dat = list(years = years, ages = ages, is_proj = rep(FALSE, 3)),
  pop = list(N = model_grid, recruitment = annual),
  rep = list(N = matrix(NA_real_, 3, 2, dimnames = list(year = years, age = ages)),
             recruitment = setNames(rep(NA_real_, 3), years)),
  obs_pred = list(), fixed_par = data.frame(par = "fixture", est = 1), random_par = list()
)
ref <- database_to_tam_ref("aligned", outputs, template = fit, age_plus_group = 3)
expected <- years * 10 + 2
stopifnot(
  identical(outputs, before),
  identical(ref$dat$years, years), identical(ref$dat$ages, ages),
  nrow(ref$pop$N) == 6L,
  all(ref$pop$N$year %in% years), all(ref$pop$N$age %in% ages),
  all(ref$pop$recruitment$est == expected),
  all(ref$pop$recruitment$unit == "thousand fish"),
  all(ref$pop$recruitment$se == 2),
  all(ref$rep$recruitment == expected),
  all(attr(ref, "native_pop")$recruitment$est == 999),
  all(attr(ref, "native_pop")$recruitment$age == 0),
  all(grepl("same calendar year", ref$pop$recruitment$notes)),
  ref$comparison_scales[["recruitment"]] == 1e-3,
  is.na(ref$pop$N$est[ref$pop$N$year == 2001 & ref$pop$N$age == 3]),
  all(is.na(ref$pop$N$se[ref$pop$N$age == 3]))
)
differences <- .assessment_percent_differences(ref)
r <- differences[differences$metric == "recruitment", ]
stopifnot(all(r$source == expected), all(r$year == years),
          all(r$comparison_status == "matched"),
          identical(tinyAM::tidy_tam(model_list = list(Accepted = ref))$pop$recruitment$est,
                    ref$pop$recruitment$est))

singleton_labels <- outputs
singleton_labels$age_group <- as.character(singleton_labels$age)
labelled_ref <- database_to_tam_ref("aligned", singleton_labels, template = fit)
stopifnot(all(labelled_ref$pop$recruitment$est == expected))

# Missing accepted N at the requested recruitment age stays missing, never shifted.
missing <- outputs[!(outputs$year == 2001 & outputs$age == 2), ]
ref_missing <- database_to_tam_ref("aligned", missing, template = fit, age_plus_group = 3)
missing_differences <- .assessment_percent_differences(ref_missing)
missing_recruitment <- missing_differences[missing_differences$metric == "recruitment" &
                                          missing_differences$year == 2001, ]
stopifnot(is.na(ref_missing$pop$recruitment$est[ref_missing$pop$recruitment$year == 2001]),
          is.na(missing_recruitment$source), is.na(missing_recruitment$percent_difference))
no_n <- outputs[outputs$type == "recruitment", ]
ref_no_n <- database_to_tam_ref("aligned", no_n, template = fit)
stopifnot(all(is.na(ref_no_n$pop$recruitment$est)),
          all(attr(ref_no_n, "native_pop")$recruitment$est == 999))

# Reference dimensions cannot silently differ from the fitted template.
for (arguments in list(list(years = 1999:2002), list(ages = 1:3),
                       list(years = c(2000, 2002)), list(ages = c(2, 2)))) {
  error <- tryCatch(do.call(database_to_tam_ref,
    c(list(assessment_id = "aligned", outputs = outputs, template = fit), arguments)),
    error = identity)
  stopifnot(inherits(error, "error"))
}
cat("Central reference alignment tests passed.\n")
