source("analysis/comp_assessments/scripts/database/002_validate_database.R")
id <- "afsc_pollock_ebs_2024"
x <- inputs[inputs$assessment_id == id, ]
y <- outputs[outputs$assessment_id == id, ]
a <- assessments[assessments$assessment_id == id, ]
stopifnot(nrow(x) == 3826L, nrow(y) == 793L,
          a$assessment_year == 2024, a$terminal_year == 2024,
          a$inputs_status == "partial", a$outputs_status == "partial",
          a$assumptions_status == "partial")
comp <- x[x$measure == "proportion_at_age", ]
sums <- aggregate(value ~ year + type + fleet + survey,
                  transform(comp, fleet = ifelse(is.na(fleet), "", fleet),
                            survey = ifelse(is.na(survey), "", survey)), sum)
stopifnot(all(abs(sums$value - 1) < 1e-10),
          nrow(comp[comp$type == "catch", ]) == 900L,
          !any(comp$year[comp$type == "catch"] == 2024))
ats <- x[!is.na(x$survey) & x$survey == "NMFS acoustic-trawl", ]
stopifnot(all(ats$age[ats$measure == "proportion_at_age"] >= 2),
          all(ats$sampling_time == 0.5))
age1 <- x[!is.na(x$survey) & x$survey == "NMFS acoustic-trawl age-1 index", ]
stopifnot(nrow(age1) == 18L, all(age1$age == 1), !any(age1$year == 2024))
catch <- x[x$type == "catch" & x$measure == "total_biomass", ]
stopifnot(nrow(catch) == 61L, catch$value[catch$year == 2024] == 1300)
maturity <- x[x$type == "maturity", ]
stopifnot(nrow(maturity) == 15L, all(maturity$year == 1964))
n <- y[y$measure == "numbers_at_age", ]
r <- y[y$measure == "recruitment", ]
stopifnot(nrow(n) == 610L, nrow(r) == 61L,
          all(n$age_group[n$age == 10] == "10+"),
          all(abs(n$value[n$age == 1] * 1000 - r$value) <= 6),
          all(is.na(y$se)), all(is.na(y$lwr)), all(is.na(y$upr)))
message("EBS pollock composition, survey timing/exclusion and output-group checks passed.")
m <- x[x$type == "M", ]
stopifnot(nrow(m) == 15, identical(m$age, 1:15),
          identical(m$value, c(.9, .45, rep(.3, 13))), all(is.na(m$year)), all(is.na(m$year_basis)))
message("Fixed age-specific M inputs verified.")
