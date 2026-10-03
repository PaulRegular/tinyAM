source("analysis/comp_assessments/scripts/database/002_validate_database.R")
id <- "nefsc_haddock_georges_bank_2026"
x <- inputs[inputs$assessment_id == id, ]
y <- outputs[outputs$assessment_id == id, ]
a <- assessments[assessments$assessment_id == id, ]
cache <- "analysis/comp_assessments/source_cache/nefsc_haddock_georges_bank_2026"
text <- readLines(file.path(cache, "2026_assessment_text.txt"), encoding = "UTF-8", warn = FALSE)
review <- paste(readLines(file.path(cache, "2026_peer_review_text.txt"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
source_values <- function(label) {
  line <- text[startsWith(text, label)]
  stopifnot(length(line) == 1)
  line <- substring(line, nchar(label) + 1)
  values <- regmatches(line, gregexpr("[0-9][0-9,]*(\\.[0-9]+)?", line))[[1]]
  values <- as.numeric(gsub(",", "", values, fixed = TRUE))
  stopifnot(length(values) == 8)
  values
}
stopifnot(nrow(x) == 17, nrow(y) == 24, nrow(a) == 1,
          a$terminal_year == 2025, a$assessment_year == 2026,
          a$model_family == "WHAM", a$framework_year == 2022,
          a$inputs_status == "partial", a$outputs_status == "partial", a$assumptions_status == "partial")
catch <- x[x$type == "catch", ]
stopifnot(identical(catch$year, 2018:2025),
          identical(catch$value, source_values("Catch for Assessment")),
          all(catch$measure == "total_biomass"), all(catch$unit == "t"))
for (measure in c("SSB", "recruitment")) {
  z <- y[y$measure == measure, ]
  label <- if (measure == "SSB") "Spawning Stock Biomass" else "Recruits (age 1)"
  stopifnot(identical(z$year, 2018:2025), identical(z$value, source_values(label)))
}
f <- y[y$measure == "Fbar", ]
line <- text[grepl("^F", text) & grepl("0.24 0.33 0.42", text, fixed = TRUE)]
stopifnot(length(line) == 1)
values <- as.numeric(regmatches(line, gregexpr("0\\.[0-9]+", line))[[1]])
stopifnot(identical(f$value, values), all(f$age_group == "5–7"),
          all(y$age[y$measure == "recruitment"] == 1), all(y$year <= a$terminal_year),
          grepl("29,037", review, fixed = TRUE), grepl("2DAR1", review, fixed = TRUE))
m <- x[x$type == "M", ]
stopifnot(identical(m$age, 1:9), all(m$value == .2), all(m$year == 1931),
          !any(x$type == "index"), !any(x$measure == "numbers_at_age"),
          all(is.na(y$se)))
message("Georges Bank haddock published catch, SSB, recruitment, Fbar and constant M verified; missing streams remain explicit.")
