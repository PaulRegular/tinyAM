source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
source(file.path(root, "R", "run_assessment.R"))
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "database_to_tam_obs.R"))

assessment_id <- "ices_plaice_english_channel_7d_2026"
stock_id <- "ck_ices_wgnssk_platessa_ns_7d"
database <- read_database()
source <- read_assessment(assessment_id, database)
inputs <- source$inputs
outputs <- source$outputs
cache <- "analysis/comp_assessments/source_cache/ices_plaice_7d_2026"
data_env <- new.env(parent = baseenv())
load(file.path(cache, "run/data.RData"), envir = data_env)
dat <- data_env$dat

stopifnot(
  nrow(source$assessment) == 1L,
  source$assessment$stock_id == stock_id,
  source$assessment$terminal_year == 2025L,
  source$assessment$estimate_terminal_year == 2026L,
  source$assessment$assumptions_status == "complete",
  source$assessment$inputs_status == "complete",
  source$assessment$outputs_status == "partial",
  source$stock$authority_stock_id == "ple.27.7d",
  source$stock$common_name == "Eastern English Channel plaice",
  identical(as.integer(dat$years), 1980:2025),
  identical(as.integer(dat$minAgePerFleet[[1L]]:dat$maxAgePerFleet[[1L]]), 1:7)
)

expected <- data.frame(
  year = dat$aux[, "year"],
  fleet = dat$aux[, "fleet"],
  age = dat$aux[, "age"],
  value = exp(dat$logobs)
)
catch <- subset(inputs, type == "catch" & measure == "numbers_at_age")
uk_bts <- subset(inputs, type == "index" & survey == "UK BTS")
fr_gfs <- subset(inputs, type == "index" & survey == "FR GFS")
expected_catch <- expected[expected$fleet == 1L, ]
expected_uk <- expected[expected$fleet == 2L & is.finite(expected$value), ]
expected_fr <- expected[expected$fleet == 3L & is.finite(expected$value), ]
key <- function(x) paste(x$year, x$age, sep = ":")
same_surface <- function(rows, surface, ages, years) {
  rows <- rows[order(as.integer(rows$year), as.integer(rows$age)), ]
  expected <- surface[cbind(match(as.integer(rows$year), years), match(as.integer(rows$age), ages))]
  expected_grid <- expand.grid(year = years, age = ages)
  !anyDuplicated(key(rows)) && setequal(key(rows), key(expected_grid)) &&
    isTRUE(all.equal(as.numeric(rows$value), as.numeric(expected), tolerance = 1e-12))
}

stopifnot(
  nrow(catch) == 322L, nrow(uk_bts) == 216L, nrow(fr_gfs) == 215L,
  setequal(key(catch), key(expected_catch)),
  isTRUE(all.equal(as.numeric(catch$value[match(key(expected_catch), key(catch))]), expected_catch$value, tolerance = 1e-12)),
  setequal(key(uk_bts), key(expected_uk)),
  setequal(key(fr_gfs), key(expected_fr)),
  all(catch$unit == "thousand fish"),
  all(c(uk_bts$unit, fr_gfs$unit) == "survey index (unit unresolved)"),
  all(as.numeric(uk_bts$sampling_time) == dat$sampleTimes[[2L]]),
  all(as.numeric(fr_gfs$sampling_time) == dat$sampleTimes[[3L]]),
  same_surface(subset(inputs, type == "weight" & measure == "weight_at_age"),
    dat$stockMeanWeight, 1:7, 1980:2025),
  same_surface(subset(inputs, type == "catch_weight" & measure == "weight_at_age"),
    dat$catchMeanWeight[, , 1L], 1:7, 1980:2025),
  same_surface(subset(inputs, type == "catch_weight" & measure == "landings_weight_at_age"),
    dat$landMeanWeight[, , 1L], 1:7, 1980:2025),
  same_surface(subset(inputs, type == "catch_weight" & measure == "discard_weight_at_age"),
    dat$disMeanWeight[, , 1L], 1:7, 1980:2025)
)

static_checks <- list(
  maturity_at_age = dat$propMat[1L, ],
  natural_mortality_at_age = dat$natMor[1L, ],
  landings_fraction_at_age = dat$landFrac[1L, , 1L],
  fraction_F_before_spawning = dat$propF[1L, , 1L],
  fraction_M_before_spawning = dat$propM[1L, ]
)
for (measure in names(static_checks)) {
  rows <- inputs[!is.na(inputs$measure) & inputs$measure == measure, , drop = FALSE]
  rows <- rows[order(as.integer(rows$age)), ]
  stopifnot(
    nrow(rows) == 7L,
    all(is.na(rows$year)),
    isTRUE(all.equal(as.numeric(rows$value), as.numeric(static_checks[[measure]]), tolerance = 1e-12))
  )
}

sag <- outputs[outputs$source_type == "official_machine_readable", ]
stopifnot(
  nrow(sag) == 140L,
  setequal(unique(as.integer(sag$year[sag$measure == "SSB"])), 1980:2026),
  setequal(unique(as.integer(sag$year[sag$measure == "recruitment"])), 1980:2026),
  setequal(unique(as.integer(sag$year[sag$measure == "Fbar"])), 1980:2025),
  all(sag$source_reference == "https://standardgraphs.ices.dk/ViewSourceData.aspx?key=22484"),
  all(as.numeric(sag$lwr[!is.na(sag$lwr)]) <= as.numeric(sag$value[!is.na(sag$lwr)])),
  all(as.numeric(sag$upr[!is.na(sag$upr)]) >= as.numeric(sag$value[!is.na(sag$upr)])),
  all(sag$age_group[sag$measure == "Fbar"] == "3-6"),
  !any(sag$measure %in% c("numbers_at_age", "fishing_mortality_at_age", "total_abundance"))
)

catch_rows <- expected[expected$fleet == 1L, ]
catch_weights <- subset(inputs, type == "catch_weight" & measure == "weight_at_age")
weights <- catch_weights$value[match(key(catch_rows), key(catch_weights))]
catch_tonnes <- tapply(catch_rows$value * as.numeric(weights), catch_rows$year, sum)
xml <- paste(readLines(file.path(cache, "ices_sag_2026_key_22484.xml"), warn = FALSE), collapse = "\n")
blocks <- regmatches(xml, gregexpr("(?s)<Fish_Data>.*?</Fish_Data>", xml, perl = TRUE))[[1L]]
read_field <- function(node, tag) {
  x <- regmatches(node, regexec(paste0("(?s)<", tag, ">(.*?)</", tag, ">"), node, perl = TRUE))[[1L]]
  if (length(x) < 2L || !nzchar(x[[2L]])) return(NA_real_)
  as.numeric(x[[2L]])
}
sag_catches <- do.call(rbind, lapply(blocks, function(node) {
  c(
    year = read_field(node, "Year"),
    catch = read_field(node, "Catches"),
    landings = read_field(node, "Landings"),
    discards = read_field(node, "Discards")
  )
}))
sag_catches <- as.data.frame(sag_catches)
sag_catches$year <- as.integer(sag_catches$year)
sag_catches[] <- lapply(sag_catches, function(x) as.numeric(as.character(x)))
sag_catches$year <- as.integer(sag_catches$year)
sag_catches$catch[is.na(sag_catches$catch)] <-
  sag_catches$landings[is.na(sag_catches$catch)] +
  sag_catches$discards[is.na(sag_catches$catch)]
sag_catches <- sag_catches[sag_catches$year %in% 1980:2025, ]
stopifnot(
  identical(sort(sag_catches$year), 1980:2025),
  max(abs(catch_tonnes[as.character(sag_catches$year)] - sag_catches$catch) /
    pmax(1, sag_catches$catch)) < 0.01
)

cat("Eastern English Channel plaice source data validated.\n")
