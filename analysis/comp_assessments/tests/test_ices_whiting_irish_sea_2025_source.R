source("analysis/comp_assessments/tests/helper_stock.R")
source("analysis/comp_assessments/R/read_database.R")

database <- read_database()
benchmark_id <- "ices_whiting_irish_sea_2025_benchmark"
summary_id <- "ices_whiting_irish_sea_2025_sag_summary"
benchmark <- read_assessment(benchmark_id, database)
summary <- read_assessment(summary_id, database)
inputs <- benchmark$inputs
outputs <- benchmark$outputs

rows <- function(x, type, measure) {
  x[x$type == type & x$measure == measure, , drop = FALSE]
}
series_matrix <- function(x, years, ages) {
  x <- x[order(x$age, x$year), , drop = FALSE]
  matrix(as.numeric(x$value), nrow = length(years), ncol = length(ages))
}

catch <- rows(inputs, "catch", "numbers_at_age")
index <- rows(inputs, "index", "numbers_at_age")
weights <- rows(inputs, "index", "relative_precision_weight")
n_age <- rows(outputs, "population", "numbers_at_age")
f_age <- rows(outputs, "mortality", "fishing_mortality_at_age")
pred_catch <- rows(outputs, "catch", "predicted_catch")
pred_index <- rows(outputs, "index", "predicted_index")
q <- rows(outputs, "catchability", "q")
stock_weight <- rows(inputs, "weight", "weight_at_age")
maturity <- rows(inputs, "maturity", "maturity_at_age")
landings_fraction <- rows(inputs, "catch", "landings_fraction_at_age")
m_fixed <- rows(inputs, "M", "natural_mortality_at_age")

stopifnot(
  nrow(benchmark$stocks) == 1L,
  benchmark$stocks$common_name == "Whiting",
  benchmark$stocks$authority_stock_id == "whg.27.7a",
  benchmark$stocks$area == "ICES Division 7.a (Irish Sea)",
  benchmark$stocks$ocean == "Northeast Atlantic",
  isTRUE(benchmark$assessment$is_current), isTRUE(benchmark$assessment$is_applied),
  benchmark$assessment$inputs_status == "complete",
  benchmark$assessment$outputs_status == "complete",
  benchmark$assessment$terminal_year == 2023L,
  benchmark$assessment$estimate_terminal_year == 2023L,
  identical(sort(unique(as.integer(catch$year))), 1980:2023),
  nrow(catch) == 308L, setequal(unique(catch$age), 0:6),
  all(catch$unit == "thousand fish"),
  nrow(index) == 444L,
  sum(index$survey == "NIGFS-Q1") == 192L,
  sum(index$survey == "NIGFS-Q4") == 224L,
  sum(index$survey == "NIMIK") == 28L,
  setequal(unique(index$survey), c("NIGFS-Q1", "NIGFS-Q4", "NIMIK")),
  all(index$unit == "native survey-index units"),
  setequal(index$year[index$survey == "NIMIK"], setdiff(1994:2023, c(2020, 2022))),
  all(index$sampling_time[index$survey == "NIGFS-Q1"] == 0.23),
  all(index$sampling_time[index$survey == "NIGFS-Q4"] == 0.855),
  all(index$sampling_time[index$survey == "NIMIK"] == 0.48),
  nrow(weights) == 416L,
  setequal(unique(weights$survey), c("NIGFS-Q1", "NIGFS-Q4")),
  all(weights$value > 0),
  nrow(stock_weight) == 44L * 7L,
  nrow(maturity) == 44L * 7L,
  nrow(m_fixed) == 44L * 7L,
  setequal(unique(m_fixed$value), c(2.20, 0.98, 0.75, 0.67, 0.63, 0.62, 0.61)),
  nrow(landings_fraction) == 44L * 7L,
  sum(landings_fraction$value > 1) == 7L,
  isTRUE(all.equal(max(landings_fraction$value), 1.0237, tolerance = 1e-4)),
  all(grepl("above one", landings_fraction$notes)),
  nrow(n_age) == 44L * 7L,
  nrow(f_age) == 44L * 7L,
  all(is.finite(n_age$value) & is.finite(n_age$se)),
  all(is.finite(f_age$value) & is.finite(f_age$se)),
  nrow(pred_catch) == 308L,
  nrow(pred_index) == 446L,
  nrow(q) == 14L,
  nrow(rows(outputs, "recruitment", "recruitment")) == 44L,
  nrow(rows(outputs, "biomass", "SSB")) == 44L,
  nrow(rows(outputs, "biomass", "total_biomass")) == 44L,
  nrow(rows(outputs, "mortality", "Fbar")) == 44L,
  all(rows(outputs, "mortality", "Fbar")$age_group == "1-3"),
  max(q$value[q$survey == "NIGFS-Q1" & q$age %in% 2:5]) - min(q$value[q$survey == "NIGFS-Q1" & q$age %in% 2:5]) < 1e-12,
  max(q$value[q$survey == "NIGFS-Q4" & q$age %in% 1:3]) - min(q$value[q$survey == "NIGFS-Q4" & q$age %in% 1:3]) < 1e-12,
  !anyDuplicated(paste(outputs$measure, outputs$year, outputs$age, outputs$survey)),
  isTRUE(summary$assessment$is_applied),
  !isTRUE(summary$assessment$is_current),
  summary$assessment$inputs_status == "partial",
  nrow(rows(summary$outputs, "biomass", "SSB")) == 46L,
  nrow(rows(summary$outputs, "recruitment", "recruitment")) == 46L,
  nrow(rows(summary$outputs, "mortality", "F")) == 45L,
  !2025L %in% rows(summary$outputs, "mortality", "F")$year,
  isTRUE(all.equal(
    rows(summary$outputs, "biomass", "SSB")$value[
      rows(summary$outputs, "biomass", "SSB")$year == 2023L
    ], 1534, tolerance = 1e-8
  )),
  abs(rows(outputs, "biomass", "SSB")$value[
    rows(outputs, "biomass", "SSB")$year == 2023L
  ] - 1590.003) < 1e-3
)

years <- 1980:2023
ages <- 0:6
n <- series_matrix(n_age, years, ages)
f <- series_matrix(f_age, years, ages)
stock_w <- series_matrix(stock_weight, years, ages)
mat <- series_matrix(maturity, years, ages)
ssb <- rows(outputs, "biomass", "SSB")
tsb <- rows(outputs, "biomass", "total_biomass")
fbar <- rows(outputs, "mortality", "Fbar")
stopifnot(
  isTRUE(all.equal(rowSums(n * stock_w * mat), as.numeric(ssb$value[match(years, ssb$year)]), tolerance = 1e-8)),
  isTRUE(all.equal(rowSums(n * stock_w), as.numeric(tsb$value[match(years, tsb$year)]), tolerance = 1e-8)),
  isTRUE(all.equal(rowMeans(f[, match(1:3, ages), drop = FALSE]), as.numeric(fbar$value[match(years, fbar$year)]), tolerance = 1e-8))
)

cat("Irish Sea whiting benchmark inputs, native outputs, and separate SAG summary passed validation.\n")
