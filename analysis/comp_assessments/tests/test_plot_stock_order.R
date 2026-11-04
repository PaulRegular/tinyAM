source("analysis/comp_assessments/R/plot_stock_order.R")
source("analysis/comp_assessments/R/run_assessment.R")
root <- "analysis/comp_assessments"
database <- read_database()
assessments <- database$assessments
stocks <- database$stocks
ids <- intersect(assessment_translations("curated")$assessment_id, assessments$assessment_id)
info <- assessment_plot_stocks(assessments, stocks, ids)
stopifnot(nrow(info) == length(ids), !anyDuplicated(info$assessment_id),
  identical(info$assessment_id, assessment_plot_stocks(
    assessments[nrow(assessments):1, ], stocks[nrow(stocks):1, ], rev(ids)
  )$assessment_id), !is.unsorted(as.integer(info$basin)))

# Stocks of the same species remain adjacent within each ocean.
for (basin in levels(info$basin)) {
  species <- info$scientific_name[info$basin == basin]
  stopifnot(!anyDuplicated(rle(species)$values))
}
atlantic <- info[info$basin == "Northwest Atlantic", ]
stopifnot(max(which(atlantic$scientific_name == "Clupea harengus")) <
  min(which(atlantic$scientific_name == "Gadus morhua")))
pacific <- info[info$basin == "Northeast Pacific", ]
stopifnot(setequal(pacific$stock[grepl("^dfo_", pacific$stock_id)], c(
  "Strait of Georgia Pacific herring", "British Columbia Pacific ocean perch",
  "British Columbia Canary Rockfish", "British Columbia petrale sole")))

missing <- stocks
selected_stock <- match(info$stock_id[1L], stocks$stock_id)
missing$scientific_name[selected_stock] <- "Unclassified fish"
testthat::expect_error(assessment_plot_stocks(assessments, missing, ids), "missing plot taxonomy")
missing <- stocks
missing$ocean[selected_stock] <- "North Pacific"
testthat::expect_error(assessment_plot_stocks(assessments, missing, ids), "classified ocean basin")
testthat::expect_error(assessment_plot_stocks(assessments, stocks, "unknown"), "assessment metadata")
cat("Shared ocean/species ordering, shuffled inputs and missing metadata checks passed.\n")

# Check the aggregation used by plots without presenting old local exports as
# current results or requiring a fresh scientific batch during a unit test.
differences <- data.frame(metric = "ssb", year = 2000:2003, age = NA_integer_,
  source = c(100, 90, 80, 70), tinyAM = c(95, 92, 85, 75))
differences$percent_difference <- with(differences, 100 * (tinyAM - source) / source)
runs <- list(example = list(assessment_id = "example", translation_tier = "curated",
  fit = list(dat = list(years = 2000:2003, is_proj = c(FALSE, FALSE, FALSE, TRUE))),
  diagnostics = data.frame(assessment_id = "example"), differences = differences,
  summary = .assessment_comparison_summary(differences[1:3, ], "example")))
aggregated <- .assessment_result_tables(runs)
summary <- attr(aggregated, "comparison_summary")
series <- attr(aggregated, "population_timeseries")
stopifnot(summary$n == 3L, summary$terminal_year == 2002L,
  summary$median_absolute_percent_difference == 5,
  identical(series$year, 2000:2002))
stopifnot(!anyDuplicated(paste(series$assessment_id, series$metric, series$year)))
rows <- summary[summary$metric %in% c("ssb", "recruitment", "abundance", "F_bar") &
  summary$comparison_status %in% c("matched", "approximate") & summary$n > 0, ]
for (i in seq_len(nrow(rows))) {
  x <- series[series$assessment_id == rows$assessment_id[i] & series$metric == rows$metric[i] &
    is.finite(series$source) & is.finite(series$tinyAM), ]
  stopifnot(nrow(x) == rows$n[i], max(x$year) == rows$terminal_year[i],
    isTRUE(all.equal(stats::median(abs(x$percent_difference), na.rm = TRUE),
                     rows$median_absolute_percent_difference[i], tolerance = 1e-8)))
}
cat("Annual plot aggregation preserves historical pairs, counts, years and median differences.\n")
