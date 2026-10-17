source("analysis/comp_assessments/R/plot_stock_order.R")
root <- "analysis/comp_assessments"
assessments <- read.csv(file.path(root, "database/assessments.csv"))
stocks <- read.csv(file.path(root, "database/stocks.csv"))
ids <- assessments$assessment_id[assessments$is_current]
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
stopifnot(all(pacific$stock[grepl("^dfo_", pacific$stock_id)] %in% c(
  "Strait of Georgia Pacific herring", "British Columbia Pacific ocean perch",
  "British Columbia Canary Rockfish")))

missing <- stocks
missing$scientific_name[1] <- "Unclassified fish"
testthat::expect_error(assessment_plot_stocks(assessments, missing, ids), "missing plot taxonomy")
missing <- stocks
missing$ocean[1] <- "North Pacific"
testthat::expect_error(assessment_plot_stocks(assessments, missing, ids), "classified ocean basin")
testthat::expect_error(assessment_plot_stocks(assessments, stocks, "unknown"), "assessment metadata")
cat("Shared ocean/species ordering, shuffled inputs and missing metadata checks passed.\n")

# Annual plot data must contain the pairs represented in the current summaries.
summary <- read.csv(file.path(root, "results/comparison_summary.csv"))
series <- read.csv(file.path(root, "results/population_timeseries.csv"))
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
cat("Annual plot data agree with the current comparison counts, years and median differences.\n")
