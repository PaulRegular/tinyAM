source("analysis/comp_assessments/R/read_database.R")
db <- read_database()
basins <- c("Northwest Atlantic", "Northeast Atlantic",
            "Northeast Pacific", "Northwest Pacific")
stopifnot(all(db$stocks$ocean %in% basins))
stopifnot(all(db$stocks$ocean[db$stocks$authority == "NOAA-AFSC"] ==
                "Northeast Pacific"))
stopifnot(db$stocks$ocean[db$stocks$stock_id == "ices_haddock_iceland"] ==
            "Northeast Atlantic")
stopifnot(all(db$stocks$ocean[grepl("^nefsc_", db$stocks$stock_id)] ==
                "Northwest Atlantic"))

# These figures use geography, not authority, so Canadian Pacific stocks
# cannot be assigned to the Atlantic by their identifier prefix.
plot_code <- readLines("analysis/comp_assessments/scripts/translation/trend_plots.R")
stopifnot(any(grepl("basin = ocean", plot_code, fixed = TRUE)))
heatmap_code <- readLines("analysis/comp_assessments/scripts/translation/stats_plots.R")
stopifnot(any(grepl("basin = ocean", heatmap_code, fixed = TRUE)))
stopifnot(!any(grepl('"North Pacific"', heatmap_code, fixed = TRUE)))
stopifnot(!any(grepl("list(Accepted", readLines(
  "analysis/comp_assessments/R/run_assessment.R"
), fixed = TRUE)))
cat("Regional basin and Official reporting conventions passed.\n")
