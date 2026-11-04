# Run manually from the tinyAM repository root after updating SAIO.
# This copies committed CSV files only; it never fetches or changes SAIO.
saio_repo <- "../saio" # Edit this path if the two repositories are not siblings.
source("analysis/comp_assessments/R/update_database_cache.R")
update_database_cache(saio_repo)
