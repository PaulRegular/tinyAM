#' Northern cod example data for TAM
#'
#' @description
#' A compact, tidy-ish list of inputs prepared for the Tiny Assessment Model (TAM),
#' derived from data used in the assessment of Northern cod, NAFO Divisions 2J3KL.
#' This is intended for vignettes, examples, and tests.
#'
#' @format
#' A named list of four data frames:
#'
#' - **catch**: commercial catch-at-age observations and blocking variables.
#'   Columns:
#'   - `year`: numeric calendar year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric observed catch in numbers (may include `NA`).
#'   - `sd_obs_block`: character; observation-variance block label (`"catch"`).
#'   - `F_y_block`: factor; year blocks for fishing-mortality mean structure
#'     (e.g., `"1954-1991"`, `"1992-1997"`, …).
#'   - `F_a_block`: factor; age blocks for fishing-mortality mean structure
#'     (e.g., `"0-1"`, `"2-3"`, …).
#'
#' - **index**: survey (index) observations and blocking variables.
#'   Columns:
#'   - `year`: integer year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric survey index-at-age (may include `NA`).
#'   - `samp_time`: numeric sampling time within the year (fraction of a year; e.g., `0.8`).
#'   - `q_block`: factor; age blocks for catchability.
#'   - `sd_obs_block`: character; observation-variance block label (`"index"`).
#'
#' - **weight**: weight-at-age and mortality assumption used for \eqn{M} mean structure.
#'   Columns:
#'   - `year`: integer year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric weight-at-age (same units as source; e.g., kg).
#'   - `M_assumption`: numeric assumed natural mortality level by age–year
#'     (used when constructing `M` mean via formula).
#'
#' - **maturity**: maturity-at-age (proportion mature).
#'   Columns:
#'   - `year`: numeric year.
#'   - `age`: numeric age (years).
#'   - `obs`: numeric proportion mature (0–1).
#'
#' @details
#' The object was assembled by:
#'
#' 1. Converting NCAM arrays to tidy data frames.
#' 2. Creating blocking variables for `F` by year (`F_y_block`) and age
#'    (`F_a_block`) using [cut_years()] and [cut_ages()].
#' 3. Joining weight-at-age with a natural mortality assumption column
#'    (`M_assumption`) derived from `NCAM::inputs$nm`.
#'
#' Exact dimensions may change if the source data evolve. At the time of inclusion,
#' the elements contained approximately:
#'
#' - `catch`: 1,065 rows × 6 columns
#' - `index`: 1,065 rows × 6 columns
#' - `weight`: 1,065 rows × 4 columns
#' - `maturity`: 1,215 rows × 3 columns
#'
#' @source
#' Derived from the \pkg{NCAM} example inputs
#' (`NCAM::inputs`). See the NCAM package for provenance and licensing.
#'
#' @examples
#' data(northern_cod_data)
#' names(northern_cod_data)
#'
#' # Peek at catch blocks
#' head(northern_cod_data$catch)
#' with(northern_cod_data$catch, table(F_y_block, F_a_block))
#'
#' # Build a minimal TAM dataset
#' dat <- make_dat(
#'   obs = northern_cod_data,
#'   years = 1983:2024,
#'   ages = 2:14,
#'   N_settings = list(process = "iid", init_N0 = FALSE),
#'   F_settings = list(process = "approx_rw", mu_form = NULL),
#'   M_settings = list(process = "off", mu_form = NULL,
#'                     assumption = ~ I(0.2), age_breaks = NULL),
#'   obs_settings = list(sd_form = ~ sd_obs_block, q_form = ~ q_block)
#' )
#'
#' @docType data
#' @keywords datasets
#' @name northern_cod_data
#' @usage data(northern_cod_data)
"northern_cod_data"
