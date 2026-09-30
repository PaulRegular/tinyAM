# File semantics checked against fishfollower/SAM at c6cfd035 (R/reading.R
# and R/conf.R). These readers do not execute R expressions from source files.
.sam_numbers <- function(lines, file) {
  tokens <- strsplit(paste(trimws(lines), collapse = " "), "[[:space:]]+")[[1]]
  tokens <- tokens[nzchar(tokens)]
  out <- suppressWarnings(as.numeric(tokens))
  if (any(is.na(out) & tokens != "NA") || any(is.infinite(out))) {
    cli::cli_abort("Non-numeric or infinite values in SAM file {.file {file}}.")
  }
  out
}

.sam_range <- function(x, file, allow_negative = FALSE) {
  if (length(x) != 2L || anyNA(x) || any(x != trunc(x)) ||
      x[1] > x[2] || (!allow_negative && any(x < 0))) {
    cli::cli_abort("Invalid year or age range in SAM file {.file {file}}.")
  }
  seq.int(x[1], x[2])
}

.read_sam_ices <- function(file) {
  lines <- readLines(file, warn = FALSE)
  if (length(lines) < 5L) cli::cli_abort("Incomplete SAM file {.file {file}}.")
  if (any(grepl("^\\s*#+'", lines))) {
    cli::cli_abort("SAM custom R-expression attributes are not supported in {.file {file}}; retain and inspect them separately.")
  }
  lines <- trimws(sub("#.*$", "", lines[-c(1, 2)]))
  lines <- lines[nzchar(lines)]
  if (!grepl("^[0-9]", lines[1])) return(.read_sam_surveys(lines, file))
  years <- .sam_range(.sam_numbers(lines[1], file), file)
  ages <- .sam_range(.sam_numbers(lines[2], file), file)
  code <- .sam_numbers(lines[3], file)
  if (length(code) != 1L || !code %in% c(1, 2, 3, 5)) {
    cli::cli_abort("Unsupported ICES table format in {.file {file}}; expected 1, 2, 3, or 5.")
  }
  values <- .sam_numbers(lines[-seq_len(3)], file)
  dims <- switch(as.character(code), `1` = c(length(years), length(ages)),
                 `2` = c(1L, length(ages)), `3` = c(1L, 1L),
                 `5` = c(length(years), 1L))
  if (length(values) != prod(dims) || length(lines) - 3L != dims[1]) {
    cli::cli_abort("Table dimensions do not match the header in {.file {file}}.")
  }
  if (any(lengths(strsplit(lines[-seq_len(3)], "[[:space:]]+")) != dims[2])) {
    cli::cli_abort("Inconsistent table row widths in {.file {file}}.")
  }
  x <- matrix(values, dims[1], dims[2], byrow = TRUE)
  x <- x[rep(seq_len(nrow(x)), length.out = length(years)),
         rep(seq_len(ncol(x)), length.out = length(ages)), drop = FALSE]
  dimnames(x) <- list(year = years, age = ages)
  x
}

.read_sam_surveys <- function(lines, file) {
  starts <- which(grepl("^[[:alpha:]]", lines))
  if (!length(starts) || starts[1] != 1L) {
    cli::cli_abort("Invalid survey header in {.file {file}}.")
  }
  ends <- c(starts[-1] - 1L, length(lines))
  out <- lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    if (length(block) < 5L) cli::cli_abort("Incomplete survey block in {.file {file}}.")
    years <- .sam_range(.sam_numbers(block[2], file), file)
    header <- .sam_numbers(block[3], file)
    ages <- .sam_range(.sam_numbers(block[4], file), file, allow_negative = TRUE)
    if (length(header) != 4L || anyNA(header) || header[3] < 0 ||
        header[4] > 1 || header[3] > header[4]) {
      cli::cli_abort("Invalid survey timing in {.file {file}}.")
    }
    rows <- block[-seq_len(4)]
    values <- .sam_numbers(rows, file)
    if (length(rows) != length(years) ||
        any(lengths(strsplit(rows, "[[:space:]]+")) != length(ages) + 1L)) {
      cli::cli_abort("Survey dimensions do not match the header in {.file {file}}.")
    }
    x <- matrix(values, nrow = length(years), byrow = TRUE)
    if (anyNA(x[, 1]) || any(x[, 1] <= 0)) {
      cli::cli_abort("Survey effort denominators must be positive in {.file {file}}.")
    }
    effort <- x[, 1]
    x <- x[, -1, drop = FALSE] / effort
    x[x < 0] <- NA_real_
    dimnames(x) <- list(year = years, age = ages)
    attr(x, "time") <- header[3:4]
    attr(x, "twofirst") <- header[1:2]
    attr(x, "effort") <- effort
    x
  })
  names(out) <- lines[starts]
  if (anyDuplicated(names(out))) cli::cli_abort("Duplicate survey names in {.file {file}}.")
  out
}

.read_sam_conf <- function(file) {
  lines <- trimws(sub("#.*$", "", readLines(file, warn = FALSE)))
  lines <- lines[nzchar(lines)]
  starts <- which(grepl("^\\$[[:alnum:]_]+$", lines))
  if (!length(starts) || starts[1] != 1L) cli::cli_abort("Invalid SAM configuration {.file {file}}.")
  nms <- substring(lines[starts], 2)
  if (anyDuplicated(nms)) cli::cli_abort("Duplicate configuration fields in {.file {file}}.")
  ends <- c(starts[-1] - 1L, length(lines))
  matrices <- c("keyLogFsta", "keyLogFpar", "keyQpow", "keyVarF", "keyVarObs",
                "keyCorObs", "keyParScaledYA", "predVarObsLink", "keyXtraSd",
                "keyCatchWeightMean", "keyCatchWeightObsVar")
  out <- lapply(seq_along(starts), function(i) {
    rows <- if (ends[i] > starts[i]) lines[(starts[i] + 1L):ends[i]] else character()
    if (nms[i] %in% c("obsCorStruct", "obsLikelihoodFlag")) {
      tokens <- scan(text = paste(rows, collapse = " "), what = character(), quiet = TRUE)
      allowed <- if (nms[i] == "obsCorStruct") c("ID", "AR", "US") else c("LN", "ALN")
      if (any(!tokens %in% c(allowed, "NA"))) cli::cli_abort("Unknown {nms[i]} label in {.file {file}}.")
      return(factor(replace(tokens, tokens == "NA", NA_character_), levels = allowed))
    }
    values <- .sam_numbers(rows, file)
    if (nms[i] %in% matrices) {
      widths <- lengths(strsplit(rows, "[[:space:]]+"))
      if (length(unique(widths)) > 1L) cli::cli_abort("Unequal configuration row widths in {.file {file}}.")
      return(if (!length(rows)) matrix(numeric(), 0L, 0L) else
        matrix(values, nrow = length(rows), byrow = TRUE))
    }
    values
  })
  stats::setNames(out, nms)
}

#' Read standard SAM assessment files without running SAM
#'
#' Read catch, survey, and biological inputs into a transparent source object
#' before converting them with [sam_to_tam_obs()]. No installed
#' \pkg{stockassessment} package is needed.
#'
#' @param path Directory containing standard SAM `.dat` files.
#' @param conf Configuration list, or path to a SAM `$field` configuration file.
#'   Relative configuration paths are resolved against `path`. `NULL` reads
#'   `model.cfg` if present, otherwise retains an empty configuration.
#' @return A list of `data`, `fleets`, `conf`, and `files`. `data` retains separate
#'   catch matrices, survey matrices, and biological inputs. File paths and MD5
#'   checksums are recorded for provenance.
#' @details
#' Supported ICES table codes are 1 (full), 2 (age row), 3 (scalar), and 5
#' (year column). Surveys use the standard effort-denominator format: observations
#' are divided by each row's denominator, negative values become missing, and
#' zeros are preserved. Timing endpoints and effort remain matrix attributes.
#'
#' Standard `cn.dat` or `cn_00001.dat` fleets are retained separately. Summed catch
#' fleets and biomass surveys are retained for inspection but their conversion
#' is unsupported. Custom executable R-expression attributes are rejected.
#' Configuration omissions are not filled with version-dependent defaults;
#' the assumption audit marks missing information as `not_checked`.
#' @seealso [sam_tam_assumptions()], [sam_reference()]
#' @export
read_sam_files <- function(path, conf = NULL) {
  if (!dir.exists(path)) cli::cli_abort("SAM directory {.file {path}} does not exist.")
  files <- list.files(path, pattern = "\\.dat$", full.names = TRUE)
  catch_files <- files[grepl("^cn(_[0-9]{5})?\\.dat$", basename(files))]
  if (!length(catch_files)) cli::cli_abort("No standard SAM catch files found in {.file {path}}.")
  read_many <- function(prefix) {
    fs <- files[grepl(paste0("^", prefix, "(_[0-9]{5})?\\.dat$"), basename(files))]
    stats::setNames(lapply(fs, .read_sam_ices), basename(fs))
  }
  data <- list(catch = read_many("cn"), surveys = .read_sam_ices(file.path(path, "survey.dat")))
  if (!is.list(data$surveys)) cli::cli_abort("survey.dat must contain survey blocks.")
  for (prefix in c("sw", "mo", "nm", "pf", "pm", "cw", "dw", "lw", "lf")) {
    fs <- files[grepl(paste0("^", prefix, "(_[0-9]{5})?\\.dat$"), basename(files))]
    if (length(fs)) data[[prefix]] <- if (length(fs) == 1L) .read_sam_ices(fs) else read_many(prefix)
  }
  data$sum_catch <- read_many("cn_sum")
  mats <- c(data$catch, data$surveys, data$sum_catch)
  types <- c(rep(0L, length(data$catch)),
             vapply(data$surveys, function(x) if (min(as.integer(colnames(x))) < 0) 3L else 2L, integer(1)),
             rep(7L, length(data$sum_catch)))
  times <- lapply(mats, function(x) if (is.null(attr(x, "time"))) c(0, 0) else attr(x, "time"))
  fleets <- data.frame(fleet_id = seq_along(mats), fleet_name = names(mats), fleet_type = types,
                       min_age = vapply(mats, function(x) min(as.integer(colnames(x))), integer(1)),
                       max_age = vapply(mats, function(x) max(as.integer(colnames(x))), integer(1)),
                       time_start = vapply(times, `[`, numeric(1), 1),
                       time_end = vapply(times, `[`, numeric(1), 2))
  fleets$samp_time <- (fleets$time_start + fleets$time_end) / 2
  if (is.null(conf)) conf <- if (file.exists(file.path(path, "model.cfg"))) "model.cfg" else list()
  if (is.character(conf) && length(conf) == 1L) {
    cfg_file <- if (file.exists(conf)) conf else file.path(path, conf)
    conf <- .read_sam_conf(cfg_file)
    files <- c(files, cfg_file)
  }
  if (!is.list(conf)) cli::cli_abort("conf must be a configuration list or file path.")
  list(data = data, fleets = fleets, conf = conf,
       files = data.frame(path = normalizePath(files, winslash = "/"), md5 = unname(tools::md5sum(files))))
}
