read_committed_database <- function() {
  repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  git <- function(args, stdout = TRUE) {
    system2("git", c("-c", paste0("safe.directory=", repo), args),
            stdout = stdout, stderr = TRUE)
  }
  commit <- git(c("rev-parse", "--verify", "HEAD^{commit}"))
  if (!is.null(attr(commit, "status"))) {
    stop("Could not resolve the committed database revision.", call. = FALSE)
  }
  commit <- commit[[1]]

  read_table <- function(name) {
    path <- file.path("analysis", "comp_assessments", "database", name)
    tempfile_path <- tempfile(pattern = "tinyAM_committed_database_",
                              fileext = ".csv")
    on.exit(unlink(tempfile_path), add = TRUE)
    status <- git(c("show", paste0(commit, ":", path)), stdout = tempfile_path)
    exit_status <- if (is.numeric(status) && length(status) == 1L) status else
      attr(status, "status")
    if (!is.null(exit_status) && exit_status != 0L ||
        !file.exists(tempfile_path) || file.info(tempfile_path)$size == 0) {
      stop("Could not read committed database table: ", name, call. = FALSE)
    }
    read.csv(tempfile_path, stringsAsFactors = FALSE,
             na.strings = c("", "NA"), check.names = FALSE)
  }

  tables <- lapply(c("stocks.csv", "assessments.csv", "assumptions.csv",
                     "inputs.csv", "outputs.csv"), read_table)
  names(tables) <- c("stocks", "assessments", "assumptions", "inputs", "outputs")
  tables$commit <- commit
  tables
}

read_committed_assessment <- function(assessment_id, database = NULL) {
  if (length(assessment_id) != 1L || is.na(assessment_id) ||
      !nzchar(assessment_id)) {
    stop("assessment_id must be one non-empty value.", call. = FALSE)
  }
  if (is.null(database)) database <- read_committed_database()
  required <- c("stocks", "assessments", "assumptions", "inputs", "outputs", "commit")
  if (!all(required %in% names(database))) {
    stop("database must come from read_committed_database().", call. = FALSE)
  }
  tables <- database
  assessment <- tables$assessments[
    !is.na(tables$assessments$assessment_id) &
      tables$assessments$assessment_id == assessment_id, , drop = FALSE]
  if (nrow(assessment) != 1L) {
    stop("Committed database must contain exactly one assessment_id: ",
         assessment_id, call. = FALSE)
  }
  stock_id <- assessment$stock_id[[1]]
  tables$stocks <- tables$stocks[
    tables$stocks$stock_id == stock_id, , drop = FALSE]
  tables$assumptions <- tables$assumptions[
    !is.na(tables$assumptions$assessment_id) &
      tables$assumptions$assessment_id == assessment_id, , drop = FALSE]
  tables$inputs <- tables$inputs[
    !is.na(tables$inputs$assessment_id) &
      tables$inputs$assessment_id == assessment_id, , drop = FALSE]
  tables$outputs <- tables$outputs[
    !is.na(tables$outputs$assessment_id) &
      tables$outputs$assessment_id == assessment_id, , drop = FALSE]
  tables$assessment <- assessment
  tables
}