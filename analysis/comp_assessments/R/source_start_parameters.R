# Called only by the conditional retry in fit_default_first().
source_start_parameters <- function(source, dat) {
  parameters <- tinyAM::make_par(dat)
  years <- dat$years[!dat$is_proj]
  ages <- dat$ages
  surface <- function(type, measure, numbers = FALSE) {
    rows <- source$outputs[source$outputs$type == type &
      source$outputs$measure == measure & source$outputs$year %in% years &
      source$outputs$age %in% ages, , drop = FALSE]
    if (anyDuplicated(paste(rows$year, rows$age))) {
      cli::cli_abort("Source starting values require a unique year-age surface; pool source dimensions explicitly first.")
    }
    value <- rows$value
    if (numbers) value <- value * tinyAM:::.unit_scale(rows$unit, "numbers") / 1000
    answer <- matrix(NA_real_, length(years), length(ages),
      dimnames = list(as.character(years), as.character(ages)))
    answer[cbind(match(rows$year, years), match(rows$age, ages))] <- value
    answer[!is.finite(answer) | answer <= 0] <- NA_real_
    answer
  }
  abundance <- surface("population", "numbers_at_age", numbers = TRUE)
  mortality <- surface("mortality", "fishing_mortality_at_age")
  changed <- FALSE
  update <- function(name, values) {
    keep <- is.finite(values)
    if (any(keep)) {
      parameters[[name]][keep] <<- log(values[keep])
      changed <<- TRUE
    }
  }
  update("log_r0", abundance[1L, 1L])
  if (!is.null(parameters$log_n0)) {
    update("log_n0", abundance[1L, names(parameters$log_n0)])
  }
  for (name in c("log_r", "log_r_init")) {
    if (!is.null(parameters[[name]])) {
      update(name, abundance[names(parameters[[name]]), 1L])
    }
  }
  if (!is.null(parameters$log_n)) {
    update("log_n", abundance[rownames(parameters$log_n),
      colnames(parameters$log_n), drop = FALSE])
  }
  if (!is.null(parameters$log_f)) {
    if (!all(colnames(parameters$log_f) %in% colnames(mortality))) {
      cli::cli_abort("Source F starts need an explicit mapping for blocked or fleet-specific states.")
    }
    update("log_f", mortality[rownames(parameters$log_f),
      colnames(parameters$log_f), drop = FALSE])
  }
  if (!changed) cli::cli_abort("No positive, aligned source N or F starting values are available.")
  parameters
}
