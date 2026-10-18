pkgload::load_all(".", quiet = TRUE, helpers = FALSE)

.test_stock <- function(source_record, ..., script = NULL) {
  e <- new.env(parent = globalenv())
  e$source <- source_record
  e$do_fit <- FALSE
  e$silent <- TRUE
  list2env(list(...), e)
  path <- file.path(
    "analysis/comp_assessments/scripts/translation/stocks",
    paste0(source_record$assessment$assessment_id[[1]], ".R")
  )
  if (!is.null(script)) path <- script
  sys.source(path, e)
  calls <- list()
  walk <- function(x) {
    if (missing(x) || !is.call(x)) {
      return(invisible(NULL))
    }
    if (identical(x[[1]], quote(tinyAM::fit_tam))) calls[[length(calls) + 1L]] <<- x
    for (z in as.list(x)[-1]) walk(z)
  }
  for (x in parse(path)) walk(x)
  call <- calls[[length(calls)]]
  call[[1]] <- quote(tinyAM::prepare_tam)
  call <- as.call(Filter(function(x) !identical(x, quote(...)), as.list(call)))
  call$silent <- call$start_par <- NULL
  e$dat <- eval(call, e)
  out <- as.list(e)
  out$fit_call <- calls[[length(calls)]]
  out$warm_call <- if (length(calls) > 1L) calls[[1L]] else NULL
  settings <- as.list(call)[intersect(
    names(call),
    c("N_settings", "rec_settings", "F_settings", "M_settings", "catch_settings", "index_settings")
  )]
  out$settings <- lapply(settings, eval, envir = e)
  out
}
