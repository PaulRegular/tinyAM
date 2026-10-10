.test_stock <- function(source_record, ...) {
  e <- new.env(parent = globalenv())
  e$source <- source_record
  e$do_fit <- FALSE
  e$silent <- TRUE
  list2env(list(...), e)
  path <- file.path('analysis/comp_assessments/scripts/translation/stocks',
                    paste0(source_record$assessment$assessment_id[[1]], '.R'))
  sys.source(path, e)
  calls <- list()
  walk <- function(x) {
    if (missing(x) || !is.call(x)) return(invisible(NULL))
    if (identical(x[[1]], quote(tinyAM::fit_tam))) calls[[length(calls)+1L]] <<- x
    for (z in as.list(x)[-1]) walk(z)
  }
  for (x in parse(path)) walk(x)
  call <- calls[[length(calls)]]
  call[[1]] <- quote(tinyAM::prepare_tam)
  call$silent <- call$start_par <- NULL
  e$dat <- eval(call, e)
  out <- as.list(e)
  out$fit_call <- calls[[length(calls)]]
  out$warm_call <- if (length(calls) > 1L) calls[[1L]] else NULL
  out$settings <- out$dat[c('N_settings','F_settings','M_settings','catch_settings','index_settings')]
  out
}


