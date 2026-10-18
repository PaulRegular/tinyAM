# Run explicitly to rebuild the cheatsheet examples; PDF rendering only reads them.
# Rscript cheatsheet/examples.R
# rmarkdown::render("cheatsheet/tinyAM_cheatsheet.Rmd")
# The PDF needs XeLaTeX and Bahnschrift. Cache files stay local; share the PDF
# together with tinyAM_dashboard.html when the interactive example is useful.
script <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(script)) sub("^--file=", "", script[1]) else sys.frame(1)$ofile
script <- normalizePath(script, winslash = "/", mustWork = TRUE)
repo <- normalizePath(dirname(dirname(script)), winslash = "/", mustWork = TRUE)
pkgload::load_all(repo, quiet = TRUE)

cache <- file.path(repo, "cheatsheet", "cache")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
revision <- system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)
lines <- readLines(script, warn = FALSE)
markers <- grep("^# ---- .* ----$", lines)
names(markers) <- sub("^# ---- (.*) ----$", "\\1", lines[markers])
fitting_chunks <- c("inputs", "baseline", "N_iid", "M_ar1")
model_code <- unlist(lapply(fitting_chunks, function(name) {
  from <- markers[[name]]
  to <- markers[which(markers > from)[1]]
  lines[seq.int(from + 1L, to - 1L)]
}))
writeLines(model_code, file.path(cache, "model-code.txt"))
provenance <- list(revision = revision, version = as.character(packageVersion("tinyAM")),
                   example_md5 = unname(tools::md5sum(file.path(cache, "model-code.txt"))),
                   generated = Sys.time())
stage <- Sys.getenv("TINYAM_CHEATSHEET_STAGE", "all")

# Replace these source tables with your own data. The included cod data make
# this complete example reproducible without external assessment files.
source <- lapply(cod_obs, function(x) {
  x <- x[x$year %in% 1983:2024 & x$age %in% 2:20, , drop = FALSE]
  keys <- intersect(c("year", "age", "survey"), names(x))
  x[do.call(order, x[keys]), , drop = FALSE]
})
catch_data <- source$catch
index_data <- source$index
weight_data <- source$weight
maturity_data <- source$maturity
population_data <- source$pop_covariates

# ---- inputs ----
dat <- make_tam_data(
  catch = catch_data, index = index_data,
  weight = weight_data, maturity = maturity_data,
  pop_covariates = population_data
)
# A plain list also works: check_tam_data(dat)
# Optional class conversion: as_tam_data(dat)
# ---- end ----
saveRDS(dat, file.path(cache, "data.rds"))
provenance$data_md5 <- unname(tools::md5sum(file.path(cache, "data.rds")))
saveRDS(list(dat = dat, provenance = provenance), file.path(cache, "examples.rds"))

records <- list()
cached <- function(name, code) {
  file <- file.path(cache, paste0(name, ".rds"))
  if (file.exists(file)) {
    record <- readRDS(file)
    same <- identical(record$provenance[c("revision", "example_md5", "data_md5")],
                      provenance[c("revision", "example_md5", "data_md5")])
    if (same) {
      records[[name]] <<- record
      cli::cli_inform("{name}: reused cached result")
      return(record$value)
    }
  }
  cli::cli_inform("{name}: started {format(Sys.time(), '%H:%M:%S')}")
  warnings <- character()
  error <- NULL
  elapsed <- system.time(value <- withCallingHandlers(
    tryCatch(code(), error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }))[["elapsed"]]
  record <- list(value = value, warnings = warnings, error = error,
                 elapsed_seconds = elapsed, provenance = provenance)
  saveRDS(record, file)
  records[[name]] <<- record
  suffix <- if (is.null(error)) "" else paste0("; error: ", error)
  cli::cli_inform("{name}: finished in {round(elapsed, 1)} s{suffix}")
  value
}

if (stage != "data") {
  base <- cached("fit_Base", function() {
# ---- baseline ----
base <- fit_tam(
  dat, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "rw"),
  F_settings = list(process = "ar1", mu_form = ~ 1),
  M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(sd_form = ~ 1, q_form = ~ mono(q_block))
)
# ---- end ----
    base
  })
  N_iid <- cached("fit_N_iid", function() {
# ---- N_iid ----
N_iid <- update(base,
  N_settings = list(process = "iid"))
# ---- end ----
    N_iid
  })
  M_ar1 <- cached("fit_M_ar1", function() {
# ---- M_ar1 ----
M_ar1 <- update(base,
  M_settings = list(process = "ar1"))
# ---- end ----
    M_ar1
  })
# ---- models ----
fits <- list(Base = base, N_iid = N_iid, M_ar1 = M_ar1)
# ---- end ----
  model_names <- names(fits)
  usable <- vapply(fits, inherits, logical(1), "tam_fit")
  fits <- fits[usable]
  checks <- lapply(fits, check_tam)
  attempts <- lapply(records, function(x) x[setdiff(names(x), "value")])
  saveRDS(list(dat = dat, fits = fits, checks = checks, records = attempts,
               provenance = provenance), file.path(cache, "examples.rds"))
}

if (stage == "all") {
  # Cache each model separately so a long validation run can be resumed.
  future::plan(future::sequential)
  hind <- lapply(model_names, function(name) {
    if (is.null(fits[[name]])) return(NULL)
    cached(paste0("hind_", name), function() {
      fit_hindcast(fits[[name]], folds = 20, progress = FALSE)
    })
  })
  names(hind) <- model_names
# ---- hindcast ----
# hind <- lapply(fits, fit_hindcast, folds = 20)
# rmse <- vapply(hind, `[[`, numeric(1), "hindcast_rmse")
# ---- end ----

  observations <- stack_list(base$obs_pred[c("catch", "index")], label = "type")
  if (!"survey" %in% names(observations)) observations$survey <- NA_character_
  keys <- c("type", "survey", "year", "age")
  key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
  observed <- observations[observations$year %in% 2005:2024 &
                            is.finite(observations$obs) & observations$obs > 0, , drop = FALSE]
  observed$is_proj <- FALSE
  expected <- key(observed)
  scores <- lapply(names(hind), function(name) {
    result <- hind[[name]]
    fit_converged <- isTRUE(checks[[name]]$is_converged)
    forecasts <- if (is.null(result)) NULL else {
      d <- stack_list(result$obs_pred[c("catch", "index")], label = "type")
      history <- unique(d[!d$is_proj & is.finite(d$obs) & d$obs > 0, keys, drop = FALSE])
      x <- d[d$is_proj & d$year %in% 2005:2024 & is.finite(d$pred) & d$pred > 0, , drop = FALSE]
      x[key(x) %in% key(history), , drop = FALSE]
    }
    if (is.null(forecasts) || !nrow(forecasts)) {
      return(data.frame(model = name, log_rmse = NA_real_, n_scored_folds = 0L,
        n_scored_observations = 0L, complete = FALSE, fit_converged = fit_converged,
        notes = "No usable forecasts; not ranked"))
    }
    forecasts <- forecasts[key(forecasts) %in% expected, , drop = FALSE]
    complete <- setequal(key(forecasts), expected) && !anyDuplicated(key(forecasts))
    data.frame(model = name, log_rmse = result$hindcast_rmse,
      n_scored_folds = length(unique(forecasts$year)),
      n_scored_observations = nrow(forecasts), complete = complete,
      fit_converged = fit_converged,
      notes = if (!fit_converged) "Experimental main fit failed; not ranked" else {
        if (complete) "Same observations and 20 held-out years" else "Incomplete; not ranked"
      })
  })
  scores <- do.call(rbind, scores)
  scores$rank <- NA_integer_
  ranked <- scores$complete & scores$fit_converged
  scores$rank[ranked] <- rank(scores$log_rmse[ranked], ties.method = "min")
  scores <- scores[order(scores$rank, scores$model, na.last = TRUE), , drop = FALSE]
  utils::write.csv(scores, file.path(cache, "forecast_scores.csv"), row.names = FALSE)
# ---- scores ----
# Lower log-RMSE is better; rank only models scoring the same observations.
# See cache/forecast_scores.csv for coverage and numerical failures.
# ---- end ----

  if (is.null(hind[["Base"]]$pop$ssb)) {
    cli::cli_abort("The baseline hindcast retained no usable SSB estimates; inspect the cached warnings before rendering the cheatsheet.")
  }
  if (!is.null(hind[["Base"]]$pop$ssb)) {
    d <- hind[["Base"]]$pop$ssb
    terminal <- max(base$dat$years[!base$dat$is_proj])
    d <- d[!d$is_proj & d$year >= terminal - 19 & d$year <= terminal, ]
    folds <- sort(unique(d$fold))
    colours <- tam_pal(length(folds))
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(file.path(cache, "hindcast_ssb.pdf"), width = 90 / 25.4,
                           height = 28 / 25.4, family = "Bahnschrift")
    } else {
      grDevices::pdf(file.path(cache, "hindcast_ssb.pdf"), width = 90 / 25.4, height = 28 / 25.4,
                     family = "Helvetica", useDingbats = FALSE)
    }
    graphics::par(mar = c(2.5, 3.3, .5, .4), mgp = c(1.6, .5, 0),
                   tcl = -.2, cex = .7, bty = "l")
    graphics::plot(range(d$year), range(c(0, d$est), finite = TRUE), type = "n",
                    xlab = "Year", ylab = "SSB (tonnes)", yaxs = "i")
    for (i in seq_along(folds)) {
      x <- d[d$fold == folds[i], , drop = FALSE]
      x <- x[order(x$year), , drop = FALSE]
      graphics::lines(x$year[!x$is_proj], x$est[!x$is_proj], col = colours[i],
                       lwd = if (i == length(folds)) 1.5 else .65)
    }
    grDevices::dev.off()
  }

  dashboard <- file.path(repo, "cheatsheet", "tinyAM_dashboard.html")
  dashboard_files <- c(dashboard, file.path(repo, "cheatsheet",
    paste0("tinyAM_", names(fits), ".html")))
  manifest <- file.path(cache, "dashboard-source.rds")
  dashboard_key <- provenance[c("revision", "example_md5", "data_md5")]
  if (!all(file.exists(dashboard_files)) || !file.exists(manifest) ||
      !identical(readRDS(manifest), dashboard_key)) {
    local({
      old <- setwd(file.path(repo, "cheatsheet"))
      on.exit(setwd(old))
# ---- dashboard ----
vis_tam(model_list = fits, output_file = "tinyAM_dashboard.html",
        open_file = FALSE)
# ---- end ----
    })
    for (name in names(fits)) {
      vis_tam(model_list = fits[name],
              output_file = file.path(repo, "cheatsheet", paste0("tinyAM_", name, ".html")),
              open_file = FALSE)
    }
    saveRDS(dashboard_key, manifest)
  }
  attempts <- lapply(records, function(x) x[setdiff(names(x), "value")])
  # Keep the render index small; individual cache files retain every fitted object.
  files <- setNames(file.path("cache", paste0(names(records), ".rds")), names(records))
  saveRDS(list(dat = dat, scores = scores, checks = checks, files = files,
               records = attempts, provenance = provenance, dashboard = dashboard),
          file.path(cache, "examples.rds"))
}
