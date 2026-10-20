# Rscript cheatsheet/examples.R prepare   # examples, hindcasts and dashboards
# Rscript cheatsheet/examples.R render    # PDF only; never refits
# Rscript cheatsheet/examples.R dashboard # main-fit dashboard only (website build)
# Generated files live in cache/ and output/. PDF rendering needs XeLaTeX,
# rsvg and Bahnschrift; the website uses the tracked PDF instead.

fingerprint <- function(x) {
  file <- tempfile()
  on.exit(unlink(file))
  saveRDS(x, file, compress = FALSE, version = 2)
  unname(tools::md5sum(file))
}

file_hashes <- function(paths) {
  stats::setNames(unname(tools::md5sum(paths)), basename(paths))
}

example_code <- function(path) {
  lines <- readLines(path, warn = FALSE)
  markers <- grep("^# ---- .* ----$", lines)
  names(markers) <- sub("^# ---- (.*) ----$", "\\1", lines[markers])
  unlist(lapply(c("inputs", "baseline", "N_iid", "M_ar1"), function(name) {
    first <- markers[[name]]
    last <- markers[which(markers > first)[1L]]
    lines[seq.int(first + 1L, last - 1L)]
  }), use.names = FALSE)
}

hindcast <- function(fit) fit_hindcast(fit, folds = 20, progress = FALSE)

render_cheatsheet <- function(repo) {
  root <- file.path(repo, "cheatsheet")
  cache <- file.path(root, "cache")
  required <- file.path(cache, c("examples.rds", "hindcast_ssb.pdf"))
  if (!all(file.exists(required))) cli::cli_abort("Run examples.R prepare before rendering the cheatsheet.")
  workspace <- tempfile("render-", tmpdir = cache)
  dir.create(workspace)
  on.exit(unlink(workspace, recursive = TRUE), add = TRUE)
  dir.create(file.path(workspace, "cache"))
  file.copy(file.path(root, c("tinyAM_cheatsheet.Rmd", "header.tex", "examples.R")), workspace)
  file.copy(required, file.path(workspace, "cache"))
  tryCatch(rmarkdown::render(file.path(workspace, "tinyAM_cheatsheet.Rmd"),
    output_file = "tinyAM_cheatsheet.pdf", output_dir = workspace,
    intermediates_dir = workspace, params = list(repo = repo), quiet = TRUE),
    error = function(e) {
      logs <- list.files(workspace, pattern = "[.]log$", full.names = TRUE)
      writeLines(c(conditionMessage(e), unlist(lapply(logs, readLines, warn = FALSE))), file.path(cache, "render.log"))
      cli::cli_abort("PDF rendering failed; see cheatsheet/cache/render.log.", parent = e)
    })
  final <- file.path(root, "tinyAM_cheatsheet.pdf")
  if (!file.copy(file.path(workspace, "tinyAM_cheatsheet.pdf"), final, overwrite = TRUE)) {
    file.copy(file.path(workspace, "tinyAM_cheatsheet.pdf"), file.path(cache, "pending.pdf"), overwrite = TRUE)
    cli::cli_abort("Close the existing PDF and rerun render. The new PDF is retained in cache/pending.pdf.")
  }
  unlink(file.path(cache, c("pending.pdf", "render.log")))
  cli::cli_inform("Updated {.file {final}}.")
  invisible(final)
}

main <- function(mode = "prepare") {
  if (length(mode) != 1L || !mode %in% c("prepare", "render", "dashboard")) cli::cli_abort("Use prepare, render, or dashboard.")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
  script <- normalizePath(script, winslash = "/", mustWork = TRUE)
  repo <- dirname(dirname(script))
  if (mode == "render") return(render_cheatsheet(repo))
  pkgload::load_all(repo, quiet = TRUE, helpers = FALSE)
  root <- file.path(repo, "cheatsheet")
  cache <- file.path(root, "cache")
  output <- file.path(root, "output")
  dir.create(cache, showWarnings = FALSE)
  dir.create(output, showWarnings = FALSE)
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  sources <- file_hashes(c(list.files(file.path(repo, "R"), "[.]R$", full.names = TRUE), file.path(repo, "NAMESPACE")))
  provenance <- list(revision = system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE),
    sources = sources, example_md5 = fingerprint(example_code(script)),
    runtime = c(R = as.character(getRversion()), RTMB = as.character(packageVersion("RTMB"))),
    generated = Sys.time())

  # Replace these tables with your own. The cod example is self-contained.
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
  catch = catch_data,
  index = index_data,
  weight = weight_data,
  maturity = maturity_data,
  pop_covariates = population_data
)
check_tam_data(dat)
# A plain list plus check_tam_data() also works
# Optional class conversion: as_tam_data(dat)
# ---- end ----
  fit_key <- list(sources = sources, example = provenance$example_md5,
                  data = fingerprint(dat), runtime = provenance$runtime)
  hind_key <- c(fit_key, list(call = deparse(body(hindcast))))
  records <- list()
  cached <- function(name, key, code) {
    file <- file.path(cache, paste0(name, ".rds"))
    if (file.exists(file)) {
      record <- readRDS(file)
      if (identical(record$key, key)) {
        records[[name]] <<- record[setdiff(names(record), "value")]
        cli::cli_inform("{name}: reused cached result")
        return(record$value)
      }
    }
    warnings <- character(); error <- NULL
    elapsed <- system.time(value <- withCallingHandlers(
      tryCatch(code(), error = function(e) { error <<- conditionMessage(e); NULL }),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
    ))[["elapsed"]]
    record <- list(value = value, key = key, provenance = provenance, warnings = warnings,
                   error = error, elapsed_seconds = elapsed)
    saveRDS(record, file)
    records[[name]] <<- record[setdiff(names(record), "value")]
    cli::cli_inform("{name}: finished in {round(elapsed, 1)} s{if (!is.null(error)) paste0('; ', error) else ''}")
    value
  }
  base <- cached("fit_Base", fit_key, function() {
# ---- baseline ----
base <- fit_tam(
  dat,
  years = 1983:2024,
  ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "rw"),
  F_settings = list(process = "ar1", mu_form = ~ 1),
  M_settings = list(process = "off", mu_supplied = ~ I(0.3)),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(
    sd_form = ~ 1,
    q_form = ~ mono(q_block)
  )
)
# ---- end ----
    base
  })
  if (is.null(base)) cli::cli_abort("The baseline could not be fitted; inspect cache/fit_Base.rds.")
  N_iid <- cached("fit_N_iid", fit_key, function() {
# ---- N_iid ----
N_iid <- update(
  base,
  N_settings = list(process = "iid")
)
# ---- end ----
    N_iid
  })
  M_ar1 <- cached("fit_M_ar1", fit_key, function() {
# ---- M_ar1 ----
M_ar1 <- update(
  base,
  M_settings = list(
    process = "ar1",
    age_breaks = c(2, 4, 6, 8, 14)
  )
)
# ---- end ----
    M_ar1
  })
# ---- models ----
fits <- list(Base = base, N_iid = N_iid, M_ar1 = M_ar1)
# ---- end ----
  fits <- Filter(function(x) inherits(x, "tam_fit"), fits)
  checks <- lapply(fits, check_tam)
  templates <- file_hashes(list.files(file.path(repo, "inst"), pattern = "[.](Rmd|js|css|svg|png)$", recursive = TRUE, full.names = TRUE))
  dashboard_file <- file.path(cache, "dashboards.rds")
  dashboards <- if (file.exists(dashboard_file)) readRDS(dashboard_file) else list()
  dashboard_key <- list(fit = fit_key, templates = templates)
  if (!file.exists(file.path(output, "tinyAM_dashboard.html")) || !identical(dashboards$combined, dashboard_key)) {
# ---- dashboard ----
vis_tam(model_list = fits,
        output_file = "output/tinyAM_dashboard.html",
        open_file = FALSE)
# ---- end ----
    dashboards$combined <- dashboard_key
    saveRDS(dashboards, dashboard_file)
  }
  if (mode == "dashboard") return(invisible(fits))
  future::plan(future::sequential)
  hind <- lapply(names(fits), function(name) cached(paste0("hind_", name), hind_key, function() hindcast(fits[[name]])))
  names(hind) <- names(fits)
# ---- hindcast ----
# hind <- lapply(fits, fit_hindcast, folds = 20)
# ---- end ----
  observations <- base$obs_pred$index
  keys <- c("survey", "year", "age")
  key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
  expected <- key(observations[observations$year %in% 2005:2024 & is.finite(observations$obs) & observations$obs > 0, ])
  scores <- do.call(rbind, lapply(names(hind), function(name) {
    result <- hind[[name]]
    d <- result$obs_pred$index
    forecasts <- if (is.null(d)) NULL else d[d$is_proj & d$year %in% 2005:2024 & is.finite(d$pred) & d$pred > 0, ]
    if (!is.null(forecasts)) forecasts <- forecasts[key(forecasts) %in% expected, ]
    complete <- !is.null(forecasts) && setequal(key(forecasts), expected) && !anyDuplicated(key(forecasts))
    data.frame(model = name, log_rmse = if (is.null(result)) NA_real_ else result$hindcast_rmse,
      n_scored_folds = length(unique(forecasts$year)), n_scored_observations = if (is.null(forecasts)) 0L else nrow(forecasts),
      complete = complete, fit_converged = isTRUE(checks[[name]]$is_converged))
  }))
  scores$rank <- NA_integer_
  ranked <- scores$complete & scores$fit_converged
  scores$rank[ranked] <- rank(scores$log_rmse[ranked], ties.method = "min")
  scores <- scores[order(scores$rank, scores$model, na.last = TRUE), ]
  utils::write.csv(scores, file.path(cache, "forecast_scores.csv"), row.names = FALSE)
  best <- scores$model[which(scores$rank == 1L)[1L]]
  if (is.na(best)) cli::cli_abort("No model has complete, numerically acceptable forecast coverage; inspect the cached diagnostics.")
  d <- hind[[best]]$pop$ssb
  terminal <- max(base$dat$years[!base$dat$is_proj])
  d <- d[!d$is_proj & d$year >= terminal - 19 & d$year <= terminal, ]
  local({
    grDevices::cairo_pdf(file.path(cache, "hindcast_ssb.pdf"), width = 90 / 25.4, height = 53 / 25.4,
                         family = if (.Platform$OS.type == "windows") "Bahnschrift" else "sans")
    on.exit(grDevices::dev.off())
    folds <- sort(unique(d$fold)); colours <- tam_pal(length(folds))
    graphics::par(mar = c(2.5, 3.3, .5, .4), mgp = c(1.6, .5, 0), tcl = -.2, cex = .7, bty = "l")
    graphics::plot(range(d$year), range(c(0, d$est), finite = TRUE), type = "n", xlab = "Year", ylab = "SSB (tonnes)", yaxs = "i")
    for (i in seq_along(folds)) {
      x <- d[d$fold == folds[i], ]; x <- x[order(x$year), ]
      graphics::lines(x$year, x$est, col = colours[i], lwd = if (i == length(folds)) 1.5 else .65)
    }
  })
  files <- stats::setNames(paste0("cache/", names(records), ".rds"), names(records))
  saveRDS(list(dat = dat, scores = scores, checks = checks, files = files, records = records, provenance = provenance), file.path(cache, "examples.rds"))
  for (name in names(hind)) {
    models <- hind[[name]]$fits
    file <- file.path(output, paste0("tinyAM_", name, ".html"))
    html_key <- list(hind = hind_key, templates = templates, self_contained = FALSE)
    if (!length(models)) { unlink(file); next }
    if (!file.exists(file) || !identical(dashboards[[name]]$key, html_key)) {
      diagnostics <- hind[[name]]$fold_diagnostics
      background <- c(paste0("### ", name, " hindcast fits"),
        paste(sum(diagnostics$is_converged), "of", nrow(diagnostics), "requested fits passed numerical checks. Only successful fits are shown; labels identify the terminal historical year."))
      error <- NULL
      tryCatch(vis_tam(model_list = models, output_file = file, open_file = FALSE,
        background = background, render_args = list(output_options = list(self_contained = FALSE,
          lib_dir = file.path(output, "assets")), intermediates_dir = tempdir())),
        error = function(e) { error <<- conditionMessage(e); unlink(file) })
      dashboards[[name]] <- list(key = if (is.null(error)) html_key else NULL, error = error)
      saveRDS(dashboards, dashboard_file)
      if (!is.null(error)) cli::cli_warn("The {name} hindcast dashboard failed: {error}")
    }
  }
  invisible(scores)
}

if (sys.nframe() == 0L) main(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1L] else "prepare")
