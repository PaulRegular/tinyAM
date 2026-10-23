# Rscript cheatsheet/examples.R prepare   # fits, two-fold hindcast, report figure and dashboard
# Rscript cheatsheet/examples.R render    # PDF only; never refits
# Rscript cheatsheet/examples.R dashboard # main-fit dashboard only (website build)
# Generated files live in cache/ and output/. PDF rendering needs XeLaTeX,
# ggplot2, rsvg and Bahnschrift; the website uses the tracked PDF instead.

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
  unlist(lapply(c("inputs", "scaa", "sam", "ncam"), function(name) {
    first <- markers[[name]]
    last <- markers[which(markers > first)[1L]]
    lines[seq.int(first + 1L, last - 1L)]
  }), use.names = FALSE)
}
render_cheatsheet <- function(repo) {
  root <- file.path(repo, "cheatsheet")
  cache <- file.path(root, "cache")
  required <- file.path(cache, c("examples.rds", "ssb_comparison.pdf"))
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
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  script <- normalizePath(if (length(script)) script[1L] else "cheatsheet/examples.R", winslash = "/", mustWork = TRUE)
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
    case_study_md5 = unname(tools::md5sum(file.path(repo, "analysis/case_study/001_models.R"))),
    runtime = c(R = as.character(getRversion()), RTMB = as.character(packageVersion("RTMB"))), generated = Sys.time())
  source <- lapply(cod_obs[c("catch", "index", "weight", "maturity")], function(x) {
    x <- x[x$year %in% 1983:2025 & x$age %in% 2:20, , drop = FALSE]
    keys <- intersect(c("year", "age", "survey"), names(x))
    x[do.call(order, x[keys]), , drop = FALSE]
  })
  catch_data <- source$catch
  index_data <- source$index
  weight_data <- source$weight
  maturity_data <- source$maturity
# ---- inputs ----
dat <- make_tam_data(
  catch = catch_data,
  index = index_data,
  weight = weight_data,
  maturity = maturity_data
)
# ---- end ----
  fit_key <- list(sources = sources, example = provenance$example_md5, data = fingerprint(dat), runtime = provenance$runtime)
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
    record <- list(value = value, key = key, provenance = provenance, warnings = warnings, error = error, elapsed_seconds = elapsed)
    saveRDS(record, file)
    records[[name]] <<- record[setdiff(names(record), "value")]
    cli::cli_inform("{name}: finished in {round(elapsed, 1)} s{if (!is.null(error)) paste0('; ', error) else ''}")
    value
  }
  scaa <- cached("fit_scaa", fit_key, function() {
# ---- scaa ----
scaa <- fit_tam(
  dat, years = 1983:2025, ages = 2:14,
  N_settings = list(process = "off", init = "exp"),
  rec_settings = list(process = "off", mu_form = ~ factor(year)),
  F_settings = list(process = "off",
    mu_form = ~ factor(year) + factor(age), mean_ages = 5:14),
  M_settings = list(process = "off",
    mu_supplied = ~ I(0.3), mean_ages = 5:14),
  catch_settings = list(sd_form = ~ 1),
  index_settings = list(sd_form = ~ 1, q_form = ~ mono(q_block))
)
# ---- end ----
    scaa
  })
  if (!inherits(scaa, "tam_fit")) cli::cli_abort("SCAA could not be fitted; inspect cache/fit_scaa.rds.")
  sam <- cached("fit_sam", fit_key, function() {
# ---- sam ----
sam <- update(scaa,
  N_settings = list(process = "iid"),
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "cor_rw", mu_form = NULL)
)
# ---- end ----
    sam
  })
  ncam <- cached("fit_ncam", fit_key, function() {
# ---- ncam ----
ncam <- update(scaa,
  rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "ar1", mu_form = ~ 1),
  M_settings = list(process = "ar1", age_breaks = c(2, 5, 8, 14))
)
# ---- end ----
    ncam
  })
  fits <- list(scaa = scaa, sam = sam, ncam = ncam)
  if (!all(vapply(fits, inherits, logical(1), "tam_fit"))) cli::cli_abort("An example fit failed; inspect cache/fit_*.rds.")
  checks <- lapply(fits, check_tam)
  failed <- names(checks)[!vapply(checks, function(x) isTRUE(x$is_converged), logical(1))]
  if (length(failed)) cli::cli_warn("Numerical checks did not pass for: {paste(failed, collapse = ', ')}. Inspect cached diagnostics.")
  templates <- file_hashes(list.files(file.path(repo, "inst"), pattern = "[.](Rmd|js|css|svg|png)$", recursive = TRUE, full.names = TRUE))
  dashboard_file <- file.path(cache, "dashboard.rds")
  background_file <- file.path(root, "background.md")
  background <- readLines(background_file, warn = FALSE, encoding = "UTF-8")
  dashboard_key <- list(fit = fit_key, templates = templates, background = file_hashes(background_file))
  prior <- if (file.exists(dashboard_file)) readRDS(dashboard_file) else NULL
# ---- dashboard ----
# vis_tam(scaa, sam, ncam)
# ---- end ----
  if (!file.exists(file.path(output, "tinyAM_dashboard.html")) || !identical(prior, dashboard_key)) {
    vis_tam(scaa, sam, ncam, background = background,
      output_file = file.path(output, "tinyAM_dashboard.html"), open_file = FALSE)
    saveRDS(dashboard_key, dashboard_file)
  }
  if (mode == "dashboard") return(invisible(fits))
# ---- tidy ----
out <- tidy_tam(scaa, sam, ncam)
ssb <- out$pop$ssb
index <- out$obs_pred$index
# ---- end ----
  ssb <- ssb[!ssb$is_proj, ]
  ssb$model <- factor(ssb$model, levels = names(fits), labels = toupper(names(fits)))
  ssb[c("est", "lwr", "upr")] <- lapply(ssb[c("est", "lwr", "upr")], function(x) x / 1000)
  colours <- stats::setNames(tam_pal(length(fits)), toupper(names(fits)))
  # ggplot2 is a cheatsheet-build prerequisite, not a package runtime dependency.
  figure <- ggplot2::ggplot(ssb, ggplot2::aes(x = year, y = est, colour = model, fill = model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = .12, colour = NA) +
    ggplot2::geom_line(linewidth = .45) +
    ggplot2::scale_colour_manual(values = colours) + ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, .04))) +
    ggplot2::labs(x = "Year", y = "SSB (kt)", colour = NULL, fill = NULL) +
    ggplot2::theme_classic(base_size = 8.5, base_family = "sans") +
    ggplot2::theme(legend.position = "top", legend.key.width = grid::unit(8, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0), legend.spacing.x = grid::unit(1, "mm"),
      axis.text = ggplot2::element_text(colour = "black"), plot.margin = ggplot2::margin(2, 2, 2, 2))
  ggplot2::ggsave(file.path(cache, "ssb_comparison.pdf"), figure,
    device = grDevices::cairo_pdf, width = 77, height = 61, units = "mm")
  # Rendering needs only main-fit summaries; keep it independent of evaluation.
  saveRDS(list(dat = dat, checks = checks, tidy = out, records = records,
    provenance = provenance), file.path(cache, "examples.rds"))
# ---- hindcast ----
# hind <- fit_hindcast(scaa, sam, ncam)
# ---- end ----
  future::plan(future::sequential)
  hind <- cached("hindcast", c(fit_key, list(folds = 2L, catch = "observed")),
    function() fit_hindcast(scaa, sam, ncam, progress = FALSE))
  if (is.null(hind)) cli::cli_abort("The two-fold hindcast failed; inspect cache/hindcast.rds.")
  saveRDS(list(dat = dat, checks = checks, tidy = out, hindcast = hind,
    records = records, provenance = provenance), file.path(cache, "examples.rds"))
  print(hind)
  invisible(fits)
}

if (sys.nframe() == 0L) main(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1L] else "prepare")
