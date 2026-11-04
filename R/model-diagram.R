#' View the structure of fitted assessment models
#'
#' See how recruitment, survival, mortality and observations are connected,
#' and which components have annual process variation. The diagram and its
#' accompanying cards use the completed settings stored in each fitted model.
#'
#' @param ... One or more `tam_fit` models or `tam_ref` reporting references,
#'   or one named list of these objects.
#' @param model_list A named list of models, instead of `...`.
#' @details
#' Multiple models have a selector, initially showing the first tinyAM fit.
#' A reporting reference displays an unavailable message because it does not
#' contain a fitted tinyAM specification. The same viewer appears on the
#' dashboard's Model structure page.
#'
#' An annual process set to `"off"` is labeled deterministic: it removes
#' process deviations, not the underlying recruitment, survival or mortality.
#' Initial states are shown separately from annual processes. Stock-recruit
#' relationships add a parent-SSB feedback arrow. Formula cards retain the
#' actual formulas, fleets, surveys, state blocks and projection controls.
#' Short formula labels appear in the diagram; dashed input arrows identify
#' the catch, index or population table supplying active predictors.
#'
#' This is a model overview, not a diagnostic of identifiability or numerical
#' convergence. It reads stored settings without evaluating the fitted
#' objective, sampling uncertainty, calculating residuals or changing the fit.
#' No additional diagramming package is required.
#'
#' @return A browsable HTML container with a responsive SVG diagram and model
#'   settings. Display it directly, or save it with [htmltools::save_html()].
#' @examples
#' \dontrun{
#' plot_model(fit)
#' plot_model(Base = fit, Alternative = other_fit)
#' }
#' @seealso [vis_tam()], [fit_tam()], [tinyAM-model]
#' @export
plot_model <- function(..., model_list = NULL) {
  info <- .dots_or_list(list(...), as.list(substitute(list(...)))[-1L],
    model_list, "model_list")
  fits <- info$fits
  first <- which(vapply(fits, inherits, logical(1), "tam_fit"))
  first <- if (length(first)) first[1L] else 1L
  panes <- lapply(seq_along(fits), function(i) {
    content <- if (inherits(fits[[i]], "tam_fit")) {
      .model_structure_view(fits[[i]]$dat, names(fits)[i])
    } else {
      htmltools::tags$p(role = "status", class = "tam-model-unavailable",
        "Model structure is unavailable for a reporting reference (tam_ref).")
    }
    htmltools::tags$div(`data-model-pane` = as.character(i),
      style = if (i != first) "display:none" else NULL, content)
  })
  selector <- if (length(fits) > 1L) htmltools::tags$label(
    "Model ", htmltools::tags$select(`aria-label` = "Model structure model",
      onchange = paste0("var c=this.closest('.tam-model-structure');",
        "var v=this.value;c.querySelectorAll('[data-model-pane]').forEach(",
        "function(p){p.style.display=p.getAttribute('data-model-pane')===v?'':'none';});"),
      lapply(seq_along(fits), function(i) htmltools::tags$option(
        value = as.character(i), selected = if (i == first) "selected" else NULL,
        names(fits)[i]))))
  htmltools::browsable(htmltools::tags$div(class = "tam-model-structure",
    htmltools::tags$style(htmltools::HTML(paste0(
      ".tam-model-structure{font-family:Arial,sans-serif;color:#1f2933;background:#fff;padding:8px;overflow:auto}",
      ".tam-model-structure select{margin:0 0 12px 8px;padding:5px;max-width:80%}",
      ".tam-model-structure svg{width:100%;height:auto;min-width:630px;display:block}",
      ".tam-model-cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(235px,1fr));gap:12px}",
      ".tam-model-card{border:1px solid #d6dce0;border-top:3px solid #007187;padding:10px}",
      ".tam-model-card h4{font-size:16px;color:#003e56;margin:0 0 6px}",
      ".tam-model-card p{font-size:13px;margin:4px 0;overflow-wrap:anywhere}",
      ".tam-model-card code{font-size:12px;white-space:pre-wrap}",
      ".tam-model-legend,.tam-model-note{font-size:13px;margin:8px 0 12px}",
      ".tam-model-unavailable{padding:25px;border:1px solid #d6dce0}"))),
    selector, panes))
}

.model_structure_spec <- function(dat) {
  form <- function(x) if (is.null(x)) "None" else paste(deparse(x, width.cutoff = 90L), collapse = " ")
  process <- function(x) if (identical(x, "off")) "Deterministic (process off)" else
    paste("Stochastic process:", if (is.null(x)) "unspecified" else toupper(x))
  block_text <- function(x) {
    if (is.null(x)) return("Not specified")
    if (is.list(x)) return(paste(vapply(names(x), function(nm)
      paste0(nm, ": ", block_text(x[[nm]])), character(1)), collapse = "; "))
    if (is.factor(x)) return(paste(levels(x), collapse = ", "))
    paste(x, collapse = ", ")
  }
  fleets <- if (isTRUE(dat$has_fleet)) dat$fleets else "Aggregate fishery"
  surveys <- unique(as.character(dat$obs$index$survey))
  curve <- dat$rec$curve
  predictors <- function(forms, table, data, exclude = character()) {
    forms <- Filter(function(x) inherits(x, "formula"), forms)
    variables <- unique(unlist(lapply(forms, .formula_input_variables), use.names = FALSE))
    # Population year/age coordinates exist even without pop_covariates.
    if (table == "pop_covariates") exclude <- c(exclude, "year", "age")
    variables <- intersect(setdiff(variables, c("obs", exclude)), names(data))
    if (length(variables)) paste0(table, ": ", paste(variables, collapse = ", ")) else ""
  }
  correlation <- if (inherits(dat$F_settings$cor_form, "formula"))
    list(dat$F_settings$cor_form) else dat$F_settings$cor_form
  inputs <- c(
    recruitment = predictors(list(dat$rec_settings$mu_form), "pop_covariates",
      dat$pop_data, exclude = "ssb"),
    population = predictors(list(dat$N_settings$mu_form,
      if (dat$N_settings$process != "off") dat$N_settings$sd_form),
      "pop_covariates", dat$pop_data),
    fishing = predictors(c(list(dat$F_settings$mu_form),
      if (dat$F_settings$process != "off") list(dat$F_settings$sd_form),
      if (dat$F_settings$process %in% c("ar1", "cor_rw")) correlation),
      "catch", dat$F_data),
    natural = predictors(list(dat$M_settings$mu_form, dat$M_settings$mu_supplied,
      if (dat$M_settings$process != "off") dat$M_settings$sd_form),
      "pop_covariates", dat$pop_data),
    catch = predictors(list(dat$catch_settings$sd_form, dat$catch_settings$sd_supplied),
      "catch", dat$obs$catch),
    survey = predictors(list(dat$index_settings$q_form, dat$index_settings$sd_form,
      dat$index_settings$sd_supplied), "index", dat$obs$index, exclude = "N"))
  formula_labels <- c(initial = paste0("N0: ", dat$N_settings$init,
      "; R0: ", dat$rec_settings$init),
    recruitment = if (is.null(dat$rec_settings$mu_form)) "Mean: previous year" else
      paste("Mean", form(dat$rec_settings$mu_form)),
    population = if (is.null(dat$N_settings$mu_form)) "Mean: survivors" else
      paste("Mean", form(dat$N_settings$mu_form)),
    fishing = if (is.null(dat$F_settings$mu_form)) "Mean: previous year" else
      paste("Mean", form(dat$F_settings$mu_form)),
    natural = if (is.null(dat$M_settings$mu_form)) paste("Baseline", form(dat$M_settings$mu_supplied)) else
      paste("Mean effects", form(dat$M_settings$mu_form)),
    catch = paste("SD", form(dat$catch_settings$sd_form)),
    survey = paste("q", form(dat$index_settings$q_form)),
    ssb = "SSB = sum of mature biomass")
  card <- function(title, lines, help) list(title = title, lines = lines, help = help)
  cards <- list(
    initial = card("Initial states", c(
      paste("Older-age initialization:", dat$N_settings$init),
      paste("Recruitment boundary:", dat$rec_settings$init),
      paste("Initial year:", min(dat$years)),
      "Initial parameters are separate from annual process deviations."), "recruitment_formulas"),
    recruitment = card("Recruitment", c(process(dat$rec_settings$process),
      paste("Mean:", form(dat$rec_settings$mu_form)),
      if (!is.null(curve)) paste("Parent-SSB relationship:", curve$type, "; lag", curve$lag),
      "Covariates: pop_covariates."), "recruitment_formulas"),
    population = card("Abundance and cohort survival", c(process(dat$N_settings$process),
      paste("Mean:", form(dat$N_settings$mu_form)),
      if (!identical(dat$N_settings$process, "off")) paste("SD:", form(dat$N_settings$sd_form)),
      paste("Modeled ages:", paste(range(dat$ages), collapse = " to "), "; oldest age is the plus group."),
      "Apparent gains/losses do not identify migration separately from mortality."), "tinyAM-model"),
    fishing = card("Fishing mortality", c(process(dat$F_settings$process),
      paste("Mean:", form(dat$F_settings$mu_form)),
      if (!identical(dat$F_settings$process, "off")) paste("SD:", form(dat$F_settings$sd_form)),
      if (dat$F_settings$process %in% c("ar1", "cor_rw")) paste("Correlation groups:", form(dat$F_settings$cor_form)),
      paste("Fleets:", paste(fleets, collapse = ", ")),
      paste("Age blocks:", block_text(if (isTRUE(dat$has_fleet)) dat$F_blocks else dat$F_settings$age_blocks)),
      if (any(dat$F_settings$state_index == 0L)) "F is exactly zero outside the declared fishing-age blocks.",
      "Total F is the sum of fleet F. Shared age blocks share actual F."), "for_fleet"),
    natural = card("Natural mortality", c(process(dat$M_settings$process),
      paste("Supplied mean:", form(dat$M_settings$mu_supplied)),
      paste("Mean effects:", form(dat$M_settings$mu_form)),
      if (!identical(dat$M_settings$process, "off")) paste("SD:", form(dat$M_settings$sd_form)),
      paste("Age blocks:", block_text(dat$M_settings$age_blocks)),
      if (!identical(dat$M_settings$process, "off")) paste("First process year:", dat$M_settings$first_dev_year),
      "Shared age blocks share actual M; other cells retain their mean."), "tinyAM-model"),
    catch = card("Catch observations", c("Lognormal numbers at age, conditional on population and F.",
      paste("Observation SD:", form(dat$catch_settings$sd_form)),
      paste("Supplied SD:", form(dat$catch_settings$sd_supplied)),
      "Zero and missing observations are excluded from the lognormal likelihood."), "fit_tam"),
    survey = card("Survey indices and catchability", c(
      paste("Surveys:", paste(surveys, collapse = ", ")),
      paste("q formula:", form(dat$index_settings$q_form)),
      paste("q link:", dat$index_settings$q_link),
      paste("Gaussian q terms:", {
        types <- vapply(dat$q_terms, `[[`, character(1), "type")
        types <- unique(types[types %in% c("iid", "rw", "ar1", "smooth")])
        if (length(types)) paste(types, collapse = ", ") else "None"
      }),
      paste("Observation SD:", form(dat$index_settings$sd_form)),
      paste("Supplied SD:", form(dat$index_settings$sd_supplied)),
      "Predictions use each survey's sampling time; q may contain Gaussian effects."), "for_fleet"),
    ssb = card("Spawning biomass and projections", c(
      "SSB uses stock weights, maturity and spawning timing, retaining plus-group biology.",
      if (is.null(dat$proj_settings) || !isTRUE(dat$proj_settings$n_proj > 0)) "Projections: none." else
        paste("Projection settings:", form(dat$proj_settings))), "tinyAM-model"))
  list(cards = cards, feedback = !is.null(curve),
    feedback_label = if (!is.null(curve)) paste0("Parent SSB: ", curve$type,
      ", lag ", curve$lag) else "", inputs = inputs, formula_labels = formula_labels,
    statuses = c(recruitment = dat$rec_settings$process,
      population = dat$N_settings$process, fishing = dat$F_settings$process,
      natural = dat$M_settings$process),
    fleets = fleets, surveys = surveys)
}

.model_structure_view <- function(dat, label) {
  spec <- .model_structure_spec(dat)
  escape <- function(x) as.character(htmltools::htmlEscape(x))
  shorten <- function(x, width = 34L) if (nchar(x) > width)
    paste0(substr(x, 1L, width - 3L), "...") else x
  positions <- list(initial = c(30, 35), recruitment = c(30, 200),
    fishing = c(30, 370), natural = c(310, 370), population = c(590, 200),
    catch = c(920, 370), survey = c(920, 35), ssb = c(590, 425))
  subtitles <- c(initial = "Boundary parameters", recruitment = "New youngest-age fish",
    fishing = "Total F = sum of fleet F", natural = "Shared population M",
    population = "Survivors + recruitment", catch = "Fleet catch predictions",
    survey = "q, sampling time and abundance", ssb = "Mature stock biomass")
  titles <- c(initial = "Initial states", recruitment = "Recruitment",
    fishing = "Fishing mortality", natural = "Natural mortality",
    population = "Population abundance", catch = "Catch observations",
    survey = "Survey observations", ssb = "Spawning biomass")
  arrow <- function(points) {
    n <- nrow(points)
    direction <- points[n, ] - points[n - 1L, ]
    direction <- direction / sqrt(sum(direction^2))
    tip <- points[n, ]; base <- tip - 9 * direction
    normal <- c(-direction[2L], direction[1L]) * 4
    paste0('<polyline fill="none" stroke="#526a75" stroke-width="2" points="',
      paste(apply(points, 1L, paste, collapse = ","), collapse = " "), '"/>',
      '<polygon fill="#526a75" points="', paste(c(
        paste(tip, collapse = ","), paste(base + normal, collapse = ","),
        paste(base - normal, collapse = ",")), collapse = " "), '"/>')
  }
  edges <- list(rbind(c(270, 90), c(550, 90), c(550, 215), c(590, 215)),
    rbind(c(270, 255), c(590, 255)),
    rbind(c(270, 425), c(295, 425), c(295, 335), c(640, 335), c(640, 310)),
    rbind(c(550, 425), c(570, 425), c(570, 290), c(590, 290)),
    rbind(c(830, 230), c(875, 230), c(875, 90), c(920, 90)),
    rbind(c(830, 280), c(875, 280), c(875, 425), c(920, 425)),
    rbind(c(150, 480), c(150, 548), c(1040, 548), c(1040, 480)),
    rbind(c(710, 310), c(710, 425)))
  if (spec$feedback) edges[[length(edges) + 1L]] <-
    rbind(c(590, 480), c(285, 480), c(285, 275), c(270, 275))
  nodes <- vapply(names(positions), function(id) {
    p <- positions[[id]]
    status <- spec$statuses[id]
    stochastic <- length(status) && !is.na(status) && status != "off"
    fill <- if (stochastic) "#d9efed" else "#f0f3f5"
    annotation <- if (is.na(status)) "" else if (stochastic)
      paste("Stochastic", toupper(status)) else "Deterministic: process off"
    title <- paste(spec$cards[[id]]$lines, collapse = "\n")
    paste0('<g data-component="', id, '"><title>', escape(title), '</title>',
      '<rect x="', p[1L], '" y="', p[2L], '" width="240" height="110" rx="8" fill="',
      fill, '" stroke="#007187"/>',
      '<text x="', p[1L] + 12, '" y="', p[2L] + 24,
      '" font-size="16" font-weight="bold">', escape(titles[id]), '</text>',
      '<text x="', p[1L] + 12, '" y="', p[2L] + 46,
      '" font-size="13">', escape(subtitles[id]), '</text>',
      '<text x="', p[1L] + 12, '" y="', p[2L] + 68,
      '" font-size="12">', escape(annotation), '</text>',
      '<text data-formula="', id, '" x="', p[1L] + 12, '" y="', p[2L] + 91,
      '" font-size="11">', escape(shorten(spec$formula_labels[id])), '</text></g>')
  }, character(1))
  input_arrows <- vapply(names(spec$inputs)[nzchar(spec$inputs)], function(id) {
    p <- positions[[id]]
    paste0('<g data-input-component="', id, '"><title>', escape(spec$inputs[id]), '</title>',
      '<rect x="', p[1L], '" y="', p[2L] - 28, '" width="240" height="18" rx="4" ',
      'fill="#fff8df" stroke="#af8a27" stroke-dasharray="3,3"/>',
      '<text x="', p[1L] + 8, '" y="', p[2L] - 15, '" font-size="11">',
      escape(shorten(spec$inputs[id], 36L)), '</text>',
      '<line x1="', p[1L] + 200, '" y1="', p[2L] - 10,
      '" x2="', p[1L] + 200, '" y2="', p[2L],
      '" stroke="#af8a27" stroke-dasharray="3,3"/>',
      '<polygon fill="#af8a27" points="', p[1L] + 200, ',', p[2L], ' ',
      p[1L] + 196, ',', p[2L] - 6, ' ', p[1L] + 204, ',', p[2L] - 6, '"/></g>')
  }, character(1))
  svg <- paste0('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 560" font-family="Arial,sans-serif" ',
    'role="img" aria-label="Model structure: ', as.character(htmltools::htmlEscape(label, attribute = TRUE)), '">',
    '<title>Model structure: ', escape(label), '</title><rect width="1200" height="560" fill="white"/>',
    paste(vapply(edges, arrow, character(1)), collapse = ""), paste(nodes, collapse = ""),
    paste(input_arrows, collapse = ""),
    if (spec$feedback) paste0('<text x="320" y="472" font-size="13">',
      escape(spec$feedback_label), '</text>') else "", '</svg>')
  cards <- lapply(spec$cards, function(card) htmltools::tags$section(class = "tam-model-card",
    htmltools::tags$h4(card$title), lapply(card$lines, htmltools::tags$p),
    htmltools::tags$a(href = paste0("https://paulregular.github.io/tinyAM/reference/", card$help, ".html"),
      "Equations and settings")))
  htmltools::tagList(htmltools::tags$p(class = "tam-model-note",
    "Stored specification for ", htmltools::tags$strong(label),
    ". Arrows show biological and observation relationships, not a computation order."),
    htmltools::HTML(svg), htmltools::tags$p(class = "tam-model-legend",
      "Teal: annual stochastic process. Gray: no annual process deviations. ",
      "Dashed arrows: active input-table predictors. Process off retains the biological component. ",
      "Initial parameters and observations are separate."),
    htmltools::tags$details(htmltools::tags$summary("Formulas and settings"),
      htmltools::tags$div(class = "tam-model-cards", cards)))
}
