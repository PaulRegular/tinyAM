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
#' @return A browsable HTML container with responsive HTML nodes, SVG connectors and model
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
      `data-model-selector` = "",
      lapply(seq_along(fits), function(i) htmltools::tags$option(
        value = as.character(i), selected = if (i == first) "selected" else NULL,
        names(fits)[i]))))
  htmltools::browsable(htmltools::tags$div(class = "tam-model-structure",
    htmltools::tags$style(htmltools::HTML(paste0(
      ".tam-model-structure{font:inherit;color:#1f2933;background:#fff;padding:8px;overflow:auto}",
      ".tam-model-structure select{margin:0 0 12px 8px;padding:5px;max-width:80%}",
      ".tam-model-graph{position:relative;display:grid;grid-template-columns:repeat(4,minmax(0,1fr));",
      "grid-template-areas:'initial . . survey' 'recruitment . population .' 'fishing natural ssb catch';",
      "gap:2.5em 2em;padding:1.5em;isolation:isolate}",
      ".tam-model-connectors{position:absolute;inset:0;width:100%;height:100%;pointer-events:none;z-index:-1}",
      ".tam-model-node-wrap{min-width:0;align-self:center}",
      ".tam-model-node{font:inherit;background:#f1f3f5;border:1px solid #c5ced3;border-radius:8px;padding:.75em}",
      ".tam-model-node[data-process='stochastic']{background:#d9e2e6;border-color:#8caab6}",
      ".tam-model-node[data-process='deterministic']{background:#fde7df;border-color:#e3ad98}",
      ".tam-model-node h4{font:inherit;font-weight:bold;margin:0 0 .4em}",
      ".tam-model-node p{font:inherit;font-size:.875em;margin:.35em 0;overflow-wrap:anywhere}",
      ".tam-model-node code{font-size:.85em;white-space:normal;overflow-wrap:anywhere;color:inherit;background:transparent;padding:0}",
      ".tam-model-input{font-size:.8em;overflow-wrap:anywhere;background:#fff8df;border:1px dashed #af8a27;",
      "border-radius:4px;padding:.35em .5em;margin-bottom:max(1em,16px)}",
      ".tam-model-connections{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0)}",
      ".tam-model-graph[data-layout='compact']{grid-template-columns:repeat(2,minmax(0,1fr));",
      "grid-template-areas:'initial recruitment' 'fishing natural' 'population population' 'survey catch' 'ssb ssb'}",
      ".tam-model-graph[data-layout='stacked']{grid-template-columns:minmax(0,1fr);",
      "grid-template-areas:'initial' 'recruitment' 'fishing' 'natural' 'population' 'survey' 'catch' 'ssb'}",
      ".tam-model-cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(235px,1fr));gap:12px}",
      ".tam-model-card{border:1px solid #d6dce0;border-top:3px solid #007187;padding:10px}",
      ".tam-model-card h4{font:inherit;font-weight:bold;color:#003e56;margin:0 0 6px}",
      ".tam-model-card p{font-size:.875em;margin:4px 0;overflow-wrap:anywhere}",
      ".tam-model-card code{font-size:.85em;white-space:pre-wrap}",
      ".tam-model-legend,.tam-model-note{font-size:.875em;margin:8px 0 12px}",
      ".tam-model-unavailable{padding:25px;border:1px solid #d6dce0}"))),
    selector, panes,
    htmltools::tags$script(htmltools::HTML(paste(readLines(
      system.file("js", "model-structure.js", package = "tinyAM"), warn = FALSE), collapse = "\n")))))
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
  shorten <- function(x, width = 65L) if (nchar(x) > width)
    paste0(substr(x, 1L, width - 3L), "...") else x
  subtitles <- c(initial = "Boundary parameters", recruitment = "New youngest-age fish",
    fishing = "Total F = sum of fleet F", natural = "Shared population M",
    population = "Survivors + recruitment", catch = "Fleet catch predictions",
    survey = "q, sampling time and abundance", ssb = "Mature stock biomass")
  titles <- c(initial = "Initial states", recruitment = "Recruitment",
    fishing = "Fishing mortality", natural = "Natural mortality",
    population = "Population abundance", catch = "Catch observations",
    survey = "Survey observations", ssb = "Spawning biomass")
  nodes <- lapply(names(titles), function(id) {
    status <- spec$statuses[id]
    stochastic <- length(status) && !is.na(status) && status != "off"
    state <- if (is.na(status)) "neutral" else if (stochastic) "stochastic" else "deterministic"
    annotation <- if (is.na(status)) "" else if (stochastic)
      paste("Stochastic", toupper(status)) else "Deterministic: process off"
    input <- spec$inputs[id]
    htmltools::tags$div(class = "tam-model-node-wrap", style = paste0("grid-area:", id),
      if (length(input) && !is.na(input) && nzchar(input)) htmltools::tags$div(
        class = "tam-model-input", `data-input-component` = id, title = input, input),
      htmltools::tags$section(class = "tam-model-node", `data-component` = id,
        `data-process` = state, title = paste(spec$cards[[id]]$lines, collapse = "\n"),
        htmltools::tags$h4(titles[id]), htmltools::tags$p(subtitles[id]),
        if (nzchar(annotation)) htmltools::tags$p(annotation),
        htmltools::tags$code(`data-formula` = id, spec$formula_labels[id] |> shorten()),
        if (id == "ssb" && spec$feedback) htmltools::tags$p(spec$feedback_label)))
  })
  edges <- list(c("initial", "population"), c("recruitment", "population"),
    c("fishing", "population"), c("natural", "population"),
    c("population", "survey"), c("population", "catch"),
    c("fishing", "catch"), c("population", "ssb"))
  if (spec$feedback) edges[[length(edges) + 1L]] <- c("ssb", "recruitment")
  graph <- htmltools::tags$div(class = "tam-model-graph", role = "group",
    `aria-label` = paste("Model structure:", label), `data-feedback` = tolower(spec$feedback),
    htmltools::HTML('<svg class="tam-model-connectors" xmlns="http://www.w3.org/2000/svg" aria-hidden="true"></svg>'),
    nodes,
    htmltools::tags$ul(class = "tam-model-connections", lapply(edges, function(edge)
      htmltools::tags$li(paste(titles[edge], collapse = " to ")))))
  cards <- lapply(spec$cards, function(card) htmltools::tags$section(class = "tam-model-card",
    htmltools::tags$h4(card$title), lapply(card$lines, htmltools::tags$p),
    htmltools::tags$a(href = paste0("https://paulregular.github.io/tinyAM/reference/", card$help, ".html"),
      "Equations and settings")))
  htmltools::tagList(graph, htmltools::tags$p(class = "tam-model-legend",
      "Navy: annual stochastic process. Orange: no annual process deviations. ",
      "Dashed arrows: active input-table predictors. Process off retains the biological component. ",
      "Initial parameters and observations are separate."),
    htmltools::tags$p(class = "tam-model-note",
      "Arrows show biological and observation relationships, not a computation order."),
    htmltools::tags$details(htmltools::tags$summary("Formulas and settings"),
      htmltools::tags$div(class = "tam-model-cards", cards)))
}
