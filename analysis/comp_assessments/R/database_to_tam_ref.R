.tam_template_na <- function(x) {
  if (is.data.frame(x)) {
    keys <- c("year", "age", "age_group", "age_block", "fleet", "survey",
              "sex", "region", "season", "samp_time", "q_block", "q_key",
              "par", "coef", "is_proj")
    for (name in setdiff(names(x), keys)) x[[name]][] <- NA
    return(x)
  }
  if (is.list(x)) {
    x[] <- lapply(x, .tam_template_na)
    return(x)
  }
  if (is.atomic(x)) {
    x[] <- NA
    return(x)
  }
  NA
}

.tam_table_key <- function(x, columns) {
  do.call(paste, c(lapply(x[columns], function(value) {
    value <- as.character(value)
    value[is.na(value)] <- "<NA>"
    value
  }), sep = "\034"))
}

.fill_tam_table <- function(template, source) {
  if (!nrow(template) || !nrow(source)) return(template)
  keys <- intersect(c("year", "age", "age_group", "age_block", "fleet", "survey",
                      "sex", "region", "season", "samp_time", "q_block", "q_key",
                      "is_proj"), intersect(names(template), names(source)))
  for (name in setdiff(names(source), names(template))) {
    template[[name]] <- source[[name]][rep(NA_integer_, nrow(template))]
  }
  if (!length(keys)) return(template)
  source_key <- .tam_table_key(source, keys)
  if (anyDuplicated(source_key)) {
    cli::cli_abort("Assessment output rows do not map uniquely to the fitted template groups.")
  }
  index <- match(.tam_table_key(template, keys), source_key)
  matched <- which(!is.na(index))
  for (name in setdiff(names(source), keys)) {
    template[[name]][matched] <- source[[name]][index[matched]]
  }
  template
}

.fill_tam_report_array <- function(x, table, years, ages) {
  if (!is.data.frame(table) || !all(c("year", "est") %in% names(table))) return(x)
  dimensions <- dim(x)
  if (length(dimensions) == 2L) {
    if (!"age" %in% names(table)) return(x)
    dim_names <- dimnames(x)
    year_dim <- match("year", names(dim_names))
    age_dim <- match("age", names(dim_names))
    if (is.na(year_dim)) year_dim <- which(dimensions == length(years))[1L]
    if (is.na(age_dim)) age_dim <- which(dimensions == length(ages))[1L]
    if (is.na(year_dim) || is.na(age_dim) || year_dim == age_dim) return(x)
    year_names <- dim_names[[year_dim]]
    age_names <- dim_names[[age_dim]]
    if (is.null(year_names)) year_names <- as.character(years)
    if (is.null(age_names)) age_names <- as.character(ages)
    row_index <- match(as.character(table$year), year_names)
    column_index <- match(as.character(table$age), age_names)
    take <- which(!is.na(row_index) & !is.na(column_index))
    if (length(take)) {
      coordinates <- matrix(NA_integer_, nrow = length(take), ncol = 2L)
      coordinates[, year_dim] <- row_index[take]
      coordinates[, age_dim] <- column_index[take]
      x[coordinates] <- table$est[take]
    }
    return(x)
  }
  if (is.null(dimensions) && length(x) == length(years)) {
    year_names <- names(x)
    if (is.null(year_names)) year_names <- as.character(years)
    index <- match(as.character(table$year), year_names)
    take <- which(!is.na(index))
    if (length(take)) x[index[take]] <- table$est[take]
  }
  x
}

database_to_tam_ref <- function(assessment_id, outputs, obs = NULL, years = NULL,
                                  ages = NULL, terminal_year = NULL,
                                  age_plus_group = NULL,
                                  comparison_scales = NULL, template = NULL,
                                  assumptions = NULL,
                                  comparison_aggregates = character(),
                                  comparison_age_groups = list(),
                                  comparison_definitions = list()) {
  required <- c("assessment_id", "type", "measure", "year", "age", "age_group",
                "value", "se", "lwr", "upr", "unit")
  missing <- setdiff(required, names(outputs))
  if (length(missing)) {
    cli::cli_abort("outputs is missing required columns: {paste(missing, collapse = ', ')}")
  }
  if (length(assessment_id) != 1L || is.na(assessment_id) || !nzchar(assessment_id)) {
    cli::cli_abort("assessment_id must be one non-empty value.")
  }
  if (!is.null(template) &&
      (!is.list(template) ||
       !all(c("dat", "pop", "rep", "obs_pred", "fixed_par", "random_par") %in% names(template)))) {
    cli::cli_abort("template must be a fitted tinyAM object with its reporting components.")
  }

  source <- outputs[!is.na(outputs$assessment_id) &
                      outputs$assessment_id == assessment_id, , drop = FALSE]
  if (!nrow(source)) cli::cli_abort("No output rows found for assessment_id {.val {assessment_id}}.")
  source$year <- suppressWarnings(as.integer(as.character(source$year)))
  source$age <- suppressWarnings(as.integer(as.character(source$age)))
  source$value <- suppressWarnings(as.numeric(as.character(source$value)))
  source$se <- suppressWarnings(as.numeric(as.character(source$se)))
  source$lwr <- suppressWarnings(as.numeric(as.character(source$lwr)))
  source$upr <- suppressWarnings(as.numeric(as.character(source$upr)))

  if (is.null(years)) {
    years <- if (is.null(template)) sort(unique(source$year[!is.na(source$year)])) else template$dat$years
  }
  if (is.null(ages)) {
    ages <- if (is.null(template)) sort(unique(source$age[!is.na(source$age)])) else template$dat$ages
    if (!length(ages) && !is.null(obs$weight$age)) ages <- sort(unique(obs$weight$age))
  }
  for (name in c("years", "ages")) {
    value <- get(name)
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
        any(value != as.integer(value)) || anyDuplicated(value) ||
        (length(value) > 1L && any(diff(value) != 1L))) {
      cli::cli_abort("{name} must be a consecutive sequence of whole numbers.")
    }
  }
  years <- as.integer(years)
  ages <- as.integer(ages)
  if (!is.null(template) &&
      (!identical(years, as.integer(template$dat$years)) ||
       !identical(ages, as.integer(template$dat$ages)))) {
    cli::cli_abort("Requested years and ages must match the fitted template.")
  }
  if (!is.null(age_plus_group) &&
      (length(age_plus_group) != 1L || !is.numeric(age_plus_group) ||
       !is.finite(age_plus_group) || age_plus_group != as.integer(age_plus_group) ||
       !age_plus_group %in% ages)) {
    cli::cli_abort("age_plus_group must be one of the requested model ages.")
  }
  if (is.null(terminal_year)) terminal_year <- if (length(years)) max(years) else NA_integer_
  if (!is.null(comparison_scales)) {
    if (!is.numeric(comparison_scales) || is.null(names(comparison_scales)) ||
        any(!nzchar(names(comparison_scales))) || anyDuplicated(names(comparison_scales)) ||
        any(!is.finite(comparison_scales)) || any(comparison_scales <= 0)) {
      cli::cli_abort("comparison_scales must be NULL or a named vector of positive finite numbers.")
    }
  }

  age_label <- function(rows) {
    age <- as.character(rows$age)
    grouped <- (is.na(age) | !nzchar(age)) &
      !is.na(rows$age_group) & nzchar(as.character(rows$age_group))
    age[grouped] <- as.character(rows$age_group[grouped])
    if (all(is.na(age) | grepl("^[0-9]+$", age))) {
      suppressWarnings(as.integer(age))
    } else {
      age
    }
  }
  reporting_table <- function(rows) {
    if (!nrow(rows)) return(NULL)
    out <- data.frame(
      year = rows$year,
      age = age_label(rows),
      est = rows$value,
      se = rows$se,
      se_scale = NA_character_,
      lwr = rows$lwr,
      upr = rows$upr,
      unit = rows$unit,
      source_type = rows$source_type,
      source_reference = rows$source_reference,
      notes = rows$notes,
      is_proj = !is.na(rows$year) & !is.na(terminal_year) & rows$year > terminal_year,
      stringsAsFactors = FALSE
    )
    if (all(is.na(out$age))) out$age <- NULL
    dimensions <- intersect(c("age_group", "fleet", "survey", "sex", "region", "season"), names(rows))
    for (name in dimensions) {
      value <- as.character(rows[[name]])
      if (any(!is.na(value) & nzchar(value))) out[[name]] <- value
    }
    out
  }

  measure_map <- c(
    numbers_at_age = "N",
    fishing_mortality_at_age = "F",
    natural_mortality_at_age = "M",
    SSB = "ssb",
    recruitment = "recruitment",
    total_biomass = "biomass",
    total_numbers = "abundance",
    biomass_at_age = "biomass_at_age",
    biomass_by_age_group = "biomass_at_age",
    mature_biomass_at_age = "ssb_mat",
    total_mortality_at_age = "Z",
    Fbar = "F_bar",
    Mbar = "M_bar",
    Zbar = "Z_bar"
  )
  type_map <- c(
    numbers_at_age = "population",
    fishing_mortality_at_age = "mortality",
    natural_mortality_at_age = "mortality",
    SSB = "biomass",
    recruitment = "recruitment",
    total_biomass = "biomass",
    total_numbers = "population",
    biomass_at_age = "biomass",
    biomass_by_age_group = "biomass",
    mature_biomass_at_age = "biomass",
    total_mortality_at_age = "mortality",
    Fbar = "mortality",
    Mbar = "mortality",
    Zbar = "mortality"
  )
  pop <- list()
  for (measure in names(measure_map)) {
    rows <- source[source$measure == measure &
                      source$type == unname(type_map[[measure]]), , drop = FALSE]
    if (nrow(rows)) pop[[unname(measure_map[[measure]])]] <- reporting_table(rows)
  }
  native_pop <- pop

  # Recruitment on the model grid is N at its first age in the same year.
  if (length(ages)) {
    recruitment_age <- unique(stats::na.omit(pop$recruitment$age))
    if (!length(recruitment_age) && !is.null(assumptions)) {
      value <- unique(assumptions$value[assumptions$setting == "recruitment_age"])
      if (length(value) == 1L && grepl("^[0-9]+$", value)) recruitment_age <- as.integer(value)
    }
    if (length(recruitment_age) != 1L || recruitment_age != min(ages)) {
      n <- pop$N
      if (!is.null(n)) {
        at_age <- suppressWarnings(as.integer(as.character(n$age))) == min(ages)
        if ("age_group" %in% names(n)) {
          grouped <- !is.na(n$age_group) & nzchar(n$age_group) &
            n$age_group != as.character(n$age)
          at_age <- at_age & !grouped
        }
        n <- n[which(at_age), , drop = FALSE]
      }
      if (!is.null(n) && nrow(n)) {
        n$notes <- paste(n$notes, "Comparable recruitment: accepted beginning-of-year N at age",
                         min(ages), "in the same calendar year; native recruitment is retained separately.")
        pop$recruitment <- n
        units <- unique(n$unit)
        if (length(units) == 1L) {
          comparison_scales["recruitment"] <- 1 / .translation_number_multiplier(units)
        }
      } else {
        pop$recruitment <- NULL
      }
    } else if (is.null(pop$recruitment$age)) {
      pop$recruitment$age <- rep(min(ages), nrow(pop$recruitment))
    }
  }

  source_m_status <- attr(obs, "translation")$M$status
  if (is.null(pop$M) && !is.null(obs$weight$M_assumption) &&
      identical(source_m_status, "fixed_numerical_input")) {
    m <- unique(obs$weight[c("year", "age", "M_assumption")])
    names(m)[[3L]] <- "est"
    m$year <- as.integer(m$year)
    m$se <- NA_real_
    m$se_scale <- NA_character_
    m$lwr <- NA_real_
    m$upr <- NA_real_
    m$unit <- "per year"
    m$source_type <- "translated_source_input"
    provenance <- attr(obs, "translation")$source_provenance
    m$source_reference <- if (is.null(provenance)) NA_character_ else
      paste(unique(provenance$source_reference[provenance$component == "M"]), collapse = "; ")
    m$notes <- "Fixed natural mortality input; not an estimated output."
    m$is_proj <- FALSE
    pop$M <- m[c("year", "age", "est", "se", "se_scale", "lwr", "upr",
                 "unit", "source_type", "source_reference", "notes", "is_proj")]
  }

  if (!is.null(age_plus_group)) {
    source_n <- pop$N
    for (metric in intersect(c("N", "F", "M", "Z", "biomass_at_age", "ssb_mat"), names(pop))) {
      x <- pop[[metric]]
      age_number <- suppressWarnings(as.integer(as.character(x$age)))
      in_plus <- !is.na(age_number) & age_number >= age_plus_group
      if (!any(in_plus)) next
      lower <- x[!in_plus, , drop = FALSE]
      plus <- x[in_plus, , drop = FALSE]
      expected_ages <- unique(age_number[in_plus])
      by_year <- split(plus, plus$year)
      rows <- lapply(by_year, function(group) {
        row <- group[1L, , drop = FALSE]
        values <- group$est
        collapsed <- nrow(group) > 1L
        complete <- all(is.finite(values)) && age_plus_group %in% group$age &&
          (!collapsed || setequal(as.integer(as.character(group$age)), expected_ages))
        if (!complete) {
          row$est <- NA_real_
        } else if (!collapsed) {
          row$est <- values[[1L]]
        } else if (metric %in% c("F", "M", "Z")) {
          weights <- source_n$est[match(
            paste(group$year, as.integer(as.character(group$age))),
            paste(source_n$year, source_n$age)
          )]
          weights <- suppressWarnings(as.numeric(weights))
          if (length(weights) != length(values) || any(!is.finite(weights)) ||
              any(weights < 0) || sum(weights) <= 0) {
            row$est <- NA_real_
          } else {
            row$est <- stats::weighted.mean(values, weights)
          }
        } else {
          row$est <- sum(values)
        }
        row$age <- age_plus_group
        if ("age_group" %in% names(row)) row$age_group <- paste0(age_plus_group, "+")
        if (collapsed || !complete) {
          row$se <- row$lwr <- row$upr <- NA_real_
          if ("source_reference" %in% names(row)) {
            row$source_reference <- paste(unique(group$source_reference), collapse = "; ")
          }
        }
        source_notes <- if ("notes" %in% names(group)) {
          unique(group$notes[!is.na(group$notes) & nzchar(group$notes)])
        } else character()
        transformation_note <- if (collapsed) {
          paste0("Collapsed ", nrow(group), " source ages for comparison with the tinyAM plus age; ",
                 if (metric %in% c("F", "M", "Z")) {
                   "Mortality is N-weighted and uncertainty is not combined."
                 } else {
                   "values are summed and uncertainty is not combined."
                 })
        } else {
          if (complete) paste0("One source row at age ", age_plus_group,
                               " was retained without aggregation.") else
            "Incomplete source ages or values; the plus-group estimate is unavailable."
        }
        if ("notes" %in% names(row)) {
          row$notes <- paste(c(source_notes, transformation_note), collapse = " ")
        }
        row
      })
      pop[[metric]] <- do.call(rbind, c(list(lower), rows))
      rownames(pop[[metric]]) <- NULL
    }
  }

  # Preserve source groups for comparison construction, then fill model grids.
  source_pop <- lapply(pop, function(x) x[x$year %in% years, , drop = FALSE])
  age_metrics <- c("N", "F", "M", "Z", "biomass_at_age", "ssb_mat")
  for (name in names(pop)) {
    x <- source_pop[[name]]
    age_specific <- name %in% age_metrics
    if (age_specific) {
      if ("age_group" %in% names(x)) {
        grouped <- !is.na(x$age_group) & grepl("^[0-9]+\\+$", x$age_group)
        start <- suppressWarnings(as.integer(sub("\\+$", "", x$age_group)))
        x <- x[!grouped | start == max(ages), , drop = FALSE]
      }
      x$age <- suppressWarnings(as.integer(as.character(x$age)))
      x <- x[x$age %in% ages, , drop = FALSE]
      grid <- expand.grid(year = years, age = ages, KEEP.OUT.ATTRS = FALSE)
    } else {
      grid <- data.frame(year = years)
      if (name == "recruitment" && length(ages)) grid$age <- min(ages)
      # Native aggregate age labels are metadata, not individual fish ages.
      if (name != "recruitment") x$age <- NULL
    }
    grid$est <- NA_real_
    grid$se <- grid$lwr <- grid$upr <- NA_real_
    grid$se_scale <- NA_character_
    grid$is_proj <- grid$year > terminal_year
    pop[[name]] <- .fill_tam_table(grid, x)
  }

  if (is.null(template)) {
    obs_pred <- lapply(obs[c("catch", "index")], function(x) {
      if (is.null(x)) return(NULL)
      x$pred <- NA_real_
      x$sd <- NA_real_
      x$std_res <- NA_real_
      x$osa_res <- NA_real_
      if ("survey" %in% names(x)) x$q <- NA_real_
      x
    })
    fixed_par <- data.frame(par = character(), coef = character(), est = numeric(),
                            se = numeric(), se_scale = character(), lwr = numeric(),
                            upr = numeric(), stringsAsFactors = FALSE)
    out <- list(
      call = match.call(),
      dat = list(obs = obs, years = years, ages = ages,
                 is_proj = rep(FALSE, length(years))),
      pop = pop, obs_pred = obs_pred, fixed_par = fixed_par, random_par = list()
    )
  } else {
    out <- .tam_template_na(template)
    out$call <- match.call()
    out$dat$obs <- if (is.null(obs)) template$dat$obs else obs
    out$dat$years <- template$dat$years
    out$dat$ages <- template$dat$ages
    out$dat$is_proj <- template$dat$is_proj

    for (name in names(template$obs_pred)) {
      input <- template$obs_pred[[name]]
      target <- out$obs_pred[[name]]
      if (!is.data.frame(input) || !is.data.frame(target)) next
      copy <- intersect(c("year", "age", "fleet", "survey", "samp_time", "q_block",
                          "q_key", "is_proj", "obs"), names(input))
      target[copy] <- input[copy]
      out$obs_pred[[name]] <- target
    }

    for (name in names(pop)) {
      source_table <- pop[[name]]
      if (is.data.frame(source_table) && "age_group" %in% names(source_table)) {
        group <- as.character(source_table$age_group)
        is_group <- grepl("^[0-9]+\\+$", group)
        group_start <- suppressWarnings(as.integer(sub("\\+$", "", group)))
        model_plus <- if (length(template$dat$ages)) max(template$dat$ages) else NA_integer_
        source_table <- source_table[!is_group | group_start == model_plus, , drop = FALSE]
      }
      if (is.data.frame(source_table) && !nrow(source_table)) source_table <- NULL
      if (is.data.frame(out$pop[[name]]) && is.data.frame(source_table)) {
        out$pop[[name]] <- .fill_tam_table(out$pop[[name]], source_table)
      } else {
        out$pop[[name]] <- source_table
      }
    }
    for (name in intersect(names(out$rep), names(out$pop))) {
      out$rep[[name]] <- .fill_tam_report_array(
        out$rep[[name]], out$pop[[name]], years = template$dat$years,
        ages = template$dat$ages
      )
    }
  }

  out$comparison_scales <- comparison_scales
  attr(out, "native_pop") <- native_pop
  attr(out, "source_pop") <- source_pop
  if (!is.null(template)) {
    common <- .tam_reference_comparisons(
      template, out, assumptions = assumptions,
      comparison_aggregates = comparison_aggregates,
      comparison_age_groups = comparison_age_groups,
      comparison_definitions = comparison_definitions
    )
    out$comparisons <- common
    out$comparison_scales <- attr(common, "scales")
    for (metric in intersect(c("abundance", "biomass", "ssb", "F_bar", "M_bar"),
                             names(out$pop))) {
      rows <- common[common$metric == metric & !is.na(common$year), , drop = FALSE]
      tab <- out$pop[[metric]]
      if (!is.data.frame(tab) || !"year" %in% names(tab)) next
      index <- match(tab$year, rows$year)
      tab$est <- rows$source[index]
      tab$unit <- rows$unit[index]
      tab$se <- tab$lwr <- tab$upr <- NA_real_
      tab$definition <- rows$definition[index]
      tab$comparison_status <- rows$comparison_status[index]
      tab$notes <- rows$reason[index]
      tab$source_type <- "common_definition_comparison"
      out$pop[[metric]] <- tab
    }
    for (name in intersect(names(out$rep), names(out$pop))) {
      out$rep[[name]] <- .fill_tam_report_array(out$rep[[name]], out$pop[[name]], years, ages)
    }
  }
  class(out) <- c("tam_ref", "list")
  out
}
