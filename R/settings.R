# Settings are merged by entry name; formulas and vectors remain whole values.
.settings_names <- c("N_settings", "rec_settings", "F_settings", "M_settings",
                     "catch_settings", "index_settings", "proj_settings")

.settings_defaults <- function(name) {
  if (name == "proj_settings") return(list(n_proj = NULL, n_mean = 3, F_mult = 1))
  eval(formals(prepare_tam)[[name]], envir = environment(prepare_tam))
}

.merge_settings <- function(base, changes, name) {
  if (name == "proj_settings" && is.null(changes)) return(NULL)
  if (!is.list(changes)) cli::cli_abort("{.arg {name}} must be a named list.")
  if (!length(changes)) return(base)
  keys <- names(changes)
  allowed <- names(.settings_defaults(name))
  if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    cli::cli_abort("Entries in {.arg {name}} must have unique, nonempty names.")
  }
  if (name == "N_settings" && "init_N0" %in% keys) {
    cli::cli_abort("N_settings$init_N0 has been retired. Use N_settings$init instead.")
  }
  if (name == "N_settings" && "rec_form" %in% keys) {
    cli::cli_abort("N_settings$rec_form has moved. Use rec_settings with separate process and mu_form entries.")
  }
  unknown <- setdiff(keys, allowed)
  if (length(unknown)) cli::cli_abort("Unknown {.arg {name}} setting(s): {paste(unknown, collapse = ', ')}.")
  if (is.null(base)) base <- .settings_defaults(name)
  base[keys] <- changes
  required <- switch(name,
    N_settings = c("process", "init"), rec_settings = "process", F_settings = "process", M_settings = "process",
    catch_settings = c("sd_form"),
    index_settings = c("sd_form", "q_form", "q_link"),
    proj_settings = c("n_proj", "n_mean", "F_mult"))
  missing <- required[vapply(base[required], is.null, logical(1))]
  if (length(missing)) cli::cli_abort("{.arg {name}} requires non-NULL {paste(missing, collapse = ', ')}.")
  if (name %in% c("N_settings", "F_settings", "M_settings") && is.null(base$sd_form)) base$sd_form <- ~ 1
  if (name == "proj_settings") {
    for (key in c("n_proj", "n_mean")) {
      value <- base[[key]]
      lower <- if (key == "n_proj") 0 else 1
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value != trunc(value) || value < lower) {
        cli::cli_abort("{.arg proj_settings}${key} must be one integer at least {lower}.")
      }
    }
  }
  base
}

.complete_settings <- function(args) {
  for (name in .settings_names) {
    base <- if (name == "proj_settings") NULL else .settings_defaults(name)
    args[name] <- list(.merge_settings(base, args[[name]], name))
  }
  args
}

# Keep Usage synchronized with function formals, including optional projections.
.settings_usage <- function(name) {
  defaults <- formals(get(name))
  entries <- vapply(names(defaults), function(key) {
    if (identical(defaults[[key]], quote(expr = ))) return(key)
    value <- trimws(deparse(defaults[[key]], width.cutoff = 70L), which = "right")
    paste0(key, " = ", paste(value, collapse = "\n    "))
  }, character(1))
  paste(c("@usage", paste0(name, "("),
    paste0("  ", entries, c(rep(",", length(entries) - 1L), "")), ")", "",
    "# Optional projections:",
    "# proj_settings = list(n_proj = 3, n_mean = 3, F_mult = 1)"), collapse = "\n")
}
