# Units are labels and scale metadata, never conversions of model inputs.
.unit_scale <- function(unit, kind) {
  key <- tolower(trimws(unit))
  key <- gsub("_", " ", key, fixed = TRUE)
  key <- gsub("thousands? of fish", "thousand fish", key)
  key <- gsub("millions? of fish", "million fish", key)
  key <- gsub("billions? of fish", "billion fish", key)
  scales <- switch(kind,
    numbers = c(fish = 1, number = 1, numbers = 1, individuals = 1,
      "thousand fish" = 1e3, "million fish" = 1e6, "billion fish" = 1e9),
    weight = c(g = 1e-3, "g/fish" = 1e-3, "g per fish" = 1e-3,
      kg = 1, "kg/fish" = 1, "kg per fish" = 1,
      t = 1e3, "t/fish" = 1e3, "tonnes per fish" = 1e3),
    mass = c(g = 1e-3, kg = 1, t = 1e3, tonne = 1e3, tonnes = 1e3,
      mt = 1e3, kt = 1e6, "thousand tonnes" = 1e6, "thousand t" = 1e6,
      "million tonnes" = 1e9, "million t" = 1e9),
    rate = c("per year" = 1, "per_year" = 1, "yr-1" = 1, "year^-1" = 1)
  )
  unname(scales[key])
}

.table_unit <- function(d, label = "Data") {
  if (is.null(d$unit)) return(NA_character_)
  if (!is.character(d$unit) || anyNA(d$unit) || any(!nzchar(trimws(d$unit)))) {
    cli::cli_abort("{label}$unit must contain nonempty character labels, without missing values.")
  }
  values <- unique(d$unit)
  if (length(values) > 1L) cli::cli_abort("{label} has inconsistent units. Convert values explicitly before combining them.")
  if (length(values)) values else NA_character_
}

# Count rates retain their effort denominator; they are not area-expanded totals.
.count_unit <- function(unit) {
  key <- tolower(trimws(unit))
  key <- gsub("_", " ", key, fixed = TRUE)
  key <- sub("^mean (numbers?|counts?)", "fish", key)
  parts <- strsplit(key, "\\s+per\\s+|\\s*/\\s*")[[1L]]
  scale <- .unit_scale(parts[1L], "numbers")
  list(scale = scale, effort = if (length(parts) == 2L) parts[2L] else "")
}

.q_unit <- function(index, numbers) {
  if (is.na(index) || is.na(numbers)) return(NA_character_)
  a <- .count_unit(index)
  b <- .count_unit(numbers)
  if (!is.finite(a$scale) || !is.finite(b$scale)) return(NA_character_)
  if (a$scale != b$scale) return(paste0("(", index, ") / (", numbers, ")"))
  if (nzchar(a$effort)) paste("per", a$effort) else "unitless"
}

.convert_report_units <- function(data, to, kind) {
  if (!is.data.frame(data) || is.null(data$unit)) return(data)
  factor <- .unit_scale(data$unit, kind) / .unit_scale(to, kind)
  known <- is.finite(factor)
  for (nm in intersect(c("est", "lwr", "upr", "obs", "pred"), names(data))) {
    data[[nm]][known] <- data[[nm]][known] * factor[known]
  }
  if ("se" %in% names(data)) {
    natural <- if ("se_scale" %in% names(data)) data$se_scale == "reported" else rep(TRUE, nrow(data))
    i <- which(known & natural)
    data$se[i] <- data$se[i] * factor[i]
    if ("se_scale" %in% names(data)) {
      data$se[which(known & factor != 1 & is.na(data$se_scale))] <- NA_real_
    }
  }
  data$unit[known] <- to
  data
}

.data_units <- function(data) {
  numbers <- .table_unit(data$catch, "catch")
  weight <- .table_unit(data$weight, "weight")
  index <- vapply(if (is.null(data$index)) list() else split(data$index, as.character(data$index$survey)),
    .table_unit, character(1), label = "index (within survey)")
  scale <- .unit_scale(numbers, "numbers") * .unit_scale(weight, "weight")
  mass <- c("g", "kg", "tonnes", "thousand tonnes", "million tonnes")[
    match(scale, c(1e-3, 1, 1e3, 1e6, 1e9))]
  list(numbers = numbers, weight = weight, mass = mass, index = index)
}

.population_unit <- function(dat, metric) {
  units <- .data_units(dat$obs)
  if (metric == "W") return(units$weight)
  if (metric == "P") return("proportion")
  if (metric == "S") return("dimensionless")
  if (metric %in% c("N", "N_plus", "abundance", "recruitment", "total_catch", "total_catch_pred", "n0")) return(units$numbers)
  if (metric %in% c("biomass", "biomass_mat", "ssb", "ssb_mat", "total_yield", "total_yield_pred")) return(units$mass)
  if (metric %in% c("F", "M", "Z", "mu_F", "mu_M", "F_bar", "M_bar", "F_full")) return("per year")
  NA_character_
}

.unit_label <- function(data, label) {
  if (!is.data.frame(data) || is.null(data$unit)) return(label)
  values <- intersect(c("est", "obs", "pred"), names(data))
  if (length(values)) data <- data[Reduce(`|`, lapply(data[values], is.finite)), , drop = FALSE]
  units <- unique(data$unit[!is.na(data$unit) & nzchar(data$unit)])
  signatures <- vapply(units, function(u) {
    for (kind in c("numbers", "mass", "rate")) {
      scale <- .unit_scale(u, kind)
      if (is.finite(scale)) return(paste(kind, scale))
    }
    u
  }, character(1))
  if (length(unique(signatures)) > 1L) cli::cli_abort("Plot data contain incompatible units: {paste(units, collapse = ', ')}. Convert them explicitly before overlaying models.")
  units <- utils::head(units, 1L)
  if (length(units) && grepl("native|unresolved|survey[ _-]*index|relative index", units, ignore.case = TRUE)) return(label)
  if (length(units)) paste0(label, if (nzchar(label)) " (" else "", units,
    if (nzchar(label)) ")" else "") else label
}

.check_unit_conversion <- function(from, to, scale) {
  kinds <- c("numbers", "mass", "rate")
  a <- vapply(kinds, function(kind) .unit_scale(from, kind), numeric(1))
  b <- vapply(kinds, function(kind) .unit_scale(to, kind), numeric(1))
  common <- is.finite(a) & is.finite(b)
  if (any(common)) {
    expected <- unname((a / b)[common])
    if (!isTRUE(all.equal(scale, expected, tolerance = 1e-10))) {
      cli::cli_abort("comparison_scales uses {scale} for {.val {from}} to {.val {to}}; these units require {expected}. Correct the scale explicitly.")
    }
  } else if (any(is.finite(a)) && any(is.finite(b))) {
    cli::cli_abort("Cannot compare incompatible units {.val {from}} and {.val {to}}.")
  }
  invisible(TRUE)
}
