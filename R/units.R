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
      mt = 1e3, kt = 1e6, "thousand tonnes" = 1e6, "million tonnes" = 1e9),
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
  units <- head(units, 1L)
  if (length(units)) paste0(label, if (nzchar(label)) " (" else "", units,
    if (nzchar(label)) ")" else "") else label
}
