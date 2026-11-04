# Pool explicitly selected sexes/ages; never use these states as observations.
pool_source_states <- function(outputs, plus_age = NULL) {
  measures <- c("numbers_at_age", "fishing_mortality_at_age", "natural_mortality_at_age")
  n <- outputs[outputs$type == "population" & outputs$measure == measures[1], ]
  if (!nrow(n)) cli::cli_abort("Pooling source states requires reported N-at-age.")
  for (dimension in c("region", "season")) {
    if (length(unique(na.omit(n[[dimension]]))) > 1L) {
      cli::cli_abort("Select one source {dimension} before pooling sex/age states.")
    }
  }
  key <- function(x) paste(x$year, x$age, x$sex)
  if (anyDuplicated(key(n)) || length(unique(n$unit)) != 1L ||
    any(!is.finite(n$value) | n$value < 0)) {
    cli::cli_abort("Source N needs unique year-age-sex rows, one unit and finite non-negative values.")
  }
  sexes <- unique(as.character(n$sex))
  if (length(sexes) > 1L && !setequal(sexes, c("female", "male"))) {
    cli::cli_abort("Mixed combined and sex-specific N must not be added together.")
  }
  if (length(sexes) == 2L && any(table(paste(n$year, n$age)) != 2L)) {
    cli::cli_abort("Both source sexes must be present before summing a year-age cell.")
  }
  if (!is.null(plus_age) && (length(plus_age) != 1L || !is.finite(plus_age) ||
    plus_age != as.integer(plus_age) || plus_age < min(n$age))) {
    cli::cli_abort("plus_age must be one whole source age or older.")
  }
  group <- paste(n$year, if (is.null(plus_age)) n$age else pmin(n$age, plus_age))
  groups <- split(seq_len(nrow(n)), group)
  template <- n[vapply(groups, `[`, integer(1), 1L), ]
  if (!is.null(plus_age)) template$age <- pmin(template$age, plus_age)
  template$sex <- "combined"
  template$fleet <- NA_character_
  template$se <- template$lwr <- template$upr <- NA_real_
  pooled <- lapply(seq_along(measures), function(i) {
    result <- template
    result$type <- if (i == 1L) "population" else "mortality"
    result$measure <- measures[i]
    if (i == 1L) {
      result$value <- vapply(groups, function(k) sum(n$value[k]), numeric(1))
    } else {
      rate <- outputs[outputs$type == "mortality" & outputs$measure == measures[i], ]
      if (!nrow(rate)) return(NULL)
      if (length(unique(rate$unit)) != 1L) cli::cli_abort("Source rates need one unit before pooling.")
      if (any(rate$value < 0, na.rm = TRUE) || any(is.infinite(rate$value))) {
        cli::cli_abort("Source mortality must be non-negative and finite where reported.")
      }
      if (i == 2L) {
        labelled <- !is.na(rate$fleet) & nzchar(rate$fleet)
        if (any(labelled) && any(!labelled)) {
          cli::cli_abort("Select source total F or component fleets before pooling; do not add both.")
        }
        if (anyDuplicated(paste(key(rate), rate$fleet))) cli::cli_abort("Duplicate source fleet F rows.")
        rates <- tapply(rate$value, key(rate), sum)
        # An omitted fleet rate is unknown, not an assumed zero contribution.
        counts <- table(key(rate))
        rates[names(counts)[counts != length(unique(rate$fleet))]] <- NA_real_
      } else {
        if (anyDuplicated(key(rate))) cli::cli_abort("Duplicate source M rows.")
        rates <- setNames(rate$value, key(rate))
      }
      value <- unname(rates[key(n)])
      contribution <- ifelse(n$value == 0, 0, n$value * value)
      result$value <- vapply(groups, function(k) {
        if (sum(n$value[k]) == 0) NA_real_ else sum(contribution[k]) / sum(n$value[k])
      }, numeric(1))
      result$unit <- unique(rate$unit)
    }
    result$notes <- if (i == 1L) "Source N summed across selected sexes and terminal ages; covariance unavailable."
      else if (i == 2L) "Component fleet F summed first, then N-weighted across selected sexes and terminal ages; covariance unavailable."
      else "Source M N-weighted across selected sexes and terminal ages; covariance unavailable."
    result
  })
  rbind(outputs[!outputs$measure %in% measures, ], do.call(rbind, pooled))
}
