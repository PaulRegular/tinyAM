audit_assumptions <- function(assessment_id, assumptions) {
  required <- c("assessment_id", "component", "setting", "value", "source_reference", "notes")
  missing <- setdiff(required, names(assumptions))
  if (length(missing)) {
    stop("assumptions is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(assessment_id) != 1L || is.na(assessment_id) || !nzchar(assessment_id)) {
    stop("assessment_id must be one non-empty value.", call. = FALSE)
  }
  x <- assumptions[assumptions$assessment_id == assessment_id, , drop = FALSE]
  if (!nrow(x)) stop("No assumptions are recorded for assessment_id: ", assessment_id, call. = FALSE)

  x$tinyam_support <- "not_checked"
  x$audit_notes <- "Review the source description against the current tinyAM model before fitting."
  key <- paste(x$component, x$setting, sep = "|")

  age_range <- key == "population|age_range"
  x$tinyam_support[age_range] <- "supported"
  x$audit_notes[age_range] <- paste(
    "tinyAM accepts integer ages; matching the range does not establish matching population dynamics."
  )

  plus_group <- key == "population|plus_group"
  x$tinyam_support[plus_group] <- "partially_supported"
  x$audit_notes[plus_group] <- paste(
    "tinyAM can represent an aggregated plus group; compare its construction with the source assessment."
  )

  for (component in c("catch", "index")) {
    likelihood <- key == paste(component, "likelihood", sep = "|") &
      grepl("lognormal", x$value, ignore.case = TRUE)
    x$tinyam_support[likelihood] <- "partially_supported"
    x$audit_notes[likelihood] <- paste(
      "tinyAM has a lognormal observation model; weighting, censoring, and variance details still need review."
    )
  }

  x[c("component", "setting", "value", "tinyam_support", "source_reference", "audit_notes")]
}
