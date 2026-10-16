assessment_review_flags <- function(summary, m_mode) {
  core <- c("ssb", "recruitment", "abundance", "N", "F", "M")
  x <- summary[summary$metric %in% core &
    summary$comparison_status %in% c("matched", "approximate") & summary$n > 0, ]
  x$M_treatment <- unname(m_mode[x$assessment_id])
  x <- x[x$metric != "M" | x$M_treatment %in% c("estimated", "fixed approximation"), ]
  flag <- function(value, predicate) is.finite(value) & predicate(value)
  keep <- flag(x$mean_absolute_percent_difference, function(z) z > 20) |
    flag(x$terminal_mean_percent_difference, function(z) abs(z) > 20) |
    flag(x$trend_correlation, function(z) z < .90)
  x[keep, , drop = FALSE]
}
