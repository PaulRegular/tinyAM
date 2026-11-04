.log_linear_buttons <- list(
  list(
    type = "buttons",
    y = 1, x = 0, pad = list(r = 60), xanchor = "right", yanchor = "top",
    buttons = list(
      list(method = "relayout", args = list("yaxis.type", "linear"), label = "Linear"),
      list(method = "relayout", args = list("yaxis.type", "log"), label = "Log")
    )
  )
)

.forecast_hover <- function(data) {
  if (!all(c("fold", "year") %in% names(data))) return(rep("", nrow(data)))
  horizon <- if ("horizon" %in% names(data)) data$horizon else pmax(0, data$year - data$fold)
  paste0("<br>Origin: ", data$fold, "<br>",
    ifelse(horizon > 0, paste0("Forecast; horizon: ", horizon), "Historical"))
}

.projection_boundary <- function(data) {
  if (!"is_proj" %in% names(data) || !any(data$is_proj, na.rm = TRUE)) return(NA_real_)
  if ("fold" %in% names(data)) {
    origins <- unique(data$fold[!is.na(data$fold)])
    if (length(origins) != 1L) return(NA_real_)
    return(as.numeric(origins))
  }
  suppressWarnings(max(data$year[!data$is_proj], na.rm = TRUE))
}

.forecast_plot_args <- function(data, args) {
  if ("fold" %in% names(data) && is.null(args$text) && is.null(args$hovertemplate)) {
    args$text <- ~paste0("Origin: ", fold, "<br>",
      ifelse(year > fold, paste0("Forecast; horizon: ", year - fold), "Historical"))
  }
  args
}

#' @name plot_tam
#' @rdname plot_tam
#' @export
add_plot_buttons <- function(p) {
  pb <- plotly::plotly_build(p)
  data_traces <- pb$x$data
  if (is.null(data_traces)) {
    data_traces <- list()
  }

  # find ribbon traces (CI)
  rib_idx <- which(vapply(
    data_traces,
    function(tr) identical(tr$type, "scatter") &&
      ((!is.null(tr$fill) && nzchar(tr$fill)) || (is.list(tr$meta) && isTRUE(tr$meta$tam_ci))),
    logical(1)
  ))

  menus <- .log_linear_buttons

  if (length(rib_idx)) {
    idx0 <- as.integer(rib_idx - 1L)
    ci_menu <- list(
      type = "buttons",
      y = 1, x = 0, pad = list(r = 60, t = 75), xanchor = "right", yanchor = "top",
      buttons = list(
        list(
          label = "Show CI",
          method = "restyle",
          args = list(list(visible = rep(TRUE, length(idx0))), idx0)
        ),
        list(
          label = "Hide CI",
          method = "restyle",
          args = list(list(visible = rep(FALSE, length(idx0))), idx0)
        )
      )
    )
    menus <- c(menus, list(ci_menu))
  }

  if (is.null(p$x$layout)) {
    p$x$layout <- list()
  }
  p$x$layout$updatemenus <- menus
  if (!is.null(pb$x$frames)) {
    p$x$frames <- pb$x$frames
  }
  p
}


