
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(patchwork)

pkgload::load_all(".", quiet = TRUE)

root <- "analysis/comp_assessments"

plot_agreement_heatmap <- function(root) {

  # Read comparison results and stock metadata.
  comparison <- read_csv(
    file.path(root, "results", "comparison_summary.csv"),
    show_col_types = FALSE
  )

  assessments <- read_csv(
    file.path(root, "database", "assessments.csv"),
    show_col_types = FALSE
  )

  stocks <- read_csv(
    file.path(root, "database", "stocks.csv"),
    show_col_types = FALSE
  )

  # Population metrics.
  metrics <- c("ssb", "recruitment", "abundance", "F_bar")

  metric_labels <- c(
    ssb = "SSB",
    recruitment = "Recruitment",
    abundance = "Abundance",
    F_bar = "Average F"
  )

  comparison <- comparison |>
    filter(metric %in% metrics)

  # Stock names and ocean-basin classification.
  stock_info <- assessments |>
    filter(assessment_id %in% comparison$assessment_id) |>
    select(assessment_id, stock_id) |>
    left_join(
      stocks |>
        select(stock_id, common_name, ocean),
      by = "stock_id"
    ) |>
    mutate(
      stock = coalesce(common_name, assessment_id),
      stock = str_remove(stock, " \\(formerly.*\\)$"),
      basin = ocean,
      basin = factor(
        basin,
        levels = c(
          "Northwest Atlantic",
          "Northeast Atlantic",
          "Northeast Pacific",
          "Northwest Pacific"
        )
      )
    ) |>
    arrange(basin, stock, assessment_id)

  if (anyNA(stock_info$basin)) {
    stop("One or more stocks have an unclassified ocean basin.")
  }

  stock_order <- stock_info$assessment_id

  stock_labels <- setNames(
    stock_info$stock,
    stock_info$assessment_id
  )

  # Complete stock-by-metric grid.
  plot_data <- expand_grid(
    assessment_id = stock_order,
    metric = metrics
  ) |>
    left_join(
      comparison |>
        select(
          assessment_id,
          metric,
          comparison_status,
          median_absolute_percent_difference,
          trend_correlation
        ),
      by = c("assessment_id", "metric")
    ) |>
    mutate(
      assessment_id = factor(
        assessment_id,
        levels = rev(stock_order)
      ),
      metric = factor(
        metric,
        levels = metrics,
        labels = unname(metric_labels[metrics])
      ),
      difference = if_else(
        comparison_status %in% c("matched", "approximate"),
        median_absolute_percent_difference,
        NA_real_
      ),
      correlation = if_else(
        comparison_status %in% c("matched", "approximate"),
        trend_correlation,
        NA_real_
      )
    )

  # Positions of basin brackets, matching heatmap rows.
  basin_ranges <- stock_info |>
    mutate(y = rev(seq_len(n()))) |>
    group_by(basin) |>
    summarise(
      ymin = min(y) - 0.45,
      ymax = max(y) + 0.45,
      ymid = (ymin + ymax) / 2,
      .groups = "drop"
    )

  # Common colour palette: green indicates better agreement.
  green <- tam_cols()[3]
  missing_colour <- "#D9D9D9"

  diff_max <- ceiling(
    max(plot_data$difference, na.rm = TRUE) / 10
  ) * 10

  corr_min <- floor(
    min(plot_data$correlation, na.rm = TRUE) * 10
  ) / 10

  # Shared plot theme.
  heatmap_theme <- theme_minimal(
    base_family = "Arial",
    base_size = 10
  ) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(
        size = 9,
        colour = "black",
        margin = margin(b = 6)
      ),
      axis.text.y = element_text(
        size = 8,
        colour = "black",
        hjust = 1
      ),
      plot.title = element_text(
        size = 13,
        face = "plain",
        hjust = 0.5,
        margin = margin(b = 6)
      ),
      legend.position = "bottom",
      legend.location = "panel",
      legend.justification.bottom = "left",
      legend.box.spacing = grid::unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.title = element_text(
        size = 9,
        hjust = 0.5
      )
    )

  # Shared heatmap function.
  make_heatmap <- function(
    value,
    title,
    limits,
    show_stocks = TRUE
  ) {

    is_difference <- value == "difference"

    p <- ggplot(
      plot_data,
      aes(x = metric, y = assessment_id)
    ) +
      geom_tile(
        aes(fill = .data[[value]]),
        width = 0.94,
        height = 0.88,
        colour = "white",
        linewidth = 0.4
      ) +
      scale_x_discrete(
        position = "top",
        expand = expansion(add = 0.5)
      ) +
      scale_y_discrete(
        limits = rev(stock_order),
        labels = if (show_stocks) {
          function(x) unname(stock_labels[x])
        } else {
          NULL
        },
        expand = expansion(add = 0.5)
      ) +
      scale_fill_gradient(
        low = if (is_difference) green else "white",
        high = if (is_difference) "white" else green,
        limits = limits,
        breaks = if (is_difference) {
          pretty(limits, n = 5)
        } else {
          seq(limits[1], limits[2], length.out = 5)
        },
        labels = if (is_difference) {
          scales::label_number(accuracy = 1)
        } else {
          scales::label_number(accuracy = 0.01)
        },
        na.value = missing_colour,
        oob = scales::squish,
        name = if (is_difference) {
          "Median absolute difference (%)"
        } else {
          "Pearson correlation (r)"
        }
      ) +
      guides(
        fill = guide_colourbar(
          direction = "horizontal",
          theme = theme(
            legend.key.width = grid::unit(1, "null"),
            legend.key.height = grid::unit(0.3, "cm"),
            legend.title.position = "top",
            legend.margin = margin(5, 5, 5, 5)
          )
        )
      ) +
      labs(title = title) +
      heatmap_theme

    if (!show_stocks) {
      p <- p +
        theme(axis.text.y = element_blank())
    }

    p
  }

  # Magnitude difference: 0% is green.
  p_diff <- make_heatmap(
    value = "difference",
    title = "Magnitude difference",
    limits = c(0, 100)
  ) +
    theme(
      plot.margin = margin(20, 20, 10, 2)
    )

  # Trajectory correlation: r = 1 is green.
  p_corr <- make_heatmap(
    value = "correlation",
    title = "Trajectory correlation",
    limits = c(0, 1),
    show_stocks = FALSE
  ) +
    theme(
      plot.margin = margin(20, 2, 10, 2)
    )

  # Narrow ocean-basin bracket panel, placed on the RIGHT.
  p_basin <- ggplot(basin_ranges) +
    geom_segment(
      aes(
        x = 0.35,
        xend = 0.35,
        y = ymin,
        yend = ymax
      ),
      linewidth = 0.5,
      colour = "grey40"
    ) +
    geom_segment(
      aes(
        x = 0.12,
        xend = 0.35,
        y = ymin,
        yend = ymin
      ),
      linewidth = 0.5,
      colour = "grey40"
    ) +
    geom_segment(
      aes(
        x = 0.12,
        xend = 0.35,
        y = ymax,
        yend = ymax
      ),
      linewidth = 0.5,
      colour = "grey40"
    ) +
    geom_text(
      aes(
        x = 0.75,
        y = ymid,
        label = basin
      ),
      angle = -90,
      size = 3,
      colour = "grey35"
    ) +
    scale_x_continuous(
      position = "top",
      breaks = 0.5,
      labels = " "
    ) +
    scale_y_continuous(
      limits = c(0.5, nrow(stock_info) + 0.5),
      expand = expansion(mult = 0),
      breaks = NULL
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      clip = "off"
    ) +
    labs(title = " ") +
    heatmap_theme +
    theme(
      axis.text.x = element_text(
        size = 9,
        colour = "transparent",
        margin = margin(b = 6)
      ),
      axis.text.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(20, 0, 10, 2)
    )

  # Combine: magnitude | correlation | ocean brackets.
  p_diff + p_corr + p_basin +
    plot_layout(
      widths = c(1, 1, 0.10)
    )
}

# Generate plot.
p_agreement <- plot_agreement_heatmap(root)
p_agreement

# Export.
figure_dir <- file.path(root, "results", "figures")

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(figure_dir, "agreement_heatmap.png"),
  plot = p_agreement,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

