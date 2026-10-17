
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)

# Run from the tinyAM repository root.
pkgload::load_all(".", quiet = TRUE)

root <- "analysis/comp_assessments"
source(file.path(root, "R", "plot_stock_order.R"))

# Read annual comparisons and stock metadata.
ts <- read_csv(
  file.path(root, "results/population_timeseries.csv"),
  show_col_types = FALSE
)

assessments <- read_csv(
  file.path(root, "database/assessments.csv"),
  show_col_types = FALSE
)

stocks <- read_csv(
  file.path(root, "database/stocks.csv"),
  show_col_types = FALSE
)

# Match the heatmap: ocean, phylogenetic group, species and stock.
comparison <- read_csv(
  file.path(root, "results/comparison_summary.csv"),
  show_col_types = FALSE
)
stock_info <- assessment_plot_stocks(assessments, stocks,
  comparison$assessment_id[comparison$metric %in% c("ssb", "recruitment", "abundance", "F_bar")])

# Colours from the tinyAM palette.
cols <- c(
  Official = tam_cols()[1],
  tinyAM = tail(tam_cols(), 1)
)

plot_comparison <- function(metric_name) {

  # Common units for each metric.
  unit_factors <- list(
    ssb = c("kg" = 1e-6, "t" = 1e-3, "kt" = 1),
    recruitment = c(
      "fish" = 1e-6,
      "thousand fish" = 1e-3,
      "million fish" = 1,
      "billion fish" = 1e3
    ),
    abundance = c(
      "fish" = 1e-6,
      "thousand fish" = 1e-3,
      "million fish" = 1,
      "billion fish" = 1e3
    ),
    F_bar = c("per year" = 1)
  )

  axis_labels <- c(
    ssb = "Spawning-stock biomass (kt)",
    recruitment = "Recruitment (million fish)",
    abundance = "Total abundance (million fish)",
    F_bar = "Average fishing mortality (yr\u207b\u00b9)"
  )

  if (!metric_name %in% names(axis_labels)) {
    stop("Unknown metric: ", metric_name)
  }

  # Select metric and convert to common units.
  x <- ts |>
    filter(metric == metric_name)

  conversion <- unname(unit_factors[[metric_name]][x$unit])

  if (anyNA(conversion)) {
    stop(
      "Unsupported units: ",
      paste(unique(x$unit[is.na(conversion)]), collapse = ", ")
    )
  }

  x <- x |>
    mutate(
      source = source * conversion,
      tinyAM = tinyAM * conversion
    )

  # Retain all stocks in a consistent order.
  panels <- stock_info

  status <- x |>
    group_by(assessment_id) |>
    summarise(
      has_pairs = any(
        comparison_status %in% c("matched", "approximate") &
          is.finite(source) & is.finite(tinyAM)
      ),
      .groups = "drop"
    )

  panels <- panels |>
    left_join(status, by = "assessment_id") |>
    mutate(
      has_pairs = coalesce(has_pairs, FALSE),
      panel = factor(assessment_id, levels = assessment_id),
      strip = paste0(
        gsub("\n", "<br>", str_wrap(stock, width = 22)),
        "<br><span style='color:#888888; font-size:7pt;'>",
        basin, "</span>"
      )
    )

  panel_labels <- setNames(
    panels$strip,
    panels$assessment_id
  )

  # Paired annual estimates.
  paired <- x |>
    filter(
      comparison_status %in% c("matched", "approximate"),
      is.finite(year),
      is.finite(source),
      is.finite(tinyAM)
    ) |>
    left_join(
      panels |> select(assessment_id, panel),
      by = "assessment_id"
    ) |>
    pivot_longer(
      cols = c(source, tinyAM),
      names_to = "model",
      values_to = "estimate"
    ) |>
    mutate(
      model = recode(
        model,
        source = "Official",
        tinyAM = "tinyAM"
      ),
      model = factor(
        model,
        levels = c("Official", "tinyAM")
      )
    )

  # Identify panels without comparable estimates.
  unavailable <- panels |>
    filter(!has_pairs) |>
    mutate(
      label = case_when(
        metric_name == "F_bar" &
          assessment_id == "afsc_pollock_ebs_2024" ~
          "Fbar definition\nunresolved",

        metric_name == "recruitment" &
          assessment_id == "afsc_cod_goa_2026" ~
          "Different\nrecruitment ages",

        metric_name == "abundance" &
          assessment_id ==
          "ices_norway_pout_north_sea_2026_benchmark" ~
          "Different\nreference timing",

        TRUE ~ "Official estimates\nunavailable"
      )
    )

  # Extract tinyAM trajectories for unavailable comparisons.
  # Model output units are kg, fish, fish and per year.
  model_scale <- c(
    ssb = 1e-6,
    recruitment = 1e-6,
    abundance = 1e-6,
    F_bar = 1
  )[[metric_name]]

  model_rows <- lapply(unavailable$assessment_id, function(id) {

    file <- file.path(
      root, "results", "cache", id, "fit.rds"
    )

    if (!file.exists(file)) return(NULL)

    fit <- readRDS(file)
    pop <- fit$pop[[metric_name]]

    if (!is.data.frame(pop) ||
        !all(c("year", "est") %in% names(pop))) {
      return(NULL)
    }

    historical <- fit$dat$years[!fit$dat$is_proj]

    pop |>
      filter(
        year %in% historical,
        is.finite(year),
        is.finite(est)
      ) |>
      transmute(
        assessment_id = id,
        year,
        estimate = est * model_scale
      )
  })

  model_only <- bind_rows(
    tibble(
      assessment_id = character(),
      year = integer(),
      estimate = numeric()
    ),
    model_rows
  ) |>
    left_join(
      panels |> select(assessment_id, panel),
      by = "assessment_id"
    )

  # Centre annotations on each tinyAM-only trajectory.
  # Join the original stock-specific labels rather than
  # replacing them with a generic message.
  annotations <- model_only |>
    group_by(assessment_id, panel) |>
    summarise(
      year = mean(range(year)),
      estimate = mean(range(estimate)),
      .groups = "drop"
    ) |>
    left_join(
      unavailable |> select(assessment_id, label),
      by = "assessment_id"
    )

  # Handle panels without cached tinyAM trajectories.
  no_series <- unavailable |>
    filter(!assessment_id %in% model_only$assessment_id) |>
    mutate(
      year = median(ts$year[is.finite(ts$year)]),
      estimate = 0,
      label = paste0(label, "\n(tinyAM fit not cached)")
    ) |>
    select(assessment_id, panel, year, estimate, label)

  annotations <- bind_rows(annotations, no_series)

  # Plot.
  ggplot() +
    geom_line(
      data = paired,
      aes(
        x = year,
        y = estimate,
        colour = model,
        linetype = model,
        group = model
      ),
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    geom_line(
      data = model_only,
      aes(
        x = year,
        y = estimate,
        group = panel
      ),
      colour = cols[["tinyAM"]],
      linetype = "11",
      linewidth = 0.7,
      alpha = 0.3,
      na.rm = TRUE
    ) +
    geom_label(
      data = annotations,
      aes(
        x = year,
        y = estimate,
        label = label
      ),
      inherit.aes = FALSE,
      colour = "grey45",
      fill = scales::alpha("white", 0.85),
      label.size = 0,
      size = 2.7,
      lineheight = 0.9
    ) +
    facet_wrap(
      ~panel,
      ncol = 7,
      nrow = 4,
      scales = "free",
      drop = FALSE,
      labeller = labeller(
        panel = as_labeller(panel_labels)
      )
    ) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_linetype_manual(
      values = c(
        Official = "solid",
        tinyAM = "11"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 3),
      labels = function(x) sprintf("%.0f", x),
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      labels = scales::label_number(big.mark = ","),
      limits = c(0, NA),
      expand = expansion(mult = c(0.03, 0.05))
    ) +
    labs(
      x = NULL,
      y = axis_labels[[metric_name]],
      colour = NULL,
      linetype = NULL
    ) +
    guides(
      colour = guide_legend(ncol = 1),
      linetype = guide_legend(ncol = 1)
    ) +
    theme_minimal(
      base_family = "Arial",
      base_size = 10
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.3
      ),
      panel.spacing.x = grid::unit(0.4, "lines"),
      panel.spacing.y = grid::unit(0.45, "lines"),
      strip.text.x = ggtext::element_markdown(
        size = 8,
        face = "plain",
        hjust = 0,
        lineheight = 1.1,
        margin = margin(b = 3)
      ),
      axis.text = element_text(size = 7),
      axis.title.y = element_text(
        size = 11,
        margin = margin(r = 8)
      ),
      legend.position = "right",
      legend.justification = "top",
      legend.direction = "vertical",
      legend.text = element_text(size = 9),
      legend.key.width = grid::unit(1.3, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(3, 5, 3, 5)
    )
}

# Generate the four figures.
p_ssb <- plot_comparison("ssb")
p_rec <- plot_comparison("recruitment")
p_abundance <- plot_comparison("abundance")
p_f <- plot_comparison("F_bar")

# Display a figure.
p_ssb

# Save publication-sized PNG files.
plots <- list(
  ssb = p_ssb,
  recruitment = p_rec,
  abundance = p_abundance,
  F_bar = p_f
)

dir.create(
  file.path(root, "results", "figures"),
  showWarnings = FALSE
)

for (metric in names(plots)) {
  ggsave(
    filename = file.path(root, "results", "figures", paste0(metric, ".png")),
    plot = plots[[metric]],
    width = 14,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}
