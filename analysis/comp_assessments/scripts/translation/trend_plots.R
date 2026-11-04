
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)

# Run from the tinyAM repository root.
pkgload::load_all(".", quiet = TRUE)

root <- "analysis/comp_assessments"
source(file.path(root, "R", "read_database.R"))
source(file.path(root, "R", "translation_registry.R"))
if (!exists("translation_tier", inherits = FALSE)) translation_tier <- "curated"
results_dir <- .assessment_results(translation_tier)
plot_database <- .saio_metadata(tier = translation_tier)
.assessment_check_results_revision(results_dir, plot_database)
translation_ids <- assessment_translations(translation_tier)$assessment_id
source(file.path(root, "R", "plot_stock_order.R"))

# Read annual comparisons and stock metadata.
ts <- read_csv(
  file.path(results_dir, "population_timeseries.csv"),
  show_col_types = FALSE
)
batch_runs <- .assessment_cached_runs(translation_tier, plot_database)

assessments <- .saio_read_table("assessments.csv", revision = plot_database$commit)

stocks <- .saio_read_table("stocks.csv", revision = plot_database$commit)

ts <- ts |> filter(assessment_id %in% translation_ids)

# Match the heatmap: ocean, phylogenetic group, species and stock.
comparison <- read_csv(
  file.path(results_dir, "comparison_summary.csv"),
  show_col_types = FALSE
)
comparison <- comparison |> filter(assessment_id %in% translation_ids)
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
    ssb = c("kg" = 1e-3, "t" = 1, "tonnes" = 1, "kt" = 1000),
    recruitment = c(
      "fish" = 1e-3,
      "thousand fish" = 1,
      "million fish" = 1000,
      "billion fish" = 1e6
    ),
    abundance = c(
      "fish" = 1e-3,
      "thousand fish" = 1,
      "million fish" = 1000,
      "billion fish" = 1e6
    ),
    F_bar = c("per year" = 1)
  )

  axis_labels <- c(
    ssb = "SSB (tonnes)",
    recruitment = "Recruitment (thousand fish)",
    abundance = "Total abundance (thousand fish)",
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
  target_unit <- if (metric_name == "ssb") "tonnes" else if (metric_name == "F_bar") "per year" else "thousand fish"
  kind <- if (metric_name == "ssb") "mass" else if (metric_name == "F_bar") "rate" else "numbers"

  model_rows <- lapply(unavailable$assessment_id, function(id) {

    fit <- if (!is.null(batch_runs)) batch_runs[[id]]$fit else
      .assessment_cached_fit(id, translation_tier, plot_database)
    if (!inherits(fit, "tam_fit")) return(NULL)
    pop <- fit$pop[[metric_name]]

    if (!is.data.frame(pop) ||
        !all(c("year", "est") %in% names(pop))) {
      return(NULL)
    }
    if (is.null(pop$unit) || !any(is.finite(tinyAM:::.unit_scale(pop$unit, kind)))) return(NULL)
    pop <- tinyAM:::.convert_report_units(pop, target_unit, kind)

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
        estimate = est
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
      label = paste0(label, "\n(no current tinyAM fit)")
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
      nrow = max(4, ceiling(nrow(panels) / 7)),
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
  file.path(results_dir, "figures"),
  showWarnings = FALSE
)

for (metric in names(plots)) {
  ggsave(
    filename = file.path(results_dir, "figures", paste0(metric, ".png")),
    plot = plots[[metric]],
    width = 14,
    height = max(8, 2 * ceiling(nrow(stock_info) / 7)),
    dpi = 300,
    bg = "white"
  )
}
