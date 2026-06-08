# 11_resource_plots.R
# Resource usage ranking figures (2×2 Iso-Seq / ONT × CPU / RAM).
# Depends on: 00_figure_config.R, 02_themes.R, patchwork, ggpattern, RColorConesa

#' Join & Call / Call & Join bar fills (resource per-tool panels).
PAPER_STRATEGY_FILL <- c("#15918A", "#F9A856")

#' Tool order and colors for resource ranking plots.
PAPER_RESOURCE_RANKING_TOOLS <- c(
  "FLAIR", "IsoQuant", "Bambu", "TALON", "Mandalorion", "IsoSeq"
)

#' Bar width for ranking panels (6 tools; narrower than default geom_col 0.9).
PAPER_RESOURCE_RANKING_BAR_WIDTH <- 0.55

#' Parse Nextflow-style elapsed/cpu time (`0-05:40:51`) to hours.
parse_time_to_hours <- function(time_str) {
  parts <- strsplit(time_str, "-", fixed = TRUE)[[1L]]
  if (length(parts) == 2L) {
    days <- as.numeric(parts[[1L]])
    time_part <- parts[[2L]]
    time_parts <- strsplit(time_part, ":", fixed = TRUE)[[1L]]
    hours <- as.numeric(time_parts[[1L]])
    minutes <- as.numeric(time_parts[[2L]])
    seconds <- as.numeric(time_parts[[3L]])
    return(days * 24 + hours + minutes / 60 + seconds / 3600)
  }
  0
}

#' Ranking bar chart: tools on x, CPU hours or max RAM on y; strategy via ggpattern.
create_resource_ranking_plot <- function(
    metric = c("cpu", "ram"),
    data_type = c("PacBio", "ONT"),
    tissue = c("brain", "kidney"),
    data,
    include_context = FALSE) {
  metric <- match.arg(metric)
  data_type <- match.arg(data_type)
  tissue <- match.arg(tissue)

  data_type_map <- c("PacBio" = "IsoSeq", "ONT" = "ONT")
  dt <- data_type_map[[data_type]]
  ti <- tissue

  df <- data %>%
    dplyr::filter(.data$tissue == !!ti, .data$data_type == !!dt) %>%
    dplyr::mutate(cpu_hours = vapply(.data$cpu_time, parse_time_to_hours, numeric(1L)))

  if (nrow(df) == 0L) {
    message("No resource data for ", dt, " ", ti)
    return(NULL)
  }

  y_col <- if (metric == "cpu") "cpu_hours" else "memory_GB"
  y_lab <- if (metric == "cpu") "Hours" else "GB"
  plot_title_base <- if (metric == "cpu") "CPU Time" else "Max. RAM"
  plot_title <- if (isTRUE(include_context)) {
    paste(plot_title_base, data_type, tissue, sep = "; ")
  } else {
    plot_title_base
  }

  df <- df %>%
    dplyr::arrange(dplyr::desc(.data[[y_col]])) %>%
    dplyr::mutate(x_pos = dplyr::row_number())

  sc <- RColorConesa::scale_fill_conesa(palette = "complete")
  base5 <- sc$palette(5)
  my6 <- c(base5, "#FDC659")
  tool_cols <- stats::setNames(
    my6[seq_along(PAPER_RESOURCE_RANKING_TOOLS)],
    PAPER_RESOURCE_RANKING_TOOLS
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x_pos, y = .data[[y_col]])) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$tool),
      width = PAPER_RESOURCE_RANKING_BAR_WIDTH,
      color = NA,
      show.legend = TRUE
    ) +
    ggpattern::geom_col_pattern(
      data = df,
      ggplot2::aes(
        fill = .data$tool,
        pattern = factor(
          .data$strategy,
          levels = c("Join&Call", "Call&Join")
        )
      ),
      width = PAPER_RESOURCE_RANKING_BAR_WIDTH,
      pattern_colour = NA,
      pattern_fill = "white",
      pattern_alpha = 0.25,
      pattern_angle = 45,
      pattern_spacing = 0.04,
      pattern_density = 0.4,
      pattern_key_scale_factor = 0.6,
      show.legend = c(fill = FALSE, pattern = TRUE)
    ) +
    ggplot2::labs(title = plot_title, y = y_lab, x = NULL) +
    ggplot2::scale_fill_manual(
      values = tool_cols,
      name = "Tool",
      breaks = PAPER_RESOURCE_RANKING_TOOLS
    ) +
    ggpattern::scale_pattern_manual(
      name = "Strategy",
      values = c("Join&Call" = "none", "Call&Join" = "stripe"),
      breaks = c("Join&Call", "Call&Join")
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_x_continuous(
      limits = c(0.5, max(df$x_pos) + 0.5),
      expand = ggplot2::expansion(mult = c(0, 0)),
      breaks = df$x_pos,
      labels = NULL
    ) +
    paper_panel_theme() +
    ggplot2::theme(
      plot.title.position = "plot",
      legend.position = "right",
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 5, unit = "pt")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1),
      pattern = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          fill = "grey50",
          pattern_fill = "white",
          pattern_alpha = 0.25
        )
      )
    )
}

#' Strip y-axis title on right-hand panels (shared ylab on left column).
resource_ranking_strip_ylab <- function(p) {
  if (inherits(p, "ggplot")) {
    p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  } else {
    p
  }
}

#' 2×2 resource ranking grid: PacBio CPU|RAM / spacer / ONT CPU|RAM + shared legend.
assemble_resource_ranking_grid <- function(
    pb_cpu,
    pb_ram,
    ont_cpu,
    ont_ram,
    title) {
  row_pacbio <- pb_cpu | pb_ram
  row_ont <- (ont_cpu + ggplot2::guides(fill = "none")) |
    (ont_ram + ggplot2::guides(fill = "none"))

  figure <- row_pacbio / patchwork::plot_spacer() / row_ont +
    patchwork::plot_layout(
      guides = "collect",
      heights = c(1, 0.08, 1)
    )

  ann <- paper_figure_annotation(title = title, x_label = NULL)
  if (!is.null(ann)) {
    figure <- figure + ann
  }

  figure &
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(
        size = .paper_font("legend_title"),
        face = "bold"
      ),
      legend.text = ggplot2::element_text(size = .paper_font("legend_text"))
    )
}

#' Per-tool CPU / RAM bars (legacy 3×4 grid; uses [paper_panel_theme]).
create_resource_tool_plots <- function(
    data_type,
    tissue,
    tool,
    data,
    cpu_ymax,
    ram_ymax,
    show_x_labels = FALSE) {
  subset_data <- data %>%
    dplyr::filter(
      .data$data_type == !!data_type,
      .data$tissue == !!tissue,
      .data$tool == !!tool
    )

  if (nrow(subset_data) == 0L) {
    message("No resource data for ", data_type, " ", tissue, " ", tool)
    return(NULL)
  }

  subset_data <- subset_data %>%
    dplyr::mutate(
      elapsed_hours = vapply(.data$elapsed_time, parse_time_to_hours, numeric(1L)),
      cpu_hours = vapply(.data$cpu_time, parse_time_to_hours, numeric(1L)),
      strategy = factor(.data$strategy, levels = sort(unique(.data$strategy)))
    )

  num_strats <- nlevels(subset_data$strategy)
  conesa_cols <- rep(PAPER_STRATEGY_FILL, length.out = num_strats)

  panel_theme <- function(p, show_x) {
    p <- p +
      paper_panel_theme() +
      ggplot2::theme(
        legend.position = "none",
        plot.title.position = "plot",
        axis.title.x = ggplot2::element_blank()
      )
    if (!isTRUE(show_x)) {
      p <- p +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    } else {
      p <- p + ggplot2::theme(axis.text.x = paper_axis_text_x(0))
    }
    p
  }

  p_cpu <- ggplot2::ggplot(
    subset_data,
    ggplot2::aes(x = .data$strategy, y = .data$cpu_hours, fill = .data$strategy)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(title = "CPU", x = NULL, y = "Hours") +
    ggplot2::scale_y_continuous(
      limits = c(0, cpu_ymax),
      breaks = seq(0, min(250, cpu_ymax), by = 50)
    ) +
    ggplot2::scale_fill_manual(values = conesa_cols) +
    panel_theme(show_x_labels)

  p_ram <- ggplot2::ggplot(
    subset_data,
    ggplot2::aes(x = .data$strategy, y = .data$memory_GB, fill = .data$strategy)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(title = "Max. RAM", x = NULL, y = "GB") +
    ggplot2::scale_y_continuous(
      limits = c(0, ram_ymax),
      breaks = seq(0, min(200, ram_ymax), by = 50)
    ) +
    ggplot2::scale_fill_manual(values = conesa_cols) +
    panel_theme(show_x_labels)

  title_grob <- grid::textGrob(
    label = tool,
    x = 0.5,
    hjust = 0.5,
    gp = grid::gpar(
      fontsize = .paper_font("panel") / ggplot2::.pt,
      fontface = "bold"
    )
  )
  title_el <- patchwork::wrap_elements(full = title_grob)

  title_el / (p_cpu + p_ram) +
    patchwork::plot_layout(heights = c(0.08, 1))
}
