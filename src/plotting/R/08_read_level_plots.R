# 08_read_level_plots.R
# Shared read-level figure styling (mouse NIH + ONT R10 SY5Y).
# Depends on: 00_figure_config.R, 02_themes.R

read_level_structural_cat_palette <- c(
  "FSM" = "#6BAED6",
  "ISM" = "#FC8D59",
  "NIC" = "#78C679",
  "NNC" = "#EE6A50",
  "Genic\nGenomic" = "#969696",
  "Antisense" = "#66C2A4",
  "Fusion" = "goldenrod1",
  "Intergenic" = "darksalmon",
  "Genic\nIntron" = "#41B6C4"
)

read_level_structural_category_map <- c(
  "full-splice_match" = "FSM",
  "incomplete-splice_match" = "ISM",
  "novel_in_catalog" = "NIC",
  "novel_not_in_catalog" = "NNC",
  "genic" = "Genic\nGenomic",
  "antisense" = "Antisense",
  "fusion" = "Fusion",
  "intergenic" = "Intergenic",
  "genic_intron" = "Genic\nIntron"
)

read_level_mouse_display_levels <- c(
  "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion",
  "Intergenic", "Genic\nIntron", "Unaligned"
)

read_level_sqanti_structural_legend_labels <- c(
  "Full\nSplice Match",
  "Incomplete\nSplice Match",
  "Novel\nIn Catalog",
  "Novel Not\nIn Catalog",
  "Genic\nGenomic",
  "Antisense",
  "Fusion",
  "Intergenic",
  "Genic\nIntron"
)

read_level_sqanti_legend_breaks <- function(display_levels, include_unstranded = FALSE) {
  extra <- if (isTRUE(include_unstranded)) {
    c("Unaligned", "Unstranded")
  } else {
    "Unaligned"
  }
  c(extra, setdiff(display_levels, extra))
}

read_level_sqanti_legend_labels <- function(include_unstranded = FALSE) {
  extra <- if (isTRUE(include_unstranded)) {
    c("Unaligned", "Unstranded")
  } else {
    "Unaligned"
  }
  c(extra, read_level_sqanti_structural_legend_labels)
}

read_level_sqanti_fill_values <- function(
    palette = read_level_structural_cat_palette,
    include_unstranded = FALSE) {
  out <- c(palette, "Unaligned" = "white")
  if (isTRUE(include_unstranded)) {
    out <- c(palette, "Unstranded" = "grey85", "Unaligned" = "white")
  }
  out
}

#' Faceted SQANTI-reads stacked bars (mouse + SY5Y; shared colours and Unaligned outline).
plot_read_level_sqanti_stacked <- function(
    counts_df,
    display_levels,
    title = NULL,
    include_unstranded = FALSE,
    show_technology_axis = TRUE,
    technology_labels = c(pacbio = "PacBio", ont = "ONT"),
    facet_col = sample) {
  facet_sym <- rlang::ensym(facet_col)
  fill_values <- read_level_sqanti_fill_values(
    include_unstranded = include_unstranded
  )
  paper_plot_sqanti_faceted(
    counts_df,
    fill_values = fill_values,
    legend_breaks = read_level_sqanti_legend_breaks(
      display_levels,
      include_unstranded = include_unstranded
    ),
    legend_labels = read_level_sqanti_legend_labels(
      include_unstranded = include_unstranded
    ),
    title = title,
    technology_labels = technology_labels,
    show_technology_axis = show_technology_axis
  ) +
    ggplot2::facet_grid(cols = ggplot2::vars(!!facet_sym), switch = "x")
}

read_level_mouse_readnum_palette <- function() {
  stats::setNames(
    RColorConesa::colorConesa(n = 2L, palette = "complete"),
    c("Aligned", "Unaligned")
  )
}

read_level_sy5y_readnum_palette <- function() {
  classified_fill <- RColorConesa::colorConesa(n = 2L, palette = "complete")[1L]
  c(
    "Classified" = classified_fill,
    "Unstranded" = "grey85",
    "Unaligned" = "white"
  )
}

#' Stacked read-number bars; white Unaligned segment only gets a black outline.
paper_plot_readnum_stacked <- function(
    tidy_counts,
    fill_values,
    legend_breaks,
    legend_labels,
    facet_col,
    title = NULL,
    technology_labels = c(pacbio = "PacBio", ont = "ONT"),
    show_technology_axis = TRUE,
    stack_order = NULL,
    unaligned_label = "Unaligned",
    bar_width = 0.7) {
  facet_sym <- rlang::ensym(facet_col)
  outline_values <- paper_sqanti_stack_outline_values(
    fill_values,
    unaligned_label = unaligned_label
  )

  if (is.null(stack_order)) {
    tidy_counts <- dplyr::mutate(
      tidy_counts,
      .stack_order = as.numeric(
        as.character(.data$assign_category) == unaligned_label
      )
    )
  } else {
    tidy_counts <- dplyr::mutate(tidy_counts, .stack_order = stack_order)
  }

  p <- ggplot2::ggplot(
    tidy_counts,
    ggplot2::aes(
      x = .data$technology,
      y = .data$num_reads,
      fill = .data$assign_category,
      colour = .data$assign_category,
      order = .data$.stack_order
    )
  ) +
    ggplot2::geom_col(
      width = bar_width,
      position = ggplot2::position_stack(reverse = TRUE),
      linewidth = 0.15
    ) +
    ggplot2::scale_colour_manual(values = outline_values, guide = "none") +
    paper_read_count_y_scale() +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = legend_breaks,
      labels = legend_labels
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(!!facet_sym), switch = "x") +
    ggplot2::labs(x = NULL, y = "Number of reads", fill = "Category", title = title) +
    paper_read_level_theme()

  if (isTRUE(show_technology_axis)) {
    p <- p + ggplot2::scale_x_discrete(labels = technology_labels)
  } else {
    p <- p +
      ggplot2::scale_x_discrete(labels = NULL) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  p
}

#' Mouse NIH: PacBio + ONT read-number stack per replicate.
plot_mouse_readnum_faceted <- function(tidy_counts, tissue) {
  assign_palette <- read_level_mouse_readnum_palette()
  paper_plot_readnum_stacked(
    tidy_counts,
    fill_values = assign_palette,
    legend_breaks = c("Unaligned", "Aligned"),
    legend_labels = c("Unaligned", "Aligned"),
    facet_col = sample_label,
    title = tissue,
    stack_order = as.numeric(tidy_counts$assign_category == "Unaligned")
  )
}

#' BAM read-length violins (99.5% trim); mouse keeps PacBio/ONT x labels, SY5Y ONT-only hides them.
paper_plot_lengths_violin <- function(
    lengths_df,
    sample_levels = NULL,
    tissue = NULL,
    title = NULL,
    technology_labels = c(pacbio = "PacBio", ont = "ONT"),
    show_technology_axis = TRUE,
    trim_quantile = 0.995) {
  df <- lengths_df
  if (!is.null(tissue)) {
    df <- dplyr::filter(df, grepl(tissue, .data$sample))
    if (is.null(sample_levels)) {
      sample_levels <- paste(tissue, 1:5)
    }
  }
  df <- df %>%
    dplyr::group_by(.data$sample, .data$technology) %>%
    dplyr::filter(
      .data$length <= stats::quantile(.data$length, trim_quantile, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  if (!is.null(sample_levels)) {
    df <- dplyr::mutate(df, sample = factor(.data$sample, levels = sample_levels))
  }

  violin_fill <- RColorConesa::colorConesa(n = 3L, palette = "complete")[3L]
  y_max <- max(df$length, na.rm = TRUE)
  breaks_1k <- seq(0, ceiling(y_max / 1000) * 1000, by = 1000)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$technology,
      y = .data$length,
      fill = .data$technology
    )
  ) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$sample), switch = "x") +
    ggplot2::scale_y_continuous(breaks = breaks_1k) +
    ggplot2::scale_fill_manual(
      values = c(pacbio = violin_fill, ont = violin_fill)
    ) +
    ggplot2::labs(x = NULL, y = "Read length (nt)", fill = "Technology", title = title) +
    paper_read_level_theme() +
    ggplot2::theme(legend.position = "none")

  if (isTRUE(show_technology_axis)) {
    # Angle the PacBio/ONT labels so adjacent narrow-facet labels don't merge.
    p <- p +
      ggplot2::scale_x_discrete(labels = technology_labels) +
      ggplot2::theme(axis.text.x = paper_axis_text_x(45))
  } else {
    p <- p +
      ggplot2::scale_x_discrete(labels = NULL) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }

  p
}
