# 02_themes.R
# Color palettes, structural category labels, and reusable ggplot themes.
# Depends on: 00_figure_config.R (PAPER_FONTS) when sourced after it.

xaxislevelsF1 <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron")
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion", "Intergenic", "Genic\nIntron")
cat.palette <- c("FSM" = "#6BAED6", "ISM" = "#FC8D59", "NIC" = "#78C679", "NNC" = "#EE6A50", "Genic\nGenomic" = "#969696", "Antisense" = "#66C2A4", "Fusion" = "goldenrod1", "Intergenic" = "darksalmon", "Genic\nIntron" = "#41B6C4")

.paper_font <- function(name) {
  if (!exists("PAPER_FONTS", inherits = TRUE)) {
    stop("Source 00_figure_config.R before 02_themes.R.", call. = FALSE)
  }
  scale <- if (exists("PAPER_SCALE", inherits = TRUE)) PAPER_SCALE else 1
  PAPER_FONTS[[name]] * scale
}

#' Publication font family. Nature requires Helvetica/Arial. Set here once so all
#' paper themes and text geoms inherit it; register on the system if missing
#' (e.g. systemfonts::register_font or ttf-mscorefonts-installer). Falls back to
#' the Cairo default sans only if Arial is unavailable.
PAPER_FONT_FAMILY <- "Arial"

#' Panel and device background (replaces `theme_minimal()` gray panel fill).
PAPER_BG_WHITE <- "white"

#' Light guideline colors (theme_minimal-style on white panels).
PAPER_GRID_MAJOR <- "#EBEBEB"
PAPER_GRID_MINOR <- "#F5F5F5"

#' Alternating metric-row fills on dumbbell plots (white / very light gray).
PAPER_BG_STRIPE_ALT <- "#F2F2F2"
PAPER_DUMBBELL_STRIPE_COLORS <- c(PAPER_BG_WHITE, PAPER_BG_STRIPE_ALT)

#' White panel/plot/strip backgrounds only (grid set separately).
paper_white_background_theme <- function() {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = PAPER_BG_WHITE, color = NA),
    plot.background  = ggplot2::element_rect(fill = PAPER_BG_WHITE, color = NA),
    strip.background = ggplot2::element_rect(fill = PAPER_BG_WHITE, color = NA)
  )
}

#' Grey major/minor guidelines on both axes (pre-refactor paper look).
paper_guideline_grid_theme <- function(minor = TRUE) {
  th <- ggplot2::theme(
    panel.grid.major = ggplot2::element_line(
      colour = PAPER_GRID_MAJOR,
      linewidth = 0.4
    ),
    panel.grid.minor = ggplot2::element_line(
      colour = PAPER_GRID_MINOR,
      linewidth = 0.25
    )
  )
  if (!isTRUE(minor)) {
    th <- th + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }
  th
}

#' Striped metric rows for dumbbell plots (white / light gray lanes).
paper_geom_stripped_rows <- function(
    mapping = NULL,
    data = NULL,
    ...,
    width = 1,
    nudge_y = 0) {
  ggstats::geom_stripped_rows(
    mapping = mapping,
    data = data,
    ...,
    width = width,
    nudge_y = nudge_y,
    odd = PAPER_DUMBBELL_STRIPE_COLORS[[1L]],
    even = PAPER_DUMBBELL_STRIPE_COLORS[[2L]]
  )
}

#' Default x breaks for dumbbell metric panels (0–100%).
PAPER_DUMBBELL_X_BREAKS <- seq(0, 100, by = 20)

#' Vertical guidelines for dumbbell plots.
#' Striped-row geoms paint over `panel.grid`, so guides are drawn as a geom on top.
paper_geom_dumbbell_vguides <- function(breaks = PAPER_DUMBBELL_X_BREAKS) {
  ggplot2::geom_vline(
    xintercept = breaks,
    colour = PAPER_GRID_MAJOR,
    linewidth = 0.4
  )
}

#' Shared x scale for dumbbell metric panels.
paper_dumbbell_x_scale <- function(
    limits = c(0, 100),
    breaks = PAPER_DUMBBELL_X_BREAKS) {
  ggplot2::scale_x_continuous(
    limits = limits,
    breaks = breaks,
    expand = ggplot2::expansion(mult = c(0.03, 0.03))
  )
}

#' Vertical x-axis grid via theme (backup; usually hidden under striped rows).
paper_dumbbell_grid_theme <- function() {
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(
      colour = PAPER_GRID_MAJOR,
      linewidth = 0.4
    ),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  )
}

#' Shared theme for single-panel dumbbell metric plots.
paper_dumbbell_theme <- function(
    legend_position = "right",
    legend_box = NULL,
    show_legend = TRUE) {
  f <- .paper_font
  th <- paper_theme() +
    paper_dumbbell_grid_theme() +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = f("panel")
      ),
      plot.margin        = ggplot2::margin(t = 5, r = 20, b = 5, l = 20, unit = "pt"),
      axis.text.y        = ggplot2::element_text(face = "bold", size = f("tick")),
      legend.title       = ggplot2::element_text(size = f("legend_title"), face = "bold"),
      legend.text        = ggplot2::element_text(size = f("legend_text")),
      legend.position    = if (isTRUE(show_legend)) legend_position else "none"
    )
  if (!is.null(legend_box)) {
    th <- th + ggplot2::theme(legend.box = legend_box)
  }
  th
}

#' Default x-axis tick angle for dense categorical labels (sample names, tools).
PAPER_AXIS_TICK_ANGLE <- 90L

#' Axis tick text (use `angle = 0` for horizontal labels, e.g. read-level tech facets).
paper_axis_text_x <- function(angle = PAPER_AXIS_TICK_ANGLE) {
  f <- .paper_font
  if (angle == 0) {
    return(ggplot2::element_text(size = f("tick"), angle = 0, hjust = 0.5, vjust = 0.5))
  }
  ggplot2::element_text(size = f("tick"), angle = angle, hjust = 1, vjust = 0.5)
}

paper_axis_text_y <- function() {
  ggplot2::element_text(size = .paper_font("tick"))
}

#' Base minimal theme with paper typography.
paper_theme <- function(base_size = 13, facet_dense = FALSE) {
  f <- .paper_font
  ggplot2::theme_minimal(base_size = base_size, base_family = PAPER_FONT_FAMILY) +
    paper_white_background_theme() +
    paper_guideline_grid_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = f("panel"), face = "bold"
      ),
      axis.text.x = paper_axis_text_x(0),
      axis.text.y = paper_axis_text_y(),
      axis.title = ggplot2::element_text(size = f("axis")),
      strip.text = ggplot2::element_text(size = f("strip"), face = "bold"),
      legend.text = ggplot2::element_text(size = f("legend_text")),
      legend.title = ggplot2::element_text(size = f("legend_title"))
    )
}

#' Per-panel typography (no coord limits).
paper_panel_theme <- function() {
  f <- .paper_font
  ggplot2::theme_minimal(base_family = PAPER_FONT_FAMILY) +
    paper_white_background_theme() +
    paper_guideline_grid_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = f("panel"), face = "bold"),
      axis.text.x = paper_axis_text_x(0),
      axis.text.y = paper_axis_text_y(),
      axis.title = ggplot2::element_text(size = f("axis")),
      legend.text = ggplot2::element_text(size = f("legend_text")),
      legend.title = ggplot2::element_text(size = f("legend_title")),
      strip.text = ggplot2::element_text(size = f("strip"), face = "bold")
    )
}

#' Patchwork figure-level title/subtitle.
paper_figure_title_theme <- function() {
  f <- .paper_font
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = f("figure"), face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = f("subtitle"), hjust = 0.5),
    plot.tag = ggplot2::element_text(size = f("tag"), face = "bold")
  )
}

#' Read-level figures: paper fonts + bordered panels/strips.
paper_read_level_theme <- function() {
  paper_panel_theme() +
    ggplot2::theme(
      # panel.border is drawn above geoms; fill must be transparent (ggplot2 docs).
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.5
      ),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(
        colour = "black", fill = PAPER_BG_WHITE, linewidth = 0.5
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )
}

#' Outline colours for stacked SQANTI bars (black border on Unaligned only).
paper_sqanti_stack_outline_values <- function(
    fill_values,
    unaligned_label = "Unaligned") {
  outline <- fill_values
  if (unaligned_label %in% names(outline)) {
    outline[[unaligned_label]] <- "black"
  }
  outline
}

#' Faceted SQANTI-reads structural-category stacked bars (mouse NIH + SY5Y).
#'
#' Only the white Unaligned segment gets a black outline; other categories use an
#' outline colour matching fill so segment borders are not visible.
paper_plot_sqanti_faceted <- function(
    counts_df,
    fill_values,
    legend_breaks,
    legend_labels,
    title = NULL,
    technology_labels = c(pacbio = "PacBio", ont = "ONT"),
    unaligned_label = "Unaligned",
    show_technology_axis = TRUE) {
  if (nrow(counts_df) == 0L || sum(counts_df$num_reads, na.rm = TRUE) <= 0) {
    stop("No read counts to plot (empty category_counts).", call. = FALSE)
  }
  if (any(is.na(counts_df$category_label))) {
    stop(
      sum(is.na(counts_df$category_label)),
      " rows have NA category_label (factor level mismatch).",
      call. = FALSE
    )
  }

  outline_values <- paper_sqanti_stack_outline_values(
    fill_values,
    unaligned_label = unaligned_label
  )

  p <- ggplot2::ggplot(
    counts_df,
    ggplot2::aes(
      x = .data$technology,
      y = .data$num_reads,
      fill = .data$category_label,
      colour = .data$category_label
    )
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      width = 0.8,
      position = ggplot2::position_stack(reverse = TRUE),
      linewidth = 0.15
    ) +
    ggplot2::scale_colour_manual(values = outline_values, guide = "none") +
    ggplot2::scale_fill_manual(
      name = "Structural Category",
      values = fill_values,
      breaks = legend_breaks,
      labels = legend_labels,
      drop = FALSE
    ) +
    paper_read_count_y_scale() +
    ggplot2::labs(x = NULL, y = "Number of reads", title = title) +
    paper_read_level_theme() +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical"
    )

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

#' Read-level count y-axis: scientific notation; limits and breaks are per-plot.
paper_read_count_y_scale <- function() {
  sci <- scales::label_scientific(digits = 1)
  ggplot2::scale_y_continuous(
    labels = function(x) ifelse(x == 0, "0", sci(x)),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  )
}

#' TPM survival line panels.
tpm_line_theme <- function() {
  paper_panel_theme() +
    ggplot2::theme(legend.position = "none")
}

#' TPM panel wrapper (SY5Y 2x2 and mouse 10-panel line tools).
tpm_panel_theme <- function(
    x_tick_angle = 0,
    x_tick_hjust = NULL,
    log_x = FALSE) {
  x_text <- paper_axis_text_x(x_tick_angle)
  x_theme <- if (isTRUE(log_x)) {
    ggplot2::theme(
      axis.text.x = x_text,
      axis.text.x.bottom = x_text,
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank()
    )
  } else {
    ggplot2::theme(
      axis.text.x = x_text,
      axis.text.x.bottom = x_text,
      axis.ticks.x = ggplot2::element_line(),
      axis.ticks.x.bottom = ggplot2::element_line()
    )
  }
  bottom_margin <- if (isTRUE(log_x)) 12 else 10
  tpm_line_theme() +
    x_theme +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        t = 5.5, r = 5.5, b = bottom_margin, l = 5.5, unit = "pt"
      )
    )
}

#' Shared x-axis label as a patchwork caption.
#' Figure-level titles are dropped for Nature Communications (the manuscript
#' caption carries the title); the `title` argument is accepted but ignored so
#' existing call sites need not change. Only the shared x-axis label is kept.
paper_figure_annotation <- function(title = NULL, x_label = NULL) {
  has_x <- length(x_label) && nzchar(x_label)
  if (!has_x) {
    return(NULL)
  }
  th <- ggplot2::theme(
    plot.caption = ggplot2::element_text(
      size = .paper_font("caption"),
      hjust = 0.5,
      margin = ggplot2::margin(t = 6)
    )
  )
  patchwork::plot_annotation(caption = x_label, theme = th)
}

#' Single centered x-axis title shared across all columns (patchwork caption).
#' Always returns a plot_annotation (so it can terminate a `... & theme` chain);
#' draws no caption when `x_label` is NULL/empty. The per-row y-axis titles are
#' kept on the panels themselves; this supplies the one shared x-axis title.
paper_figure_x_caption <- function(x_label = NULL) {
  patchwork::plot_annotation(
    caption = x_label,
    theme = ggplot2::theme(
      plot.caption = ggplot2::element_text(
        size = .paper_font("caption"),
        hjust = 0.5,
        margin = ggplot2::margin(t = 6)
      )
    )
  )
}

#' Patchwork design for mouse 10-panel figures (`all_plots.Rmd`).
PAPER_MOUSE_10PANEL_DESIGN_IJKM <- "
        ABCD
        EFGG
        HHHH
        IJKM
"

PAPER_MOUSE_10PANEL_DESIGN_IJKL <- "
        ABCD
        EFGG
        HHHH
        IJKL
"

#' @param heights Row heights (top → bottom), including spacer row.
#' @param design One of [PAPER_MOUSE_10PANEL_DESIGN_IJKM] or `_IJKL`.
paper_mouse_10panel_plot_layout <- function(
    heights,
    design = PAPER_MOUSE_10PANEL_DESIGN_IJKM,
    axes = "collect_y",
    axis_titles = "collect") {
  patchwork::plot_layout(
    design = design,
    heights = heights,
    axes = axes,
    axis_titles = axis_titles
  )
}

#' Theme applied to assembled patchwork figures (legend off; shared x title centered).
paper_figure_patchwork_theme <- function() {
  ggplot2::theme(
    legend.position = "none",
    axis.title.x = ggplot2::element_text(hjust = 0.5)
  )
}

#' UpSet / boxplot subpanels (replaces theme_Publication for paper figures).
paper_upset_theme <- function() {
  paper_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}

#' ComplexUpset sub-panel themes (matrix set labels at `matrix_label_size` or PAPER_FONTS$tick).
paper_upset_modify_themes <- function(matrix_label_size = NULL, extra = list()) {
  matrix_text <- if (is.null(matrix_label_size)) {
    paper_axis_text_y()
  } else {
    # matrix_label_size is a reference-design pt; scale to print like other fonts.
    ggplot2::element_text(size = matrix_label_size * PAPER_SCALE)
  }
  base <- list(
    # Zero the touching margins (bar bottom / matrix top) so the two stacked
    # sub-panels sit flush — default ggplot plot.margin otherwise leaves a
    # visible gap between the intersection-size bars and the matrix below.
    intersections_matrix = ggplot2::theme(
      axis.text.y = matrix_text,
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 1, b = 1, l = 1)
    ),
    `Intersection size` = ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 1, r = 1, b = 0, l = 1)
    )
  )
  if (length(extra)) {
    base <- c(base, extra)
  }
  ComplexUpset::upset_modify_themes(base)
}

#' @rdname paper_upset_theme
theme_Publication <- function(base_size = 14, base_family = "sans") {
  paper_upset_theme()
}

bar_theme <- function(plot, title, ylims) {
  plot +
    ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    ggplot2::ggtitle(title) +
    paper_panel_theme() +
    ggplot2::labs(x = NULL) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(ylim = ylims)
}

bar_theme_noaxis <- function(plot, title, ylims) {
  bar_theme(plot, title, ylims) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}

compare_theme <- function(plot, title) {
  plot +
    ggplot2::ggtitle(title) +
    paper_panel_theme() +
    ggplot2::theme(
      axis.text.x = paper_axis_text_x(),
      legend.position = "none"
    ) +
    ggplot2::labs(x = "Strategy")
}

compare_theme_noaxis <- function(plot, title) {
  compare_theme(plot, title) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}

#' NPC top for panel-tag placement. Figure-level titles are no longer drawn
#' (dropped for Nature), so tags sit near the very top regardless of the legacy
#' `has_figure_title` argument (kept for call-site compatibility).
paper_tag_content_top <- function(has_figure_title = TRUE) {
  0.97
}

#' Place a panel tag (a, b, …) at an arbitrary figure-relative position via inset.
#' Single unified tag mechanism: generalizes `paper_inset_tag_at()` (row-based)
#' to any (x, y) so side-by-side / rearranged panels can be tagged too. Replaces
#' the former `labs(tag=)` mechanism. Size comes from `.paper_font("tag")` (scaled).
paper_inset_tag <- function(plot, label, x = 0.0, y = 0.98,
                            box_w = 0.06, box_h = 0.06) {
  plot +
    patchwork::inset_element(
      paper_panel_tag_grob(label),
      left = x,
      bottom = y - box_h,
      right = x + box_w,
      top = y,
      align_to = "full"
    )
}

#' Panel tag grob for non-ggplot panels (ComplexUpset, etc.) via inset_element.
paper_panel_tag_grob <- function(label) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 1,
      label = label,
      hjust = 0,
      vjust = 1,
      size = .paper_font("tag") / ggplot2::.pt,
      fontface = "bold",
      family = PAPER_FONT_FAMILY
    ) +
    ggplot2::theme_void(base_family = PAPER_FONT_FAMILY)
}

#' Place one inset tag in a small top-left box at `y_top` (NPC, 0 = bottom).
paper_inset_tag_at <- function(
    plot,
    label,
    y_top,
    has_figure_title = TRUE) {
  y_cap <- paper_tag_content_top(has_figure_title)
  y_top <- min(y_top, y_cap)
  tag_h <- 0.075
  plot +
    patchwork::inset_element(
      paper_panel_tag_grob(label),
      left = -0.02,
      bottom = y_top - tag_h,
      right = 0.055,
      top = y_top
    )
}

#' Top edge (NPC) of row `row_from_top` in a vertical patchwork stack (`heights` top → bottom).
paper_row_top_npc <- function(row_from_top, heights, has_figure_title = TRUE) {
  h <- heights / sum(heights)
  y_max <- paper_tag_content_top(has_figure_title)
  if (row_from_top <= 1L) {
    return(y_max)
  }
  y_max * (1 - sum(h[seq_len(row_from_top - 1L)]))
}

#' Inset tags on specific rows of a multi-row figure (row 1 = top).
paper_inset_tags_rows <- function(
    plot,
    tag_by_row,
    heights,
    has_figure_title = TRUE) {
  out <- plot
  for (row_idx in names(tag_by_row)) {
    tag <- tag_by_row[[row_idx]]
    if (!nzchar(tag)) {
      next
    }
    y_top <- paper_row_top_npc(
      as.integer(row_idx),
      heights,
      has_figure_title = has_figure_title
    )
    out <- paper_inset_tag_at(
      out,
      tag,
      y_top,
      has_figure_title = FALSE
    )
  }
  out
}

#' Two stacked panels: `tags[1]` on top row, `tags[2]` on bottom row.
paper_inset_panel_tags <- function(
    plot,
    tags = c("a", "b"),
    heights = c(1, 1),
    has_figure_title = TRUE) {
  if (length(tags) < 1L) {
    return(plot)
  }
  tag_map <- c("1" = tags[1L])
  if (length(tags) >= 2L) {
    tag_map <- c(tag_map, setNames(tags[2L], as.character(length(heights))))
  }
  paper_inset_tags_rows(
    plot,
    tag_by_row = tag_map,
    heights = heights,
    has_figure_title = has_figure_title
  )
}
