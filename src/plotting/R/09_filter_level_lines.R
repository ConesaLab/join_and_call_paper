# 09_filter_level_lines.R
# TPM-normalized continuous filter (survival) curves for J&C vs C&J.
# Depends on: 00_figure_config.R, 02_themes.R (cat.palette), 04_data_loading.R, patchwork, tidyverse

simplified_cat_palette <- c(
  "FSM" = "#6BAED6", "ISM" = "#FC8D59", "NIC" = "#78C679",
  "NNC" = "#EE6A50", "Other" = "#969696"
)

simplify_category <- function(cat) {
  ifelse(cat %in% c("FSM", "ISM", "NIC", "NNC"), cat, "Other")
}

build_tpm_curve_data <- function(df_lists, tissue) {

  df_list_key <- if (tissue == "brain") "Bclass_df_list" else "Kclass_df_list"
  class_df_list <- df_lists[[df_list_key]]

  jc_raw <- class_df_list[["Join&Call"]]
  cj_raw <- class_df_list[["Call&Join"]]

  process_strategy <- function(df, strategy_label) {
    df <- df %>% filter(!is.na(FL) & FL > 0)
    total_fl <- sum(df$FL)
    df %>%
      mutate(
        TPM = (FL / total_fl) * 1e6,
        strategy = strategy_label,
        full_category = factor(as.character(structural_category),
                               levels = names(cat.palette)),
        simplified_category = factor(
          simplify_category(as.character(structural_category)),
          levels = names(simplified_cat_palette))
      ) %>%
      select(strategy, FL, TPM, full_category, simplified_category)
  }

  bind_rows(
    process_strategy(jc_raw, "J&C"),
    process_strategy(cj_raw, "C&J")
  )
}

compute_survival_curve <- function(tpm_data, thresholds, simplify = TRUE) {

  cat_col <- if (simplify) "simplified_category" else "full_category"
  cat_levels <- if (simplify) names(simplified_cat_palette) else names(cat.palette)

  groups <- tpm_data %>%
    group_by(strategy, .data[[cat_col]]) %>%
    summarise(tpm_sorted = list(sort(TPM)), .groups = "drop") %>%
    rename(category = !!sym(cat_col))

  groups %>%
    crossing(tpm_threshold = thresholds) %>%
    mutate(n_transcripts = purrr::map2_int(
      tpm_sorted, tpm_threshold, ~ sum(.x >= .y)
    )) %>%
    mutate(category = factor(category, levels = cat_levels)) %>%
    select(strategy, category, tpm_threshold, n_transcripts)
}

pick_tpm_quantiles <- function(tpm_data) {
  probs <- c(0.01, 0.10, 0.25, 0.50, 0.75)
  quantile(tpm_data$TPM, probs = probs)
}

tpm_x_axis_label <- function(log_x) {
  if (log_x) "TPM threshold (\u2265, log scale)" else "TPM threshold (\u2265)"
}

#' Major log-axis breaks at decades within the visible range.
tpm_log10_decade_breaks <- function(x) {
  rng <- range(x, na.rm = TRUE, finite = TRUE)
  if (any(!is.finite(rng)) || any(rng <= 0)) {
    return(numeric(0))
  }
  10^seq(ceiling(log10(rng[1])), floor(log10(rng[2])))
}

#' Minor log-axis breaks (2–9 between decades) for grid lines and ticks.
tpm_log10_minor_breaks <- function(xlim = NULL) {
  br <- as.numeric(outer(1:9, 10^(-3:2)))
  if (!is.null(xlim) && length(xlim) == 2L) {
    br <- br[br >= xlim[1] & br <= xlim[2]]
  }
  br
}

tpm_log10_axis_labels <- function(x) {
  formatC(x, drop0trailing = TRUE, format = "fg")
}

create_tpm_line_plot <- function(curve_data, title, log_x = FALSE,
                                 ylim_max = NULL, xlim = NULL) {

  jc_data <- curve_data %>% filter(strategy == "J&C")
  cj_data <- curve_data %>% filter(strategy == "C&J")

  p <- ggplot(curve_data,
         aes(x = tpm_threshold, y = n_transcripts,
             color = category,
             group = interaction(category, strategy))) +
    geom_line(data = jc_data, linewidth = 0.8 * PAPER_SCALE) +
    geom_line(data = cj_data, linewidth = 0.8 * PAPER_SCALE, linetype = "dashed") +
    scale_color_manual(values = simplified_cat_palette, name = "Category") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",",
                                                    scientific = FALSE)) +
    labs(x = tpm_x_axis_label(log_x), y = "# transcripts", title = title) +
    tpm_line_theme()

  if (log_x) {
    log_scale_args <- list(
      breaks       = tpm_log10_decade_breaks,
      labels       = tpm_log10_axis_labels,
      minor_breaks = tpm_log10_minor_breaks(xlim)
    )
    if (!is.null(xlim)) {
      log_scale_args$limits <- xlim
    }
    p <- p +
      do.call(ggplot2::scale_x_log10, log_scale_args) +
      ggplot2::annotation_logticks(sides = "b")
  }

  coord_args <- list()
  if (!is.null(xlim) && !isTRUE(log_x)) {
    coord_args$xlim <- xlim
  }
  if (!is.null(ylim_max)) {
    coord_args$ylim <- c(0, ylim_max)
  }
  if (length(coord_args) > 0) {
    p <- p + do.call(ggplot2::coord_cartesian, coord_args)
  }
  p
}

strategy_label_map <- c("J&C" = "Join & Call", "C&J" = "Call & Join")

category_label_map <- c(
  "FSM"   = "Full Splice Match",
  "ISM"   = "Incomplete Splice Match",
  "NIC"   = "Novel in Catalog",
  "NNC"   = "Novel Not in Catalog",
  "Other" = "Other"
)

.legend_font_sizes <- function(text_size = NULL, title_size = NULL) {
  if (is.null(text_size) && exists("PAPER_FONTS", inherits = TRUE)) {
    text_size <- PAPER_FONTS$legend_text
  }
  if (is.null(title_size) && exists("PAPER_FONTS", inherits = TRUE)) {
    title_size <- PAPER_FONTS$legend_title
  }
  list(
    text_size = if (is.null(text_size)) 13L else text_size,
    title_size = if (is.null(title_size)) 14L else title_size
  )
}

legend_only_theme <- function(text_size = NULL, title_size = NULL) {
  fs <- .legend_font_sizes(text_size, title_size)
  use_inside <- utils::packageVersion("ggplot2") >= "3.5.0"
  base <- theme_void(base_size = fs$text_size) +
    theme(
      legend.justification  = c(0.5, 0.5),
      legend.title          = element_text(size = fs$title_size, face = "bold"),
      legend.text           = element_text(size = fs$text_size),
      legend.key.width      = unit(2.5, "lines"),
      legend.key.height     = unit(1.5, "lines"),
      legend.background     = element_blank(),
      legend.box.background = element_blank(),
      plot.margin           = margin(4, 4, 4, 4)
    )
  if (use_inside) {
    base + theme(legend.position = "inside",
                 legend.position.inside = c(0.5, 0.5))
  } else {
    base + theme(legend.position = c(0.5, 0.5))
  }
}

build_tpm_strategy_legend <- function(text_size = NULL, title_size = NULL) {
  fs <- .legend_font_sizes(text_size, title_size)
  strategy_df <- data.frame(
    x = c(1, 2, 1, 2),
    y = NA_real_,
    strategy = factor(c("J&C", "J&C", "C&J", "C&J"),
                      levels = c("J&C", "C&J"))
  )
  ggplot(strategy_df, aes(x, y, linetype = strategy)) +
    geom_line(linewidth = 0.9 * PAPER_SCALE, na.rm = TRUE) +
    scale_linetype_manual(
      values = c("J&C" = "solid", "C&J" = "dashed"),
      labels = strategy_label_map,
      name   = "Strategy",
      drop   = FALSE
    ) +
    guides(linetype = guide_legend(
      nrow = 2, title.position = "top",
      override.aes = list(linewidth = 1.1 * PAPER_SCALE)
    )) +
    legend_only_theme(text_size = fs$text_size, title_size = fs$title_size)
}

build_tpm_category_legend <- function(nrow = 2, ncol = 3,
                                      text_size = NULL, title_size = NULL) {
  fs <- .legend_font_sizes(text_size, title_size)
  cat_levels <- names(simplified_cat_palette)
  cat_df <- data.frame(
    x = seq_along(cat_levels),
    y = NA_real_,
    category = factor(cat_levels, levels = cat_levels)
  )
  ggplot(cat_df, aes(x, y, fill = category)) +
    geom_tile(na.rm = TRUE) +
    scale_fill_manual(
      values = simplified_cat_palette,
      labels = category_label_map[cat_levels],
      name   = "Structural Category",
      drop   = FALSE
    ) +
    guides(fill = guide_legend(
      nrow = nrow, ncol = ncol, title.position = "top",
      override.aes = list(color = NA)
    )) +
    legend_only_theme(text_size = fs$text_size, title_size = fs$title_size)
}

build_tpm_line_legend_figure <- function(category_layout = c(1, 5),
                                         text_size = NULL, title_size = NULL) {
  cat_layout <- as.integer(category_layout)
  strategy_legend <- build_tpm_strategy_legend(
    text_size = text_size, title_size = title_size
  )
  category_legend <- build_tpm_category_legend(
    nrow = cat_layout[1], ncol = cat_layout[2],
    text_size = text_size, title_size = title_size
  )

  (strategy_legend | category_legend) + plot_layout(widths = c(1, 2))
}

filter_level_strip_ylab <- function(p) {
  if (inherits(p, "ggplot")) {
    p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  } else {
    p
  }
}

filter_level_mouse_figure_layout <- function() {
  patchwork::plot_layout(
    design = "
      ABCD
      EFGG
      HHHH
      IJKL
    ",
    heights = c(1, 1, 0.4, 1),
    guides = "keep",
    axes = "keep",
    axis_titles = "keep"
  )
}

assemble_tpm_line_figure <- function(all_line_plots, title, log_x = FALSE) {

  line_1_1 <- all_line_plots$IsoSeq$IsoQuant
  line_1_2 <- filter_level_strip_ylab(all_line_plots$IsoSeq$FLAIR)
  line_1_3 <- filter_level_strip_ylab(all_line_plots$IsoSeq$Bambu)
  line_1_4 <- filter_level_strip_ylab(all_line_plots$IsoSeq$TALON)
  line_2_1 <- all_line_plots$IsoSeq$Mandalorion
  line_2_2 <- filter_level_strip_ylab(all_line_plots$IsoSeq$isoseq_sqanti)

  line_3_1 <- all_line_plots$ONT$IsoQuant
  line_3_2 <- filter_level_strip_ylab(all_line_plots$ONT$FLAIR)
  line_3_3 <- filter_level_strip_ylab(all_line_plots$ONT$Bambu)
  line_3_4 <- filter_level_strip_ylab(all_line_plots$ONT$TALON)

  shared_legend <- wrap_elements(full = build_tpm_line_legend_figure(
    category_layout = c(2, 3)
  ))

  no_x_title <- theme(axis.title.x = element_blank())

  figure <- (
    (line_1_1 + no_x_title) + (line_1_2 + no_x_title) +
    (line_1_3 + no_x_title) + (line_1_4 + no_x_title) +
    (line_2_1 + no_x_title) + (line_2_2 + no_x_title) +
    shared_legend +
    plot_spacer() +
    (line_3_1 + no_x_title) + (line_3_2 + no_x_title) +
    (line_3_3 + no_x_title) + (line_3_4 + no_x_title)
  ) +
    filter_level_mouse_figure_layout()

  ann <- paper_figure_annotation(
    title = title,
    x_label = tpm_x_axis_label(log_x)
  )
  if (!is.null(ann)) {
    figure <- figure + ann
  }

  paper_inset_tags_rows(
    wrap_elements(full = figure),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.4, 1)
  )
}
