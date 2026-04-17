# 09_filter_level_lines.R
# Line plots comparing J&C vs C&J transcript counts across FL filter levels
# Depends on: 01_config.R, 02_themes.R (cat.palette), 04_data_loading.R, patchwork, tidyverse

simplified_cat_palette <- c(
  "FSM" = "#6BAED6", "ISM" = "#FC8D59", "NIC" = "#78C679",
  "NNC" = "#EE6A50", "Other" = "#969696"
)

simplify_category <- function(cat) {
  ifelse(cat %in% c("FSM", "ISM", "NIC", "NNC"), cat, "Other")
}


build_filter_level_data <- function(df_lists, tissue,
                                    fl_levels = c(1, 3, 5, 10, 20, 30, 40, 50),
                                    simplify = TRUE) {

  df_list_key <- if (tissue == "brain") "Bclass_df_list" else "Kclass_df_list"
  class_df_list <- df_lists[[df_list_key]]

  jc_raw <- class_df_list[["Join&Call"]]
  cj_raw <- class_df_list[["Call&Join"]]

  cat_levels <- if (simplify) names(simplified_cat_palette) else names(cat.palette)

  map_dfr(fl_levels, function(fl) {
    jc_df <- jc_raw %>% filter(!is.na(FL) & FL >= fl)
    cj_df <- cj_raw %>% filter(!is.na(FL) & FL >= fl)

    bind_rows(
      jc_df %>% mutate(strategy = "J&C"),
      cj_df %>% mutate(strategy = "C&J")
    ) %>%
      mutate(category = if (simplify) {
        factor(simplify_category(as.character(structural_category)),
               levels = cat_levels)
      } else {
        factor(as.character(structural_category), levels = cat_levels)
      }) %>%
      count(strategy, category, name = "n_transcripts") %>%
      mutate(fl_filter = fl)
  })
}


create_filter_level_line_plot <- function(plot_data, title, ylim_max = NULL) {

  jc_data <- plot_data %>% filter(strategy == "J&C")
  cj_data <- plot_data %>% filter(strategy == "C&J")

  p <- ggplot(plot_data,
         aes(x = factor(fl_filter), y = n_transcripts,
             color = category,
             group = interaction(category, strategy))) +
    geom_line(data = jc_data, linewidth = 0.8) +
    geom_line(data = cj_data, linewidth = 0.8, linetype = "dashed") +
    geom_point(data = jc_data, size = 2, shape = 16) +
    geom_point(data = cj_data, size = 2, shape = 17) +
    scale_color_manual(values = simplified_cat_palette, name = "Category") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    labs(x = "FL filter (\u2265)", y = "# transcripts", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )

  if (!is.null(ylim_max)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_max))
  }
  p
}


create_filter_level_bar_plot <- function(plot_data, title) {

  plot_data <- plot_data %>%
    mutate(fl_label = factor(paste0("\u2265", fl_filter),
                             levels = paste0("\u2265", sort(unique(fl_filter)))))

  ggplot(plot_data, aes(x = strategy, y = n_transcripts, fill = category)) +
    geom_col(width = 0.7) +
    facet_wrap(~ fl_label, nrow = 1) +
    scale_x_discrete(limits = c("J&C", "C&J")) +
    scale_fill_manual(values = cat.palette, name = "Structural Category") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    labs(x = "Strategy", y = "# transcripts", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      axis.title   = element_text(size = 10),
      strip.text   = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
}


create_filter_level_pct_bar_plot <- function(plot_data, title) {

  plot_data <- plot_data %>%
    mutate(fl_label = factor(paste0("\u2265", fl_filter),
                             levels = paste0("\u2265", sort(unique(fl_filter)))))

  ggplot(plot_data, aes(x = strategy, y = n_transcripts, fill = category)) +
    geom_col(position = "fill", width = 0.7) +
    facet_wrap(~ fl_label, nrow = 1) +
    scale_x_discrete(limits = c("J&C", "C&J")) +
    scale_fill_manual(values = cat.palette, name = "Structural Category") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Strategy", y = "% transcripts", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      axis.title   = element_text(size = 10),
      strip.text   = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
}


# --- TPM-normalized continuous filter functions ---

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


create_tpm_line_plot <- function(curve_data, title, log_x = FALSE,
                                 ylim_max = NULL, xlim = NULL) {

  jc_data <- curve_data %>% filter(strategy == "J&C")
  cj_data <- curve_data %>% filter(strategy == "C&J")

  p <- ggplot(curve_data,
         aes(x = tpm_threshold, y = n_transcripts,
             color = category,
             group = interaction(category, strategy))) +
    geom_line(data = jc_data, linewidth = 0.8) +
    geom_line(data = cj_data, linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(values = simplified_cat_palette, name = "Category") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",",
                                                    scientific = FALSE)) +
    labs(x = "TPM threshold (\u2265)", y = "# transcripts", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )

  if (log_x) p <- p + scale_x_log10()

  coord_args <- list()
  if (!is.null(xlim)) coord_args$xlim <- xlim
  if (!is.null(ylim_max)) coord_args$ylim <- c(0, ylim_max)
  if (length(coord_args) > 0) p <- p + do.call(coord_cartesian, coord_args)
  p
}


create_tpm_bar_plot <- function(curve_data, title, quantile_names = NULL) {

  thresholds_sorted <- sort(unique(curve_data$tpm_threshold))
  if (is.null(quantile_names)) {
    quantile_names <- paste0("Q", seq_along(thresholds_sorted))
  }

  curve_data <- curve_data %>%
    mutate(tpm_label = factor(
      quantile_names[match(tpm_threshold, thresholds_sorted)],
      levels = quantile_names
    ))

  ggplot(curve_data, aes(x = strategy, y = n_transcripts, fill = category)) +
    geom_col(width = 0.7) +
    facet_wrap(~ tpm_label, nrow = 1) +
    scale_x_discrete(limits = c("J&C", "C&J")) +
    scale_fill_manual(values = cat.palette, name = "Structural Category") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",",
                                                    scientific = FALSE)) +
    labs(x = "Strategy", y = "# transcripts", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      axis.title   = element_text(size = 10),
      strip.text   = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
}


create_tpm_pct_bar_plot <- function(curve_data, title, quantile_names = NULL) {

  thresholds_sorted <- sort(unique(curve_data$tpm_threshold))
  if (is.null(quantile_names)) {
    quantile_names <- paste0("Q", seq_along(thresholds_sorted))
  }

  curve_data <- curve_data %>%
    mutate(tpm_label = factor(
      quantile_names[match(tpm_threshold, thresholds_sorted)],
      levels = quantile_names
    ))

  ggplot(curve_data, aes(x = strategy, y = n_transcripts, fill = category)) +
    geom_col(position = "fill", width = 0.7) +
    facet_wrap(~ tpm_label, nrow = 1) +
    scale_x_discrete(limits = c("J&C", "C&J")) +
    scale_fill_manual(values = cat.palette, name = "Structural Category") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Strategy", y = "% transcripts", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      axis.title   = element_text(size = 10),
      strip.text   = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
}


assemble_tpm_line_figure <- function(all_line_plots, title) {

  line_1_1 <- all_line_plots$IsoSeq$IsoQuant
  line_1_2 <- all_line_plots$IsoSeq$FLAIR
  line_1_3 <- all_line_plots$IsoSeq$Bambu
  line_1_4 <- all_line_plots$IsoSeq$TALON
  line_2_1 <- all_line_plots$IsoSeq$Mandalorion
  line_2_2 <- all_line_plots$IsoSeq$isoseq_sqanti

  line_3_1 <- all_line_plots$ONT$IsoQuant
  line_3_2 <- all_line_plots$ONT$FLAIR
  line_3_3 <- all_line_plots$ONT$Bambu
  line_3_4 <- all_line_plots$ONT$TALON

  color_legend_plot <- line_1_1 +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 2, ncol = 3,
                                title = "Structural Category",
                                title.position = "top"))
  color_legend <- cowplot::get_legend(color_legend_plot)

  strategy_df <- data.frame(
    x = c(1, 2, 1, 2), y = c(1, 1, 2, 2),
    strategy = factor(c("J&C", "C&J", "J&C", "C&J"),
                      levels = c("J&C", "C&J"))
  )
  strategy_legend_plot <- ggplot(strategy_df,
      aes(x, y, linetype = strategy)) +
    geom_line() +
    scale_linetype_manual(values = c("J&C" = "solid", "C&J" = "dashed"),
                          name = "Strategy") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(linetype = guide_legend(nrow = 2, title.position = "top"))
  linetype_legend <- cowplot::get_legend(strategy_legend_plot)

  combined_legend <- (wrap_elements(full = linetype_legend) |
                      wrap_elements(full = color_legend)) +
    plot_layout(widths = c(1, 2))
  shared_legend <- wrap_elements(full = combined_legend)

  no_x_title <- theme(axis.title.x = element_blank())
  x_title_grob <- wrap_elements(full = grid::textGrob(
    "TPM threshold (\u2265)", gp = grid::gpar(fontsize = 12)
  ))

  figure <- (
    (line_1_1 + no_x_title) + (line_1_2 + no_x_title) +
    (line_1_3 + no_x_title) + (line_1_4 + no_x_title) +
    (line_2_1 + no_x_title) + (line_2_2 + no_x_title) +
    shared_legend +
    plot_spacer() +
    (line_3_1 + no_x_title) + (line_3_2 + no_x_title) +
    (line_3_3 + no_x_title) + (line_3_4 + no_x_title) +
    x_title_grob
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKL
        MMMM
      ",
      heights = c(1, 1, 0.4, 1, 0.05),
      axes = "collect_y"
    ) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
      )
    )

  label_a <- ggplot() +
    annotate("text", x = 1, y = 1, label = "a", hjust = 0.5, vjust = 1,
             size = 10, fontface = "bold") +
    theme_void()

  label_b <- ggplot() +
    annotate("text", x = 1, y = 1, label = "b", hjust = 0.5, vjust = 1,
             size = 10, fontface = "bold") +
    theme_void()

  final_plot <- wrap_elements(full = figure)
  final_plot <- final_plot +
    inset_element(label_a, left = -0.02, bottom = 0.85, right = 0.05, top = 0.95) +
    inset_element(label_b, left = -0.02, bottom = 0.28, right = 0.05, top = 0.43)

  final_plot
}


assemble_filter_level_bar_figure <- function(all_bar_plots, title) {

  bar_1_1 <- all_bar_plots$IsoSeq$IsoQuant
  bar_1_2 <- all_bar_plots$IsoSeq$FLAIR
  bar_1_3 <- all_bar_plots$IsoSeq$Bambu
  bar_1_4 <- all_bar_plots$IsoSeq$TALON
  bar_2_1 <- all_bar_plots$IsoSeq$Mandalorion
  bar_2_2 <- all_bar_plots$IsoSeq$isoseq_sqanti

  bar_3_1 <- all_bar_plots$ONT$IsoQuant
  bar_3_2 <- all_bar_plots$ONT$FLAIR
  bar_3_3 <- all_bar_plots$ONT$Bambu
  bar_3_4 <- all_bar_plots$ONT$TALON

  fill_legend_plot <- bar_1_1 +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 3, ncol = 3,
                               title = "Structural Category",
                               title.position = "top"))
  fill_legend <- cowplot::get_legend(fill_legend_plot)

  strategy_df <- data.frame(
    x = factor(c("J&C", "C&J"), levels = c("J&C", "C&J")),
    y = c(1, 1)
  )
  strategy_legend_plot <- ggplot(strategy_df, aes(x = x, y = y, fill = x)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("J&C" = "grey50", "C&J" = "grey80"),
                      name = "Strategy") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2, title.position = "top"))
  strategy_legend <- cowplot::get_legend(strategy_legend_plot)

  combined_legend <- (wrap_elements(full = strategy_legend) |
                      wrap_elements(full = fill_legend)) +
    plot_layout(widths = c(1, 2))
  shared_legend <- wrap_elements(full = combined_legend)

  no_x_title <- theme(axis.title.x = element_blank())
  x_title_grob <- wrap_elements(full = grid::textGrob(
    "Strategy", gp = grid::gpar(fontsize = 12)
  ))

  figure <- (
    (bar_1_1 + no_x_title) + (bar_1_2 + no_x_title) +
    (bar_1_3 + no_x_title) + (bar_1_4 + no_x_title) +
    (bar_2_1 + no_x_title) + (bar_2_2 + no_x_title) +
    shared_legend +
    plot_spacer() +
    (bar_3_1 + no_x_title) + (bar_3_2 + no_x_title) +
    (bar_3_3 + no_x_title) + (bar_3_4 + no_x_title) +
    x_title_grob
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKL
        MMMM
      ",
      heights = c(1, 1, 0.4, 1, 0.05)
    ) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
      )
    )

  label_a <- ggplot() +
    annotate("text", x = 1, y = 1, label = "a", hjust = 0.5, vjust = 1,
             size = 10, fontface = "bold") +
    theme_void()

  label_b <- ggplot() +
    annotate("text", x = 1, y = 1, label = "b", hjust = 0.5, vjust = 1,
             size = 10, fontface = "bold") +
    theme_void()

  final_plot <- wrap_elements(full = figure)
  final_plot <- final_plot +
    inset_element(label_a, left = -0.02, bottom = 0.85, right = 0.05, top = 0.95) +
    inset_element(label_b, left = -0.02, bottom = 0.28, right = 0.05, top = 0.43)

  final_plot
}


assemble_filter_level_figure <- function(all_line_plots, title) {

  line_1_1 <- all_line_plots$IsoSeq$IsoQuant
  line_1_2 <- all_line_plots$IsoSeq$FLAIR
  line_1_3 <- all_line_plots$IsoSeq$Bambu
  line_1_4 <- all_line_plots$IsoSeq$TALON
  line_2_1 <- all_line_plots$IsoSeq$Mandalorion
  line_2_2 <- all_line_plots$IsoSeq$isoseq_sqanti

  line_3_1 <- all_line_plots$ONT$IsoQuant
  line_3_2 <- all_line_plots$ONT$FLAIR
  line_3_3 <- all_line_plots$ONT$Bambu
  line_3_4 <- all_line_plots$ONT$TALON

  color_legend_plot <- line_1_1 +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 2, ncol = 3,
                                title = "Structural Category",
                                title.position = "top"))
  color_legend <- cowplot::get_legend(color_legend_plot)

  strategy_df <- data.frame(
    x = c(1, 2, 1, 2), y = c(1, 1, 2, 2),
    strategy = factor(c("J&C", "C&J", "J&C", "C&J"),
                      levels = c("J&C", "C&J"))
  )
  strategy_legend_plot <- ggplot(strategy_df,
      aes(x, y, linetype = strategy, shape = strategy)) +
    geom_line() +
    geom_point(size = 2) +
    scale_linetype_manual(values = c("J&C" = "solid", "C&J" = "dashed"),
                          name = "Strategy") +
    scale_shape_manual(values = c("J&C" = 16, "C&J" = 17),
                       name = "Strategy") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(linetype = guide_legend(nrow = 2, title.position = "top"),
           shape    = guide_legend(nrow = 2, title.position = "top"))
  linetype_legend <- cowplot::get_legend(strategy_legend_plot)

  combined_legend <- (wrap_elements(full = linetype_legend) |
                      wrap_elements(full = color_legend)) +
    plot_layout(widths = c(1, 2))
  shared_legend <- wrap_elements(full = combined_legend)

  no_x_title <- theme(axis.title.x = element_blank())
  x_title_grob <- wrap_elements(full = grid::textGrob(
    "FL filter level (\u2265)", gp = grid::gpar(fontsize = 12)
  ))

  figure <- (
    (line_1_1 + no_x_title) + (line_1_2 + no_x_title) +
    (line_1_3 + no_x_title) + (line_1_4 + no_x_title) +
    (line_2_1 + no_x_title) + (line_2_2 + no_x_title) +
    shared_legend +
    plot_spacer() +
    (line_3_1 + no_x_title) + (line_3_2 + no_x_title) +
    (line_3_3 + no_x_title) + (line_3_4 + no_x_title) +
    x_title_grob
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKL
        MMMM
      ",
      heights = c(1, 1, 0.4, 1, 0.05),
      axes = "collect_y"
    ) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
      )
    )

  label_a <- ggplot() +
    annotate("text", x = 1, y = 1, label = "a", hjust = 0.5, vjust = 1,
             size = 10, fontface = "bold") +
    theme_void()

  label_b <- ggplot() +
    annotate("text", x = 1, y = 1, label = "b", hjust = 0.5, vjust = 1,
             size = 10, fontface = "bold") +
    theme_void()

  final_plot <- wrap_elements(full = figure)
  final_plot <- final_plot +
    inset_element(label_a, left = -0.02, bottom = 0.85, right = 0.05, top = 0.95) +
    inset_element(label_b, left = -0.02, bottom = 0.28, right = 0.05, top = 0.43)

  final_plot
}
