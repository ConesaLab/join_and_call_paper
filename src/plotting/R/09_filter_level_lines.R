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

  legend_plot <- bar_1_1 +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1))
  shared_legend <- wrap_elements(full = cowplot::get_legend(legend_plot))

  no_x_title <- theme(axis.title.x = element_blank())
  figure <- (
    (bar_1_1 + no_x_title) + (bar_1_2 + no_x_title) +
    (bar_1_3 + no_x_title) + (bar_1_4 + no_x_title) +
    (bar_2_1 + no_x_title) + (bar_2_2 + no_x_title) +
    shared_legend +
    plot_spacer() +
    bar_3_1 + bar_3_2 + bar_3_3 + bar_3_4
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKM
      ",
      heights = c(1, 1, 0.15, 1)
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
    inset_element(label_b, left = -0.02, bottom = 0.30, right = 0.05, top = 0.45)

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
    guides(color = guide_legend(order = 1, nrow = 1))
  color_legend <- cowplot::get_legend(color_legend_plot)

  strategy_df <- data.frame(
    x = c(1, 2, 1, 2), y = c(1, 1, 2, 2),
    strategy = rep(c("J&C", "C&J"), each = 2)
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
    guides(linetype = guide_legend(order = 2, nrow = 1),
           shape    = guide_legend(order = 2, nrow = 1))
  linetype_legend <- cowplot::get_legend(strategy_legend_plot)

  combined_legend <- (wrap_elements(full = color_legend) |
                      wrap_elements(full = linetype_legend)) +
    plot_layout(widths = c(3, 1))
  shared_legend <- wrap_elements(full = combined_legend)

  no_x_title <- theme(axis.title.x = element_blank())
  figure <- (
    (line_1_1 + no_x_title) + (line_1_2 + no_x_title) +
    (line_1_3 + no_x_title) + (line_1_4 + no_x_title) +
    (line_2_1 + no_x_title) + (line_2_2 + no_x_title) +
    shared_legend +
    plot_spacer() +
    line_3_1 + line_3_2 + line_3_3 + line_3_4
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKM
      ",
      heights = c(1, 1, 0.2, 1),
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
    inset_element(label_b, left = -0.02, bottom = 0.30, right = 0.05, top = 0.45)

  final_plot
}
