# 08_cross_tool_upset.R
# Cross-tool UJC upset/bar plots: compare UJC overlap between tools within a platform
# Depends on: 02_themes.R (cat.palette), ComplexUpset, patchwork

tool_display_names <- c(
  IsoQuant      = "IsoQuant",
  FLAIR         = "FLAIR",
  TALON         = "TALON",
  Bambu         = "Bambu",
  Mandalorion   = "Mandalorion",
  isoseq_sqanti = "IsoSeq + SQANTI3"
)

cross_tool_category_colors <- c(
  "FSM" = "#6BAED6", "ISM" = "#FC8D59", "NIC" = "#78C679", "NNC" = "#EE6A50",
  "Genic\nGenomic" = "#969696", "Antisense" = "#66C2A4", "Fusion" = "goldenrod1",
  "Intergenic" = "darksalmon", "Genic\nIntron" = "#41B6C4"
)


build_cross_tool_ujc_matrix <- function(all_results, platform, tissue,
                                        fl_filter_level = "1") {

  fl_filter_level <- as.character(fl_filter_level)
  df_list_key <- if (tissue == "brain") "b_class_df_list" else "k_class_df_list"
  tools <- names(all_results[[platform]])

  per_tool_ujcs <- list()
  for (tool in tools) {
    class_df_list <- all_results[[platform]][[tool]][[fl_filter_level]][[df_list_key]]
    individual_dfs <- class_df_list[3:7]

    tool_ujcs <- bind_rows(individual_dfs) %>%
      filter(!grepl("NA", UJC)) %>%
      select(structural_category, UJC) %>%
      distinct()

    display_name <- tool_display_names[[tool]]
    if (is.na(display_name)) display_name <- tool
    per_tool_ujcs[[display_name]] <- tool_ujcs
  }

  all_ujcs <- bind_rows(per_tool_ujcs) %>%
    distinct(UJC, .keep_all = TRUE)

  display_names <- names(per_tool_ujcs)
  for (dn in display_names) {
    all_ujcs[[dn]] <- as.integer(all_ujcs$UJC %in% per_tool_ujcs[[dn]]$UJC)
  }

  list(all_ujcs = all_ujcs, display_names = display_names)
}


create_cross_tool_upset <- function(all_results, platform, tissue,
                                    fl_filter_level = "1",
                                    n_intersections = 30) {

  mat <- build_cross_tool_ujc_matrix(all_results, platform, tissue, fl_filter_level)
  all_ujcs      <- mat$all_ujcs
  display_names <- mat$display_names

  intersection_size_annotation <- intersection_size(
    counts = FALSE,
    mapping = aes(fill = structural_category)
  ) +
    scale_fill_manual(values = cross_tool_category_colors, name = "Structural Category") +
    ggtitle(platform) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_text(size = 10),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )

  upset_plot <- ComplexUpset::upset(
    all_ujcs,
    intersect       = rev(display_names),
    sort_sets        = FALSE,
    name             = element_blank(),
    base_annotations = list(
      "Intersection size" = intersection_size_annotation
    ),
    n_intersections = n_intersections,
    width_ratio     = 0.1,
    set_sizes       = FALSE
  )

  upset_plot
}


create_cross_tool_bar <- function(all_results, platform, tissue,
                                  fl_filter_level = "1") {

  mat <- build_cross_tool_ujc_matrix(all_results, platform, tissue, fl_filter_level)
  all_ujcs      <- mat$all_ujcs
  display_names <- mat$display_names

  n_tools_total <- length(display_names)

  bar_data <- all_ujcs %>%
    mutate(n_tools = rowSums(across(all_of(display_names)))) %>%
    count(n_tools, structural_category, name = "n_ujc") %>%
    mutate(n_tools = factor(n_tools, levels = seq_len(n_tools_total)))

  ggplot(bar_data, aes(x = n_tools, y = n_ujc, fill = structural_category)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = cross_tool_category_colors, name = "Structural Category") +
    scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    labs(x = "# tools detecting UJC", y = "# UJCs", title = platform) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
}


assemble_cross_tool_figure <- function(isoseq_upset, ont_upset, title) {

  combined <- (
    wrap_elements(full = isoseq_upset) /
    wrap_elements(full = ont_upset)
  ) +
    plot_layout(heights = c(1.2, 1)) +
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

  final_plot <- wrap_elements(full = combined)
  final_plot <- final_plot +
    inset_element(label_a, left = -0.02, bottom = 0.94, right = 0.05, top = 1) +
    inset_element(label_b, left = -0.02, bottom = 0.35, right = 0.05, top = 0.45)

  final_plot
}


assemble_cross_tool_bar_figure <- function(isoseq_bar, ont_bar, title) {

  legend_plot <- isoseq_bar +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1))
  legend_grob <- cowplot::get_legend(legend_plot)
  shared_legend <- wrap_elements(full = legend_grob)

  combined <- (
    isoseq_bar + ont_bar + shared_legend
  ) +
    plot_layout(
      design = "A#B\nCCC",
      heights = c(1, 0.15),
      widths = c(1, 0.05, 1),
      axis_titles = "collect"
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

  final_plot <- wrap_elements(full = combined)
  final_plot <- final_plot +
    inset_element(label_a, left = -0.02, bottom = 0.88, right = 0.05, top = 0.98) +
    inset_element(label_b, left = 0.48,  bottom = 0.88, right = 0.55, top = 0.98)

  final_plot
}
