# 08_cross_tool_upset.R
# Cross-tool UJC upset/bar plots: compare UJC overlap between tools within a platform
# Depends on: 01_config.R, 04_data_loading.R, 02_themes.R, ComplexUpset, patchwork

#' Load filtered class lists for cross-tool figures (no process_and_plot).
load_cross_tool_all_results <- function(src_dirs, fl_filter_level = 1L) {
  fl <- as.character(fl_filter_level)
  all_results <- list()

  for (platform in names(src_dirs)) {
    all_results[[platform]] <- list()
    for (method in names(src_dirs[[platform]])) {
      all_results[[platform]][[method]] <- list()
      cat(platform, method, "\n")
      src_dir <- src_dirs[[platform]][[method]]
      paths <- get_paths(src_dir)
      df_lists <- load_df_lists(paths)
      filtered <- process_classification_only(
        df_lists$Bclass_df_list,
        df_lists$Kclass_df_list,
        fl_threshold = fl_filter_level
      )
      all_results[[platform]][[method]][[fl]] <- list(
        b_class_df_list = filtered$Bclass_df_list,
        k_class_df_list = filtered$Kclass_df_list
      )
    }
  }

  all_results
}

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
    paper_panel_theme() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none"
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
    paper_theme() +
    theme(legend.position = "none")
}


assemble_cross_tool_figure <- function(isoseq_upset, ont_upset, title) {

  combined <- (
    wrap_elements(full = isoseq_upset) /
    wrap_elements(full = ont_upset)
  ) +
    plot_layout(heights = c(1.2, 1)) +
    plot_annotation(
      title = title,
      theme = paper_figure_title_theme()
    )

  paper_inset_panel_tags(
    wrap_elements(full = combined),
    tags = c("a", "b"),
    heights = c(1.2, 1)
  )
}


assemble_cross_tool_bar_figure <- function(isoseq_bar, ont_bar, title) {

  legend_src <- isoseq_bar +
    ggplot2::theme(legend.position = "right") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        ncol = 1,
        title = "Structural Category"
      )
    )
  legend_grob <- cowplot::get_legend(legend_src)
  if (is.null(legend_grob)) {
    stop("Could not extract structural-category legend for cross-tool bar figure.", call. = FALSE)
  }

  isoseq_panel <- paper_tag_panel(
    isoseq_bar + ggplot2::theme(legend.position = "none"),
    "a"
  )
  ont_panel <- paper_tag_panel(
    ont_bar + ggplot2::theme(legend.position = "none"),
    "b"
  )

  (isoseq_panel | ont_panel | legend_grob) +
    patchwork::plot_layout(
      widths = c(1, 1, 0.32),
      guides = "keep"
    ) +
    patchwork::plot_annotation(
      title = title,
      theme = paper_figure_title_theme()
    )
}
