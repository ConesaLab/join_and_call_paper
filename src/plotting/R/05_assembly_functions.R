# 05_assembly_functions.R
# Functions that assemble individual plots into multi-panel figures
# NOTE: all_results$ references from the original code have been fixed to use the results parameter

assemble_bar_plots <- function(results, tissue, fl_filter_level, title, subtitle,
                               pb_ylims = c(0, 325000), ont_ylims = c(0, 125000), 
                               plot_name_suffix = "_plots", use_scientific = FALSE) {
  
  fl_filter_level = as.character(fl_filter_level)
  
  plot_name = paste(tissue, plot_name_suffix, sep="")
  
  bar_1_1 <- bar_theme(
    results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_name]][[1]],
    "IsoQuant",
    pb_ylims
  )
    
  bar_1_2 <- bar_theme(
    results$IsoSeq$FLAIR[[fl_filter_level]][[plot_name]][[1]],
    "FLAIR",
    pb_ylims
  )
  
  bar_1_3 <- bar_theme(
    results$IsoSeq$Bambu[[fl_filter_level]][[plot_name]][[1]],
    "Bambu",
    pb_ylims
  )
  
  bar_1_4 <- bar_theme(
    results$IsoSeq$TALON[[fl_filter_level]][[plot_name]][[1]],
    "TALON",
    pb_ylims
  )
    
  bar_2_1 <- bar_theme(
    results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_name]][[1]],
    "Mandalorion",
    pb_ylims
  )
  
  bar_2_2 <- bar_theme(
    results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_name]][[1]],
    "IsoSeq + SQANTI3",
    pb_ylims
  )
  
  bar_3_1 <- bar_theme(
    results$ONT$IsoQuant[[fl_filter_level]][[plot_name]][[1]],
    "IsoQuant",
    ont_ylims
  ) +
    ylab(" # isoforms")
  
  bar_3_2 <- bar_theme(
    results$ONT$FLAIR[[fl_filter_level]][[plot_name]][[1]],
    "FLAIR",
    ont_ylims
  ) +
    ylab(" # isoforms")
  
  bar_3_3 <- bar_theme(
    results$ONT$Bambu[[fl_filter_level]][[plot_name]][[1]],
    "Bambu",
    ont_ylims
  ) +
    ylab(" # isoforms")
  
  bar_3_4 <- bar_theme(
    results$ONT$TALON[[fl_filter_level]][[plot_name]][[1]],
    "TALON",
    ont_ylims
  ) +
    ylab(" # isoforms")
  
  if (use_scientific) {
    sci_format <- scales::scientific_format()
    bar_1_1 <- bar_1_1 + scale_y_continuous(labels = sci_format)
    bar_1_2 <- bar_1_2 + scale_y_continuous(labels = sci_format)
    bar_1_3 <- bar_1_3 + scale_y_continuous(labels = sci_format)
    bar_1_4 <- bar_1_4 + scale_y_continuous(labels = sci_format)
    bar_2_1 <- bar_2_1 + scale_y_continuous(labels = sci_format)
    bar_2_2 <- bar_2_2 + scale_y_continuous(labels = sci_format)
    bar_3_1 <- bar_3_1 + scale_y_continuous(labels = sci_format)
    bar_3_2 <- bar_3_2 + scale_y_continuous(labels = sci_format)
    bar_3_3 <- bar_3_3 + scale_y_continuous(labels = sci_format)
    bar_3_4 <- bar_3_4 + scale_y_continuous(labels = sci_format)
  }
  
  bar_plot <- (
    bar_1_1 + bar_1_2 + bar_1_3 + bar_1_4 +
    bar_2_1 + bar_2_2 + plot_spacer() +
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
      heights = c(1, 1, 0.2, 1),
      axes = "collect_y"
    ) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = paper_figure_title_theme()
    ) &
    ggplot2::theme(legend.position = "none")
  
  paper_inset_tags_rows(
    wrap_elements(full = bar_plot),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.2, 1)
  )
}


assemble_expr_bar_plots <- function(results, tissue, fl_filter_level, title,
                               pb_ylims = c(0, 25000000), ont_ylims = c(0, 45000000), 
                               plot_name_suffix = "_plots", use_scientific = FALSE) {
  
  fl_filter_level = as.character(fl_filter_level)
  
  plot_name = paste(tissue, plot_name_suffix, sep="")
  
  bar_1_1 <- bar_theme(
    results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_name]][[1]],
    "IsoQuant",
    pb_ylims
  )
    
  bar_1_2 <- bar_theme(
    results$IsoSeq$FLAIR[[fl_filter_level]][[plot_name]][[1]],
    "FLAIR",
    pb_ylims
  )
  
  bar_1_3 <- bar_theme(
    results$IsoSeq$Bambu[[fl_filter_level]][[plot_name]][[1]],
    "Bambu",
    pb_ylims
  )
  
  bar_1_4 <- bar_theme(
    results$IsoSeq$TALON[[fl_filter_level]][[plot_name]][[1]],
    "TALON",
    pb_ylims
  )
    
  bar_2_1 <- bar_theme(
    results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_name]][[1]],
    "Mandalorion",
    pb_ylims
  )
  
  bar_2_2 <- bar_theme(
    results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_name]][[1]],
    "IsoSeq + SQANTI3",
    pb_ylims
  )
  
  bar_3_1 <- bar_theme(
    results$ONT$IsoQuant[[fl_filter_level]][[plot_name]][[1]],
    "IsoQuant",
    ont_ylims
  ) +
    ylab(" # reads")
  
  bar_3_2 <- bar_theme(
    results$ONT$FLAIR[[fl_filter_level]][[plot_name]][[1]],
    "FLAIR",
    ont_ylims
  ) +
    ylab(" # reads")
  
  bar_3_3 <- bar_theme(
    results$ONT$Bambu[[fl_filter_level]][[plot_name]][[1]],
    "Bambu",
    ont_ylims
  ) +
    ylab(" # reads")
  
  bar_3_4 <- bar_theme(
    results$ONT$TALON[[fl_filter_level]][[plot_name]][[1]],
    "TALON",
    ont_ylims
  ) +
    ylab(" # reads")
  
  if (use_scientific) {
    millions_format <- function(x) {
      paste0(x / 1e6, "M")
    }
    bar_1_1 <- bar_1_1 + scale_y_continuous(labels = millions_format)
    bar_1_2 <- bar_1_2 + scale_y_continuous(labels = millions_format)
    bar_1_3 <- bar_1_3 + scale_y_continuous(labels = millions_format)
    bar_1_4 <- bar_1_4 + scale_y_continuous(labels = millions_format)
    bar_2_1 <- bar_2_1 + scale_y_continuous(labels = millions_format)
    bar_2_2 <- bar_2_2 + scale_y_continuous(labels = millions_format)
    bar_3_1 <- bar_3_1 + scale_y_continuous(labels = millions_format)
    bar_3_2 <- bar_3_2 + scale_y_continuous(labels = millions_format)
    bar_3_3 <- bar_3_3 + scale_y_continuous(labels = millions_format)
    bar_3_4 <- bar_3_4 + scale_y_continuous(labels = millions_format)
  }
  
  bar_plot <- (
    bar_1_1 + bar_1_2 + bar_1_3 + bar_1_4 +
    bar_2_1 + bar_2_2 + plot_spacer() +
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
      heights = c(1, 1, 0.2, 1),
      axes = "collect_y"
    ) +
    plot_annotation(
      title = title,
      theme = paper_figure_title_theme()
    ) &
    ggplot2::theme(legend.position = "none")
  
  paper_inset_tags_rows(
    wrap_elements(full = bar_plot),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.2, 1)
  )
  
}


assemble_bar_sirv_plots <- function(results, tissue, fl_filter_level, title, subtitle, pb_ylims, ont_ylims, pb_talon_ylims, pb_isoseq_ylims, ont_talon_ylims) {
  
  fl_filter_level = as.character(fl_filter_level)
  
  plot_name = paste(tissue, "_plots", sep="")

  bar_1_1 <- bar_theme(
    results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_name]][[1]],
    "IsoQuant",
    pb_ylims
  ) + geom_hline(yintercept = 69, color = "red")
    
  bar_1_2 <- bar_theme(
    results$IsoSeq$FLAIR[[fl_filter_level]][[plot_name]][[1]],
    "FLAIR",
    pb_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  bar_1_3 <- bar_theme(
    results$IsoSeq$Bambu[[fl_filter_level]][[plot_name]][[1]],
    "Bambu",
    pb_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  bar_1_4 <- bar_theme(
    results$IsoSeq$TALON[[fl_filter_level]][[plot_name]][[1]],
    "TALON",
    pb_talon_ylims
  ) + geom_hline(yintercept = 69, color = "red")
    
  bar_2_1 <- bar_theme(
    results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_name]][[1]],
    "Mandalorion",
    pb_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  bar_2_2 <- bar_theme(
    results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_name]][[1]],
    "IsoSeq + SQANTI3",
    pb_isoseq_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  bar_3_1 <- bar_theme(
    results$ONT$IsoQuant[[fl_filter_level]][[plot_name]][[1]],
    "IsoQuant",
    ont_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  bar_3_2 <- bar_theme(
    results$ONT$FLAIR[[fl_filter_level]][[plot_name]][[1]],
    "FLAIR",
    ont_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  bar_3_3 <- bar_theme(
    results$ONT$Bambu[[fl_filter_level]][[plot_name]][[1]],
    "Bambu",
    ont_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  bar_3_4 <- bar_theme(
    results$ONT$TALON[[fl_filter_level]][[plot_name]][[1]],
    "TALON",
    ont_talon_ylims
  ) + geom_hline(yintercept = 69, color = "red")
  
  
  bar_plot <-   bar_1_1 + bar_1_2 + bar_1_3 + bar_1_4 +
                bar_2_1 + bar_2_2 + plot_spacer() +
                plot_spacer() +
                bar_3_1 + bar_3_2 + bar_3_3 + bar_3_4 +
                plot_layout(design = "
                ABCD
                EFGG
                HHHH
                IJKL
                ",
                heights = c(1, 1, 0.1, 1),
                axes = "collect_y") +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = PAPER_FONTS$caption, hjust = 0.5)
      )
    ) &
    ggplot2::theme(legend.position = "none")
  
  paper_inset_tags_rows(
    wrap_elements(full = bar_plot),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.1, 1)
  )
}


assemble_upset_plots <- function(results, tissue, fl_filter_level, title, plot_name_suffix = "_upset_plot") {
  
  fl_filter_level <- as.character(fl_filter_level)
  plot_name <- paste(tissue, plot_name_suffix, sep = "")
  
  upset_1_1 <- wrap_elements(full = results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_name]])
  upset_1_2 <- wrap_elements(full = results$IsoSeq$FLAIR[[fl_filter_level]][[plot_name]])
  upset_1_3 <- wrap_elements(full = results$IsoSeq$Bambu[[fl_filter_level]][[plot_name]])
  upset_1_4 <- wrap_elements(full = results$IsoSeq$TALON[[fl_filter_level]][[plot_name]])
  
  upset_2_1 <- wrap_elements(full = results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_name]])
  upset_2_2 <- wrap_elements(full = results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_name]])
  
  upset_3_1 <- wrap_elements(full = results$ONT$IsoQuant[[fl_filter_level]][[plot_name]])
  upset_3_2 <- wrap_elements(full = results$ONT$FLAIR[[fl_filter_level]][[plot_name]])
  upset_3_3 <- wrap_elements(full = results$ONT$Bambu[[fl_filter_level]][[plot_name]])
  upset_3_4 <- wrap_elements(full = results$ONT$TALON[[fl_filter_level]][[plot_name]])
  
  upset_plot <- (
    upset_1_1 + upset_1_2 + upset_1_3 + upset_1_4 +
    upset_2_1 + upset_2_2 + plot_spacer() +
    plot_spacer() +
    upset_3_1 + upset_3_2 + upset_3_3 + upset_3_4
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKM
      ",
      heights = c(1, 1, 0.1, 1)
    ) +
    plot_annotation(
      title = title,
      theme = paper_figure_title_theme()
    )
  
  paper_inset_tags_rows(
    wrap_elements(full = upset_plot),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.1, 1)
  )
}


assemble_comb_plots <- function(results, tissue, fl_filter_level, 
                                title, pb_ylims = c(0, 50), ont_ylims = c(0, 50), 
                                plot_name_suffix = "_comb_plot", y_label = NULL, x_label = NULL) {
  
  fl_filter_level <- as.character(fl_filter_level)
  plot_name <- paste(tissue, plot_name_suffix, sep = "")
  
  comb_1_1 <- results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_name]] +
    labs(title = "IsoQuant") +
    scale_y_continuous(limits = pb_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_1_2 <- results$IsoSeq$FLAIR[[fl_filter_level]][[plot_name]] +
    labs(title = "FLAIR") +
    scale_y_continuous(limits = pb_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_1_3 <- results$IsoSeq$Bambu[[fl_filter_level]][[plot_name]] +
    labs(title = "Bambu") +
    scale_y_continuous(limits = pb_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_1_4 <- results$IsoSeq$TALON[[fl_filter_level]][[plot_name]] +
    labs(title = "TALON") +
    scale_y_continuous(limits = pb_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_2_1 <- results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_name]] +
    labs(title = "Mandalorion") +
    scale_y_continuous(limits = pb_ylims)+
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_2_2 <- results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_name]] +
    labs(title = "IsoSeq + SQANTI3") +
    scale_y_continuous(limits = pb_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_3_1 <- results$ONT$IsoQuant[[fl_filter_level]][[plot_name]] +
    labs(title = "IsoQuant") +
    scale_y_continuous(limits = ont_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_3_2 <- results$ONT$FLAIR[[fl_filter_level]][[plot_name]] +
    labs(title = "FLAIR") +
    scale_y_continuous(limits = ont_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_3_3 <- results$ONT$Bambu[[fl_filter_level]][[plot_name]] +
    labs(title = "Bambu") +
    scale_y_continuous(limits = ont_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_3_4 <- results$ONT$TALON[[fl_filter_level]][[plot_name]] +
    labs(title = "TALON") +
    scale_y_continuous(limits = ont_ylims) +
    paper_panel_theme() +
    theme(legend.position = "none")
  
  comb_plot <- (
    comb_1_1 + comb_1_2 + comb_1_3 + comb_1_4 +
    comb_2_1 + comb_2_2 + plot_spacer() +
    plot_spacer() +
    comb_3_1 + comb_3_2 + comb_3_3 + comb_3_4
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
      theme = paper_figure_title_theme()
    ) &
    ggplot2::theme(legend.position = "none")
  
  if (!is.null(y_label)) {
    comb_plot <- comb_plot &
      ylab(y_label)
  }
  if(!is.null(x_label)) {
    comb_plot <- comb_plot &
      xlab(x_label)
  }
  
  paper_inset_tags_rows(
    wrap_elements(full = comb_plot),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.2, 1)
  )
}


assemble_comb_stack_plots <- function(results, tissue, fl_filter_level, 
                                title, pb_ylims = c(0, 50), ont_ylims = c(0, 50), 
                                plot_name_suffix = "_comb_plot", y_label = NULL, x_label = NULL) {
  
  fl_filter_level <- as.character(fl_filter_level)
  plot_name <- paste(tissue, plot_name_suffix, sep = "")
  
  comb_1_1 <- results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_name]] +
    labs(title = "IsoQuant") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = pb_ylims)
  
  comb_1_2 <- results$IsoSeq$FLAIR[[fl_filter_level]][[plot_name]] +
    labs(title = "FLAIR") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = pb_ylims)
  
  comb_1_3 <- results$IsoSeq$Bambu[[fl_filter_level]][[plot_name]] +
    labs(title = "Bambu") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = pb_ylims)
  
  comb_1_4 <- results$IsoSeq$TALON[[fl_filter_level]][[plot_name]] +
    labs(title = "TALON") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = pb_ylims)
  
  comb_2_1 <- results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_name]] +
    labs(title = "Mandalorion") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = pb_ylims)
  
  comb_2_2 <- results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_name]] +
    labs(title = "IsoSeq + SQANTI3") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = pb_ylims)
  
  comb_3_1 <- results$ONT$IsoQuant[[fl_filter_level]][[plot_name]] +
    labs(title = "IsoQuant") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = ont_ylims)
  
  comb_3_2 <- results$ONT$FLAIR[[fl_filter_level]][[plot_name]] +
    labs(title = "FLAIR") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = ont_ylims)
  
  comb_3_3 <- results$ONT$Bambu[[fl_filter_level]][[plot_name]] +
    labs(title = "Bambu") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = ont_ylims)
  
  comb_3_4 <- results$ONT$TALON[[fl_filter_level]][[plot_name]] +
    labs(title = "TALON") +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none") +
    scale_y_continuous(limits = ont_ylims)
  
  comb_plot <- (
    comb_1_1 + comb_1_2 + comb_1_3 + comb_1_4 +
    comb_2_1 + comb_2_2 + plot_spacer() +
    plot_spacer() +
    comb_3_1 + comb_3_2 + comb_3_3 + comb_3_4
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKL
      ",
      heights = c(1, 1, 0.2, 1),
      axes = "collect_y"
    ) +
    plot_annotation(
      title = title,
      theme = paper_figure_title_theme()
    ) &
    ggplot2::theme(legend.position = "none") 
  
  if (!is.null(y_label)) {
    comb_plot <- comb_plot &
      ylab(y_label)
  }
  if(!is.null(x_label)) {
    comb_plot <- comb_plot &
      xlab(x_label)
  }
  
  paper_inset_tags_rows(
    wrap_elements(full = comb_plot),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.2, 1)
  )
}


assemble_tama_st_upset_plots <- function(results, tissue, fl_filter_level, title) {

  fl_filter_level <- as.character(fl_filter_level)
  plot_name <- paste0(tissue, "_tama_st_upset")

  upset_1_1 <- wrap_elements(full = results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_name]])
  upset_1_2 <- wrap_elements(full = results$IsoSeq$FLAIR[[fl_filter_level]][[plot_name]])
  upset_1_3 <- wrap_elements(full = results$IsoSeq$Bambu[[fl_filter_level]][[plot_name]])
  upset_1_4 <- wrap_elements(full = results$IsoSeq$TALON[[fl_filter_level]][[plot_name]])

  upset_2_1 <- wrap_elements(full = results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_name]])
  upset_2_2 <- wrap_elements(full = results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_name]])

  upset_3_1 <- wrap_elements(full = results$ONT$IsoQuant[[fl_filter_level]][[plot_name]])
  upset_3_2 <- wrap_elements(full = results$ONT$FLAIR[[fl_filter_level]][[plot_name]])
  upset_3_3 <- wrap_elements(full = results$ONT$Bambu[[fl_filter_level]][[plot_name]])
  upset_3_4 <- wrap_elements(full = results$ONT$TALON[[fl_filter_level]][[plot_name]])

  upset_plot <- (
    upset_1_1 + upset_1_2 + upset_1_3 + upset_1_4 +
    upset_2_1 + upset_2_2 + plot_spacer() +
    plot_spacer() +
    upset_3_1 + upset_3_2 + upset_3_3 + upset_3_4
  ) +
    plot_layout(
      design = "
        ABCD
        EFGG
        HHHH
        IJKM
      ",
      heights = c(1, 1, 0.1, 1)
    ) +
    plot_annotation(
      title = title,
      theme = paper_figure_title_theme()
    )

  paper_inset_tags_rows(
    wrap_elements(full = upset_plot),
    tag_by_row = c("1" = "a", "4" = "b"),
    heights = c(1, 1, 0.1, 1)
  )
}


assemble_compare_plots <- function(results, tissue, fl_filter_level) {
  
  fl_filter_level = as.character(fl_filter_level)
  plot_base_name = paste0(tissue, "_compare_plots")
  
  subplots = c(
    "count_transcript",
    "total_isoforms_transcript",
    "count_ujc",
    "total_isoforms_ujc"
  )
  subplot_titles <- c(
    paste0("Redundancy of transcripts (# of transcripts); ", tissue),
    paste0("Redundancy of transcripts (# of isoforms); ", tissue),
    paste0("Redundancy of UJCs (# of UJCs); ", tissue),
    paste0("Redundancy of UJCs (# of isoforms); ", tissue)
  )
  names(subplot_titles) <- subplots

  final_plots = list()

  for (subplot_name in subplots) {
    subplot_title <- subplot_titles[[subplot_name]]
    
    bar_1_1_base <- results$IsoSeq$IsoQuant[[fl_filter_level]][[plot_base_name]][[subplot_name]]
    
    bar_1_1 <- compare_theme(
      bar_1_1_base,
      "IsoQuant"
    )
      
    bar_1_2 <- compare_theme(
      results$IsoSeq$FLAIR[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "FLAIR"
    )
    
    bar_1_3 <- compare_theme(
      results$IsoSeq$Bambu[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "Bambu"
    )
    
    bar_1_4 <- compare_theme(
      results$IsoSeq$TALON[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "TALON"
    )
      
    bar_2_1 <- compare_theme(
      results$IsoSeq$Mandalorion[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "Mandalorion"
    )
    
    bar_2_2 <- compare_theme(
      results$IsoSeq$isoseq_sqanti[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "IsoSeq + SQANTI3"
    )
    
    bar_3_1 <- compare_theme(
      results$ONT$IsoQuant[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "IsoQuant"
    )
    
    bar_3_2 <- compare_theme(
      results$ONT$FLAIR[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "FLAIR"
    )
    
    bar_3_3 <- compare_theme(
      results$ONT$Bambu[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "Bambu"
    )
    
    bar_3_4 <- compare_theme(
      results$ONT$TALON[[fl_filter_level]][[plot_base_name]][[subplot_name]],
      "TALON"
    )
    
    bar_plot <- (
      bar_1_1 + bar_1_2 + bar_1_3 + bar_1_4 +
      bar_2_1 + bar_2_2 + plot_spacer() +
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
        heights = c(1, 1, 0.2, 1),
        axes = "collect_y",
        axis_titles = "collect"
      ) +
      plot_annotation(
        title = paste0(subplot_title, "; ", tissue),
        theme = paper_figure_title_theme()
      ) &
      ggplot2::theme(legend.position = "none")
    
    final_plots[[subplot_name]] <- paper_inset_tags_rows(
      wrap_elements(full = bar_plot),
      tag_by_row = c("1" = "a", "4" = "b"),
      heights = c(1, 1, 0.2, 1)
    )
  }
  
  return(final_plots)
}
