# 07_ont_r10_sy5y.R
# ONT R10 SY5Y (single condition) plotting helpers — mouse scripts unchanged.
#
# FL repair: SQANTI3 *_classification.txt for concat/TAMA often has NA in `FL` but
# per-sample columns `FL.<suffix>` (e.g. FL.SRR31732191, FL.SRR31732191_primary_aln_sorted).
# `ensure_fl_sy5y()` sums every column whose name starts with `FL` except the main `FL` column.
# FOFN sort order (prepare_report.sh): 1_concat, 2_TAMA, 3_SRR* lexicographic → S1..S6 in generate_report_rdata.R.

SY5Y_TOOL_ORDER <- c("IsoQuant", "FLAIR", "Bambu")

SY5Y_UJC_PANEL_NAMES <- c(
  "upset", "comb", "comb_bar", "comb_fl_bar",
  "ujc_curve", "ujc_stack", "ujc_fl_stack"
)

SY5Y_COMPARE_PANEL_NAMES <- c(
  "count_transcript", "total_isoforms_transcript",
  "count_ujc", "total_isoforms_ujc"
)

# Mouse `all_plots.Rmd`: `brain_comb_bar.pdf` (# UJCs) and `brain_comb_fl_bar.pdf` (read support).
SY5Y_UJC_OCCURRENCE_TRANSCRIPT <- "comb_bar"
SY5Y_UJC_OCCURRENCE_READS <- "comb_fl_bar"

#' Output directory for ONT R10 SY5Y figures.
sy5y_plot_creation_dir <- function() {
  file.path(ont_r10_repo_root, "plot_creation_ont_r10")
}

#' Short x-axis labels (FOFN order: concat, TAMA, four replicates).
sy5y_sample_labels <- function() {
  c("J&C", "C&J", "191", "198", "209", "210")
}

#' Replicate-only axis labels (positions 3–6; excludes J&C / C&J strategies).
sy5y_replicate_sample_labels <- function() {
  ont_r10_sy5y_default_sample_labels
}

#' Internal `sample` factor levels from a classification list.
sy5y_sample_levels <- function(class_df_list) {
  lev <- names(class_df_list)
  if (!length(lev)) {
    lev <- paste0("S", seq_along(class_df_list))
  }
  lev
}

#' @description Sum per-sample FL columns when main `FL` is mostly NA (concat / TAMA).
ensure_fl_sy5y <- function(class_df_list) {
  lapply(class_df_list, function(df) {
    fl_aux <- grep("^FL", names(df), value = TRUE, ignore.case = TRUE)
    fl_aux <- setdiff(fl_aux, "FL")
    if (!"FL" %in% names(df)) {
      df$FL <- NA_real_
    }
    fl_ok <- !all(is.na(df$FL)) && mean(is.na(df$FL)) <= 0.95
    if (fl_ok) {
      return(df)
    }
    if (length(fl_aux)) {
      mat <- as.matrix(df[, fl_aux, drop = FALSE])
      mode(mat) <- "numeric"
      df$FL <- rowSums(replace(mat, is.na(mat), 0), na.rm = TRUE)
    }
    df
  })
}

load_sy5y_bclass_list <- function(report_dir) {
  primary <- file.path(report_dir, "Bclass_df_list.RData")
  if (file.exists(primary)) {
    e <- new.env(parent = emptyenv())
    load(primary, envir = e)
    return(e$Bclass_df_list)
  }
  sy5y_primary <- file.path(report_dir, "SY5Y_class_df_list.RData")
  if (file.exists(sy5y_primary)) {
    e <- new.env(parent = emptyenv())
    load(sy5y_primary, envir = e)
    if (!exists("SY5Y_class_df_list", envir = e)) {
      stop(
        "Found SY5Y_class_df_list.RData but object SY5Y_class_df_list is missing: ",
        sy5y_primary
      )
    }
    return(e$SY5Y_class_df_list)
  }
  fallback <- file.path(report_dir, "all_class_data.RData")
  if (file.exists(fallback)) {
    e <- new.env(parent = emptyenv())
    load(fallback, envir = e)
    if (!exists("all_class_df_lists", envir = e) || !"SY5Y" %in% names(e$all_class_df_lists)) {
      stop("Found all_class_data.RData but no SY5Y entry — check prepare_report output.")
    }
    warning(
      "Using SY5Y from all_class_data.RData; prefer SY5Y_class_df_list.RData or Bclass_df_list.RData in the report directory.",
      call. = FALSE
    )
    return(e$all_class_df_lists[["SY5Y"]])
  }
  stop(
    "No classification RData in: ", report_dir,
    "\nExpected one of: Bclass_df_list.RData, SY5Y_class_df_list.RData, all_class_data.RData",
    "\n(Re-run prepare_report / generate_report_rdata for this tool.)"
  )
}

filter_fl_threshold <- function(class_df_list, thr) {
  nm <- names(class_df_list)
  out <- lapply(seq_along(class_df_list), function(i) {
    df <- class_df_list[[i]]
    if ("FL" %in% names(df) && mean(is.na(df$FL)) > 0.95) {
      stop(
        "FL column still >95% NA for list element ", i,
        " after ensure_fl_sy5y() — inspect classification column names."
      )
    }
    ok <- !is.na(df$FL) & df$FL >= thr
    df[ok, , drop = FALSE]
  })
  if (length(nm)) names(out) <- nm
  out
}

sy5y_combined_classification_df <- function(class_df_list, fl_threshold = 1) {
  df_list <- filter_fl_threshold(ensure_fl_sy5y(class_df_list), fl_threshold)
  lev <- sy5y_sample_levels(df_list)
  combined <- do.call(
    rbind,
    lapply(df_list, function(class_df) class_df[, c("sample", "structural_category", "chrom", "FL")])
  )
  combined$sample <- factor(combined$sample, levels = lev)
  combined
}

#' Structural category counts only (same pipeline as mouse `all_plots.Rmd` count panel).
sy5y_classification_count_panel <- function(class_df_list, sample_labels, fl_threshold = 1) {
  combined <- sy5y_combined_classification_df(class_df_list, fl_threshold)
  lev <- levels(combined$sample)
  lbl <- sample_labels
  if (length(lev) != length(lbl)) {
    warning(
      "SY5Y structural plot: ", length(lev), " classification table(s) vs ", length(lbl),
      " axis labels (expected 6). Check this tool's classification FOFN and SQANTI FL columns.",
      call. = FALSE
    )
    lbl <- lbl[seq_along(lev)]
  }
  ps <- plot_classification_data(
    combined,
    lbl,
    x_discrete_limits     = lev,
    x_drop_missing_levels = FALSE
  )
  list(count = ps$count)
}

sy5y_classification_two_panel <- sy5y_classification_count_panel

style_sy5y_count_y <- function(p) {
  p + ggplot2::scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    expand = ggplot2::expansion(mult = c(0, 0.03))
  )
}

style_sy5y_expr_y <- function(p) {
  p + ggplot2::scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    expand = ggplot2::expansion(mult = c(0, 0.03))
  )
}

style_sy5y_axes <- function(p) {
  p + ggplot2::theme(
    axis.text = ggplot2::element_text(size = 7.5),
    axis.text.x = ggplot2::element_text(
      size = 7,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = ggplot2::element_text(size = 8)
  )
}

strip_legend <- function(p) {
  p + ggplot2::theme(legend.position = "none")
}

#' Per-panel tool title — same as mouse `assemble_comb_plots` / `bar_theme` (`labs(title = tool)`).
sy5y_panel_tool_title <- function(plot, tool_name) {
  plot +
    ggplot2::labs(title = tool_name) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12)
    )
}

sy5y_structural_panel <- function(panels, tool_name) {
  ylims <- ont_r10_sy5y_ylims$structural
  p <- style_sy5y_axes(
    style_sy5y_count_y(strip_legend(panels$count)) +
      ggplot2::labs(y = "# isoforms")
  )
  if (!is.null(ylims) && length(ylims) == 2L && !any(is.na(ylims))) {
    bar_theme(p, tool_name, ylims)
  } else {
    sy5y_panel_tool_title(p, tool_name)
  }
}

sy5y_comb_panel_theme <- function(plot, tool_name, ylims = NULL, ylab = NULL) {
  p <- strip_legend(plot)
  p <- sy5y_panel_tool_title(p, tool_name)
  p <- sy5y_apply_ylims(p, ylims)
  if (!is.null(ylab) && nzchar(ylab)) {
    p <- p + ggplot2::labs(y = ylab)
  }
  style_sy5y_axes(p)
}

sy5y_legend_from_plot <- function(legend_plot) {
  cowplot::get_legend(
    legend_plot +
      ggplot2::theme(
        legend.position = "right",
        legend.key.size = grid::unit(0.4, "cm"),
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 9)
      )
  )
}

sy5y_legend_column <- function(legend_grob) {
  patchwork::wrap_elements(full = legend_grob, clip = FALSE)
}

#' Drop duplicate y-axis titles (keep ylab on top-left panel only; mouse uses `collect_y`).
sy5y_strip_ylab <- function(p) {
  p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
}

#' Core 2×2 layout: IsoQuant | FLAIR / Bambu | legend.
assemble_ont_r10_three_tool_grid <- function(
    iso_col,
    flair_col,
    bambu_col,
    legend_grob,
    title = NULL) {
  leg_col <- sy5y_legend_column(legend_grob)
  core <- (iso_col | sy5y_strip_ylab(flair_col)) / (sy5y_strip_ylab(bambu_col) | leg_col) +
    patchwork::plot_layout(
      widths = c(1, 1),
      heights = c(1, 1),
      guides = "keep",
      axes = "collect_y"
    )
  if (length(title) && nzchar(title)) {
    core <- core +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13)
        )
      )
  }
  core
}

#' 2×2 structural category count grid.
assemble_ont_r10_tool_grid <- function(iso_panels, flair_panels, bambu_panels, title = NULL) {
  leg <- sy5y_legend_from_plot(iso_panels$count)
  assemble_ont_r10_three_tool_grid(
    sy5y_structural_panel(iso_panels, "IsoQuant"),
    sy5y_structural_panel(flair_panels, "FLAIR"),
    sy5y_structural_panel(bambu_panels, "Bambu"),
    leg,
    title = title
  )
}

sy5y_apply_ylims <- function(p, ylims) {
  if (is.null(ylims) || length(ylims) != 2L || any(is.na(ylims))) {
    return(p)
  }
  p + ggplot2::coord_cartesian(ylim = ylims)
}

#' Save one multi-tool 2×2 UJC panel (used for mouse-style occurrence / curve / stack figures).
sy5y_ggsave_ujc_multitool_panel <- function(
    all_results,
    panel_name,
    title,
    out_dir,
    fl_threshold = 1L,
    ylab = NULL,
    ylims = NULL,
    width = 14,
    height = 11) {
  if (is.null(ylims) && exists("ont_r10_sy5y_ylims", inherits = TRUE)) {
    ylims <- ont_r10_sy5y_ylims[[panel_name]]
    if (is.null(ylims)) {
      ylims <- ont_r10_sy5y_ylims$comb_bar
    }
  }
  p <- assemble_ont_r10_ujc_ggplot_grid(
    all_results,
    panel_name = panel_name,
    title = title,
    ylims = ylims,
    ylab = ylab
  )
  out_pdf <- file.path(
    out_dir,
    sprintf("sy5y_%s_fl%s.pdf", panel_name, fl_threshold)
  )
  ggplot2::ggsave(
    out_pdf,
    plot = p,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    units = "in"
  )
  invisible(out_pdf)
}

#' 2×2 grid for one UJC ggplot panel type across three tools.
assemble_ont_r10_ujc_ggplot_grid <- function(
    all_results,
    panel_name,
    title = NULL,
    ylims = NULL,
    ylab = NULL,
    legend_panel = "comb_bar") {
  iso_p <- all_results$IsoQuant$ujc[[panel_name]]
  fl_p <- all_results$FLAIR$ujc[[panel_name]]
  bm_p <- all_results$Bambu$ujc[[panel_name]]
  leg_plot <- all_results$IsoQuant$ujc[[legend_panel]]
  if (is.null(leg_plot)) {
    leg_plot <- iso_p
  }
  leg <- sy5y_legend_from_plot(leg_plot)
  assemble_ont_r10_three_tool_grid(
    sy5y_comb_panel_theme(iso_p, "IsoQuant", ylims, ylab),
    sy5y_comb_panel_theme(fl_p, "FLAIR", ylims, ylab),
    sy5y_comb_panel_theme(bm_p, "Bambu", ylims, ylab),
    leg,
    title = title
  )
}

#' 2×2 grid for UpSet panels (ComplexUpset grobs; title from `create_upset_plot()` only).
assemble_ont_r10_upset_grid <- function(all_results, title = NULL) {
  upset_col <- function(tool_name) {
    patchwork::wrap_elements(
      full = all_results[[tool_name]]$ujc$upset,
      clip = FALSE
    )
  }
  leg_plot <- all_results$IsoQuant$ujc$comb_bar
  leg <- sy5y_legend_from_plot(leg_plot)
  assemble_ont_r10_three_tool_grid(
    upset_col("IsoQuant"),
    sy5y_strip_ylab(upset_col("FLAIR")),
    sy5y_strip_ylab(upset_col("Bambu")),
    leg,
    title = title
  )
}

#' 2×2 grid for one compare panel across tools.
assemble_ont_r10_compare_grid <- function(
    all_results,
    panel_name,
    title = NULL) {
  iso_base <- all_results$IsoQuant$compare[[panel_name]]
  leg <- cowplot::get_legend(
    iso_base +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3, nrow = 2)) +
      ggplot2::theme(
        legend.position = "right",
        legend.key.size = grid::unit(0.35, "cm"),
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 9)
      )
  )
  cmp_col <- function(tool_name) {
    style_sy5y_axes(
      compare_theme(all_results[[tool_name]]$compare[[panel_name]], tool_name)
    )
  }
  assemble_ont_r10_three_tool_grid(
    cmp_col("IsoQuant"),
    cmp_col("FLAIR"),
    cmp_col("Bambu"),
    leg,
    title = title
  )
}

sy5y_jc_compare_inputs <- function(class_df_list, fl_threshold = 1) {
  df_list <- filter_fl_threshold(ensure_fl_sy5y(class_df_list), fl_threshold)
  if (length(df_list) < 2L) {
    stop("Join&Call vs Call&Join compare needs at least two list elements.")
  }
  nm <- names(df_list)
  if (!is.null(nm) && all(c("Join&Call", "Call&Join") %in% nm)) {
    return(list(
      `Join&Call` = df_list[["Join&Call"]],
      `Call&Join` = df_list[["Call&Join"]]
    ))
  }
  list(
    `Join&Call` = df_list[[1L]],
    `Call&Join` = df_list[[2L]]
  )
}

sy5y_compare_isoform_bundle <- function(class_df_list, fl_threshold = 1,
                                        condition_label = "SY5Y") {
  compare_isoform_plots(
    sy5y_jc_compare_inputs(class_df_list, fl_threshold),
    include_percentage_plots = FALSE,
    condition_label = condition_label
  )
}

sy5y_build_single_tool_results <- function(
    report_dir,
    tool_name,
    fl_threshold = 1,
    include_expression = FALSE) {
  raw <- load_sy5y_bclass_list(report_dir)
  filtered <- filter_fl_threshold(ensure_fl_sy5y(raw), fl_threshold)
  out <- list(
    raw        = raw,
    filtered   = filtered,
    structural = sy5y_classification_count_panel(raw, sy5y_sample_labels(), fl_threshold),
    ujc        = sy5y_ujc_bundle(filtered, tool_name),
    compare    = sy5y_compare_isoform_bundle(raw, fl_threshold, condition_label = "SY5Y")
  )
  if (isTRUE(include_expression)) {
    out$expression <- sy5y_expression_panel(raw, fl_threshold = fl_threshold)
  }
  out
}

#' Load IsoQuant, FLAIR, Bambu SY5Y results for one FL threshold.
sy5y_build_all_tools_results <- function(fl_threshold = 1, include_expression = FALSE) {
  out <- list()
  for (tool in SY5Y_TOOL_ORDER) {
    out[[tool]] <- sy5y_build_single_tool_results(
      ont_r10_report_dirs[[tool]],
      tool,
      fl_threshold,
      include_expression = include_expression
    )
  }
  out
}

#' Read totals for expression bars (J&C/C&J = sum of replicate FASTQ; replicates from TSV).
sy5y_read_numbers_for_expression <- function(
    class_df_list,
    read_numbers_path = file.path(ont_r10_read_qc_dir, "read_numbers_joint.tsv"),
    sample_ids = ont_r10_sy5y_default_sample_ids,
    read_col = "ont_fastq") {
  if (!file.exists(read_numbers_path)) {
    stop("Missing read numbers file for expression plots: ", read_numbers_path)
  }
  rn <- utils::read.table(
    read_numbers_path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!read_col %in% names(rn)) {
    stop("Column ", read_col, " not found in ", read_numbers_path)
  }
  rep_vals <- stats::setNames(as.numeric(rn[[read_col]]), rn$sample)
  lev <- sy5y_sample_levels(class_df_list)
  if (length(lev) < 2L + length(sample_ids)) {
    stop(
      "Expression read mapping needs ", 2L + length(sample_ids),
      " classification tables; got ", length(lev), "."
    )
  }
  jc_total <- sum(rep_vals[sample_ids], na.rm = TRUE)
  totals <- stats::setNames(rep(NA_real_, length(lev)), lev)
  totals[lev[1L]] <- jc_total
  totals[lev[2L]] <- jc_total
  for (i in seq_along(sample_ids)) {
    sid <- sample_ids[[i]]
    if (!sid %in% names(rep_vals)) {
      stop("Sample ", sid, " not found in ", read_numbers_path)
    }
    totals[lev[2L + i]] <- rep_vals[[sid]]
  }
  totals
}

sy5y_expression_panel <- function(
    class_df_list,
    fl_threshold = 1,
    read_numbers_path = file.path(ont_r10_read_qc_dir, "read_numbers_joint.tsv")) {
  combined <- sy5y_combined_classification_df(class_df_list, fl_threshold)
  lev <- levels(combined$sample)
  lbl <- sy5y_sample_labels()
  if (length(lev) != length(lbl)) {
    lbl <- lbl[seq_along(lev)]
  }
  totals <- sy5y_read_numbers_for_expression(class_df_list, read_numbers_path)
  ps <- plot_classification_expression_data(
    combined,
    totals,
    sample_codes  = lev,
    sample_labels = lbl
  )
  list(count = ps[[1L]])
}

assemble_ont_r10_expression_grid <- function(all_results, title = NULL) {
  ylims <- ont_r10_sy5y_ylims$expression
  expr_col <- function(tool_name) {
    p <- all_results[[tool_name]]$expression$count
    p <- strip_legend(style_sy5y_axes(style_sy5y_expr_y(p)))
    p <- sy5y_apply_ylims(p, ylims)
    p <- sy5y_panel_tool_title(p, tool_name)
    p + ggplot2::labs(y = "# reads")
  }
  leg <- sy5y_legend_from_plot(all_results$IsoQuant$expression$count)
  assemble_ont_r10_three_tool_grid(
    expr_col("IsoQuant"),
    expr_col("FLAIR"),
    expr_col("Bambu"),
    leg,
    title = title
  )
}

#' Columns for UJC combination / curve / stack plots (four biological replicates only).
#' UpSet still uses all six tracks via [sy5y_sample_labels()]; J&C and C&J are excluded here
#' so "number of samples present" reflects replicates, matching mouse brain/kidney (n = 5).
sy5y_ujc_combination_columns <- function() {
  sy5y_replicate_sample_labels()
}

sy5y_ujc_bundle <- function(class_df_list, method_title, upset_n_intersections = 12L) {
  create_upset_plot(
    class_df_list,
    sy5y_sample_labels(),
    method_title,
    n = upset_n_intersections,
    combination_columns = sy5y_ujc_combination_columns(),
    include_percentage_plots = FALSE
  )
}

sy5y_save_compare_pdfs <- function(bundle, tool_slug, out_dir,
                                 width = 12,
                                 height = 8) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  for (nm in names(bundle)) {
    p <- bundle[[nm]]
    if (is.null(p)) next
    ggplot2::ggsave(
      file.path(out_dir, sprintf("sy5y_%s_compare_%s.pdf", tool_slug, nm)),
      plot = p,
      device = grDevices::cairo_pdf,
      width = width,
      height = height,
      units = "in"
    )
  }
  invisible(out_dir)
}

sy5y_compare_patchwork <- function(bundle) {
  (bundle$count_transcript | bundle$total_isoforms_transcript) /
    (bundle$count_ujc | bundle$total_isoforms_ujc)
}

sy5y_save_ujc_pdfs <- function(bundle, tool_slug, out_dir,
                             width = 12,
                             height = 8,
                             upset_width = 16,
                             upset_height = 9) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  nm_dims <- list(
    upset         = c(upset_width, upset_height),
    comb          = c(width, height),
    comb_bar      = c(width, height),
    comb_fl_bar   = c(width, height),
    ujc_curve     = c(width, height * 0.85),
    ujc_stack     = c(width, height),
    ujc_fl_stack  = c(width, height)
  )
  for (nm in names(nm_dims)) {
    p <- bundle[[nm]]
    if (is.null(p)) next
    wd <- nm_dims[[nm]][1]
    ht <- nm_dims[[nm]][2]
    out_pdf <- file.path(out_dir, sprintf("sy5y_%s_ujc_%s.pdf", tool_slug, nm))
    ggplot2::ggsave(
      out_pdf,
      plot = p,
      device = grDevices::cairo_pdf,
      width = wd,
      height = ht,
      units = "in"
    )
  }
  invisible(out_dir)
}

sy5y_ujc_overview_patchwork <- function(bundle, tool_title) {
  top <- patchwork::wrap_elements(full = bundle$upset)
  p_comb <- bundle$comb + ggplot2::labs(subtitle = "UJCs per intersection size")
  p_bar <- bundle$comb_bar + ggplot2::labs(subtitle = "By structural category")
  p_fl <- bundle$comb_fl_bar + ggplot2::labs(subtitle = "Read support (FL)")
  p_cv <- bundle$ujc_curve + ggplot2::labs(subtitle = "Discovery curve")
  p_st <- bundle$ujc_stack + ggplot2::labs(subtitle = "UJC counts by reproducibility")
  p_fst <- bundle$ujc_fl_stack + ggplot2::labs(subtitle = "Read support by reproducibility")
  mid <- (p_comb | p_bar) / (p_fl | p_cv)
  bot <- (p_st | p_fst)
  (top / mid / bot) +
    patchwork::plot_layout(heights = c(1.15, 1, 0.95)) +
    patchwork::plot_annotation(
      title = tool_title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14))
    )
}

#' TPM curve data for J&C vs C&J (same logic as mouse build_tpm_curve_data).
build_sy5y_tpm_curve_data <- function(class_df_list, fl_threshold = 1) {
  jc_inputs <- sy5y_jc_compare_inputs(class_df_list, fl_threshold)
  if (!exists("simplify_category", mode = "function")) {
    stop("Source 09_filter_level_lines.R before build_sy5y_tpm_curve_data().")
  }
  process_strategy <- function(df, strategy_label) {
    df <- df %>% dplyr::filter(!is.na(.data$FL) & .data$FL > 0)
    total_fl <- sum(df$FL)
    df %>%
      dplyr::mutate(
        TPM = (.data$FL / total_fl) * 1e6,
        strategy = strategy_label,
        full_category = factor(
          as.character(.data$structural_category),
          levels = names(cat.palette)
        ),
        simplified_category = factor(
          simplify_category(as.character(.data$structural_category)),
          levels = names(simplified_cat_palette)
        )
      ) %>%
      dplyr::select(
        dplyr::all_of(c("strategy", "FL", "TPM", "full_category", "simplified_category"))
      )
  }
  dplyr::bind_rows(
    process_strategy(jc_inputs[["Join&Call"]], "J&C"),
    process_strategy(jc_inputs[["Call&Join"]], "C&J")
  )
}

#' 2×2 TPM line-plot figure (log or linear x).
assemble_ont_r10_tpm_line_figure <- function(
    line_plots,
    title,
    log_x = FALSE) {
  no_x <- ggplot2::theme(axis.title.x = ggplot2::element_blank())
  x_grob <- patchwork::wrap_elements(full = grid::textGrob(
    tpm_x_axis_label(log_x),
    gp = grid::gpar(fontsize = 11)
  ))
  line_col <- function(p) {
    p + no_x
  }
  shared_legend <- patchwork::wrap_elements(
    full = build_tpm_line_legend_figure(category_layout = c(2, 3)),
    clip = FALSE
  )
  leg_col <- shared_legend
  core <- (
    line_col(line_plots$IsoQuant) |
      sy5y_strip_ylab(line_col(line_plots$FLAIR))
  ) / (
    sy5y_strip_ylab(line_col(line_plots$Bambu)) |
      leg_col
  ) / x_grob +
    patchwork::plot_layout(
      widths = c(1, 1),
      heights = c(1, 1, 0.04),
      axes = "collect_y"
    ) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  core
}

#' 2×2 TPM quantile bar figure (count or percent).
assemble_ont_r10_tpm_bar_figure <- function(bar_plots, title) {
  no_x <- ggplot2::theme(axis.title.x = ggplot2::element_blank())
  bar_col <- function(p) {
    p + no_x
  }
  fill_legend_plot <- bar_plots$IsoQuant +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 3, ncol = 3,
        title = "Structural Category",
        title.position = "top"
      )
    )
  fill_legend <- cowplot::get_legend(fill_legend_plot)
  strategy_df <- data.frame(
    x = factor(c("J&C", "C&J"), levels = c("J&C", "C&J")),
    y = c(1, 1)
  )
  strategy_legend_plot <- ggplot2::ggplot(strategy_df, ggplot2::aes(x = x, y = y, fill = x)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(
      values = c("J&C" = "grey50", "C&J" = "grey80"),
      name = "Strategy"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, title.position = "top"))
  strategy_legend <- cowplot::get_legend(strategy_legend_plot)
  combined_legend <- (
    patchwork::wrap_elements(full = strategy_legend) |
      patchwork::wrap_elements(full = fill_legend)
  ) +
    patchwork::plot_layout(widths = c(1, 2))
  leg_col <- patchwork::wrap_elements(full = combined_legend)
  x_grob <- patchwork::wrap_elements(full = grid::textGrob(
    "Strategy", gp = grid::gpar(fontsize = 11)
  ))
  (
    bar_col(bar_plots$IsoQuant) |
      sy5y_strip_ylab(bar_col(bar_plots$FLAIR))
  ) / (
    sy5y_strip_ylab(bar_col(bar_plots$Bambu)) |
      leg_col
  ) / x_grob +
    patchwork::plot_layout(
      widths = c(1, 1),
      heights = c(1, 1, 0.04),
      axes = "collect_y"
    ) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
}

#' Print max y values per tool/panel to tune ont_r10_sy5y_ylims.
sy5y_panel_ymax <- function(p) {
  tryCatch(
    max(ggplot2::ggplot_build(p)$data[[1]]$y, na.rm = TRUE),
    error = function(e) NA_real_
  )
}

sy5y_print_ylim_diagnostics <- function(all_results) {
  cat("=== SY5Y y-limit diagnostics (FL filter used in all_results) ===\n")
  for (tool in SY5Y_TOOL_ORDER) {
    tr <- all_results[[tool]]
    n_iso <- sy5y_panel_ymax(tr$structural$count)
    cat(sprintf(
      "%s structural max # isoforms: %s\n",
      tool,
      if (is.na(n_iso)) "n/a" else format(n_iso, big.mark = ",")
    ))
    for (pn in c("comb_bar", "comb_fl_bar", "ujc_curve")) {
      ymax <- sy5y_panel_ymax(tr$ujc[[pn]])
      cat(sprintf(
        "  %s max y: %s\n",
        pn,
        if (is.na(ymax)) "n/a" else format(ymax, big.mark = ",")
      ))
    }
  }
  invisible(NULL)
}
