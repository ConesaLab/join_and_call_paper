# 01_config_ont_r10.R
# Paths for ONT R10 SY5Y human reports (does not modify mouse 01_config.R).

if (!exists("paper_repo_root", mode = "function")) {
  if (exists("r_source_dir", inherits = TRUE)) {
    source(file.path(r_source_dir, "00_figure_config.R"))
  }
}

if (!exists("ont_r10_repo_root", inherits = TRUE)) {
  if (exists("paper_repo_root", mode = "function")) {
    ont_r10_repo_root <- paper_repo_root()
  } else {
    ev <- Sys.getenv("JOIN_AND_CALL_REPO", unset = "")
    if (nzchar(ev)) {
      ont_r10_repo_root <- normalizePath(ev, mustWork = TRUE)
    } else if (exists("r_source_dir", inherits = TRUE)) {
      ont_r10_repo_root <- normalizePath(
        file.path(r_source_dir, "..", "..", ".."),
        mustWork = TRUE
      )
    } else {
      stop("Define `r_source_dir` before sourcing 01_config_ont_r10.R, or set JOIN_AND_CALL_REPO.")
    }
  }
}

ont_r10_report_dirs <- list(
  IsoQuant = file.path(ont_r10_repo_root, "reports", "ont_r10", "isoquant_report"),
  FLAIR    = file.path(ont_r10_repo_root, "reports", "ont_r10", "flair_report"),
  Bambu    = file.path(ont_r10_repo_root, "reports", "ont_r10", "bambu_report")
)

# Local data layout (e.g. after rsync from cluster) — read_level_ont_r10_sy5y_plots.Rmd
ont_r10_data_root <- file.path(ont_r10_repo_root, "data", "ont_r10")
ont_r10_read_qc_dir <- file.path(ont_r10_data_root, "read_qc")
ont_r10_sqanti_reads_root <- file.path(ont_r10_data_root, "sqanti_reads")

ont_r10_sy5y_default_sample_ids <- c(
  "SRR31732191",
  "SRR31732198",
  "SRR31732209",
  "SRR31732210"
)
ont_r10_sy5y_default_sample_labels <- c("191", "198", "209", "210")

ont_r10_sy5y_fl_filter_levels <- c(1, 3, 5, 10, 20)

# Condition suffix for figure titles (mouse: `...; brain` / `...; kidney`).
ont_r10_sy5y_condition_label <- "SY5Y"

#' Mouse-style assembled figure title: `main; condition` (see `all_plots.Rmd`).
sy5y_plot_title <- function(
    main,
    condition = ont_r10_sy5y_condition_label) {
  paste(main, condition, sep = "; ")
}

#' Per-panel theme for 2×2 data panels (fonts/grid; legend stripped — see corner grob).
sy5y_mouse_panel_theme <- function() {
  if (!exists("paper_panel_theme", mode = "function")) {
    source(file.path(
      if (exists("r_source_dir", inherits = TRUE)) r_source_dir else "R",
      "02_themes.R"
    ))
  }
  paper_panel_theme() +
    ggplot2::theme(legend.position = "none")
}

#' Patchwork figure title (alias of `paper_figure_title_theme()`).
sy5y_mouse_figure_annotation_theme <- function() {
  if (!exists("paper_figure_title_theme", mode = "function")) {
    source(file.path(
      if (exists("r_source_dir", inherits = TRUE)) r_source_dir else "R",
      "02_themes.R"
    ))
  }
  paper_figure_title_theme()
}

#' Centered figure title and/or x-axis label (alias of `paper_figure_annotation()`).
sy5y_figure_annotation <- function(title = NULL, x_label = NULL) {
  if (!exists("paper_figure_annotation", mode = "function")) {
    source(file.path(
      if (exists("r_source_dir", inherits = TRUE)) r_source_dir else "R",
      "02_themes.R"
    ))
  }
  if (length(x_label) == 1L && identical(x_label, "Sample")) {
    x_label <- NULL
  }
  paper_figure_annotation(title = title, x_label = x_label)
}

# Denominator for expression-bar "Unassigned" (see sy5y_read_numbers_for_expression).
ont_r10_sy5y_expression_total_col <- "ont_fastq"

#' PDF size classes for ONT R10 figures (see `FIG_SIZE` in `00_figure_config.R`):
#' `sy5y_2x2`, `sy5y_tpm_2x2` (12×9 in), `sy5y_legend_strip`, `sy5y_read_level_panel`, `sy5y_read_level_combo`.
ont_r10_sy5y_ylims <- list(
  structural      = NULL,
  expression      = NULL,
  comb            = NULL,
  comb_bar        = NULL,
  comb_fl_bar     = NULL,
  perc_comb_bar    = c(0, 65),
  perc_comb_fl_bar = c(0, 100),
  ujc_curve       = NULL,
  ujc_stack       = NULL,
  ujc_fl_stack    = NULL
)

check_paths_ont_r10 <- function() {
  #' Plot drivers load one of these (see [load_sy5y_bclass_list]).
  rdata_files <- c(
    "SY5Y_class_df_list.RData",
    "Bclass_df_list.RData",
    "all_class_data.RData"
  )
  for (nm in names(ont_r10_report_dirs)) {
    d <- ont_r10_report_dirs[[nm]]
    if (!dir.exists(d)) {
      cat(sprintf("[ONT R10] Missing directory [%s]: %s\n", nm, d))
      next
    }
    cat(sprintf("[ONT R10] OK [%s]: %s\n", nm, d))

    found_rdata <- NULL
    for (f in rdata_files) {
      fp <- file.path(d, f)
      if (file.exists(fp)) {
        cat(sprintf("  Found: %s\n", f))
        found_rdata <- f
        break
      }
    }
    if (is.null(found_rdata)) {
      cat(
        "  Missing plot RData (need SY5Y_class_df_list.RData from generate_report_rdata.R)\n",
        "  Note: Bclass_df_list.RData is only written for mouse condition B100K0, not SY5Y.\n",
        sep = ""
      )
    }

    d_sy5y <- file.path(d, "data", "classification_SY5Y")
    if (!dir.exists(d_sy5y)) {
      cat("  Missing raw data: data/classification_SY5Y/\n")
    } else {
      n_class <- length(list.files(
        d_sy5y,
        pattern = "_classification\\.txt$",
        full.names = FALSE
      ))
      cat(sprintf(
        "  Raw SQANTI classifications: %d file(s) under data/classification_SY5Y/ (expect 6)\n",
        n_class
      ))
      if (n_class < 6L) {
        cat("  Incomplete replicate set — re-run prepare_report.sh for this tool.\n")
      } else if (is.null(found_rdata)) {
        cat(
          "  Raw inputs look complete; run prepare_report.sh (R step) to create SY5Y_class_df_list.RData.\n"
        )
      }
    }

    fofn <- file.path(d, "data", "classification_SY5Y.fofn")
    if (file.exists(fofn)) {
      paths <- readLines(fofn, warn = FALSE)
      paths <- paths[nzchar(paths)]
      n_missing <- sum(!file.exists(paths))
      if (n_missing > 0L) {
        cat(sprintf(
          "  Warning: %d path(s) in classification_SY5Y.fofn do not exist here (cluster paths?).\n",
          n_missing
        ))
        cat("  Regenerate .fofn locally before running generate_report_rdata.R.\n")
      }
    }
  }

  rn <- file.path(ont_r10_read_qc_dir, "read_numbers_joint.tsv")
  if (file.exists(rn)) {
    cat(sprintf("[ONT R10 read-level] OK: %s\n", rn))
  } else {
    cat(sprintf(
      "[ONT R10 read-level] Missing: %s (needed for expression + read_level Rmds)\n",
      rn
    ))
  }

  invisible(ont_r10_report_dirs)
}

#' Check inputs for [read_level_ont_r10_sy5y_plots.Rmd].
check_paths_ont_r10_read_level <- function(
    sample_ids = ont_r10_sy5y_default_sample_ids,
    read_qc_dir = ont_r10_read_qc_dir,
    sqanti_reads_root = ont_r10_sqanti_reads_root) {
  rn <- file.path(read_qc_dir, "read_numbers_joint.tsv")
  ok <- file.exists(rn)
  cat(sprintf("[ONT R10 read-level] read_numbers_joint.tsv: %s\n  %s\n", if (ok) "OK" else "MISSING", rn))

  len_dir <- file.path(read_qc_dir, "lengths_seq_only", "ont")
  cat("[ONT R10 read-level] BAM read lengths dir: ", len_dir, "\n", sep = "")
  for (id in sample_ids) {
    lp <- file.path(len_dir, paste0(id, ".bam.readlen.txt"))
    cat(sprintf("  %s: %s\n", id, if (file.exists(lp)) "OK" else "MISSING"))
  }

  cat("[ONT R10 read-level] SQANTI-reads root: ", sqanti_reads_root, "\n", sep = "")
  for (id in sample_ids) {
    sp <- file.path(sqanti_reads_root, id, paste0(id, "_reads_classification.txt"))
    ok <- file.exists(sp)
    cat(sprintf("  %s: %s\n    %s\n", id, if (ok) "OK" else "MISSING", sp))
  }
  invisible(list(
    read_numbers_tsv = rn,
    lengths_ont_dir = len_dir,
    sqanti_reads_root = sqanti_reads_root
  ))
}
