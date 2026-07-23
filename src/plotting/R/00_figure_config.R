# 00_figure_config.R
# Repository paths, canonical figure dimensions, typography constants, paper_ggsave().

#' TRUE when R runs on native Windows (e.g. RStudio on Windows, not WSL R).
paper_is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

#' Map WSL paths (`/mnt/c/...`) to Windows drive paths (`C:/...`) when needed.
paper_wsl_to_native_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    return(path)
  }
  if (!paper_is_windows()) {
    return(path)
  }
  m <- regmatches(path, regexec("^/mnt/([a-zA-Z])/(.*)$", path))[[1]]
  if (length(m) < 3L) {
    return(path)
  }
  sprintf("%s:/%s", toupper(m[2]), m[3])
}

#' Normalize a path for the current OS; optional `mustWork` check.
paper_normalize_path <- function(path, must_work = FALSE) {
  path <- paper_wsl_to_native_path(path[1L])
  if (isTRUE(must_work)) {
    return(normalizePath(path, mustWork = TRUE, winslash = "/"))
  }
  if (dir.exists(path) || file.exists(path)) {
    return(normalizePath(path, mustWork = FALSE, winslash = "/"))
  }
  path
}

#' NIH analysis data root (`C:/data/nih` on Windows, `/mnt/c/data/nih` in WSL).
#' Override with env var `JOIN_AND_CALL_NIH_DATA`.
paper_nih_data_root <- function() {
  ev <- Sys.getenv("JOIN_AND_CALL_NIH_DATA", unset = "")
  if (nzchar(ev)) {
    return(paper_normalize_path(ev, must_work = FALSE))
  }
  default <- if (paper_is_windows()) "C:/data/nih" else "/mnt/c/data/nih"
  paper_normalize_path(default, must_work = FALSE)
}

#' SQANTI-reads outputs for mouse read-level figures (`.../sqanti_reads`).
paper_sqanti_reads_root <- function() {
  file.path(paper_nih_data_root(), "sqanti_reads")
}

#' Read-length sidecar files (`.../sqanti_reads/lengths_seq_only`).
paper_sqanti_lengths_root <- function() {
  file.path(paper_sqanti_reads_root(), "lengths_seq_only")
}

#' Mouse reference + SIRV GTF (`mm39.ncbiRefSeq_SIRV.gtf` under NIH data root).
#' Override with env var `JOIN_AND_CALL_MOUSE_GTF`.
paper_mouse_ref_sirv_gtf <- function() {
  ev <- Sys.getenv("JOIN_AND_CALL_MOUSE_GTF", unset = "")
  if (nzchar(ev)) {
    return(paper_normalize_path(ev, must_work = FALSE))
  }
  file.path(paper_nih_data_root(), "mm39.ncbiRefSeq_SIRV.gtf")
}

#' Exonic transcript lengths TSV (from `get_transcript_lengths.py` on the mouse GTF).
#' Override with env var `JOIN_AND_CALL_MOUSE_TRANSCRIPT_LENGTHS`.
paper_mouse_transcript_lengths_tsv <- function() {
  ev <- Sys.getenv("JOIN_AND_CALL_MOUSE_TRANSCRIPT_LENGTHS", unset = "")
  if (nzchar(ev)) {
    return(paper_normalize_path(ev, must_work = FALSE))
  }
  paste0(tools::file_path_sans_ext(paper_mouse_ref_sirv_gtf()), ".transcript_lengths.tsv")
}

#' Transcript ID to gene symbol TSV (built from GTF `transcript` rows; optional cache).
paper_mouse_transcript_gene_index_tsv <- function() {
  ev <- Sys.getenv("JOIN_AND_CALL_MOUSE_TRANSCRIPT_GENES", unset = "")
  if (nzchar(ev)) {
    return(paper_normalize_path(ev, must_work = FALSE))
  }
  paste0(tools::file_path_sans_ext(paper_mouse_ref_sirv_gtf()), ".transcript_genes.tsv")
}

#' Mouse reference GTF without SIRV spike-ins (`mm39.ncbiRefSeq.gtf`).
paper_mouse_reference_gtf <- function() {
  ev <- Sys.getenv("JOIN_AND_CALL_MOUSE_REF_GTF", unset = "")
  if (nzchar(ev)) {
    return(paper_normalize_path(ev, must_work = FALSE))
  }
  file.path(paper_nih_data_root(), "mm39.ncbiRefSeq.gtf")
}

#' documenting_NIH fabian tree (reports, stringtie merge, tusco, etc.).
#' Override with env var `JOIN_AND_CALL_NIH_REPORTS`.
paper_documenting_nih_root <- function() {
  ev <- Sys.getenv("JOIN_AND_CALL_NIH_REPORTS", unset = "")
  if (nzchar(ev)) {
    return(paper_normalize_path(ev, must_work = FALSE))
  }
  default <- if (paper_is_windows()) {
    "C:/Users/jetzi/other_repos/documenting_NIH/fabian"
  } else {
    "/mnt/c/Users/jetzi/other_repos/documenting_NIH/fabian"
  }
  paper_normalize_path(default, must_work = FALSE)
}

#' Locate `src/plotting/R` regardless of knitr working directory.
paper_find_r_source_dir <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  candidates <- unique(c(
    file.path(wd, "src", "plotting", "R"),
    file.path(wd, "R"),
    if (basename(wd) == "plotting") file.path(wd, "R") else NA_character_,
    if (basename(wd) == "src") file.path(wd, "plotting", "R") else NA_character_
  ))
  candidates <- candidates[!is.na(candidates)]
  for (cand in candidates) {
    marker <- file.path(cand, "00_figure_config.R")
    if (file.exists(marker)) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Cannot find src/plotting/R. Knit from the repo or src/plotting, or set JOIN_AND_CALL_REPO.",
    call. = FALSE
  )
}

#' Active plotting R source directory (`src/plotting/R`).
paper_r_source_dir <- function() {
  if (exists("r_source_dir", inherits = TRUE)) {
    cand <- normalizePath(r_source_dir, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(cand, "00_figure_config.R"))) {
      return(cand)
    }
  }
  paper_find_r_source_dir()
}

#' Repository root for plot output paths.
paper_repo_root <- function() {
  ev <- Sys.getenv("JOIN_AND_CALL_REPO", unset = "")
  if (nzchar(ev)) {
    return(paper_normalize_path(ev, must_work = TRUE))
  }
  src <- paper_r_source_dir()
  root <- normalizePath(file.path(src, "..", "..", ".."), winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "src", "plotting", "R", "00_figure_config.R"))) {
    stop(
      "Resolved repo root does not look like join_and_call_paper: ", root,
      call. = FALSE
    )
  }
  root
}

#' Mouse / main paper figure output directory (repo root).
plot_creation_dir <- function() {
  file.path(paper_repo_root(), "plot_creation")
}

#' ONT R10 SY5Y figure output directory (repo root).
plot_creation_ont_r10_dir <- function() {
  file.path(paper_repo_root(), "plot_creation_ont_r10")
}

#' TUSCO/SIRV metric `.RData` tree (`reports/tusco_results` at repo root).
paper_tusco_results_dir <- function() {
  file.path(paper_repo_root(), "reports", "tusco_results")
}

#' Design-to-print scale factor.
#' Figures are authored on a uniform 18-inch reference canvas and printed at
#' Nature Communications 2-column width (180 mm = 7.0866 in). Multiplying the
#' reference sizes (fonts, explicit geom sizes, export dimensions) by PAPER_SCALE
#' emits directly at 180 mm with text at the required 5-7 pt.
PAPER_SCALE <- 7.0866 / 18  # ~= 0.3937

#' Canonical figure widths/heights in *reference design inches* (width, height).
#' Every class uses an 18-inch reference width so a single PAPER_SCALE yields a
#' uniform 180 mm print width and uniform final font sizes. Heights preserve each
#' figure's aspect ratio; `paper_figure_size()` applies PAPER_SCALE at export.
FIG_SIZE <- list(
  mouse_10panel      = c(18, 12),
  #' ONT R10 SY5Y 2×2 (3 tools + legend): ~one mouse 10-panel cell per subplot.
  sy5y_2x2           = c(18, 13.5),
  tpm_lines_10panel  = c(18, 12),
  sy5y_tpm_2x2       = c(18, 13.5),
  legend_strip       = c(18, 2.5),
  sy5y_legend_strip  = c(18, 2.5),
  read_level_combo   = c(18, 10),
  #' ONT R10 read-level: one faceted panel (4 replicates) or two-panel combo.
  sy5y_read_level_panel = c(18, 10.8),
  sy5y_read_level_combo = c(18, 7.7),
  #' 3 violin panels rearranged to a wider/shorter layout (was portrait 8×14).
  ref_tusco_sirv     = c(18, 13),
  cross_tool_upset   = c(18, 16),
  #' Mouse UJC UpSet 10-panel: taller than the bar grid so each 7-set matrix
  #' has vertical room (dots/labels otherwise overlap).
  mouse_upset        = c(18, 17),
  cross_tool_bar     = c(18, 10.3),
  st_upset           = c(18, 14),
  presentation_wide  = c(14, 6.5),
  sirv_dumbbell      = c(18, 12),
  #' 2×2 ranking grid: ~same per-panel size as one cell in `mouse_10panel` (18×12, 4×3 layout).
  resource_ranking_2x2 = c(18, 13.5)
)

#' Canonical y-axis labels for UJC occurrence / expression panels.
PAPER_UJC_YLAB <- list(
  comb_bar           = "# of UJCs",
  comb_fl_bar        = "# of reads in UJC",
  perc_comb_bar      = "% of UJCs",
  perc_comb_fl_bar   = "% of reads in UJC",
  ujc_stack          = "Mean # of UJCs",
  perc_ujc_stack     = "Mean % of UJCs",
  ujc_fl_stack       = "Mean # of reads in UJC",
  perc_ujc_fl_stack  = "Mean % of reads in UJC"
)

#' Typography sizes in points (ggplot2 element_text).
#' `tick` is the minimum size used anywhere on a figure; all other roles are >= tick.
PAPER_FONTS <- list(
  tick         = 13,
  axis         = 13,
  panel        = 16,
  figure       = 22,
  subtitle     = 16,
  legend_text  = 13,
  legend_title = 14,
  strip        = 13,
  tag          = 14,
  caption      = 13
)

#' UpSet matrix set labels (J&C, C&J, B1–B5) in mouse `all_plots` 10-panel grids only.
PAPER_MOUSE_UPSET_MATRIX_LABEL <- 13L

#' ComplexUpset matrix dot/segment sizes, in FINAL print units (mm-ish) — the
#' UpSet is rendered directly at 180 mm, so these are NOT multiplied by
#' PAPER_SCALE. Tune here if matrix dots merge or look too small.
PAPER_UPSET_DOT_SIZE     <- 0.9
PAPER_UPSET_SEGMENT_SIZE <- 0.35

#' Resolve FIG_SIZE name or numeric c(width, height) to *print* inches.
#' Reference design sizes are multiplied by PAPER_SCALE so the exported PDF is
#' 180 mm wide with 5-7 pt text (see PAPER_SCALE). A numeric c(width, height) is
#' also treated as reference design inches and scaled.
paper_figure_size <- function(size) {
  if (is.character(size) && length(size) == 1L) {
    wh <- FIG_SIZE[[size]]
    if (is.null(wh)) {
      stop("Unknown FIG_SIZE class: ", size, call. = FALSE)
    }
    return(wh * PAPER_SCALE)
  }
  if (is.numeric(size) && length(size) == 2L) {
    return(size * PAPER_SCALE)
  }
  stop("`size` must be a FIG_SIZE name or numeric c(width, height).", call. = FALSE)
}

#' Save a figure PDF with canonical device and optional FIG_SIZE class.
paper_ggsave <- function(
    filename,
    plot,
    size = "mouse_10panel",
    width = NULL,
    height = NULL,
    device = grDevices::cairo_pdf,
    units = "in",
    ...) {
  if (!is.null(width) && !is.null(height)) {
    wh <- c(width, height)
  } else {
    wh <- paper_figure_size(size)
    if (!is.null(width)) {
      wh[1] <- width
    }
    if (!is.null(height)) {
      wh[2] <- height
    }
  }
  out_dir <- dirname(filename)
  if (nzchar(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    device = device,
    width = wh[1],
    height = wh[2],
    units = units,
    ...
  )
  invisible(filename)
}

# ---- Source Data collection -------------------------------------------------
# Nature Communications requires one Source Data workbook with one sheet per
# data-bearing figure, in the TUSCO/ncomms layout: row 1 = descriptive title,
# row 2 = blank, row 3 = column header, row 4+ = tidy data.
#
# Figures are built across several .Rmd files that are run independently, so
# each captured sheet is persisted to disk as an .rds fragment; a final
# write_source_data() call assembles the single .xlsx. Requires the 'openxlsx'
# package (used via ::, so no library() call is needed in the drivers).

#' Directory holding per-sheet Source Data fragments.
source_data_fragment_dir <- function() {
  file.path(plot_creation_dir(), "source_data_fragments")
}

#' Fragment filename for a sheet name (sanitised).
.source_data_frag_file <- function(sheet_name) {
  safe <- gsub("[^A-Za-z0-9]+", "_", sheet_name)
  file.path(source_data_fragment_dir(), paste0(safe, ".rds"))
}

#' Register one Source Data sheet, persisted as an .rds fragment.
#' @param sheet_name Excel sheet name, e.g. "Fig 2", "Supplementary Fig 11".
#' @param title Descriptive one-line title placed in row 1 of the sheet.
#' @param df Tidy data.frame of the plotted values (full precision).
add_source_sheet <- function(sheet_name, title, df) {
  stopifnot(is.character(sheet_name), length(sheet_name) == 1L, is.data.frame(df))
  dir.create(source_data_fragment_dir(), recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(sheet_name = sheet_name, title = title, df = as.data.frame(df)),
    .source_data_frag_file(sheet_name)
  )
  message("Source Data captured: ", sheet_name, " (", nrow(df), " rows)")
  invisible(sheet_name)
}

#' Delete all captured Source Data fragments (call once before a full rebuild).
source_data_reset <- function() {
  d <- source_data_fragment_dir()
  if (dir.exists(d)) {
    unlink(list.files(d, pattern = "\\.rds$", full.names = TRUE))
  }
  invisible(NULL)
}

#' Assemble captured fragments into one Source Data .xlsx (TUSCO layout).
#' @param path Output .xlsx path.
#' @param order Optional character vector of sheet names giving output order;
#'   any unlisted sheets follow in captured (alphabetical fragment) order.
write_source_data <- function(path, order = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for write_source_data().", call. = FALSE)
  }
  d <- source_data_fragment_dir()
  files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) {
    stop("No Source Data fragments found in ", d, call. = FALSE)
  }
  frags <- lapply(files, readRDS)
  names(frags) <- vapply(frags, function(x) x$sheet_name, character(1))
  if (!is.null(order)) {
    frags <- frags[c(order[order %in% names(frags)], setdiff(names(frags), order))]
  }
  wb <- openxlsx::createWorkbook()
  bold <- openxlsx::createStyle(textDecoration = "bold")
  used <- character(0)
  for (frag in frags) {
    sn <- substr(frag$sheet_name, 1, 31) # Excel sheet-name limit
    base_sn <- sn
    i <- 1L
    while (sn %in% used) {
      i <- i + 1L
      sn <- paste0(substr(base_sn, 1, 28), "_", i)
    }
    used <- c(used, sn)
    openxlsx::addWorksheet(wb, sn)
    openxlsx::writeData(wb, sn, x = frag$title, startCol = 1, startRow = 1)
    openxlsx::writeData(
      wb, sn,
      x = frag$df, startCol = 1, startRow = 3, headerStyle = bold
    )
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  message("Source Data written: ", path, " (", length(frags), " sheets)")
  invisible(path)
}
