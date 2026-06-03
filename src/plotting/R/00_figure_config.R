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

#' Canonical figure widths/heights in inches (width, height).
FIG_SIZE <- list(
  mouse_10panel      = c(18, 12),
  sy5y_2x2           = c(18, 12),
  tpm_lines_10panel  = c(18, 12),
  sy5y_tpm_2x2       = c(18, 12),
  legend_strip       = c(18, 2.5),
  read_level_combo   = c(18, 10),
  ref_tusco_sirv     = c(8, 14),
  cross_tool_upset   = c(18, 16),
  cross_tool_bar     = c(14, 8),
  st_upset           = c(18, 14),
  presentation_wide  = c(14, 6.5),
  sirv_dumbbell      = c(18, 12)
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
  tick         = 11,
  axis         = 12,
  panel        = 16,
  figure       = 22,
  subtitle     = 16,
  legend_text  = 12,
  legend_title = 14,
  strip        = 11,
  tag          = 14,
  caption      = 12
)

#' Resolve FIG_SIZE name or numeric c(width, height).
paper_figure_size <- function(size) {
  if (is.character(size) && length(size) == 1L) {
    wh <- FIG_SIZE[[size]]
    if (is.null(wh)) {
      stop("Unknown FIG_SIZE class: ", size, call. = FALSE)
    }
    return(wh)
  }
  if (is.numeric(size) && length(size) == 2L) {
    return(size)
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
