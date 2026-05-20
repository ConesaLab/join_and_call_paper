# 01_config_ont_r10.R
# Paths for ONT R10 SY5Y human reports (does not modify mouse 01_config.R).

if (!exists("ont_r10_repo_root", inherits = TRUE)) {
  ev <- Sys.getenv("JOIN_AND_CALL_REPO", unset = "")
  if (nzchar(ev)) {
    ont_r10_repo_root <- normalizePath(ev, mustWork = TRUE)
  } else if (exists("r_source_dir", inherits = TRUE)) {
    # r_source_dir is "R" or "src/plotting/R" — three parents up = repository root
    ont_r10_repo_root <- normalizePath(
      file.path(r_source_dir, "..", "..", ".."),
      mustWork = TRUE
    )
  } else {
    stop("Define `r_source_dir` before sourcing 01_config_ont_r10.R, or set JOIN_AND_CALL_REPO.")
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

ont_r10_sy5y_ylims <- list(
  structural   = NULL,
  expression   = NULL,
  comb         = NULL,
  comb_bar     = NULL,
  comb_fl_bar  = NULL,
  ujc_curve    = NULL,
  ujc_stack    = NULL,
  ujc_fl_stack = NULL
)

check_paths_ont_r10 <- function() {
  required <- c("Bclass_df_list.RData")
  for (nm in names(ont_r10_report_dirs)) {
    d <- ont_r10_report_dirs[[nm]]
    if (!dir.exists(d)) {
      cat(sprintf("[ONT R10] Missing directory [%s]: %s\n", nm, d))
      next
    }
    cat(sprintf("[ONT R10] OK [%s]: %s\n", nm, d))
    for (f in required) {
      fp <- file.path(d, f)
      if (!file.exists(fp)) {
        cat(sprintf("  Missing file: %s (re-run prepare_report after generate_report_rdata.R update)\n", f))
      } else {
        cat(sprintf("  Found: %s\n", f))
      }
    }
    d_sy5y <- file.path(d, "data", "classification_SY5Y")
    if (!dir.exists(d_sy5y)) {
      cat(sprintf("  Optional: no data/classification_SY5Y under report\n"))
    }
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
