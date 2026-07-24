# assemble_source_data.R
# ---------------------------------------------------------------------------
# Master Source Data assembly for the Nature Communications submission.
#
# Each figure-building driver (all_plots.Rmd, read_figure_plots.Rmd,
# sirv_comparison_plots.Rmd, cross_tool_upset_plots.Rmd, resource_plots.Rmd,
# filter_level_line_plots.Rmd, tusco_metrics.Rmd, and the four SY5Y drivers)
# calls add_source_sheet(), which persists one .rds fragment per sheet to
#   plot_creation/source_data_fragments/
# This script assembles those fragments into a single source_data.xlsx in the
# TUSCO/ncomms layout (title row, blank row, header row, tidy data).
#
# USAGE (from RStudio, working dir = repo root or src/plotting):
#   1. (optional, for a clean build) source this file with reset first:
#          source("src/plotting/assemble_source_data.R"); reset_source_data()
#      OR just call source_data_reset() once, then re-knit every driver.
#   2. Knit ALL figure drivers so every fragment is (re)written to disk.
#   3. source("src/plotting/assemble_source_data.R")
#      -> writes plot_creation/source_data.xlsx
#
# Fragments persist across R sessions, so the drivers and this script do not
# need to run in the same session — only before this script.
# ---------------------------------------------------------------------------

r_source_dir <- if (dir.exists("R")) {
  "R"
} else if (dir.exists("src/plotting/R")) {
  "src/plotting/R"
} else {
  stop("Cannot find the R/ helper directory; set the working directory to the ",
       "repo root or src/plotting.")
}
source(file.path(r_source_dir, "00_figure_config.R"))

# Convenience wrapper so callers can wipe stale fragments before a full rebuild.
reset_source_data <- function() {
  source_data_reset()
  message("Cleared Source Data fragments in ", source_data_fragment_dir())
}

# Desired workbook sheet order: main figures F1..F10 (F1 is a schematic with no
# data), then supplementary SF1..SF26. Heterogeneous figures are split into
# a/b sub-sheets (read-level bars vs lengths; ref/TUSCO/SIRV length vs
# expression). Any captured sheet not listed here is appended in fragment order.
source_data_sheet_order <- c(
  # --- Main figures ---
  "Fig 2",
  "Fig 3a", "Fig 3b",
  "Fig 4",
  "Fig 5",
  "Fig 6",
  "Fig 7",
  "Fig 8",
  "Fig 9",
  "Fig 10",
  # --- Supplementary figures ---
  "Supplementary Fig 1a", "Supplementary Fig 1b",
  "Supplementary Fig 2",
  "Supplementary Fig 3",
  "Supplementary Fig 4",
  "Supplementary Fig 5",
  "Supplementary Fig 6",
  "Supplementary Fig 7",
  "Supplementary Fig 8",
  "Supplementary Fig 9",
  "Supplementary Fig 10",
  "Supplementary Fig 11a", "Supplementary Fig 11b",
  "Supplementary Fig 12a", "Supplementary Fig 12b",
  "Supplementary Fig 13",
  "Supplementary Fig 14",
  "Supplementary Fig 15",
  "Supplementary Fig 16",
  "Supplementary Fig 17",
  "Supplementary Fig 18",
  # --- Supplementary results (SH-SY5Y; former SRF1-8) ---
  "Supplementary Fig 19",
  "Supplementary Fig 20",
  "Supplementary Fig 21",
  "Supplementary Fig 22a", "Supplementary Fig 22b",
  "Supplementary Fig 23",
  "Supplementary Fig 24",
  "Supplementary Fig 25",
  "Supplementary Fig 26"
)

# Report which expected sheets are present / missing (so a forgotten driver run
# is obvious rather than silently dropped from the workbook).
frag_dir <- source_data_fragment_dir()
present <- if (dir.exists(frag_dir)) {
  vapply(
    list.files(frag_dir, pattern = "\\.rds$", full.names = TRUE),
    function(f) readRDS(f)$sheet_name,
    character(1)
  )
} else {
  character(0)
}
missing <- setdiff(source_data_sheet_order, present)
extra   <- setdiff(present, source_data_sheet_order)
message("Source Data fragments found: ", length(present), " of ",
        length(source_data_sheet_order), " expected sheets.")
if (length(missing)) {
  warning("Missing sheets (driver not run since last reset?): ",
          paste(missing, collapse = ", "), call. = FALSE)
}
if (length(extra)) {
  message("Extra (unlisted) sheets, appended at end: ",
          paste(extra, collapse = ", "))
}

out_path <- file.path(plot_creation_dir(), "source_data.xlsx")
write_source_data(out_path, order = source_data_sheet_order)
message("Wrote ", out_path)
