# 07_st_merge_helpers.R
# StringTie-merge specific configuration, data loading, and combination helpers
# Depends on: 01_config.R (base_path), 02_themes.R (xaxislevelsF1, xaxislabelsF1)

isoseq_st_merge_path <- file.path(base_path, "data", "output", "stringtie_merge_results", "pacbio")
ont_st_merge_path    <- file.path(base_path, "data", "output", "stringtie_merge_results", "ont")

st_merge_dirs <- list(
  IsoSeq = list(
    IsoQuant      = file.path(isoseq_st_merge_path, "isoquant"),
    FLAIR         = file.path(isoseq_st_merge_path, "flair_ar_sr"),
    TALON         = file.path(isoseq_st_merge_path, "talon"),
    Mandalorion   = file.path(isoseq_st_merge_path, "mandalorion"),
    Bambu         = file.path(isoseq_st_merge_path, "bambu"),
    isoseq_sqanti = file.path(isoseq_st_merge_path, "isoseq_pipeline")
  ),
  ONT = list(
    IsoQuant = file.path(ont_st_merge_path, "isoquant"),
    FLAIR    = file.path(ont_st_merge_path, "flair_ar_sr"),
    TALON    = file.path(ont_st_merge_path, "talon"),
    Bambu    = file.path(ont_st_merge_path, "bambu")
  )
)

# Map raw SQANTI structural_category names to the short names used by cat.palette
sqanti_category_map <- setNames(xaxislabelsF1, xaxislevelsF1)

# Extended sample names/labels for 8-bar plots (J&C, C&J TAMA, C&J ST, 5 individuals)
ext_sample_names_Brain  <- c("Join&Call", "Call&Join", "Call&Join (ST)", paste0("Brain ", 1:5))
ext_sample_names_Kidney <- c("Join&Call", "Call&Join", "Call&Join (ST)", paste0("Kidney ", 1:5))

ext_sample_labels_Brain  <- c("J&C", "C&J (TAMA)", "C&J (ST)", paste0("B", 1:5))
ext_sample_labels_Kidney <- c("J&C", "C&J (TAMA)", "C&J (ST)", paste0("K", 1:5))


check_st_merge_paths <- function(dirs, parent = "") {
  for (name in names(dirs)) {
    path_or_sublist <- dirs[[name]]
    current <- if (parent == "") name else paste(parent, name, sep = " -> ")

    if (is.list(path_or_sublist)) {
      check_st_merge_paths(path_or_sublist, current)
    } else {
      if (!dir.exists(path_or_sublist)) {
        cat(sprintf("ST path does NOT exist: [%s] %s\n", current, path_or_sublist))
      } else {
        cat(sprintf("ST path OK: [%s]\n", current))
        for (cond in c("B100K0", "B0K100")) {
          f <- file.path(path_or_sublist, paste0(cond, "_STMERGE_classification.txt"))
          if (!file.exists(f)) {
            cat(sprintf("  Missing: %s\n", basename(f)))
          } else {
            cat(sprintf("  Found:   %s\n", basename(f)))
          }
        }
      }
    }
  }
}


load_st_classification <- function(st_dir, tissue) {
  filename <- if (tissue == "brain") {
    "B100K0_STMERGE_classification.txt"
  } else {
    "B0K100_STMERGE_classification.txt"
  }

  df <- read.delim(file.path(st_dir, filename), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)

  df$structural_category <- sqanti_category_map[df$structural_category]
  df <- df[!is.na(df$structural_category), ]

  sample_code <- if (tissue == "brain") "B_ST" else "K_ST"
  df$sample <- sample_code

  df[, c("isoform", "chrom", "structural_category", "sample")]
}


process_classification_only <- function(Bclass_df_list, Kclass_df_list, fl_threshold = 1) {
  for (df_name in names(Bclass_df_list)) {
    df <- Bclass_df_list[[df_name]]
    if ("FL" %in% names(df) && mean(is.na(df$FL)) > 0.95) {
      stop(paste("Error: Data frame", df_name, "has over 95% NA values in the 'FL' column."))
    }
    Bclass_df_list[[df_name]] <- df %>% filter(!is.na(FL) & FL >= fl_threshold)
  }

  kidney_sample_mapping <- setNames(paste0("K", 1:7), paste0("B", 1:7))
  for (df_name in names(Kclass_df_list)) {
    df <- Kclass_df_list[[df_name]]
    if ("FL" %in% names(df) && mean(is.na(df$FL)) > 0.95) {
      stop(paste("Error: K Data frame", df_name, "has over 95% NA values in the 'FL' column."))
    }
    if ("sample" %in% names(df)) {
      df$sample <- recode(df$sample, !!!kidney_sample_mapping)
    }
    Kclass_df_list[[df_name]] <- df %>% filter(!is.na(FL) & FL >= fl_threshold)
  }

  list(Bclass_df_list = Bclass_df_list, Kclass_df_list = Kclass_df_list)
}


build_extended_combined_df <- function(filtered_df_list, st_df) {
  cols <- c("sample", "structural_category", "chrom")
  existing_dfs <- lapply(filtered_df_list, function(df) df[, cols])

  original_codes <- vapply(existing_dfs, function(df) as.character(unique(df$sample)[1]),
                           character(1))
  st_code <- as.character(unique(st_df$sample)[1])
  extended_levels <- c(original_codes[1:2], st_code, original_codes[3:length(original_codes)])

  combined <- do.call(rbind, c(existing_dfs, list(st_df[, cols])))
  combined$sample <- factor(combined$sample, levels = extended_levels)
  combined
}
