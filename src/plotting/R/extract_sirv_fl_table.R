library(tidyverse)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(file.path(script_dir, "01_config.R"))
source(file.path(script_dir, "04_data_loading.R"))

results <- list()

for (platform in names(src_dirs)) {
  for (method in names(src_dirs[[platform]])) {
    cat(sprintf("Loading %s / %s...\n", platform, method))
    
    paths <- get_paths(src_dirs[[platform]][[method]])
    df_lists <- tryCatch(load_df_lists(paths), error = function(e) {
      cat(sprintf("  Skipping: %s\n", e$message))
      return(NULL)
    })
    if (is.null(df_lists)) next
    
    Bclass_df_list <- df_lists$Bclass_df_list
    
    fl_by_sample <- list()
    for (sample_name in names(Bclass_df_list)) {
      df <- Bclass_df_list[[sample_name]]
      sirv_rows <- df[grepl("^SIRV", df$chrom), ]
      
      if (nrow(sirv_rows) == 0 || !("FL" %in% names(sirv_rows))) next
      
      fl_by_sample[[sample_name]] <- sirv_rows %>%
        group_by(chrom, associated_transcript) %>%
        summarise(FL = sum(FL, na.rm = TRUE), .groups = "drop") %>%
        rename(!!sample_name := FL)
    }
    
    if (length(fl_by_sample) == 0) next
    
    merged <- fl_by_sample %>%
      reduce(full_join, by = c("chrom", "associated_transcript")) %>%
      mutate(
        platform = platform,
        method = method
      )
    
    results[[length(results) + 1]] <- merged
  }
}

sirv_fl_table <- bind_rows(results) %>%
  select(platform, method, chrom, associated_transcript, all_of(sample_names_Brain)) %>%
  arrange(platform, method, chrom, associated_transcript)

cat("\n=== Example: IsoSeq / IsoQuant / SIRV1 ===\n\n")
example <- sirv_fl_table %>%
  filter(platform == "IsoSeq", method == "IsoQuant", chrom == "SIRV1")
print(example, n = 50, width = 200)

write_csv(sirv_fl_table, "sirv_fl_per_transcript_brain.csv")
cat(sprintf("\nFull table written: sirv_fl_per_transcript_brain.csv (%d rows)\n", nrow(sirv_fl_table)))
