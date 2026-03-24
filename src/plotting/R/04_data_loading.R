# 04_data_loading.R
# Data loading, preprocessing, and summary extraction helpers

load_df_lists <- function(paths) {
  knitr::opts_knit$set(root.dir = paths$src_dir)
  
  load(paste0(paths$src_dir, "/Bclass_df_list.RData"))
  load(paste0(paths$src_dir, "/Kclass_df_list.RData"))
  
  names(Bclass_df_list) <- sample_names_Brain
  names(Kclass_df_list) <- sample_names_Kidney
  
  if (all(!sapply(paste0("B3", 1:5), function(pattern) any(grepl(pattern, names(Bclass_df_list$`Call&Join`), ignore.case = TRUE))))) {
    Bquant_cj <- read_tsv(paths$quant_cj_brain)
    Bclass_df_list$`Call&Join` <- merge(x = Bclass_df_list$`Call&Join`, y = Bquant_cj, by.x = 1, by.y = 1, all.x = TRUE, suffixes = c("", "_FL"))
  }
  if (all(!sapply(paste0("K3", 1:5), function(pattern) any(grepl(pattern, names(Kclass_df_list$`Call&Join`), ignore.case = TRUE))))) {
    Kquant_cj <- read_tsv(paths$quant_cj_kidney)
    Kclass_df_list$`Call&Join` <- merge(x = Kclass_df_list$`Call&Join`, y = Kquant_cj, by.x = 1, by.y = 1, all.x = TRUE, suffixes = c("", "_FL"))
  }
  
  if (!("FL" %in% names(Bclass_df_list$`Join&Call`)) || all(is.na(Bclass_df_list$`Join&Call`$FL))) {
    b_columns <- grep("B3[1-5]", names(Bclass_df_list$`Join&Call`), value = TRUE, ignore.case = TRUE)
    if (length(b_columns) != 5) {
      stop("Error: The number of 'B31' to 'B35' columns in J&C is not exactly 5.")
    }
    Bclass_df_list$`Join&Call`$FL <- rowSums(Bclass_df_list$`Join&Call`[, b_columns], na.rm = TRUE)
  }
  if (!("FL" %in% names(Kclass_df_list$`Join&Call`)) || all(is.na(Kclass_df_list$`Join&Call`$FL))) {
    k_columns <- grep("K3[1-5]", names(Kclass_df_list$`Join&Call`), value = TRUE, ignore.case = TRUE)
    if (length(k_columns) != 5) {
      stop("Error: The number of 'K31' to 'K35' columns in J&C is not exactly 5.")
    }
    Kclass_df_list$`Join&Call`$FL <- rowSums(Kclass_df_list$`Join&Call`[, k_columns], na.rm = TRUE)
  }
  
  if (!("FL" %in% names(Bclass_df_list$`Call&Join`)) || all(is.na(Bclass_df_list$`Call&Join`$FL))) {
    b_columns <- grep("B3[1-5]", names(Bclass_df_list$`Call&Join`), value = TRUE, ignore.case = TRUE)
    if (length(b_columns) != 5) {
      stop("Error: The number of 'B31' to 'B35' columns in C&J is not exactly 5.")
    }
    Bclass_df_list$`Call&Join`$FL <- rowSums(Bclass_df_list$`Call&Join`[, b_columns], na.rm = TRUE)
  }
  if (!("FL" %in% names(Kclass_df_list$`Call&Join`)) || all(is.na(Kclass_df_list$`Call&Join`$FL))) {
    k_columns <- grep("K3[1-5]", names(Kclass_df_list$`Call&Join`), value = TRUE, ignore.case = TRUE)
    if (length(k_columns) != 5) {
      stop("Error: The number of 'K31' to 'K35' columns in C&J is not exactly 5.")
    }
    Kclass_df_list$`Call&Join`$FL <- rowSums(Kclass_df_list$`Call&Join`[, k_columns], na.rm = TRUE)
  }
  return(list(
    Bclass_df_list = Bclass_df_list,
    Kclass_df_list = Kclass_df_list
  ))
}


process_and_plot <- function(Bclass_df_list, Kclass_df_list, technology, fl_threshold=0) {
  
  for (df_name in names(Bclass_df_list)) {
    df <- Bclass_df_list[[df_name]]
    
    if ("FL" %in% names(df) && mean(is.na(df$FL)) > 0.95) {
      stop(paste("Error: Data frame", df_name, "has over 95% NA values in the 'FL' column."))
    }
    
    Bclass_df_list[[df_name]] <- df %>% 
      filter(!is.na(FL) & FL >= fl_threshold)
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
    
    Kclass_df_list[[df_name]] <- df %>% 
      filter(!is.na(FL) & FL >= fl_threshold)
  }

  class_combined_df_Brain <- do.call(rbind, lapply(Bclass_df_list, function(class_df) class_df[, c("sample", "structural_category", "chrom", "FL")]))
  plots_Brain <- plot_classification_data(class_combined_df_Brain, sample_labels_Brain)

  class_combined_df_Kidney <- do.call(rbind, lapply(Kclass_df_list, function(class_df) class_df[, c("sample", "structural_category", "chrom", "FL")]))
  plots_Kidney <- plot_classification_data(class_combined_df_Kidney, sample_labels_Kidney)

  total_transcripts_plot <- ggarrange(
    annotate_figure(ggarrange(plots_Brain[[1]], plots_Brain[[2]], nrow = 2), top = "Brain"),
    annotate_figure(ggarrange(plots_Kidney[[1]], plots_Kidney[[2]], nrow = 2, common.legend = TRUE, legend = "right"), top = "Kidney"),
    widths = c(2, 3)
  )

  class_combined_df_Brain_Mouse <- class_combined_df_Brain[!grepl("SIRV", class_combined_df_Brain$chrom), ]
  plots_Brain_Mouse <- plot_classification_data(class_combined_df_Brain_Mouse, sample_labels_Brain)

  class_combined_df_Kidney_Mouse <- class_combined_df_Kidney[!grepl("SIRV", class_combined_df_Kidney$chrom), ]
  plots_Kidney_Mouse <- plot_classification_data(class_combined_df_Kidney_Mouse, sample_labels_Kidney)

  mouse_transcripts_plot <- ggarrange(
    annotate_figure(ggarrange(plots_Brain_Mouse[[1]], plots_Brain_Mouse[[2]], nrow = 2), top = "Brain"),
    annotate_figure(ggarrange(plots_Kidney_Mouse[[1]], plots_Kidney_Mouse[[2]], nrow = 2, common.legend = TRUE, legend = "right"), top = "Kidney"),
    widths = c(2, 3)
  )

  class_combined_df_Brain_SIRVs <- class_combined_df_Brain[grepl("SIRV", class_combined_df_Brain$chrom), ]
  plots_Brain_SIRVs <- plot_classification_data(class_combined_df_Brain_SIRVs, sample_labels_Brain)

  class_combined_df_Kidney_SIRVs <- class_combined_df_Kidney[grepl("SIRV", class_combined_df_Kidney$chrom), ]
  plots_Kidney_SIRVs <- plot_classification_data(class_combined_df_Kidney_SIRVs, sample_labels_Kidney)

  sirv_transcripts_plot <- ggarrange(
    annotate_figure(ggarrange(plots_Brain_SIRVs[[1]], plots_Brain_SIRVs[[2]], nrow = 2), top = "Brain"),
    annotate_figure(ggarrange(plots_Kidney_SIRVs[[1]], plots_Kidney_SIRVs[[2]], nrow = 2, common.legend = TRUE, legend = "right"), top = "Kidney"),
    widths = c(2, 3)
  )
  
  read_numbers_brain <- pacbio_read_numbers_brain
  read_numbers_kidney <- pacbio_read_numbers_kidney
  
  if (technology == "ONT"){
    read_numbers_brain <- ont_read_numbers_brain
    read_numbers_kidney <- ont_read_numbers_kidney
  }
  
  plots_Brain_Expression <- plot_classification_expression_data(class_combined_df_Brain, read_numbers_brain, sample_codes_brain, sample_labels_Brain)
  plots_Kidney_Expression <- plot_classification_expression_data(class_combined_df_Kidney, read_numbers_kidney, sample_codes_kidney, sample_labels_Kidney)

  UJC_counts_Brain <- UJC_count_matrix(Bclass_df_list, sample_labels_Brain)
  UJC_counts_Brain[, sample_labels_Brain[1]] <- rowSums(UJC_counts_Brain[, 2:ncol(UJC_counts_Brain)], na.rm = TRUE)
  UJC_counts_Brain[UJC_counts_Brain[, sample_labels_Brain[1]] == 0, sample_labels_Brain[1]] <- NA
  UJC_counts_Brain[, 2:ncol(UJC_counts_Brain)] <- apply(UJC_counts_Brain[, 2:ncol(UJC_counts_Brain)], 2, cpm)
  
  brain_upset_plots <- create_upset_plot(Bclass_df_list, sample_labels_Brain, method)

  UJC_counts_Kidney <- UJC_count_matrix(Kclass_df_list, sample_labels_Kidney)
  UJC_counts_Kidney[, sample_labels_Kidney[1]] <- rowSums(UJC_counts_Kidney[, 2:ncol(UJC_counts_Kidney)], na.rm = TRUE)
  UJC_counts_Kidney[UJC_counts_Kidney[, sample_labels_Kidney[1]] == 0, sample_labels_Kidney[1]] <- NA
  UJC_counts_Kidney[, 2:ncol(UJC_counts_Kidney)] <- apply(UJC_counts_Kidney[, 2:ncol(UJC_counts_Kidney)], 2, cpm)
  
  kidney_upset_plots <- create_upset_plot(Kclass_df_list, sample_labels_Kidney, method)
  
  brain_compare_plots <- compare_isoform_plots(Bclass_df_list)
  kidney_compare_plots <- compare_isoform_plots(Kclass_df_list)

  return(list(
    total_transcripts_plot = total_transcripts_plot,
    mouse_transcripts_plot = mouse_transcripts_plot,
    sirv_transcripts_plot = sirv_transcripts_plot,
    brain_plots = plots_Brain,
    kidney_plots = plots_Kidney,
    sirv_brain_plots = plots_Brain_SIRVs,
    sirv_kidney_plots = plots_Kidney_SIRVs,
    brain_upset_plot = brain_upset_plots[["upset"]],
    brain_comb_plot = brain_upset_plots[["comb"]],
    brain_perc_comb_plot = brain_upset_plots[["perc_comb"]],
    brain_comb_bar_plot = brain_upset_plots[["comb_bar"]],
    brain_perc_comb_bar_plot = brain_upset_plots[["perc_comb_bar"]],
    brain_comb_fl_bar_plot = brain_upset_plots[["comb_fl_bar"]],
    brain_perc_comb_fl_bar_plot = brain_upset_plots[["perc_comb_fl_bar"]],
    brain_ujc_curve_plot = brain_upset_plots[["ujc_curve"]],
    brain_perc_ujc_curve_plot = brain_upset_plots[["perc_ujc_curve"]],
    brain_ujc_stack_plot = brain_upset_plots[["ujc_stack"]],
    brain_perc_ujc_stack_plot = brain_upset_plots[["perc_ujc_stack"]],
    brain_ujc_fl_stack_plot = brain_upset_plots[["ujc_fl_stack"]],
    brain_perc_ujc_fl_stack_plot = brain_upset_plots[["perc_ujc_fl_stack"]],
    brain_ujc_stack_data = brain_upset_plots[["ujc_stack_data"]],
    brain_ujc_fl_stack_data = brain_upset_plots[["ujc_fl_stack_data"]],
    brain_ujc_fl_stack_pct_data = brain_upset_plots[["ujc_fl_stack_pct_data"]],
    kidney_upset_plot = kidney_upset_plots[["upset"]],
    kidney_comb_plot = kidney_upset_plots[["comb"]],
    kidney_perc_comb_plot = kidney_upset_plots[["perc_comb"]],
    kidney_comb_bar_plot = kidney_upset_plots[["comb_bar"]],
    kidney_perc_comb_bar_plot = kidney_upset_plots[["perc_comb_bar"]],
    kidney_comb_fl_bar_plot = kidney_upset_plots[["comb_fl_bar"]],
    kidney_perc_comb_fl_bar_plot = kidney_upset_plots[["perc_comb_fl_bar"]],
    kidney_ujc_curve_plot = kidney_upset_plots[["ujc_curve"]],
    kidney_perc_ujc_curve_plot = kidney_upset_plots[["perc_ujc_curve"]],
    kidney_ujc_stack_plot = kidney_upset_plots[["ujc_stack"]],
    kidney_perc_ujc_stack_plot = kidney_upset_plots[["perc_ujc_stack"]],
    kidney_ujc_fl_stack_plot = kidney_upset_plots[["ujc_fl_stack"]],
    kidney_perc_ujc_fl_stack_plot = kidney_upset_plots[["perc_ujc_fl_stack"]],
    kidney_ujc_stack_data = kidney_upset_plots[["ujc_stack_data"]],
    kidney_ujc_fl_stack_data = kidney_upset_plots[["ujc_fl_stack_data"]],
    kidney_ujc_fl_stack_pct_data = kidney_upset_plots[["ujc_fl_stack_pct_data"]],
    brain_expression_plot = plots_Brain_Expression,
    kidney_expression_plot = plots_Kidney_Expression,
    brain_compare_plots = brain_compare_plots,
    kidney_compare_plots = kidney_compare_plots,
    b_class_df_list = Bclass_df_list,
    k_class_df_list = Kclass_df_list
  ))
}


extract_ujc_data <- function(all_results, data_type, tool, tissue, fl_filter_level = "1") {
  tryCatch({
    data <- all_results[[data_type]][[tool]][[fl_filter_level]][[paste0(tissue, "_ujc_stack_data")]]
    if (!is.null(data)) {
      data %>%
        mutate(
          data_type = data_type,
          tool = tool,
          tissue = tissue
        )
    } else {
      NULL
    }
  }, error = function(e) {
    NULL
  })
}

extract_fl_data <- function(all_results, data_type, tool, tissue, fl_filter_level = "1") {
  tryCatch({
    fl_data <- all_results[[data_type]][[tool]][[fl_filter_level]][[paste0(tissue, "_ujc_fl_stack_data")]]
    pct_data <- all_results[[data_type]][[tool]][[fl_filter_level]][[paste0(tissue, "_ujc_fl_stack_pct_data")]]

    if (!is.null(fl_data) && !is.null(pct_data)) {
      fl_data %>%
        left_join(pct_data %>% select(nr_samples, samples_present, mean_FL_pct, sd_FL_pct),
                  by = c("nr_samples", "samples_present")) %>%
        mutate(
          data_type = data_type,
          tool = tool,
          tissue = tissue
        )
    } else {
      NULL
    }
  }, error = function(e) {
    NULL
  })
}
