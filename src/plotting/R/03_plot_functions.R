# 03_plot_functions.R
# Individual plot-generating functions: classification, expression, UpSet, isoform comparison

plot_classification_data <- function(class_combined_df, sample_labels,
                                     x_discrete_limits = NULL,
                                     x_drop_missing_levels = TRUE) {
  if (!is.null(x_discrete_limits)) {
    if (length(sample_labels) != length(x_discrete_limits)) {
      stop(
        "plot_classification_data: length(sample_labels) must equal length(x_discrete_limits)."
      )
    }
    class_combined_df <- class_combined_df %>%
      mutate(sample = factor(.data$sample, levels = x_discrete_limits))
    sx <- ggplot2::scale_x_discrete(
      limits   = x_discrete_limits,
      labels   = sample_labels,
      drop     = x_drop_missing_levels
    )
  } else {
    sx <- ggplot2::scale_x_discrete(labels = sample_labels)
  }

  p1 <- ggplot(class_combined_df, aes(x = sample, y = (..count..), fill = structural_category)) +
    geom_bar() +
    scale_fill_manual(values = cat.palette, labels = c("Full\nSplice Match", "Incomplete\nSplice Match", "Novel\nIn Catalog", "Novel Not\nIn Catalog", "Genic\nGenomic", "Antisense", "Fusion", "Intergenic", "Genic\nIntron")) +
    sx +
    paper_theme() +
    labs(x = NULL) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    ylab("# isoforms") +
    theme(axis.text.x = paper_axis_text_x())

  list(count = p1)
}

plot_classification_expression_data <- function(class_combined_df,
                                                     total_reads,
                                                     sample_codes,
                                                     sample_labels,
                                                     unassigned_color = "white",
                                                     border_color     = "black") {

  aggregated_df <- class_combined_df %>%
    group_by(sample, structural_category) %>%
    summarise(FL = sum(FL), .groups = "drop")

  assigned_df <- aggregated_df %>%
    group_by(sample) %>%
    summarise(assigned = sum(FL), .groups = "drop")

  unassigned_df <- tibble(
    sample = sample_codes,
    total  = unname(total_reads[sample_codes])
  ) %>%
    left_join(assigned_df, by = "sample") %>%
    mutate(
      assigned = replace_na(assigned, 0),
      FL = total - assigned,
      structural_category = "Unassigned"
    ) %>%
    select(sample, structural_category, FL)

  plot_df <- bind_rows(aggregated_df, unassigned_df) %>%
    mutate(
      sample = factor(sample, levels = sample_codes),
      structural_category = factor(structural_category, levels = c("Unassigned", names(cat.palette))),
      border_color = ifelse(structural_category == "Unassigned", "Unassigned", NA_character_)
    )

  extended_fill   <- c(Unassigned = unassigned_color, cat.palette)
  extended_border <- c(Unassigned = border_color,     cat.palette)

  base_theme <- paper_theme() +
    labs(x = NULL) +
    theme(
      legend.position = "none",
      axis.title.x    = element_blank(),
      axis.text.x     = paper_axis_text_x()
    )

  p1 <- ggplot(plot_df, aes(x = sample, y = FL)) +
    geom_bar(
      aes(fill = structural_category, color = border_color),
      stat = "identity", size = 0.3
    ) +
    scale_x_discrete(labels = sample_labels) +
    scale_fill_manual(values = extended_fill, na.value = NA) +
    scale_color_manual(values = extended_border, na.value = NA) +
    base_theme +
    ylab("# reads")

  p2 <- ggplot(plot_df, aes(x = sample, y = FL)) +
    geom_bar(
      aes(fill = structural_category, color = border_color),
      stat = "identity", position = "fill", size = 0.3
    ) +
    scale_x_discrete(labels = sample_labels) +
    scale_fill_manual(values = extended_fill, na.value = NA) +
    scale_color_manual(values = extended_border, na.value = NA) +
    base_theme +
    ylab("% of reads")


  return(list(p1, p2))
}


# Helper functions for UpSet plots
cpm <- function(counts) {
  counts / sum(counts, na.rm = TRUE) * 1e6
}

UJC_count_matrix <- function(class_df_list, sample_labels) {
  counts_list <- list()
  i <- 1
  for (class_df in class_df_list) {
    class_df <- class_df[!grepl("NA", class_df$UJC), c("UJC", "FL")]
    counts_list[[sample_labels[i]]] <- class_df %>% 
      group_by(UJC) %>% 
      summarise(counts = sum(FL))
    i <- i + 1
  }
  counts <- counts_list[[1]]
  for (i in 2:length(counts_list)) {
    counts <- merge(counts, counts_list[[i]], by = "UJC", all = TRUE)
  }
  colnames(counts) <- c("UJC", sample_labels)
  return(counts)
}

create_boxplot <- function(listUJC) {
  max_value <- max(log10(listUJC$counts)) + max(log10(listUJC$counts)) * 0.1
  min_value <- min(log10(listUJC$counts)) - min(log10(listUJC$counts)) * 0.1

  structural_categories <- c("FSM", "ISM", "NIC", "NNC")
  l_plots <- list()
  for (category in structural_categories) {
    tmp <- listUJC[listUJC$structural_category == category, ]
    p <- ggplot(tmp, aes(presence_code, log10(counts), color = structural_category)) +
      geom_boxplot() + ylim(c(min_value, max_value)) +
      scale_color_manual(values = cat.palette, name = "Structural Category") +
      scale_x_discrete(drop = FALSE) + 
      paper_upset_theme()
    l_plots[[category]] <- p
  }
  tmp <- listUJC[!listUJC$structural_category %in% structural_categories, ]
  p <- ggplot(tmp, aes(presence_code, log10(counts), color = "a")) +
    geom_boxplot() + ylim(c(min_value, max_value)) +
    scale_color_manual(values = "#969696", name = "Structural Category") +
    paper_upset_theme()
  l_plots[["rest"]] <- p
  return(l_plots)
}


create_upset_plot <- function(class_df_list, sample_labels, method, n = 10,
                              combination_columns = NULL,
                              include_percentage_plots = FALSE,
                              upset_matrix_label_size = NULL) {
  if (method == "isoseq_sqanti") {
    method <- "IsoSeq + SQANTI3"
  }

  if (length(sample_labels) != length(class_df_list)) {
    nm <- names(class_df_list)
    stop(
      "create_upset_plot: need one sample label per classification table (",
      length(sample_labels), " labels vs ", length(class_df_list), " tables). ",
      if (length(nm)) paste0("Table names: ", paste(nm, collapse = ", "), ". ") else "",
      "Regenerate SY5Y_class_df_list.RData from the full classification+junction FOFN (concat, TAMA, four replicates)."
    )
  }

  listUJC <- do.call(rbind, lapply(class_df_list, function(class_df) class_df[!grepl("NA",class_df$UJC), c("structural_category", "UJC")]))
  listUJC <- distinct(listUJC)
  
  for (i in 1:length(class_df_list)){
    sample_id <- names(class_df_list)[[i]]
    listUJC[, sample_id] <- ifelse(listUJC$UJC %in% class_df_list[[i]]$UJC, 1, 0)
  }
  colnames(listUJC) <- c("structural_category", "UJC", sample_labels)
  
  category_colors <- c(
    "FSM" = "#6BAED6",
    "ISM" = "#FC8D59",
    "NIC" = "#78C679",
    "NNC" = "#EE6A50",
    "Genic\nGenomic" = "#969696",
    "Antisense" = "#66C2A4",
    "Fusion" = "goldenrod1",
    "Intergenic" = "darksalmon",
    "Genic\nIntron" = "#41B6C4"
  )
  
  ############################## 
  ######### UPSET PLOT ######### 
  ############################## 
  
  intersection_size_annotation <- ComplexUpset::intersection_size(
      counts = FALSE,
      mapping = aes(fill = structural_category)
    ) +
    scale_fill_manual(
      values = category_colors,
      name = "Structural Category"
    ) +
    ggtitle(method) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = paper_axis_text_y(),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5, size = .paper_font("panel"), face = "bold"
      )
    )
  
  rev_sample_labels <- rev(sample_labels)
  
  upset_plot <- (
    ComplexUpset::upset(
      listUJC,
      intersect = rev_sample_labels,
      sort_sets = FALSE,
      name = element_blank(),
      base_annotations = list(
        'Intersection size' = intersection_size_annotation
      ),
      # ComplexUpset dot/segment sizes are absolute (final print units); tune via
      # PAPER_UPSET_DOT_SIZE / PAPER_UPSET_SEGMENT_SIZE so the 7-set matrix reads
      # correctly at 180 mm instead of merging into blobs.
      matrix = ComplexUpset::intersection_matrix(
        geom = ggplot2::geom_point(size = PAPER_UPSET_DOT_SIZE),
        segment = ggplot2::geom_segment(linewidth = PAPER_UPSET_SEGMENT_SIZE)
      ),
      n_intersections = n,
      width_ratio=0.1,
      set_sizes=FALSE,
      themes = paper_upset_modify_themes(matrix_label_size = upset_matrix_label_size)
    )
  )
  
  #####################################
  ######### COMBINATION PLOTs ######### 
  #####################################

  # UpSet uses all `sample_labels`. Combination / curve / stack plots use a subset of
  # listUJC columns: default = names ending in a digit (mouse B1–B5 / Brain 1–5 only).
  # Pass `combination_columns` explicitly for layouts like SY5Y (J&C, C&J, …).
  if (is.null(combination_columns)) {
    sample_cols <- names(listUJC)[str_detect(names(listUJC), "\\d+$")]
    if (length(sample_cols) != 5) {
      stop(
        "Default combination logic expects 5 sample columns whose names end in a digit ",
        "(replicates only). Pass combination_columns explicitly for other layouts."
      )
    }
    prefix <- sample_cols %>%
      str_remove("\\d+$") %>%
      unique()
    if (length(prefix) != 1) {
      stop("Failed to detect a single common prefix among the sample columns.")
    }
  } else {
    miss <- setdiff(combination_columns, names(listUJC))
    if (length(miss)) {
      stop("combination_columns not found in UJC table: ", paste(miss, collapse = ", "))
    }
    sample_cols <- combination_columns
    prefix <- "samples"
  }

  n_sc <- length(sample_cols)
  n_sc_m1 <- max(1L, n_sc - 1L)

  comb_counts <- listUJC %>%
    unite("pattern", all_of(sample_cols), sep = "") %>%
    filter(pattern != str_dup("0", length(sample_cols))) %>%
    group_by(pattern) %>%
    summarise(n_UJC = n(), .groups = "drop") %>%
    mutate(n_samples = str_count(pattern, "1")) %>%
    mutate(percent_UJC = (n_UJC / sum(n_UJC)) * 100)
  
  comb_plot <- ggplot(comb_counts, aes(x = factor(n_samples), y = n_UJC)) +
    geom_boxplot(
      data = comb_counts %>% filter(n_samples %in% seq_len(n_sc_m1)),
      aes(group = factor(n_samples)),
      width = 0.25,
      outlier.shape = NA
    ) +
    geom_point(
      data = comb_counts %>% filter(n_samples == n_sc),
      aes(group = pattern),
      position = position_dodge(width = 0.8),
      size = 2
    ) +
    scale_x_discrete(drop = FALSE) +
    labs(
      x     = "Number of samples present",
      y     = "Number of UJCs",
      title = paste0("Distribution of UJCs across ", prefix, " (n=", n_sc, ")")
    ) +
    paper_theme()
  comb_plot$full_data <- comb_counts

  perc_comb_plot <- NULL
  if (isTRUE(include_percentage_plots)) {
    perc_comb_plot <- ggplot(comb_counts, aes(x = factor(n_samples), y = percent_UJC)) +
      geom_boxplot(
        data = comb_counts %>% filter(n_samples %in% seq_len(n_sc_m1)),
        aes(group = factor(n_samples)),
        width = 0.25,
        outlier.shape = NA
      ) +
      geom_point(
        data = comb_counts %>% filter(n_samples == n_sc),
        aes(group = pattern),
        position = position_dodge(width = 0.8),
        size = 2
      ) +
      scale_x_discrete(drop = FALSE) +
      labs(
        x     = "Number of samples present",
        y     = "Percent of total UJCs",
        title = paste0("Distribution of UJCs across ", prefix, " (n=", n_sc, ")")
      ) +
      paper_theme()
  }


  #########################################
  ######### COMBINATION BAR PLOTs ######### 
  ######################################### 
  
  comb_counts_sc <- listUJC %>%
    unite("pattern", all_of(sample_cols), sep = "") %>%
    filter(pattern != str_dup("0", length(sample_cols))) %>%
    group_by(pattern, structural_category) %>%
    summarise(n_UJC = n(), .groups = "drop") %>%
    mutate(n_samples = str_count(pattern, "1")) %>%
    mutate(percent_UJC = (n_UJC / sum(n_UJC)) * 100)
  
  comb_bar_data <- comb_counts_sc %>%
    group_by(n_samples, structural_category) %>%
    summarise(
      sum_n_UJC = sum(n_UJC, na.rm = TRUE),
      sum_percent_UJC = sum(percent_UJC, na.rm = TRUE),
      .groups = "drop"
    )
  
  comb_bar_plot <- ggplot(comb_bar_data, aes(x = factor(n_samples), y = sum_n_UJC, fill = factor(structural_category))) +
    geom_col(width = 0.6) +
    labs(
      x     = "Number of samples present",
      y     = "# UJCs",
      title = paste0("UJC distribution by structural category (", prefix, ", n=", n_sc, ")")
    ) +
    paper_theme() +
    scale_fill_manual(
      values = category_colors,
      name = "Structural Category"
    )

  perc_comb_bar_plot <- NULL
  if (isTRUE(include_percentage_plots)) {
    perc_comb_bar_plot <- ggplot(comb_bar_data, aes(x = factor(n_samples), y = sum_percent_UJC, fill = factor(structural_category))) +
      geom_col(width = 0.6) +
      labs(
        x     = "Number of samples present",
        y     = "Percent of total UJCs",
        title = paste0("Average UJC distribution across ", prefix, " samples")
      ) +
      paper_theme() +
      scale_fill_manual(
        values = category_colors,
        name = "Structural Category"
      )
  }

  ############################################
  ######### COMBINATION BAR FL PLOTs ######### 
  ############################################
  
  fl_sums <- bind_rows(class_df_list, .id = "sample") %>%
    filter(!grepl("^NA", UJC)) %>%
    group_by(UJC, sample, structural_category) %>%
    summarise(FL_sum = sum(FL, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from   = sample,
      values_from  = FL_sum,
      names_prefix = "FL_"
    )
  drop_fl <- intersect(c("FL_Join&Call", "FL_Call&Join"), names(fl_sums))
  if (length(drop_fl)) {
    fl_sums <- fl_sums %>% select(-all_of(drop_fl))
  }

  cdf_names <- names(class_df_list)
  idx_map <- match(sample_cols, sample_labels)
  if (any(is.na(idx_map))) {
    stop("combination_columns must all appear in sample_labels (same order as class_df_list).")
  }
  fl_keep <- paste0("FL_", cdf_names[idx_map])
  miss_fl <- setdiff(fl_keep, names(fl_sums))
  if (length(miss_fl)) {
    stop("Missing FL columns for combination set: ", paste(miss_fl, collapse = ", "))
  }
  fl_sums <- fl_sums %>%
    select(UJC, structural_category, all_of(fl_keep)) %>%
    mutate(across(all_of(fl_keep), ~ replace_na(.x, 0)))
  
  sample_cols_fl <- grep("^FL", names(fl_sums), value = TRUE)
  
  fl_sums <- fl_sums %>%
    mutate(
      samples_present = rowSums(across(all_of(sample_cols_fl), ~ .x != 0))
    )
  
  fl_comb <- fl_sums %>%
    filter(samples_present != 0) %>%
    mutate(
      row_total = rowSums(across(all_of(sample_cols_fl)), na.rm = TRUE)
    ) %>%
    group_by(samples_present, structural_category) %>%
    summarise(
      FL_sum = sum(row_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      FL_percent = (FL_sum / sum(FL_sum)) * 100
    )
  
  comb_fl_bar_plot <- ggplot(fl_comb, aes(x = factor(samples_present), y = FL_sum, fill = factor(structural_category))) +
    geom_col(width = 0.6) +
    labs(
      x     = "Number of samples present",
      y     = PAPER_UJC_YLAB$comb_fl_bar,
      title = paste0("Total read distribution of UJCs across ", prefix, " samples")
    ) +
    paper_theme() +
    scale_fill_manual(
      values = category_colors,
      name = "Structural Category"
    )
  
  perc_comb_FL_bar_plot <- NULL
  if (isTRUE(include_percentage_plots)) {
    perc_comb_FL_bar_plot <- ggplot(fl_comb, aes(x = factor(samples_present), y = FL_percent, fill = factor(structural_category))) +
      geom_col(width = 0.6) +
      labs(
        x     = "Number of samples present",
        y     = PAPER_UJC_YLAB$perc_comb_fl_bar,
        title = paste0("Relative read distribution of UJCs across ", prefix, " samples")
      ) +
      paper_theme() +
      scale_fill_manual(
        values = category_colors,
        name = "Structural Category"
      )
  }


  #########################################
  ######### DISCOVERY CURVE PLOTs ######### 
  #########################################
  
  discovery_curve <- map_dfr(1:length(sample_cols), function(n) {
    combos <- combn(sample_cols, n, simplify = FALSE)
    
    ujc_counts <- map_int(combos, function(sample_set) {
      listUJC %>%
        select(all_of(sample_set)) %>%
        transmute(present = rowSums(.) > 0) %>%
        summarise(n_UJC = sum(present)) %>%
        pull(n_UJC)
    })
    
    tibble(
      n_samples = n,
      mean_UJCs = mean(ujc_counts),
      min_UJCs  = min(ujc_counts),
      max_UJCs  = max(ujc_counts)
    )
  })
  
  max_total <- discovery_curve %>% filter(n_samples == max(n_samples)) %>% pull(max_UJCs)

  discovery_curve <- discovery_curve %>%
    mutate(
      mean_pct = 100 * mean_UJCs / max_total,
      min_pct  = 100 * min_UJCs  / max_total,
      max_pct  = 100 * max_UJCs  / max_total
    )
  
  cumul_ujc_curve_plot <- ggplot(discovery_curve, aes(x = n_samples, y = mean_UJCs)) +
    geom_ribbon(aes(ymin = min_UJCs, ymax = max_UJCs), alpha = 0.2) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = seq_len(n_sc)) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(
      title = "Mean number of UJCs discovered across sample combinations",
      x = "Number of samples included",
      y = "Mean unique UJCs discovered"
    ) +
    paper_theme()

  perc_cumul_ujc_curve_plot <- NULL
  if (isTRUE(include_percentage_plots)) {
    perc_cumul_ujc_curve_plot <- ggplot(discovery_curve, aes(x = n_samples, y = mean_pct)) +
      geom_ribbon(aes(ymin = min_pct, ymax = max_pct), alpha = 0.2) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_x_continuous(breaks = seq_len(n_sc)) +
      scale_y_continuous(limits = c(0, NA)) +
      labs(
        title = "Mean percentage of UJCs discovered across sample combinations",
        x = "Number of samples included",
        y = "Mean percentage UJCs discovered"
      ) +
      paper_theme()
  }


  #################################################
  ######### COMBINATION STACKED BAR PLOTs ######### 
  #################################################
  
  stack_data <- map_dfr(1:length(sample_cols), function(k) {
    combos <- combn(sample_cols, k, simplify = FALSE)

    bind_rows(lapply(combos, function(sample_set) {
      ujc_data <- listUJC %>%
        mutate(n_present = rowSums(across(all_of(sample_set)))) %>%
        filter(n_present > 0) %>%
        count(n_present) %>%
        complete(n_present = 1:k, fill = list(n = 0)) %>%
        mutate(n_samples = k)
      
      return(ujc_data)
    }))
  }) %>%
    group_by(n_samples, n_present) %>%
    summarise(
      mean_count = mean(n),
      sd_count = sd(n),
      .groups = "drop"
    ) %>%
    mutate(
      percent = 100 * mean_count / max_total,
      sd_percent = 100 * sd_count / max_total
    )
  
  cumul_ujc_stack_plot <- ggplot(stack_data, aes(x = factor(n_samples), y = mean_count, fill = factor(n_present))) +
    geom_col(width = 0.8) +
    scale_fill_conesa(palette = "complete", name = "Samples present in") +
    labs(
      title = "UJC reproducibility across sample combinations",
      x = "Number of samples included",
      y = "Mean number of UJCs",
      fill = "Reproducibility"
    ) +
    paper_theme()

  perc_cumul_ujc_stack_plot <- NULL
  if (isTRUE(include_percentage_plots)) {
    perc_cumul_ujc_stack_plot <- ggplot(stack_data, aes(x = factor(n_samples), y = percent, fill = factor(n_present))) +
      geom_col(width = 0.8) +
      scale_fill_conesa(palette = "complete", name = "Samples present in") +
      labs(
        title = "UJC reproducibility across sample combinations (relative)",
        x = "Number of samples included",
        y = "Mean % of UJCs",
        fill = "Reproducibility"
      ) +
      paper_theme()
  }


  ####################################################
  ######### COMBINATION FL STACKED BAR PLOTs ######### 
  ####################################################
  
  
  sample_combs_fl <- seq_along(sample_cols_fl) %>%
    map(~ combn(sample_cols_fl, .x, simplify = FALSE)) %>%
    flatten()
  
  combos_tbl <- tibble(
    cols = sample_combs_fl
  ) %>% 
    rowwise() %>% 
    mutate(
      df = list(
        fl_sums %>% 
          select(all_of(cols)) %>%
          mutate(
            samples_present = rowSums(across(everything()) > 0)
          ) %>%
          group_by(samples_present) %>%
          summarise(
            across(
              .cols = all_of(cols),
              .fns  = ~sum(.x, na.rm = TRUE)
            ),
            .groups = "drop"
          ) %>%
          mutate(
            FL_total = rowSums(across(all_of(cols)))
          ) %>%
          filter(samples_present > 0) %>%
          select(samples_present, FL_total)
      )
    ) %>% 
    ungroup()
  
  combos_fl <- combos_tbl %>%
    mutate(nr_samples = map_int(cols, length)) %>%
    unnest(df) %>%
    group_by(nr_samples, samples_present) %>%
    summarise(
      min_FL_total = min(FL_total),
      mean_FL_total = mean(FL_total),
      sd_FL_total = sd(FL_total),
      max_FL_total = max(FL_total),
      .groups  = "drop"
    ) %>%
    arrange(nr_samples, samples_present)
  
  fl_max <- combos_fl %>%
    filter(nr_samples == max(nr_samples)) %>%
    summarise(sum_FL_total = sum(max_FL_total)) %>%
    pull(sum_FL_total)

  combos_fl_pct <- NULL
  if (isTRUE(include_percentage_plots)) {
    combos_fl_pct <- combos_fl %>%
      mutate(
        mean_FL_pct = max_FL_total / fl_max * 100,
        sd_FL_pct = sd_FL_total / fl_max * 100
      )
  }


  cumul_ujc_fl_stack_plot <- ggplot(combos_fl,
                              aes(x = factor(nr_samples),
                                  y = mean_FL_total,
                                  fill = factor(samples_present))) +
    geom_col(width = 0.8) +
    scale_fill_conesa(palette = "complete", name = "Samples present in") +
    labs(
      title = "Mean FL total across sample combinations",
      x     = "Number of samples included",
      y     = PAPER_UJC_YLAB$ujc_fl_stack
    ) +
    paper_theme()

  perc_cumul_ujc_fl_stack_plot <- NULL
  if (isTRUE(include_percentage_plots)) {
    perc_cumul_ujc_fl_stack_plot <- ggplot(combos_fl_pct,
                                 aes(x = factor(nr_samples),
                                     y = mean_FL_pct,
                                     fill = factor(samples_present))) +
      geom_col(width = 0.8) +
      scale_fill_conesa(palette = "complete", name = "Samples present in") +
      labs(
        title = "FL total as % of max across sample combinations",
        x     = "Number of samples included",
        y     = PAPER_UJC_YLAB$perc_ujc_fl_stack
      ) +
      paper_theme()
  }




  ###################################
  ######### RETURNING PLOTS ######### 
  ###################################
  
  
  upset_plots <- list()
  upset_plots[["upset"]] <- upset_plot
  upset_plots[["comb"]] <- comb_plot
  upset_plots[["perc_comb"]] <- perc_comb_plot
  upset_plots[["comb_bar"]] <- comb_bar_plot
  upset_plots[["perc_comb_bar"]] <- perc_comb_bar_plot
  upset_plots[["comb_fl_bar"]] <- comb_fl_bar_plot
  upset_plots[["perc_comb_fl_bar"]] <- perc_comb_FL_bar_plot
  upset_plots[["ujc_curve"]] <- cumul_ujc_curve_plot
  upset_plots[["perc_ujc_curve"]] <- perc_cumul_ujc_curve_plot
  upset_plots[["ujc_stack"]] <- cumul_ujc_stack_plot
  upset_plots[["perc_ujc_stack"]] <- perc_cumul_ujc_stack_plot
  upset_plots[["ujc_fl_stack"]] <- cumul_ujc_fl_stack_plot
  upset_plots[["perc_ujc_fl_stack"]] <- perc_cumul_ujc_fl_stack_plot
  
  upset_plots[["ujc_stack_data"]] <- stack_data
  upset_plots[["ujc_fl_stack_data"]] <- combos_fl
  upset_plots[["ujc_fl_stack_pct_data"]] <- combos_fl_pct
  # UJC presence/absence membership matrix (Source Data for the UpSet figures).
  upset_plots[["membership_data"]] <- listUJC

  return(upset_plots)
}


compare_isoform_plots <- function(class_df_list, include_percentage_plots = FALSE,
                                  condition_label = "brain") {

  df <- bind_rows(class_df_list, .id = "source") %>%
    mutate(source = factor(source, levels = c("Join&Call", "Call&Join"))) %>%
    filter(!grepl("NA", UJC)) %>%
    select(source, isoform, chrom, associated_gene,
           associated_transcript, UJC, structural_category) %>%
    filter(source %in% c("Join&Call", "Call&Join"))

  iso_per_transcript <- df %>%
    filter(associated_transcript != "novel") %>%
    group_by(source, associated_transcript) %>%
    summarise(n_isoforms = n_distinct(isoform), .groups = "drop") %>%
    mutate(
      iso_cat = case_when(
        n_isoforms == 1 ~ "1",
        n_isoforms == 2 ~ "2",
        n_isoforms == 3 ~ "3",
        n_isoforms == 4 ~ "4",
        n_isoforms >= 5 ~ "5+"
      ),
      iso_cat = factor(iso_cat, levels = c("1","2","3","4","5+"))
    ) %>%
      group_by(source, iso_cat) %>%
      summarise(
        n_assoc_transcript = n_distinct(associated_transcript),
        sum_n_isoforms     = sum(n_isoforms),
        .groups = "drop"
      ) %>%
        group_by(source) %>%
        mutate(
          total_assoc = sum(n_assoc_transcript),
          total_iso   = sum(sum_n_isoforms)
        ) %>%
        ungroup() %>%
        mutate(
          perc_assoc_transcript = n_assoc_transcript / max(total_assoc),
          perc_n_isoforms       = sum_n_isoforms   / max(total_iso)
        ) %>%
        select(-total_assoc, -total_iso)

  p_count_transcript <- iso_per_transcript %>%
    ggplot(aes(x = source, y = n_assoc_transcript, fill = iso_cat)) +
      geom_col() +
      scale_fill_conesa(palette = "complete", name = "# of isoforms per transcript") +
      labs(
        title = paste0("Redundancy of transcripts; ", condition_label),
        x     = "Strategy",
        y     = "# of associated transcripts"
      ) +
      paper_theme()

  p_perc_transcript <- NULL
  if (isTRUE(include_percentage_plots)) {
    p_perc_transcript <- iso_per_transcript %>%
      ggplot(aes(x = source, y = perc_assoc_transcript * 100, fill = iso_cat)) +
        geom_col() +
        scale_fill_conesa(palette = "complete", name = "# of isoforms per transcript") +
        labs(
          title = paste0("Relative redundancy of transcripts; ", condition_label),
          x     = "Strategy",
          y     = "% of associated transcripts"
        ) +
        paper_theme()
  }

  p_iso_transcript <- iso_per_transcript %>%
    ggplot(aes(x = source, y = sum_n_isoforms, fill = iso_cat)) +
      geom_bar(stat = "identity") +
      scale_fill_conesa(palette = "complete", name = "# of isoforms per transcript") +
      labs(
        title = paste0("Redundancy of transcripts; ", condition_label),
        x     = "Strategy",
        y     = "# of isoforms"
      ) +
      paper_theme()

  p_perc_iso_transcript <- NULL
  if (isTRUE(include_percentage_plots)) {
    p_perc_iso_transcript <- iso_per_transcript %>%
      ggplot(aes(x = source, y = perc_n_isoforms * 100, fill = iso_cat)) +
        geom_bar(stat = "identity") +
        scale_fill_conesa(palette = "complete", name = "# of isoforms per transcript") +
        labs(
          title = paste0("Relative redundancy of transcripts; ", condition_label),
          x     = "Strategy",
          y     = "% of isoforms"
        ) +
        paper_theme()
  }

  iso_per_ujc <- df %>%
    group_by(source, UJC) %>%
    summarise(n_isoforms = n_distinct(isoform), .groups = "drop") %>%
    mutate(
      iso_cat = case_when(
        n_isoforms == 1 ~ "1",
        n_isoforms == 2 ~ "2",
        n_isoforms == 3 ~ "3",
        n_isoforms == 4 ~ "4",
        n_isoforms >= 5 ~ "5+"
      ),
      iso_cat = factor(iso_cat, levels = c("1","2","3","4","5+"))
    ) %>%
    group_by(source, iso_cat) %>%
    summarise(
      n_UJC           = n_distinct(UJC),
      sum_n_isoforms  = sum(n_isoforms),
      .groups = "drop"
    ) %>%
    group_by(source) %>%
    mutate(
      total_UJC       = sum(n_UJC),
      total_isoforms  = sum(sum_n_isoforms)
    ) %>%
    ungroup() %>%
    mutate(
      perc_UJC          = n_UJC          / max(total_UJC),
      perc_n_isoforms   = sum_n_isoforms / max(total_isoforms)
    ) %>%
    select(-total_UJC, -total_isoforms)


  p_count_ujc <- iso_per_ujc %>%
    ggplot(aes(x = source, y = n_UJC, fill = iso_cat)) +
      scale_fill_conesa(palette = "complete", name = "# of isoforms per UJC") +
      geom_bar(stat = "identity") +
      labs(
        title = paste0("Redundancy of UJCs; ", condition_label),
        x     = "Strategy",
        y     = "# of UJCs"
      ) +
      paper_theme()
  
  p_perc_ujc <- NULL
  if (isTRUE(include_percentage_plots)) {
    p_perc_ujc <- iso_per_ujc %>%
      ggplot(aes(x = source, y = perc_UJC * 100, fill = iso_cat)) +
        scale_fill_conesa(palette = "complete", name = "# of isoforms per UJC") +
        geom_bar(stat = "identity") +
        labs(
          title = paste0("Redundancy of UJCs; ", condition_label),
          x     = "Strategy",
          y     = "% of UJCs"
        ) +
        paper_theme()
  }

  p_iso_ujc <- iso_per_ujc %>%
    ggplot(aes(x = source, y = sum_n_isoforms, fill = iso_cat)) +
      scale_fill_conesa(palette = "complete", name = "# of isoforms per UJC") +
      geom_bar(stat = "identity") +
      labs(
        title = paste0("Redundancy of UJCs; ", condition_label),
        x     = "Strategy",
        y     = "# of isoforms"
      ) +
      paper_theme()

  p_perc_iso_ujc <- NULL
  if (isTRUE(include_percentage_plots)) {
    p_perc_iso_ujc <- iso_per_ujc %>%
      ggplot(aes(x = source, y = perc_n_isoforms * 100, fill = iso_cat)) +
        scale_fill_conesa(palette = "complete", name = "# of isoforms per UJC") +
        geom_bar(stat = "identity") +
        labs(
          title = paste0("Redundancy of UJCs; ", condition_label),
          x     = "Strategy",
          y     = "% of isoforms"
        ) +
        paper_theme()
  }

  out <- list(
    count_transcript = p_count_transcript,
    total_isoforms_transcript = p_iso_transcript,
    count_ujc = p_count_ujc,
    total_isoforms_ujc = p_iso_ujc
  )
  if (isTRUE(include_percentage_plots)) {
    out <- c(
      out,
      list(
        perc_transcript = p_perc_transcript,
        perc_isoforms_transcript = p_perc_iso_transcript,
        perc_ujc = p_perc_ujc,
        perc_isoforms_ujc = p_perc_iso_ujc
      )
    )
  }
  out
}
