# 08_read_level_sy5y.R
# Read-level plots (SQANTI-reads stacked categories, FASTQ vs primary-mapped counts,
# BAM read-length violins) for ONT R10 SY5Y — mirrors mouse read_figure_plots.Rmd logic.

read_level_sy5y_cat_palette <- c(
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

read_level_sy5y_category_map <- c(
  "full-splice_match" = "FSM",
  "incomplete-splice_match" = "ISM",
  "novel_in_catalog" = "NIC",
  "novel_not_in_catalog" = "NNC",
  "genic" = "Genic\nGenomic",
  "antisense" = "Antisense",
  "fusion" = "Fusion",
  "intergenic" = "Intergenic",
  "genic_intron" = "Genic\nIntron"
)

read_level_sy5y_display_levels <- c(
  "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion",
  "Intergenic", "Genic\nIntron", "Unstranded", "Unaligned"
)

#' One SQANTI-reads classification file -> category counts (memory: structural_category only).
summarize_sqanti_classification_file <- function(path, sample_label, tech_value) {
  if (!file.exists(path)) {
    stop("Missing classification file: ", path, call. = FALSE)
  }
  df <- readr::read_tsv(path, show_col_types = FALSE, col_select = "structural_category")
  df <- dplyr::mutate(df, structural_category = as.character(.data$structural_category))
  result <- df %>%
    dplyr::transmute(
      category_label = dplyr::recode(structural_category, !!!read_level_sy5y_category_map)
    ) %>%
    dplyr::filter(!is.na(category_label)) %>%
    dplyr::count(category_label, name = "num_reads") %>%
    dplyr::mutate(sample = sample_label, technology = tech_value)
  rm(df)
  gc()
  result
}

#' Per-sample SQANTI-reads classification diagnostics (read `structural_category` only).
#'
#' Reports rows in file, how many become `NA` after [read_level_sy5y_category_map] recode
#' (dropped from stacked counts), and comparison to `read_numbers` FASTQ vs primary BAM.
#'
#' @return A tibble (one row per sample) suitable for `knitr::kable()`.
read_level_sy5y_sqanti_diagnostics <- function(read_numbers,
                                               sample_ids,
                                               sample_labels,
                                               sqanti_reads_root) {
  if (length(sample_ids) != length(sample_labels)) {
    stop("`sample_ids` and `sample_labels` must have the same length.", call. = FALSE)
  }

  purrr::map_dfr(seq_along(sample_ids), function(i) {
    id <- sample_ids[[i]]
    lab <- sample_labels[[i]]
    path <- file.path(sqanti_reads_root, id, paste0(id, "_reads_classification.txt"))
    if (!file.exists(path)) {
      return(tibble::tibble(
        sample = id,
        sample_label = lab,
        error = paste0("missing file: ", path)
      ))
    }

    sc <- as.character(
      readr::read_tsv(
        path,
        show_col_types = FALSE,
        col_select = "structural_category"
      )[[1L]]
    )

    n_row <- length(sc)
    cat_lab <- dplyr::recode(sc, !!!read_level_sy5y_category_map)
    n_na <- sum(is.na(cat_lab))
    n_used <- sum(!is.na(cat_lab))

    unk <- sc[is.na(cat_lab)]
    top_u <- if (length(unk)) {
      tb <- sort(table(unk), decreasing = TRUE)
      k <- min(5L, length(tb))
      idx <- seq_len(k)
      paste0(names(tb)[idx], " (n=", unname(tb)[idx], ")", collapse = "; ")
    } else {
      NA_character_
    }

    rn1 <- read_numbers[read_numbers$sample == id, , drop = FALSE]
    if (!nrow(rn1)) {
      stop("No read_numbers row for sample: ", id, call. = FALSE)
    }

    fq <- as.numeric(rn1$ont_fastq[[1L]])
    pr <- as.numeric(rn1$ont_prim_aln[[1L]])

    assigned_used <- n_used
    unmapped <- as.integer(pmax(0, fq - pr))
    unstranded <- as.integer(pmax(0, pr - assigned_used))

    tibble::tibble(
      sample = id,
      sample_label = lab,
      n_rows_in_classification_file = n_row,
      n_rows_na_after_recode = as.integer(n_na),
      pct_rows_na_after_recode = round(100 * n_na / pmax(1L, n_row), 3),
      n_rows_used_in_category_counts = assigned_used,
      ont_fastq = fq,
      ont_prim_aln = pr,
      n_unmapped_fastq = unmapped,
      n_unstranded_primary = unstranded,
      n_primary_mapped_minus_rows = pr - n_row,
      pct_sqanti_rows_of_primary_mapped = round(100 * n_row / pmax(1, pr), 2),
      n_fastq_minus_rows = fq - n_row,
      top_unknown_structural_categories = top_u
    )
  })
}

#' Per-replicate SQANTI-reads row counts used in stacked categories (`assigned`).
#'
#' @return Named numeric vector (`names` = `sample_ids`).
sy5y_sqanti_assigned_per_sample <- function(sample_ids, sqanti_reads_root) {
  stats::setNames(
    purrr::map_dbl(sample_ids, function(id) {
      path <- file.path(
        sqanti_reads_root,
        id,
        paste0(id, "_reads_classification.txt")
      )
      if (!file.exists(path)) {
        stop("Missing classification file: ", path, call. = FALSE)
      }
      sc <- as.character(
        readr::read_tsv(
          path,
          show_col_types = FALSE,
          col_select = "structural_category"
        )[[1L]]
      )
      sum(!is.na(dplyr::recode(sc, !!!read_level_sy5y_category_map)))
    }),
    sample_ids
  )
}

#' SQANTI structural categories + Unstranded + Unaligned (ONT only).
#'
#' With **minimap2 `-ub`** (Dorado-trimmed, not strand-oriented), decompose reads not
#' in the SQANTI-reads table into:
#'
#' - **Unaligned**: `ont_fastq - ont_prim_aln` (no primary alignment).
#' - **Unstranded**: `ont_prim_aln - assigned` (primary-mapped but no classification row;
#'   BAM-derived SQANTI input often yields fewer rows than primary alignments).
#'
#' SQANTI categories (including Antisense) sum to `assigned`.
build_sy5y_sqanti_category_counts <- function(read_numbers,
                                              sample_ids,
                                              sample_labels,
                                              sqanti_reads_root) {
  if (length(sample_ids) != length(sample_labels)) {
    stop("`sample_ids` and `sample_labels` must have the same length.", call. = FALSE)
  }

  paths <- file.path(
    sqanti_reads_root,
    sample_ids,
    paste0(sample_ids, "_reads_classification.txt")
  )

  counts <- purrr::pmap_dfr(
    list(paths, sample_labels),
    function(p, lab) summarize_sqanti_classification_file(p, lab, "ont")
  )

  rn <- read_numbers %>%
    dplyr::filter(.data$sample %in% sample_ids) %>%
    dplyr::mutate(
      sample_label = sample_labels[match(.data$sample, sample_ids)]
    )

  assigned <- counts %>%
    dplyr::group_by(.data$sample, .data$technology) %>%
    dplyr::summarise(assigned = sum(.data$num_reads), .groups = "drop")

  extra_base <- rn %>%
    dplyr::transmute(
      sample = .data$sample_label,
      technology = "ont",
      ont_fastq = .data$ont_fastq,
      ont_prim_aln = .data$ont_prim_aln
    ) %>%
    dplyr::left_join(assigned, by = c("sample", "technology")) %>%
    dplyr::mutate(assigned = dplyr::coalesce(.data$assigned, 0L))

  extra <- dplyr::bind_rows(
    extra_base %>%
      dplyr::transmute(
        category_label = "Unstranded",
        num_reads = as.integer(pmax(0L, .data$ont_prim_aln - .data$assigned)),
        sample = .data$sample,
        technology = .data$technology
      ),
    extra_base %>%
      dplyr::transmute(
        category_label = "Unaligned",
        num_reads = as.integer(pmax(0L, .data$ont_fastq - .data$ont_prim_aln)),
        sample = .data$sample,
        technology = .data$technology
      )
  )

  dplyr::bind_rows(counts, extra) %>%
    dplyr::mutate(
      category_label = factor(
        .data$category_label,
        levels = read_level_sy5y_display_levels
      ),
      sample = factor(.data$sample, levels = sample_labels),
      technology = factor(.data$technology, levels = c("pacbio", "ont"))
    ) %>%
    dplyr::arrange(dplyr::desc(.data$num_reads))
}

#' Same layout as mouse `plot_sqanti_faceted` in read_figure_plots.Rmd.
plot_read_level_sqanti_faceted <- function(counts_df, title = NULL) {
  pal <- read_level_sy5y_cat_palette
  extended_fill <- c(
    pal,
    "Unstranded" = "grey85",
    "Unaligned" = "white"
  )
  dl <- read_level_sy5y_display_levels
  legend_breaks <- c("Unaligned", "Unstranded", dl[!dl %in% c("Unaligned", "Unstranded")])
  legend_labels <- c(
    "Unaligned",
    "Unstranded",
    "Full\nSplice Match",
    "Incomplete\nSplice Match",
    "Novel\nIn Catalog",
    "Novel Not\nIn Catalog",
    "Genic\nGenomic",
    "Antisense",
    "Fusion",
    "Intergenic",
    "Genic\nIntron"
  )

  ggplot2::ggplot(
    counts_df,
    ggplot2::aes(
      x = .data$technology,
      y = .data$num_reads,
      fill = .data$category_label
    )
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      width = 0.8,
      position = ggplot2::position_stack(reverse = TRUE),
      colour = "black",
      linewidth = 0.15
    ) +
    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::scale_fill_manual(
      name = "Structural Category",
      values = extended_fill,
      breaks = legend_breaks,
      labels = legend_labels,
      drop = FALSE
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$sample), switch = "x") +
    paper_read_count_y_scale() +
    ggplot2::labs(x = NULL, y = "Number of reads", title = title) +
    paper_read_level_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right",
      legend.box = "vertical"
    )
}

#' Tidy table for stacked SQANTI-classified / unstranded / unmapped reads (one facet per replicate).
#'
#' Requires `sqanti_reads_root` to sum per-read classifications, or pass precomputed
#' `category_counts` from [build_sy5y_sqanti_category_counts].
build_sy5y_readnum_tidy <- function(read_numbers,
                                    sample_ids,
                                    sample_labels,
                                    sqanti_reads_root = NULL,
                                    category_counts = NULL) {
  req <- c("sample", "ont_fastq", "ont_prim_aln", "pb_fastq", "pb_prim_aln")
  if (!all(req %in% names(read_numbers))) {
    stop(
      "read_numbers must contain: ",
      paste(req, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(category_counts)) {
    if (is.null(sqanti_reads_root)) {
      stop("Provide `category_counts` or `sqanti_reads_root`.", call. = FALSE)
    }
    category_counts <- build_sy5y_sqanti_category_counts(
      read_numbers,
      sample_ids,
      sample_labels,
      sqanti_reads_root
    )
  }

  assigned <- category_counts %>%
    dplyr::filter(
      !.data$category_label %in% c("Unstranded", "Unaligned"),
      .data$technology == "ont"
    ) %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::summarise(classified = sum(.data$num_reads), .groups = "drop")

  extra <- category_counts %>%
    dplyr::filter(
      .data$category_label %in% c("Unstranded", "Unaligned"),
      .data$technology == "ont"
    ) %>%
    dplyr::transmute(
      sample_label = as.character(.data$sample),
      assign_category = as.character(.data$category_label),
      num_reads = .data$num_reads
    )

  classified <- assigned %>%
    dplyr::transmute(
      sample_label = as.character(.data$sample),
      assign_category = "Classified",
      num_reads = .data$classified
    )

  dplyr::bind_rows(classified, extra) %>%
    dplyr::mutate(
      sample_label = factor(.data$sample_label, levels = sample_labels),
      technology = factor("ont", levels = c("pacbio", "ont")),
      assign_category = forcats::fct_relevel(
        .data$assign_category,
        "Classified",
        "Unstranded",
        "Unaligned"
      )
    )
}

#' Faceted read-number plot (ONT only; matches mouse styling).
plot_sy5y_readnum_faceted <- function(tidy_counts, title = "SY5Y") {
  classified_fill <- RColorConesa::colorConesa(n = 2L, palette = "complete")[1L]
  assign_palette <- c(
    "Classified" = classified_fill,
    "Unstranded" = "grey85",
    "Unaligned" = "white"
  )

  ggplot2::ggplot(
    tidy_counts,
    ggplot2::aes(
      x = .data$technology,
      y = .data$num_reads,
      fill = .data$assign_category,
      order = dplyr::case_when(
        .data$assign_category == "Unaligned" ~ 3L,
        .data$assign_category == "Unstranded" ~ 2L,
        TRUE ~ 1L
      )
    )
  ) +
    ggplot2::geom_col(
      width = 0.7,
      position = ggplot2::position_stack(reverse = TRUE),
      color = "black",
      linewidth = 0.2
    ) +
    paper_read_count_y_scale() +
    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::scale_fill_manual(
      values = assign_palette,
      breaks = c("Unaligned", "Unstranded", "Classified"),
      labels = c("Unaligned", "Unstranded", "Classified")
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$sample_label), switch = "x") +
    ggplot2::labs(x = NULL, y = "Number of reads", fill = "Category", title = title) +
    paper_read_level_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}

read_lengths_sy5y_file <- function(file_path, sample_label, technology) {
  if (!file.exists(file_path)) {
    stop("Missing read-length file: ", file_path, call. = FALSE)
  }
  vals <- readr::read_tsv(file_path, col_names = FALSE, show_col_types = FALSE)[[1L]]
  tibble::tibble(
    length = as.numeric(vals),
    sample = sample_label,
    technology = technology
  )
}

#' Load `*.bam.readlen.txt` files (one length per primary alignment) under `.../lengths_seq_only/ont/`.
build_sy5y_lengths_df <- function(length_ont_dir, sample_ids, sample_labels) {
  paths <- file.path(length_ont_dir, paste0(sample_ids, ".bam.readlen.txt"))
  out <- purrr::map2_dfr(
    paths,
    sample_labels,
    ~ read_lengths_sy5y_file(.x, .y, "ont")
  )
  out %>%
    dplyr::mutate(
      sample = factor(.data$sample, levels = sample_labels),
      technology = factor(.data$technology, levels = c("pacbio", "ont"))
    )
}

#' Violin + boxplot per replicate (mouse read_figure_plots.Rmd style).
plot_sy5y_lengths_violin <- function(lengths_df, title = "SY5Y") {
  df <- lengths_df %>%
    dplyr::group_by(.data$sample, .data$technology) %>%
    dplyr::filter(
      .data$length <= stats::quantile(.data$length, 0.995, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  violin_fill <- RColorConesa::colorConesa(n = 3, palette = "complete")[3L]
  y_max <- max(df$length, na.rm = TRUE)
  breaks_1k <- seq(0, ceiling(y_max / 1000) * 1000, by = 1000)

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$technology,
    y = .data$length,
    fill = .data$technology
  )) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$sample), switch = "x") +
    ggplot2::scale_y_continuous(breaks = breaks_1k) +
    ggplot2::scale_fill_manual(values = c(pacbio = violin_fill, ont = violin_fill)) +
    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::labs(x = NULL, y = "Read length (nt)", fill = "Technology", title = title) +
    paper_read_level_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none"
    )
}
