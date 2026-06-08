# 08_read_level_sy5y.R
# ONT R10 SY5Y read-level data loaders (plots in 08_read_level_plots.R).
# Source 08_read_level_plots.R before this file.

read_level_sy5y_cat_palette <- read_level_structural_cat_palette
read_level_sy5y_category_map <- read_level_structural_category_map

read_level_sy5y_display_levels <- c(
  "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion",
  "Intergenic", "Genic\nIntron", "Unstranded", "Unaligned"
)

#' Column spec: read only `structural_category`, skip all other fields (multi-GB TSV safe).
sqanti_classification_col_spec <- function(path) {
  hdr <- readr::read_lines(path, n_max = 1L)
  cols <- strsplit(hdr, "\t", fixed = TRUE)[[1L]]
  if (!"structural_category" %in% cols) {
    stop(
      "Column structural_category not found in ",
      path,
      call. = FALSE
    )
  }
  spec <- stats::setNames(
    rep(list(readr::col_skip()), length(cols)),
    cols
  )
  spec[["structural_category"]] <- readr::col_character()
  do.call(readr::cols, spec)
}

#' Per-category row counts from a SQANTI-reads classification TSV (chunked; safe for multi-GB files).
sqanti_structural_category_tally <- function(path) {
  if (!file.exists(path)) {
    stop("Missing classification file: ", path, call. = FALSE)
  }
  tallies <- new.env(hash = TRUE, parent = emptyenv())
  add_chunk <- function(sc) {
    sc <- as.character(sc)
    tab <- table(sc, useNA = "ifany")
    for (nm in names(tab)) {
      key <- if (length(nm) == 1L && is.na(nm)) {
        "<NA>"
      } else {
        nm
      }
      prev <- tallies[[key]]
      tallies[[key]] <- (if (is.null(prev)) 0L else prev) + as.integer(tab[[nm]])
    }
  }
  cb <- readr::SideEffectChunkCallback$new(function(x, pos) {
    add_chunk(x$structural_category)
  })
  readr::read_tsv_chunked(
    file = path,
    callback = cb,
    chunk_size = 250000L,
    col_types = sqanti_classification_col_spec(path),
    show_col_types = FALSE,
    progress = FALSE
  )
  keys <- ls(tallies, all.names = TRUE)
  if (!length(keys)) {
    return(tibble::tibble(structural_category = character(), n = integer()))
  }
  sc <- keys
  sc[sc == "<NA>"] <- NA_character_
  tibble::tibble(
    structural_category = sc,
    n = vapply(keys, function(k) tallies[[k]], integer(1L))
  )
}

#' One SQANTI-reads classification file -> category counts (uses [sqanti_structural_category_tally]).
summarize_sqanti_classification_file <- function(path, sample_label, tech_value) {
  tab <- sqanti_structural_category_tally(path)
  tab %>%
    dplyr::mutate(
      category_label = dplyr::recode(
        .data$structural_category,
        !!!read_level_sy5y_category_map
      )
    ) %>%
    dplyr::filter(!is.na(.data$category_label)) %>%
    dplyr::transmute(
      category_label = .data$category_label,
      num_reads = .data$n,
      sample = sample_label,
      technology = tech_value
    )
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

    tab <- sqanti_structural_category_tally(path)
    n_row <- sum(tab$n)
    tab <- dplyr::mutate(
      tab,
      category_label = dplyr::recode(
        .data$structural_category,
        !!!read_level_sy5y_category_map
      )
    )
    n_na <- sum(tab$n[is.na(tab$category_label)])
    n_used <- sum(tab$n[!is.na(tab$category_label)])

    unk_tab <- tab[is.na(tab$category_label), , drop = FALSE]
    top_u <- if (nrow(unk_tab)) {
      unk_tab <- unk_tab[order(-unk_tab$n), , drop = FALSE]
      k <- min(5L, nrow(unk_tab))
      paste0(
        unk_tab$structural_category[seq_len(k)],
        " (n=",
        unk_tab$n[seq_len(k)],
        ")",
        collapse = "; "
      )
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
      tab <- sqanti_structural_category_tally(path)
      tab <- dplyr::mutate(
        tab,
        category_label = dplyr::recode(
          .data$structural_category,
          !!!read_level_sy5y_category_map
        )
      )
      sum(tab$n[!is.na(tab$category_label)])
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

#' SQANTI-reads stacked bars for SY5Y (ONT-only facets; includes Unstranded bucket).
plot_read_level_sqanti_faceted <- function(counts_df, title = NULL) {
  plot_read_level_sqanti_stacked(
    counts_df = counts_df,
    display_levels = read_level_sy5y_display_levels,
    title = title,
    include_unstranded = TRUE,
    show_technology_axis = FALSE,
    technology_labels = NULL,
    facet_col = sample
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

#' Faceted read-number plot (ONT only; same styling as mouse read-number stacks).
plot_sy5y_readnum_faceted <- function(tidy_counts, title = "SY5Y") {
  paper_plot_readnum_stacked(
    tidy_counts,
    fill_values = read_level_sy5y_readnum_palette(),
    legend_breaks = c("Unaligned", "Unstranded", "Classified"),
    legend_labels = c("Unaligned", "Unstranded", "Classified"),
    facet_col = sample_label,
    title = title,
    show_technology_axis = FALSE,
    stack_order = dplyr::case_when(
      tidy_counts$assign_category == "Unaligned" ~ 3L,
      tidy_counts$assign_category == "Unstranded" ~ 2L,
      TRUE ~ 1L
    )
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

#' Violin + boxplot per replicate (same as mouse; ONT-only hides technology x labels).
plot_sy5y_lengths_violin <- function(lengths_df, title = "SY5Y") {
  paper_plot_lengths_violin(
    lengths_df,
    title = title,
    show_technology_axis = FALSE
  )
}

#' Side-by-side SQANTI-reads stack + read-length violins (mouse `read_level_combo` layout).
assemble_sy5y_read_level_combined <- function(sq_reads_plot, lengths_plot) {
  combined <- patchwork::wrap_plots(
    L = sq_reads_plot +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5.5, r = 12, b = 5.5, l = 5.5, unit = "pt")
      ),
    R = lengths_plot +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5, l = 12, unit = "pt")
      ),
    design = "LR"
  )
  combined &
    ggplot2::theme(
      legend.title = ggplot2::element_text(
        size = .paper_font("legend_title"),
        face = "bold"
      ),
      legend.text = ggplot2::element_text(size = .paper_font("legend_text"))
    )
}
