# 10_ref_tusco_sirv_violins.R
# Reference / TUSCO / SIRV length + expression violins (PacBio + ONT) for mouse read figures.
# Knit via read_figure_plots.Rmd. Depends on: 00_figure_config.R, 02_themes.R

REF_TUSCO_SIRV_SET_LEVELS <- c("Reference", "TUSCO", "SIRVs")

# Session cache: GTF gene-symbol index is tissue-independent (brain/kidney share it).
.paper_mouse_annotation_cache <- new.env(parent = emptyenv())

#' Conesa palette first color (green), used for all violins in these panels.
paper_ref_tusco_sirv_fill <- function() {
  RColorConesa::colorConesa(n = 2L, palette = "complete")[1L]
}

#' Resolve tissue-specific TUSCO transcript ID list (column 2 of TUSCO TSV).
paper_tusco_transcript_tsv <- function(tissue = c("brain", "kidney")) {
  tissue <- match.arg(tissue)
  fname <- if (tissue == "brain") "tusco_mmu_brain.tsv" else "tusco_mmu_kidney.tsv"
  candidates <- c(
    file.path(paper_tusco_results_dir(), fname),
    file.path(paper_repo_root(), "reports", "tusco_results", fname),
    file.path(paper_repo_root(), "src", "nextflow", "scripts", "tusco", fname)
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    stop("TUSCO transcript list not found: ", fname, call. = FALSE)
  }
  hit
}

paper_strip_transcript_version <- function(x) {
  stringr::str_remove(as.character(x), "\\.\\d+$")
}

#' Full TUSCO tissue table (`tusco_mmu_brain.tsv` / `tusco_mmu_kidney.tsv`).
read_tusco_tissue_table <- function(tissue = c("brain", "kidney")) {
  tissue <- match.arg(tissue)
  readr::read_tsv(
    paper_tusco_transcript_tsv(tissue),
    show_col_types = FALSE,
    comment = "#",
    col_names = c(
      "ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"
    )
  )
}

#' TUSCO gene/transcript sets (same columns as `tusco_ind.R`).
load_tusco_annotation_sets <- function(tissue = c("brain", "kidney")) {
  tbl <- read_tusco_tissue_table(tissue)
  strip_ver <- paper_strip_transcript_version
  list(
    ensembl    = unique(strip_ver(tbl$ensembl)),
    transcript = unique(strip_ver(tbl$transcript)),
    gene_name  = unique(as.character(tbl$gene_name)),
    refseq     = unique(strip_ver(tbl$refseq))
  )
}

#' Extract one GTF attribute value (`gene_name "X"` or `gene "X"`).
extract_gtf_attribute <- function(attribute_field, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regexec(pattern, attribute_field, perl = TRUE)[[1L]]
  if (length(m) >= 2L) {
    return(m[[2L]])
  }
  NA_character_
}

#' Transcript ID and gene symbol from GTF `transcript` rows (or cached TSV).
load_gtf_transcript_gene_index <- function(gtf_path = NULL) {
  tsv_path <- paper_mouse_transcript_gene_index_tsv()
  if (file.exists(tsv_path)) {
    return(
      readr::read_tsv(
        tsv_path,
        col_names = c("transcript_id", "gene_name"),
        show_col_types = FALSE
      ) %>%
        dplyr::mutate(
          transcript_id = paper_strip_transcript_version(.data$transcript_id)
        )
    )
  }
  if (is.null(gtf_path)) {
    gtf_path <- paper_mouse_ref_sirv_gtf()
  }
  gene_idx <- parse_gtf_transcript_gene_index(gtf_path)
  readr::write_tsv(gene_idx, tsv_path, col_names = FALSE)
  message("Wrote transcript gene index: ", tsv_path)
  gene_idx
}

#' Transcript ID and gene symbol from GTF `transcript` rows (one pass).
parse_gtf_transcript_gene_index <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    stop("GTF not found: ", gtf_path, call. = FALSE)
  }
  rows <- list()
  con <- file(gtf_path, open = "r")
  on.exit(close(con), add = TRUE)
  while (length(line <- readLines(con, n = 1L)) > 0L) {
    if (!nzchar(line) || startsWith(line, "#")) {
      next
    }
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 9L || parts[[3L]] != "transcript") {
      next
    }
    attrs <- parts[[9L]]
    tid <- extract_gtf_attribute(attrs, "transcript_id")
    if (is.na(tid) || !nzchar(tid)) {
      next
    }
    gn <- extract_gtf_attribute(attrs, "gene_name")
    if (is.na(gn) || !nzchar(gn)) {
      gn <- extract_gtf_attribute(attrs, "gene")
    }
    if (is.na(gn) || !nzchar(gn)) {
      gn <- extract_gtf_attribute(attrs, "gene_id")
    }
    if (is.na(gn) || !nzchar(gn)) {
      next
    }
    rows[[length(rows) + 1L]] <- list(
      transcript_id = paper_strip_transcript_version(tid),
      gene_name = gn
    )
  }
  if (length(rows) == 0L) {
    stop("No transcript rows parsed from GTF: ", gtf_path, call. = FALSE)
  }
  dplyr::bind_rows(rows)
}

#' Pick one reference transcript per gene symbol (NM > NR > XM > other; then longest).
choose_primary_transcript <- function(transcript_ids, lengths_by_id) {
  if (length(transcript_ids) == 1L) {
    return(transcript_ids[[1L]])
  }
  rank_prefix <- function(tid) {
    if (stringr::str_starts(tid, "NM_")) {
      return(1L)
    }
    if (stringr::str_starts(tid, "NR_")) {
      return(2L)
    }
    if (stringr::str_starts(tid, "XM_")) {
      return(3L)
    }
    4L
  }
  ord <- order(
    vapply(transcript_ids, rank_prefix, integer(1L)),
    -lengths_by_id[transcript_ids],
    transcript_ids
  )
  transcript_ids[[ord[[1L]]]]
}

#' Length + gene-symbol index from length TSV and reference GTF transcript rows.
build_reference_annotation_lookup <- function(
    length_tbl,
    gtf_path = NULL) {
  if (is.null(gtf_path)) {
    gtf_path <- paper_mouse_ref_sirv_gtf()
  }
  tbl <- length_tbl %>%
    dplyr::mutate(
      transcript_id = paper_strip_transcript_version(.data$transcript_id)
    ) %>%
    dplyr::filter(!is.na(.data$length), .data$length > 0L)

  gene_idx <- load_gtf_transcript_gene_index(gtf_path)
  tbl <- tbl %>%
    dplyr::left_join(gene_idx, by = "transcript_id")

  by_id <- stats::setNames(tbl$length, tbl$transcript_id)

  by_gene <- tbl %>%
    dplyr::filter(!is.na(.data$gene_name), nzchar(.data$gene_name)) %>%
    dplyr::group_by(.data$gene_name) %>%
    dplyr::summarise(
      transcript_ids = list(unique(.data$transcript_id)),
      .groups = "drop"
    )

  by_gene_primary <- stats::setNames(
    vapply(
      by_gene$transcript_ids,
      choose_primary_transcript,
      character(1L),
      lengths_by_id = by_id
    ),
    by_gene$gene_name
  )

  # ncbiRefSeq IG/TCR: id-{gene} transcript IDs (supplement GTF symbol index).
  id_tx <- names(by_id)[startsWith(names(by_id), "id-")]
  id_genes <- sub("^id-", "", id_tx)
  for (i in seq_along(id_genes)) {
    g <- id_genes[[i]]
    if (!g %in% names(by_gene_primary)) {
      by_gene_primary[[g]] <- id_tx[[i]]
    }
  }

  list(
    by_id = by_id,
    by_gene_primary = by_gene_primary,
    gene_index = gene_idx
  )
}

#' Shared mouse reference lookup (lengths + GTF gene symbols); built once per R session.
#'
#' Brain and kidney figures only differ in the TUSCO TSV subset; the ncbiRefSeq+SIRV
#' annotation index is identical. Pass the return value as `lookup` to avoid rebuilding.
get_mouse_reference_annotation_lookup <- function(
    length_tbl = NULL,
    gtf_path = NULL,
    rebuild = FALSE) {
  if (isTRUE(rebuild)) {
    rm(list = ls(envir = .paper_mouse_annotation_cache), envir = .paper_mouse_annotation_cache)
  }
  if (exists("lookup", envir = .paper_mouse_annotation_cache, inherits = FALSE)) {
    return(get("lookup", envir = .paper_mouse_annotation_cache))
  }
  if (is.null(length_tbl)) {
    length_tbl <- load_mouse_transcript_lengths()
  }
  message(
    "Building reference annotation index from GTF gene symbols ",
    "(cached for this R session) ..."
  )
  lookup <- build_reference_annotation_lookup(length_tbl, gtf_path = gtf_path)
  assign("lookup", lookup, envir = .paper_mouse_annotation_cache)
  lookup
}

#' Map one TUSCO TSV row to a reference transcript ID (RefSeq accession or GTF gene symbol).
resolve_tusco_reference_transcript_id <- function(
    refseq,
    prot_refseq,
    gene_name,
    lookup) {
  rs <- paper_strip_transcript_version(refseq)
  if (!is.na(rs) && nzchar(rs) && rs %in% names(lookup$by_id)) {
    return(rs)
  }
  gn <- as.character(gene_name)
  if (!is.na(gn) && nzchar(gn) && gn %in% names(lookup$by_gene_primary)) {
    return(lookup$by_gene_primary[[gn]])
  }
  np <- paper_strip_transcript_version(prot_refseq)
  if (!is.na(np) && nzchar(np) && np %in% names(lookup$by_id)) {
    return(np)
  }
  NA_character_
}

#' Map SQANTI `associated_transcript` / `associated_gene` to a reference transcript ID.
resolve_sqanti_reference_transcript_id <- function(
    associated_transcript,
    associated_gene,
    lookup) {
  resolve_sqanti_reference_transcript_id_vec(
    associated_transcript,
    associated_gene,
    lookup
  )[[1L]]
}

#' Vectorized SQANTI transcript / gene → reference ID (for large classification tables).
resolve_sqanti_reference_transcript_id_vec <- function(
    associated_transcript,
    associated_gene,
    lookup) {
  tx <- paper_strip_transcript_version(associated_transcript)
  ag1 <- paper_strip_transcript_version(associated_gene)
  ag2 <- as.character(associated_gene)
  n <- length(tx)
  out <- rep(NA_character_, n)

  hit_tx <- !is.na(tx) & nzchar(tx) & tx %in% names(lookup$by_id)
  out[hit_tx] <- tx[hit_tx]

  miss <- is.na(out)
  hit_g1 <- miss & !is.na(ag1) & nzchar(ag1) & ag1 %in% names(lookup$by_gene_primary)
  out[hit_g1] <- unname(lookup$by_gene_primary[ag1[hit_g1]])

  miss <- is.na(out)
  hit_g2 <- miss & !is.na(ag2) & nzchar(ag2) & ag2 %in% names(lookup$by_gene_primary)
  out[hit_g2] <- unname(lookup$by_gene_primary[ag2[hit_g2]])

  out
}

#' Map one SQANTI-reads row (any category with `associated_transcript`) to a catalog ID.
resolve_expression_catalog_id <- function(
    associated_transcript,
    associated_gene,
    lookup,
    tusco_tbl,
    tusco_resolved) {
  resolve_expression_catalog_id_vec(
    associated_transcript,
    associated_gene,
    lookup,
    tusco_tbl,
    tusco_resolved
  )[[1L]]
}

#' Vectorized read → catalog transcript ID.
resolve_expression_catalog_id_vec <- function(
    associated_transcript,
    associated_gene,
    lookup,
    tusco_tbl,
    tusco_resolved) {
  out <- resolve_sqanti_reference_transcript_id_vec(
    associated_transcript,
    associated_gene,
    lookup
  )
  miss <- is.na(out)
  if (!any(miss)) {
    return(out)
  }

  tx <- paper_strip_transcript_version(associated_transcript[miss])
  tusco_tx <- paper_strip_transcript_version(tusco_tbl$transcript)
  idx <- match(tx, tusco_tx)
  fill <- tusco_resolved[idx]
  ok <- !is.na(fill)
  if (any(ok)) {
    miss_sub <- which(miss)
    out[miss_sub[ok]] <- fill[ok]
  }

  miss <- is.na(out)
  if (!any(miss)) {
    return(out)
  }

  ag <- paper_strip_transcript_version(associated_gene[miss])
  tusco_en <- paper_strip_transcript_version(tusco_tbl$ensembl)
  idx <- match(ag, tusco_en)
  fill <- tusco_resolved[idx]
  ok <- !is.na(fill)
  if (any(ok)) {
    miss_sub <- which(miss)
    out[miss_sub[ok]] <- fill[ok]
  }

  out
}

#' Resolve all TUSCO TSV rows to reference transcript IDs (with validation message).
resolve_tusco_reference_ids <- function(tusco_tbl, lookup, tissue_label) {
  resolved <- purrr::pmap_chr(
    list(tusco_tbl$refseq, tusco_tbl$prot_refseq, tusco_tbl$gene_name),
    resolve_tusco_reference_transcript_id,
    lookup = lookup
  )
  n_ok <- sum(!is.na(resolved))
  n_total <- nrow(tusco_tbl)
  message(
    "TUSCO (", tissue_label, "): ", n_ok, "/", n_total,
    " transcript IDs found in reference lengths"
  )
  if (n_ok < n_total) {
    missing <- tusco_tbl[is.na(resolved), , drop = FALSE]
    warning(
      "TUSCO ", tissue_label, ": ", nrow(missing), " ID(s) missing from reference: ",
      paste(
        paste0(missing$gene_name, " (", missing$transcript, ")"),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  resolved
}

#' TUSCO exonic lengths: subset of reference lengths for rows in the tissue TSV.
subset_tusco_lengths_from_reference <- function(
    tusco_tbl,
    lookup,
    tissue_label,
    resolved_ids = NULL) {
  if (is.null(resolved_ids)) {
    resolved_ids <- resolve_tusco_reference_ids(tusco_tbl, lookup, tissue_label)
  }
  keep <- !is.na(resolved_ids)
  tibble::tibble(
    transcript_set = "TUSCO",
    value = unname(lookup$by_id[resolved_ids[keep]])
  )
}

#' Match Reference / TUSCO / SIRVs (aligned with `tusco_ind.R` TUSCO logic).
assign_ref_tusco_sirv_set <- function(
    associated_transcript,
    associated_gene,
    chrom,
    tusco_sets) {
  at <- stringr::str_remove(as.character(associated_transcript), "\\.\\d+$")
  ag <- stringr::str_remove(as.character(associated_gene), "\\.\\d+$")
  ch <- as.character(chrom)

  dplyr::case_when(
    is.na(at) & is.na(ag) ~ NA_character_,
    stringr::str_starts(at, "SIRV") |
      stringr::str_starts(ag, "SIRV") |
      stringr::str_starts(ch, "SIRV") ~ "SIRVs",
    (!is.na(at) & at %in% tusco_sets$transcript) |
      (!is.na(ag) & ag %in% tusco_sets$ensembl) |
      (!is.na(ag) & ag %in% tusco_sets$gene_name) |
      (!is.na(ag) & ag %in% tusco_sets$refseq) ~ "TUSCO",
    TRUE ~ "Reference"
  )
}

#' Load SQANTI-reads classifications for one tissue (PacBio + ONT tables separate).
load_sqanti_reads_tissue <- function(pacbio_dir, ont_dir, sample_codes, sample_labels) {
  read_one <- function(paths, labels, technology) {
    out <- purrr::map2(
      paths,
      labels,
      function(path, sample_label) {
        if (!file.exists(path)) {
          return(tibble::tibble())
        }
        df <- readr::read_tsv(
          path,
          show_col_types = FALSE,
          col_select = c(
            "chrom",
            "structural_category",
            "associated_gene",
            "associated_transcript"
          )
        )
        df$sample <- sample_label
        df$technology <- technology
        df
      }
    )
    dplyr::bind_rows(out)
  }

  pacbio_files <- file.path(pacbio_dir, paste0(sample_codes, "_reads_classification.txt"))
  ont_files <- file.path(ont_dir, paste0(sample_codes, "_classification.txt"))

  list(
    pacbio = read_one(pacbio_files, sample_labels, "pacbio"),
    ont = read_one(ont_files, sample_labels, "ont")
  )
}

#' Stacked-bar counts from preloaded SQANTI-reads (one TSV read per tissue).
build_category_counts_from_class <- function(
    class_lists,
    prefix,
    label_base,
    read_numbers,
    category_map,
    display_levels) {
  sample_labels <- paste(label_base, 1:5)
  sample_codes <- paste0(prefix, 31:35)

  counts <- dplyr::bind_rows(class_lists$pacbio, class_lists$ont) %>%
    dplyr::mutate(
      category_label = dplyr::recode(.data$structural_category, !!!category_map)
    ) %>%
    dplyr::filter(!is.na(.data$category_label)) %>%
    dplyr::count(
      .data$category_label,
      .data$sample,
      .data$technology,
      name = "num_reads"
    )

  rn <- read_numbers %>%
    dplyr::filter(.data$sample %in% sample_codes) %>%
    dplyr::mutate(
      sample_label = paste(
        label_base,
        as.integer(stringr::str_sub(.data$sample, -1, -1))
      )
    )

  assigned <- counts %>%
    dplyr::group_by(.data$sample, .data$technology) %>%
    dplyr::summarise(assigned = sum(.data$num_reads), .groups = "drop")

  unaligned <- dplyr::bind_rows(
    rn %>%
      dplyr::transmute(
        sample = .data$sample_label,
        technology = "pacbio",
        total = .data$pb_fastq
      ),
    rn %>%
      dplyr::transmute(
        sample = .data$sample_label,
        technology = "ont",
        total = .data$ont_fastq
      )
  ) %>%
    dplyr::left_join(assigned, by = c("sample", "technology")) %>%
    dplyr::mutate(
      assigned = dplyr::coalesce(.data$assigned, 0),
      num_reads = pmax(0, .data$total - .data$assigned),
      category_label = "Unaligned"
    ) %>%
    dplyr::select(
      .data$category_label,
      .data$num_reads,
      .data$sample,
      .data$technology
    )

  dplyr::bind_rows(counts, unaligned) %>%
    dplyr::mutate(
      category_label = factor(.data$category_label, levels = display_levels),
      sample = factor(.data$sample, levels = sample_labels),
      technology = factor(.data$technology, levels = c("pacbio", "ont"))
    ) %>%
    dplyr::arrange(dplyr::desc(.data$num_reads))
}

#' Parse exonic transcript lengths from a GTF (sum of merged exon intervals per ID).
parse_gtf_transcript_lengths <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    stop("GTF not found: ", gtf_path, call. = FALSE)
  }
  intervals <- list()
  con <- file(gtf_path, open = "r")
  on.exit(close(con), add = TRUE)
  while (length(line <- readLines(con, n = 1L)) > 0L) {
    if (!nzchar(line) || startsWith(line, "#")) {
      next
    }
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 9L || parts[[3L]] != "exon") {
      next
    }
    start <- as.integer(parts[[4L]])
    end <- as.integer(parts[[5L]])
    if (is.na(start) || is.na(end)) {
      next
    }
    attrs <- sub("^[^;]*;\\s*", "", parts[[9L]])
    m <- regexec('transcript_id "([^"]+)"', attrs)
    hit <- regmatches(attrs, m)[[1L]]
    if (length(hit) < 2L) {
      next
    }
    tid <- hit[[2L]]
    if (start > end) {
      tmp <- start
      start <- end
      end <- tmp
    }
    intervals[[tid]] <- c(intervals[[tid]], start, end)
  }
  if (length(intervals) == 0L) {
    stop("No exon intervals parsed from GTF: ", gtf_path, call. = FALSE)
  }
  exonic_length <- function(x) {
    if (length(x) == 0L) {
      return(0L)
    }
    starts <- x[seq(1L, length(x), 2L)]
    ends <- x[seq(2L, length(x), 2L)]
    ord <- order(starts, ends)
    starts <- starts[ord]
    ends <- ends[ord]
    total <- 0L
    cur_s <- starts[[1L]]
    cur_e <- ends[[1L]]
    for (i in seq_along(starts)[-1L]) {
      if (starts[[i]] <= cur_e + 1L) {
        if (ends[[i]] > cur_e) {
          cur_e <- ends[[i]]
        }
      } else {
        total <- total + cur_e - cur_s + 1L
        cur_s <- starts[[i]]
        cur_e <- ends[[i]]
      }
    }
    total + cur_e - cur_s + 1L
  }
  tibble::tibble(
    transcript_id = names(intervals),
    length = vapply(intervals, exonic_length, integer(1L))
  )
}

#' Load mouse transcript lengths (TSV sidecar or GTF parse).
load_mouse_transcript_lengths <- function() {
  tsv_path <- paper_mouse_transcript_lengths_tsv()
  if (file.exists(tsv_path)) {
    return(
      readr::read_tsv(
        tsv_path,
        col_names = c("transcript_id", "length"),
        show_col_types = FALSE
      )
    )
  }
  message("Transcript lengths TSV not found; parsing GTF: ", paper_mouse_ref_sirv_gtf())
  parse_gtf_transcript_lengths(paper_mouse_ref_sirv_gtf())
}

#' Full transcript catalog per set (reference / TUSCO / SIRVs) for lengths and expression.
build_ref_tusco_sirv_transcript_catalog <- function(
    tusco_tissue = c("brain", "kidney"),
    length_tbl = NULL,
    lookup = NULL) {
  tusco_tissue <- match.arg(tusco_tissue)
  if (is.null(length_tbl)) {
    length_tbl <- load_mouse_transcript_lengths()
  }
  if (is.null(lookup)) {
    lookup <- get_mouse_reference_annotation_lookup(length_tbl = length_tbl)
  }

  tusco_tbl <- read_tusco_tissue_table(tusco_tissue)
  tusco_resolved <- resolve_tusco_reference_ids(
    tusco_tbl,
    lookup,
    tissue_label = tusco_tissue
  )
  tusco_ref_ids <- unique(tusco_resolved[!is.na(tusco_resolved)])
  all_ids <- names(lookup$by_id)
  sirv_ids <- all_ids[stringr::str_starts(all_ids, "SIRV")]
  ref_ids <- setdiff(all_ids, c(sirv_ids, tusco_ref_ids))

  tusco_sets <- load_tusco_annotation_sets(tusco_tissue)

  catalog <- dplyr::bind_rows(
    tibble::tibble(transcript_set = "Reference", catalog_id = ref_ids),
    tibble::tibble(transcript_set = "TUSCO", catalog_id = tusco_ref_ids),
    tibble::tibble(transcript_set = "SIRVs", catalog_id = sirv_ids)
  ) %>%
    dplyr::mutate(
      transcript_set = factor(.data$transcript_set, levels = REF_TUSCO_SIRV_SET_LEVELS)
    )

  list(
    catalog = catalog,
    lookup = lookup,
    tusco_tbl = tusco_tbl,
    tusco_resolved = tusco_resolved,
    tusco_ref_ids = tusco_ref_ids,
    tusco_sets = tusco_sets,
    tusco_tissue = tusco_tissue
  )
}

#' Expression catalog: Reference rows include tissue TUSCO IDs (for zero counts there).
expression_catalog_df <- function(catalog_info) {
  tusco_ref_ids <- catalog_info$tusco_ref_ids
  if (is.null(tusco_ref_ids)) {
    tusco_ref_ids <- unique(
      catalog_info$tusco_resolved[!is.na(catalog_info$tusco_resolved)]
    )
  }
  dplyr::bind_rows(
    catalog_info$catalog,
    tibble::tibble(
      transcript_set = factor("Reference", levels = REF_TUSCO_SIRV_SET_LEVELS),
      catalog_id = tusco_ref_ids
    )
  ) %>%
    dplyr::distinct(.data$transcript_set, .data$catalog_id)
}

#' GTF exonic lengths by Reference / TUSCO / SIRVs (99.5% trim on Reference only).
build_ref_tusco_sirv_length_df <- function(
    tusco_tissue = c("brain", "kidney"),
    length_tbl = NULL,
    catalog_info = NULL) {
  if (is.null(catalog_info)) {
    catalog_info <- build_ref_tusco_sirv_transcript_catalog(
      tusco_tissue = match.arg(tusco_tissue),
      length_tbl = length_tbl
    )
  }

  out <- catalog_info$catalog %>%
    dplyr::mutate(value = catalog_info$lookup$by_id[.data$catalog_id])

  ref_cap <- stats::quantile(
    out$value[out$transcript_set == "Reference"],
    0.995,
    na.rm = TRUE
  )
  out %>%
    dplyr::filter(
      .data$transcript_set != "Reference" | .data$value <= ref_cap
    ) %>%
    dplyr::transmute(
      transcript_set = .data$transcript_set,
      value = .data$value
    )
}

#' Summarize per-set read counts (including catalog transcripts with zero reads).
summarize_ref_tusco_sirv_expression <- function(expr_df) {
  expr_df %>%
    dplyr::group_by(.data$transcript_set) %>%
    dplyr::summarise(
      n_catalog = dplyr::n(),
      n_zero = sum(.data$value == 0L),
      n_positive = sum(.data$value > 0L),
      min_positive = if (any(.data$value > 0L)) {
        min(.data$value[.data$value > 0L])
      } else {
        NA_integer_
      },
      max_reads = max(.data$value),
      total_reads = sum(.data$value),
      .groups = "drop"
    )
}

#' Per-transcript read counts; catalog left-joined so missing transcripts = 0.
#'
#' Counts every SQANTI-reads row with a non-empty `associated_transcript` (FSM, ISM,
#' NIC, NNC, etc.), not FSM only. Each read is assigned to Reference / TUSCO / SIRVs
#' ([assign_ref_tusco_sirv_set()]), mapped to a catalog transcript ID, and counted
#' only when read set and catalog set agree.
#'
#' Reference rows with zero reads are dropped except tissue TUSCO transcript IDs
#' (still shown at 0 in the Reference panel). TUSCO and SIRV keep all catalog zeros.
build_ref_tusco_sirv_expression_df <- function(class_df, catalog_info) {
  tusco_ref_ids <- catalog_info$tusco_ref_ids
  if (is.null(tusco_ref_ids)) {
    tusco_ref_ids <- unique(
      catalog_info$tusco_resolved[!is.na(catalog_info$tusco_resolved)]
    )
  }
  catalog_df <- expression_catalog_df(catalog_info)
  lookup <- catalog_info$lookup
  tusco_tbl <- catalog_info$tusco_tbl
  tusco_resolved <- catalog_info$tusco_resolved
  tusco_sets <- catalog_info$tusco_sets
  if (is.null(tusco_sets)) {
    tusco_sets <- load_tusco_annotation_sets(catalog_info$tusco_tissue)
  }

  catalog_keys <- catalog_df %>%
    dplyr::select("transcript_set", "catalog_id") %>%
    dplyr::distinct()

  if (anyDuplicated(interaction(catalog_keys$transcript_set, catalog_keys$catalog_id))) {
    stop(
      "Duplicate (transcript_set, catalog_id) in expression catalog.",
      call. = FALSE
    )
  }

  with_tx <- class_df %>%
    dplyr::filter(
      !is.na(.data$associated_transcript),
      nzchar(as.character(.data$associated_transcript))
    )

  chrom_col <- if ("chrom" %in% names(with_tx)) {
    with_tx$chrom
  } else {
    rep(NA_character_, nrow(with_tx))
  }

  with_tx$read_set <- assign_ref_tusco_sirv_set(
    with_tx$associated_transcript,
    with_tx$associated_gene,
    chrom_col,
    tusco_sets
  )
  with_tx$catalog_id <- resolve_expression_catalog_id_vec(
    with_tx$associated_transcript,
    with_tx$associated_gene,
    lookup,
    tusco_tbl,
    tusco_resolved
  )

  matched <- with_tx %>%
    dplyr::filter(!is.na(.data$catalog_id), !is.na(.data$read_set)) %>%
    dplyr::inner_join(
      catalog_keys,
      by = c("catalog_id", "read_set" = "transcript_set")
    )

  n_matched <- nrow(matched)
  if (n_matched < sum(!is.na(with_tx$catalog_id))) {
    message(
      "Expression: dropped ",
      sum(!is.na(with_tx$catalog_id)) - n_matched,
      " reads with associated_transcript (not in catalog or set mismatch)"
    )
  }

  observed <- matched %>%
    dplyr::count(
      .data$read_set,
      .data$catalog_id,
      name = "expression"
    ) %>%
    dplyr::rename(transcript_set = read_set)

  out <- catalog_df %>%
    dplyr::left_join(observed, by = c("transcript_set", "catalog_id")) %>%
    dplyr::transmute(
      transcript_set = .data$transcript_set,
      catalog_id = .data$catalog_id,
      value = dplyr::coalesce(.data$expression, 0L)
    )

  out %>%
    dplyr::filter(
      .data$transcript_set != "Reference" |
        .data$value > 0L |
        .data$catalog_id %in% tusco_ref_ids
    ) %>%
    dplyr::select("transcript_set", "value")
}

plot_ref_tusco_sirv_violin <- function(
    df,
    ylab,
    log_y = FALSE,
    show_x_axis = TRUE) {
  fill_col <- paper_ref_tusco_sirv_fill()
  plot_df <- df
  if (isTRUE(log_y)) {
    # log10(reads + 1): catalog transcripts with 0 reads sit at 0 on the axis.
    plot_df <- dplyr::mutate(
      plot_df,
      value = log10(.data$value + 1L)
    )
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$transcript_set,
      y = .data$value,
      fill = .data$transcript_set
    )
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = TRUE,
      fill = fill_col,
      color = NA
    ) +
    ggplot2::geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.6,
      fill = fill_col,
      color = "grey20"
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(
        rep(fill_col, length(REF_TUSCO_SIRV_SET_LEVELS)),
        REF_TUSCO_SIRV_SET_LEVELS
      ),
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(limits = REF_TUSCO_SIRV_SET_LEVELS, drop = FALSE) +
    ggplot2::labs(x = NULL, y = ylab) +
    paper_panel_theme() +
    ggplot2::theme(legend.position = "none")

  if (isTRUE(log_y)) {
    p <- p +
      ggplot2::scale_y_continuous(
        labels = function(x) {
          counts <- pmax(0L, as.integer(round(10^x - 1L)))
          ifelse(
            counts < 1L,
            "0",
            scales::label_scientific(digits = 1)(counts)
          )
        }
      )
  }

  if (!isTRUE(show_x_axis)) {
    p <- p +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  } else {
    p <- p +
      ggplot2::theme(axis.text.x = paper_axis_text_x(0))
  }

  p
}

assemble_ref_tusco_sirv_figure <- function(
    length_plot,
    pacbio_plot,
    ont_plot,
    title) {
  # Wider/shorter layout for 2-column (180 mm) fit: length panel on top,
  # PacBio + ONT expression side-by-side below (was a tall 3-row stack).
  combined <- (length_plot / (pacbio_plot | ont_plot)) +
    patchwork::plot_layout(heights = c(1, 1))

  patchwork::wrap_elements(full = combined)
}

build_ref_tusco_sirv_tissue_figure_from_class <- function(
    class_lists,
    tissue_label,
    tusco_tissue = c("brain", "kidney"),
    lookup = NULL) {
  tusco_tissue <- match.arg(tusco_tissue)
  catalog_info <- build_ref_tusco_sirv_transcript_catalog(
    tusco_tissue = tusco_tissue,
    lookup = lookup
  )

  length_df <- build_ref_tusco_sirv_length_df(catalog_info = catalog_info)
  tb <- table(length_df$transcript_set)
  message(
    "Ref/TUSCO/SIRV length panel (", tissue_label, "): ",
    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", ")
  )

  expr_zero_frac <- function(class_df, tech_label) {
    expr_df <- build_ref_tusco_sirv_expression_df(class_df, catalog_info)
    summ <- summarize_ref_tusco_sirv_expression(expr_df)
    message(
      "Expression read counts (", tissue_label, ", ", tech_label, "): ",
      paste(
        sprintf(
          "%s catalog=%d zero=%d positive=%d (min>0=%s, max=%d, sum=%s)",
          summ$transcript_set,
          summ$n_catalog,
          summ$n_zero,
          summ$n_positive,
          ifelse(is.na(summ$min_positive), "—", format(summ$min_positive, big.mark = ",")),
          summ$max_reads,
          format(summ$total_reads, big.mark = ",")
        ),
        collapse = "; "
      )
    )
    expr_df
  }

  # Hoist the expression frames into variables so they can be exposed for
  # Source Data capture below (order preserved: PacBio then ONT).
  pacbio_expr_df <- expr_zero_frac(class_lists$pacbio, "PacBio")
  ont_expr_df    <- expr_zero_frac(class_lists$ont, "ONT")

  length_plot <- paper_inset_tag(
    plot_ref_tusco_sirv_violin(
      length_df,
      ylab = "Transcript length (nt)",
      show_x_axis = TRUE
    ),
    "a"
  )
  pacbio_plot <- paper_inset_tag(
    plot_ref_tusco_sirv_violin(
      pacbio_expr_df,
      ylab = "Expression (log10)",
      log_y = TRUE,
      show_x_axis = TRUE
    ),
    "b"
  )
  ont_plot <- paper_inset_tag(
    plot_ref_tusco_sirv_violin(
      ont_expr_df,
      ylab = "Expression (log10)",
      log_y = TRUE,
      show_x_axis = TRUE
    ),
    "c"
  )

  fig <- assemble_ref_tusco_sirv_figure(
    length_plot,
    pacbio_plot,
    ont_plot,
    title = tissue_label
  )
  # Attach the plotted frames for Source Data capture (ignored by rendering).
  # length: columns transcript_set, value(= transcript length in nt);
  # expression_*: columns transcript_set, value(= read count per catalog id).
  attr(fig, "source_data") <- list(
    length = length_df,
    expression_pacbio = pacbio_expr_df,
    expression_ont = ont_expr_df
  )
  fig
}

build_ref_tusco_sirv_tissue_figure <- function(
    pacbio_dir,
    ont_dir,
    sample_codes,
    sample_labels,
    tissue_label,
    tusco_tissue = c("brain", "kidney"),
    lookup = NULL) {
  class_lists <- load_sqanti_reads_tissue(
    pacbio_dir, ont_dir, sample_codes, sample_labels
  )
  build_ref_tusco_sirv_tissue_figure_from_class(
    class_lists,
    tissue_label = tissue_label,
    tusco_tissue = tusco_tissue,
    lookup = lookup
  )
}
