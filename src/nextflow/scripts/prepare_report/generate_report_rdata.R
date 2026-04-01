# Generate the Rdata necesary to run the report

# Load packages
# library(ggpubr)
# library(ggfortify)
# library(NOISeq)
# library(RColorConesa)
# library(grid)
# library(UpSetR)
library(tidyverse)
# library(cowplot)

args <- commandArgs(trailingOnly=TRUE)

helper_dir <- args[1]

source(file.path(helper_dir, "gc_content_and_length.R")) # BioToyBox GitHub
source(file.path(helper_dir, "sqanti_generateUJC.R")) # BioToyBox GitHub


conditions_config_file <- file.path(args[2])
quant_ind <- file.path(args[3])
quant_concat <- file.path(args[4])
class_ind <- file.path(args[5])
class_concat <- file.path(args[6])
ref_genome <- file.path(args[7])
annot_ind_quantification <- file.path(args[8])
annot_concat_quantification <- file.path(args[9])
out_dir <- file.path(args[10])
class_ind_quantification <- class_ind
class_concat_quantification <- class_concat

conditions_config <- read.table(conditions_config_file, header = FALSE, sep = "\t",
                                col.names = c("cond_name", "class_fofn", "junc_fofn"))

# Themes and colors
xaxislevelsF1 <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron")
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion", "Intergenic", "Genic\nIntron")
cat.palette <- c("FSM" = "#6BAED6", "ISM" = "#FC8D59", "NIC" = "#78C679", "NNC" = "#EE6A50", "Genic\nGenomic" = "#969696", "Antisense" = "#66C2A4", "Fusion" = "goldenrod1", "Intergenic" = "darksalmon", "Genic\nIntron" = "#41B6C4")

# Define a function to read and process classification files
read_and_process_data <- function(class_fofn, junc_fofn) {
  class_fofn_lines <- readLines(class_fofn)
  junc_fofn_lines <- readLines(junc_fofn)

  class_df_list <- list()

  for (i in 1:length(class_fofn_lines)) {
    class_file <- class_fofn_lines[i]
    junc_file <- junc_fofn_lines[i]
    print(i)
    print(junc_file)
    file_name_without_ext <- sub("\\.[^.]+$", "", basename(class_file))
    file_name_without_ext <- paste0("S", i)
    class_df <- read.table(class_file, header = TRUE, sep = "\t")
    class_df$sample <- file_name_without_ext
    class_df$structural_category <- factor(
      class_df$structural_category,
      labels = xaxislabelsF1,
      levels = xaxislevelsF1,
      ordered = TRUE
    )
    if ("filter_result" %in% colnames(class_df)){
      class_df <- class_df[class_df$filter_result != "Artifact",]
    }
    junc_df <- sqanti_generateUJC(junc_file)

    class_df_list[[file_name_without_ext]] <- merge(class_df, junc_df, by = "isoform", all.x = TRUE)
    class_df_list[[file_name_without_ext]]$UJC <- paste(class_df_list[[file_name_without_ext]]$UJC, 
                                                    class_df_list[[file_name_without_ext]]$structural_category, 
                                                    sep = "_")
  }

  return(class_df_list)
}
all_class_df_lists <- list()
all_sample_labels <- list()
all_class_combined <- list()

for (j in 1:nrow(conditions_config)) {
  cond_name <- conditions_config$cond_name[j]
  class_fofn <- conditions_config$class_fofn[j]
  junc_fofn <- conditions_config$junc_fofn[j]

  print(paste("Processing condition:", cond_name))
  cond_list <- read_and_process_data(class_fofn, junc_fofn)
  all_class_df_lists[[cond_name]] <- cond_list

  var_name <- paste0(cond_name, "_class_df_list")
  assign(var_name, cond_list)
  save(list = var_name, file = paste0(out_dir, var_name, ".RData"))

  class_combined <- do.call(rbind, lapply(cond_list, function(df) df[, c("sample", "structural_category", "chrom")]))
  all_class_combined[[cond_name]] <- class_combined

  n_replicates <- length(cond_list) - 2  # subtract concat + TAMA entries
  sample_labels <- c('Concat\nReads', 'TAMA\nisoforms', paste0(cond_name, "_", seq_len(n_replicates)))
  all_sample_labels[[cond_name]] <- sample_labels
}

save(all_class_df_lists, all_sample_labels, all_class_combined,
     file = paste0(out_dir, "all_class_data.RData"))

# Backward-compatible aliases for the mouse dataset (B100K0 = brain, B0K100 = kidney)
if ("B100K0" %in% names(all_class_df_lists)) {
  Bclass_df_list <- all_class_df_lists[["B100K0"]]
  save(Bclass_df_list, file = paste0(out_dir, "Bclass_df_list.RData"))
}
if ("B0K100" %in% names(all_class_df_lists)) {
  Kclass_df_list <- all_class_df_lists[["B0K100"]]
  save(Kclass_df_list, file = paste0(out_dir, "Kclass_df_list.RData"))
}

# Get length, GC content, chrom, exons from reconstructed transcripts
# Function can be found in BioToyBox GitHub
ind_gclength <- gc_length(g=ref_genome, 
                      a= annot_ind_quantification)
save(ind_gclength, file = paste0(out_dir, "ind_gclength.RData"))

# Get length, GC content, chrom, exons from reconstructed transcripts
# Function can be found in BioToyBox GitHub
concat_gclength <- gc_length(g=ref_genome, 
                      a= annot_concat_quantification)
save(concat_gclength, file = paste0(out_dir, "concat_gclength.RData"))