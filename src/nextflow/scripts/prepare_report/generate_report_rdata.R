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


# Bclass_fofn <- file.path("classification_brain.fofn")
# Kclass_fofn <- file.path("classification_kidney.fofn")
# Bjunc_fofn <- file.path("junctions_brain.fofn")
# Kjunc_fofn <- file.path("junctions_kidney.fofn")
# quant_ind <- file.path("mouse_isoseq_ind_BK_rules.counts.tsv")
# quant_concat <- file.path("mouse_isoseq_jc_BK_rules_jc.counts.tsv")
# class_ind <- file.path("BK_rules_classification.txt")
# class_concat <- file.path("BK_rules_jc_classification.txt")
# ref_genome <- file.path(src_dir, "mm39_SIRV.fa")
# annot_ind_quantification <- file.path("BK_rules.gtf")
# annot_concat_quantification <- file.path("BK_rules_jc.gtf")

Bclass_fofn <- file.path(args[2])
Kclass_fofn <- file.path(args[3])
Bjunc_fofn <- file.path(args[4])
Kjunc_fofn <- file.path(args[5])
quant_ind <- file.path(args[6])
quant_concat <- file.path(args[7])
class_ind <- file.path(args[8])
class_concat <- file.path(args[9])
ref_genome <- file.path(args[10])
annot_ind_quantification <- file.path(args[11])
annot_concat_quantification <- file.path(args[12])
out_dir <- file.path(args[13])
class_ind_quantification <- class_ind
class_concat_quantification <- class_concat

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
    file_name_without_ext <- paste0("B", i) # sample name is not used as real names, sample labels are always used
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
print("Processing brain samples")
# Process Brain data
Bclass_df_list <- read_and_process_data(Bclass_fofn, Bjunc_fofn)

print("Processing kidney samples")
# Process Kidney data
Kclass_df_list <- read_and_process_data(Kclass_fofn, Kjunc_fofn)

# Save processed data
save(Bclass_df_list, file = paste0(out_dir, "Bclass_df_list.RData"))
save(Kclass_df_list, file = paste0(out_dir, "Kclass_df_list.RData"))

# Combine and plot classification files for Brain samples
class_combined_df_Brain <- do.call(rbind, lapply(Bclass_df_list, function(class_df) class_df[, c("sample", "structural_category", "chrom")]))
sample_labels_Brain <- c('Concat\nReads',  'TAMA\nisoforms', paste0("B3", 1:5))
#plots_Brain <- plot_classification_data(class_combined_df_Brain, sample_labels_Brain)

# Combine and plot classification files for Kidney samples
class_combined_df_Kidney <- do.call(rbind, lapply(Kclass_df_list, function(class_df) class_df[, c("sample", "structural_category", "chrom")]))
sample_labels_Kidney <- c('Concat\nReads',  'TAMA\nisoforms', paste0("K3", 1:5))
#plots_Kidney <- plot_classification_data(class_combined_df_Kidney, sample_labels_Kidney)

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