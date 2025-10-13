
# Packages
library(tidyverse)
library(stringr)

# Input
args <- commandArgs(trailingOnly = TRUE)
quantification_fofn <- args[1]
tama_merge_file <- args[2]
out_file <- args[3]

quantification_fofn <- readLines(quantification_fofn)
tama_merge_table <- read.table(tama_merge_file, header = FALSE, sep = "\t")

# Format TAMA merge dataframe
# remove_strings <- c("_talon")
# tama_merge_table_cleaned <- tama_merge_table
# tama_merge_table_cleaned[, 4] <- stringr::str_replace_all(tama_merge_table[, 4], 
#                                                           paste(remove_strings, collapse = "|"), "")
tama_id_map <- as.data.frame(
  stringr::str_split_fixed(tama_merge_table[, 4], pattern = ";|_", 3)
)

colnames(tama_id_map) <- c("new_id", "sample_id", "old_id")

print("head of tama_id_map")
print(head(tama_id_map))

# Change old transcript id to merged transcript id
quant_mat_list <- list()
for (i in 1:length(quantification_fofn)) {
  quant_fofn_line <- strsplit(quantification_fofn[i], "\t")[[1]]
  sample_id <- quant_fofn_line[1]
  file <- quant_fofn_line[2]
  quant_file <- read.table(file, header = TRUE, sep = "\t")
  if (ncol(quant_file) == 2) {
    colnames(quant_file) <- c("pbid", sample_id)
  } else {
    colnames(quant_file)[1] <- "pbid"
  }

  print("colnames of quant file")
  print(colnames(quant_file))

  tama_name_by_sample <- tama_id_map[tama_id_map$sample_id == sample_id, ]
  match_ids <- match(quant_file$pbid, tama_name_by_sample$old_id)
  quant_file$transcript_id <- tama_name_by_sample$new_id[match_ids]
  quant_file <- subset(quant_file, select = -c(pbid))

  print("head of quant_file after transcript_id reassignment")
  print(head(quant_file))

  quant_file <- quant_file %>%
    pivot_longer(-transcript_id, names_to = "sample", values_to = "count") %>%
    group_by(transcript_id, sample) %>%
    summarise(tot_count = sum(count)) %>%
    pivot_wider(names_from = "sample", values_from = "tot_count")

  quant_mat_list[[i]] <- as.data.frame(quant_file)
}

# Merge quantification matrices
quant_merged_matrix <- quant_mat_list %>%
  reduce(full_join, by = "transcript_id") %>%
  replace(is.na(.), 0) %>%
  relocate(transcript_id)

# rename to fit sqanti FL format (will always be more than 1 sample)
names(quant_merged_matrix)[names(quant_merged_matrix) == 'transcript_id'] <- 'superPBID'

# Write output
write.table(quant_merged_matrix, file = out_file,
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
