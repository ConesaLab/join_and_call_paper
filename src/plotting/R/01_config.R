# 01_config.R
# Paths, constants, sample labels, read numbers, FL filter levels

# base_path <- "/mnt/c/Users/jetzi/other_repos/documenting_NIH/fabian"
base_path <- "C:/Users/jetzi/other_repos/documenting_NIH/fabian"

isoseq_report_path <- file.path(base_path, "reports", "isoseq")
ont_report_path    <- file.path(base_path, "reports", "ont")

get_paths <- function(src_dir) {
  data_dir <- file.path(src_dir, "data")
  list(
    src_dir = src_dir,
    data_dir = data_dir,
    quant_cj_brain = file.path(src_dir, "B100K0_tama.counts.tsv"),
    quant_cj_kidney = file.path(src_dir, "B0K100_tama.counts.tsv"),
    Bclass_fofn = file.path(data_dir, "classification_brain.fofn"),
    Kclass_fofn = file.path(data_dir, "classification_kidney.fofn"),
    Bjunc_fofn = file.path(data_dir, "junctions_brain.fofn"),
    Kjunc_fofn = file.path(data_dir, "junctions_kidney.fofn"),
    quant_ind = file.path(data_dir, "ind_tama.counts.tsv"),
    quant_concat = file.path(data_dir, "concat_tama.counts.tsv"),
    class_ind = file.path(data_dir, "ind_TAMA_classification.txt"),
    class_concat = file.path(data_dir, "concat_TAMA_classification.txt"),
    annot_ind_quantification = file.path(data_dir, "ind_TAMA.gtf"),
    annot_concat_quantification = file.path(data_dir, "concat_TAMA.gtf"),
    class_ind_quantification = file.path(data_dir, "ind_TAMA_classification.txt"),
    class_concat_quantification = file.path(data_dir, "concat_TAMA_classification.txt")
  )
}

# Directory structure for all platform/method combinations
src_dirs <- list(
  IsoSeq = list(
    IsoQuant      = file.path(isoseq_report_path, "isoquant", "run3_report"),
    FLAIR         = file.path(isoseq_report_path, "flair_ar_sr", "run7_report"),
    TALON         = file.path(isoseq_report_path, "talon", "run3_report"),
    Mandalorion   = file.path(isoseq_report_path, "mandalorion", "run1_report"),
    Bambu         = file.path(isoseq_report_path, "bambu", "run3_report"),
    isoseq_sqanti = file.path(isoseq_report_path, "isoseq_isoseq_new_filter_data", "report")
  ),
  ONT = list(
    IsoQuant = file.path(ont_report_path, "isoquant", "run2_report"),
    FLAIR    = file.path(ont_report_path, "flair_ar_sr", "run2_report"),
    TALON    = file.path(ont_report_path, "talon", "run1_report"),
    Bambu    = file.path(ont_report_path, "bambu", "run2_report")
  )
)

check_paths <- function(dirs, parent = "") {
  required_files <- c("Bclass_df_list.RData", "Kclass_df_list.RData")

  for (name in names(dirs)) {
    path_or_sublist <- dirs[[name]]
    current <- if (parent == "") name else paste(parent, name, sep = " -> ")

    if (is.list(path_or_sublist)) {
      check_paths(path_or_sublist, current)
    } else {
      if (!dir.exists(path_or_sublist)) {
        cat(sprintf("Path does NOT exist: [%s] %s\n", current, path_or_sublist))
      } else {
        cat(sprintf("Path OK: [%s]\n", current))
        for (file in required_files) {
          file_path <- file.path(path_or_sublist, file)
          if (!file.exists(file_path)) {
            cat(sprintf("  Missing file: %s\n", file))
          } else {
            cat(sprintf("  Found file: %s\n", file))
          }
        }
      }
    }
  }
}

# Sample codes
sample_codes_brain <- paste0("B", 1:7)
sample_codes_kidney <- paste0("K", 1:7)

# PacBio read numbers
fastq_reads_pb_brain <- c(0, 0, 3687280, 3857138, 4415023, 3889091, 3915568)
fastq_reads_pb_kidney <- c(0, 0, 4598058, 5138219, 3903934, 5437762, 5097532)

pacbio_read_numbers_brain <- setNames(fastq_reads_pb_brain, sample_codes_brain)
pacbio_read_numbers_kidney <- setNames(fastq_reads_pb_kidney, sample_codes_kidney)

pacbio_read_numbers_brain["B1"] <- sum(fastq_reads_pb_brain)
pacbio_read_numbers_brain["B2"] <- sum(fastq_reads_pb_brain)

pacbio_read_numbers_kidney["K1"] <- sum(fastq_reads_pb_kidney)
pacbio_read_numbers_kidney["K2"] <- sum(fastq_reads_pb_kidney)

# ONT read numbers
fastq_reads_ont_brain <- c(0, 0, 11225327, 7712637, 9969075, 9282465, 5311045)
fastq_reads_ont_kidney <- c(0, 0, 6433527, 8626753, 6705708, 4542656, 4857879)

ont_read_numbers_brain <- setNames(fastq_reads_ont_brain, sample_codes_brain)
ont_read_numbers_kidney <- setNames(fastq_reads_ont_kidney, sample_codes_kidney)

ont_read_numbers_brain["B1"] <- sum(fastq_reads_ont_brain)
ont_read_numbers_brain["B2"] <- sum(fastq_reads_ont_brain)

ont_read_numbers_kidney["K1"] <- sum(fastq_reads_ont_kidney)
ont_read_numbers_kidney["K2"] <- sum(fastq_reads_ont_kidney)

# Sample names and labels
sample_names_Brain <- c('Join&Call', 'Call&Join', paste0("Brain ", 1:5))
sample_names_Kidney <- c('Join&Call', 'Call&Join', paste0("Kidney ", 1:5))

sample_labels_Brain <- c('J&C', 'C&J', paste0("B", 1:5))
sample_labels_Kidney <- c('J&C', 'C&J', paste0("K", 1:5))

# FL filter levels
fl_filter_levels <- list(1, 3, 5, 10, 20)
