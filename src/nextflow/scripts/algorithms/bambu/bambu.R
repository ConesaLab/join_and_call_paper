library(bambu)

write_bambu <- function(bambu_output, path, prefix){
  read_manifest <- read_association(bambu_output)
  counts <- data.frame(assays(bambu_output)$counts)
  sel <- apply(counts, 1, function(x){!all(x==0)})
  counts <- counts[sel,, drop=F]
  bambu_found <- bambu_output[sel,]
  #writeToGTF(rowRanges(bambu_output), paste0(path, "/",prefix, ".gtf"))
  writeToGTF(rowRanges(bambu_found), paste0(path, "/",prefix, ".gtf"))
  write.table(read_manifest, 
              file = paste0(path, "/",prefix, "_read2trasncript.tsv"), 
              quote = F,
              row.names = F,
              sep = "\t")
  counts <- counts2SQ(counts)
  write.table(counts, 
              file = paste0(path, "/",prefix, "_counts.tsv"), 
              quote = F,
              row.names = F,
              sep = "\t")
}

# trasncripts associations are provided in different df for different samples
# and the associated transcript is the index in the se returned by bambu 
read_association <- function(bambu_output){
  process_row <- function(line) {
    # Convert idx to id and return a named vector
    formated_line <- c(
      readId = line[1],
      equalMatches = idx2id(line[2], transcripts_id),
      compatibleMatches = idx2id(line[3], transcripts_id)
    )
    return(formated_line)
  }
  transcripts_id <- rownames(bambu_output)
  # For multiple samples
  for (read_transcript in metadata(bambu_output)$readToTranscriptMaps){
    # read all the reads in the sample
    read_manifest <- apply(read_transcript, 1, process_row)
    read_manifest <- do.call(rbind, read_manifest)
  }
  return(read_manifest)
}

idx2id <- function(idx, ids){
  # No associated trasncripts were found for the read
  if (is.null(idx[[1]])){
    return("NA")
  # at least one associated transcripts
  } else{
    return(paste0(ids[idx[[1]]], collapse=","))
  }
}

# format counts the way SQANTI likes
counts2SQ <- function(counts){
  original_cols <- colnames(counts)
  if (ncol(counts) > 1){
    counts$superPBID <- rownames(counts)
    counts <- counts[, c("superPBID", original_cols)]
  } else {
    counts$pbid <- rownames(counts)
    counts <- counts[, c("pbid", original_cols)]
    colnames(counts) <- c("pbid", "count_fl")
  }
  return(counts)
}


args <- commandArgs(trailingOnly = TRUE)

genome = args[[1]]
annotation = args[[2]]
reads = args[[3]]
out_dir = args[[4]]
n_cores = args[[5]]

# Options
print("Running bambu with inputs...")
print("Genome:")
print(genome)
print("Annotation:")
print(annotation)
print("Reads:")
print(reads)
print("Out dir:")
print(out_dir)

# Prepare the path to write the outputs
path <- dirname(out_dir)
prefix <- basename(out_dir)

# Prepare the referecne annotation
annotation <- prepareAnnotations(annotation)

# Run bambu
se <- bambu(reads = reads, annotations = annotation, genome = genome, trackReads = TRUE, ncore = n_cores)

# save se
print("save se")
saveRDS(se, file = paste0(path, "/", prefix, ".RDS"))
# Write the outputs (gtf, counts and read2trasncript info)
print("writing output")
write_bambu(se, path, prefix)
