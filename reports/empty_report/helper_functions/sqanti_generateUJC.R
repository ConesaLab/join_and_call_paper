#' Generate UJC from SQANTI3 junctions file (Unique Junction Chain)
#'
#' This function reads a tab-separated junctions file, extracts and generates a unique 
#' junction chain (UJC) identifier for each isoform, and returns a dataframe containing 
#' the isoform and its corresponding UJC.
#'
#' @param junc_file A character string specifying the path to the junctions file.
#' @return A dataframe containing the isoform and its corresponding UJC.
#'
#' @examples
#' junc_file <- "junctions.txt"
#' ujc_data <- sqanti_generateUJC(junc_file)
#'
#' @import dplyr
#' @export
sqanti_generateUJC <- function(junc_file) {
  # Read the junctions file into a dataframe
  junc_df <- read.table(junc_file, header = TRUE, sep = "\t")
  
  # Group, summarize, and generate the UJC
  UJC_df <- junc_df %>%
    group_by(isoform, chrom, strand) %>%
    summarise(junction_chain = paste0(genomic_start_coord, "_", genomic_end_coord, collapse = "_")) %>%
    ungroup() %>%
    mutate(UJC = paste0(chrom, "_", strand, "_", junction_chain)) %>%
    select(c(isoform, UJC))
  
  return(UJC_df)
}
