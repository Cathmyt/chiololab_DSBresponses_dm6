#!/usr/bin/env Rscript

# ============================================================
# Identify Restriction Sites (e.g., AsiSI) in a Genome FASTA
# ============================================================
# This script:
#   1) Reads a reference genome (FASTA file).
#   2) Searches for a specified DNA sequence (restriction site).
#   3) Outputs all matching sites as a BED file.
#
# Dependencies: Biostrings, rtracklayer
# ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(rtracklayer)
  library(GenomicRanges)
})

# --------------------- Utils ---------------------
timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

# --------------------- User Config ---------------------
# setwd("~/wdir/Cathy/Methods/hg19/bed_files")

# Input reference genome (FASTA) and target motif
fasta_file <- "../example_files/ref_genome/dm6.fa"
site_seq   <- "GCGATCGC"   # Target sequence to search (AsiSI recognition site)

# Output BED file
output_bed <- "AsiSI_sites_dm6.bed"

# --------------------- Main ---------------------
timestamp_message("Loading genome FASTA: ", fasta_file)
fasta <- readDNAStringSet(fasta_file)

timestamp_message("Searching for sequence motif: ", site_seq)
bed_list <- lapply(names(fasta), function(chr_name) {
  seq <- fasta[[chr_name]]
  matches <- matchPattern(site_seq, seq)
  
  if (length(matches) > 0) {
    data.frame(
      seqnames = chr_name,
      start    = start(matches) - 1,  # BED is 0-based
      end      = end(matches),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
})

timestamp_message("Combining results...")
bed_df <- do.call(rbind, bed_list)
bed_gr <- makeGRangesFromDataFrame(bed_df, keep.extra.columns = TRUE)

timestamp_message("Writing BED file: ", output_bed)
export.bed(bed_gr, output_bed)
