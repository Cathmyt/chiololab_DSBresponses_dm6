#!/usr/bin/env Rscript

# ============================================================
# Filter and Expand BED Sites
# ============================================================
# This script:
#   1) Reads a BED file of genomic sites.
#   2) Identifies "isolated" sites separated by a given distance.
#   3) Expands the window size around each site (optional).
#   4) Optionally filters sites based on a value column and threshold.
#   5) Saves processed BED files (original and sorted by value if applicable).
#
# Dependencies: data.table, GenomicRanges, tools, BSgenome.Dmelanogaster.UCSC.dm6
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(tools)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
})

# --------------------- Utils ---------------------
timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

# ----------------- Functions ---------------------

#' Find isolated sites based on minimum distance
#'
#' @param gr A GRanges object.
#' @param isolate_distance Minimum distance to other sites to be considered isolated.
#' @return A GRanges object of isolated sites.
find_isolated_sites <- function(gr, isolate_distance) {
  keep <- rep(FALSE, length(gr))
  
  if (length(gr) > 1) {
    same_chr <- as.character(seqnames(gr))
    prev_samechr <- c(FALSE, same_chr[-length(same_chr)] == same_chr[-1])
    next_samechr <- c(same_chr[-length(same_chr)] == same_chr[-1], FALSE)
    
    prev_dist <- c(Inf, start(gr)[-1] - end(gr)[-length(gr)])
    next_dist <- c(start(gr)[-1] - end(gr)[-length(gr)], Inf)
    
    keep <- ((!prev_samechr) | (prev_dist > isolate_distance)) &
            ((!next_samechr) | (next_dist > isolate_distance))
  } else {
    keep <- TRUE  # Keep single entry
  }
  
  return(gr[keep])
}

#' Expand genomic sites by a given window
#'
#' @param gr A GRanges object.
#' @param expand_window Integer, number of bp to expand on both ends.
#' @return Expanded GRanges object.
expand_sites <- function(gr, expand_window) {
  expanded_start <- pmax(0, start(gr) - expand_window)
  expanded_end   <- end(gr) + expand_window
  
  chr_names <- as.character(seqnames(gr))
  max_ends <- sapply(chr_names, function(chr) chrom_sizes[[chr]], USE.NAMES = FALSE)
  expanded_end <- pmin(expanded_end, max_ends)
  
  ranges(gr) <- IRanges(start = expanded_start, end = expanded_end)
  return(gr)
}

#' Read reference sequence lengths from FASTA index
#'
#' @param fasta_index Path to .fai file
#'
#' @return Named vector of contig lengths
read_length <- function(fasta_index) {
  if (!file.exists(fasta_index)) stop("FASTA index file not found: ", fasta_index)
  df <- read.delim(fasta_index, header = FALSE)
  setNames(as.integer(df[[2]]), df[[1]])
}

# ----------------- User Config -------------------
# setwd("/Users/chiololab/wdir/Cathy/bed_files")

# Parameters
isolate_distance <- 100000   # Minimum distance to consider a site "isolated"
expand_window    <- 0        # Expand window size for start/end
value_filter     <- FALSE    # Apply value-based filtering?
value_threshold  <- 0        # Threshold for filtering if enabled

# Input BED file (format: chr start end [optional value])
input_bed_file <- "../example_files/bed_files/AsiSIsite_dm6.bed"

# Genome chromosome sizes
chrom_sizes <- as.list(seqlengths(BSgenome.Dmelanogaster.UCSC.dm6))
# chrom_sizes <- as.list(read_length("../example_files/ref_genome/dm6.fa.fai"))

# ----------------- File Names --------------------
file_base <- file_path_sans_ext(basename(input_bed_file))
if (value_filter) {
  output_bed_file <- paste0(file_base, "_", isolate_distance, "bpIsolated_above", 
                            value_threshold, "_+-", expand_window, "bp.bed")
  output_sorted_bed_file <- paste0(file_base, "_", isolate_distance, "bpIsolated_above", 
                                   value_threshold, "_+-", expand_window, "bp_sorted.bed")
} else {
  output_bed_file <- paste0(file_base, "_", isolate_distance, "bpIsolated_+-", 
                            expand_window, "bp.bed")
  output_sorted_bed_file <- paste0(file_base, "_", isolate_distance, "bpIsolated_+-", 
                                   expand_window, "bp_sorted.bed")
}

# ----------------- Main --------------------------

timestamp_message("Loading BED file: ", input_bed_file)
bed_dt <- fread(input_bed_file, header = FALSE, sep = "\t")

if (value_filter) {
  bed_dt <- bed_dt[, 1:4]
  setnames(bed_dt, c("chr", "start", "end", "value"))
  bed_dt <- bed_dt[value >= value_threshold]
} else {
  bed_dt <- bed_dt[, 1:3]
  setnames(bed_dt, c("chr", "start", "end"))
}

bed_gr <- sort(makeGRangesFromDataFrame(bed_dt, keep.extra.columns = TRUE))

timestamp_message("Finding isolated sites (distance: ", isolate_distance, " bp)")
bed_isolated <- find_isolated_sites(bed_gr, isolate_distance)

timestamp_message("Expanding sites by ", expand_window, " bp")
bed_expanded <- expand_sites(bed_isolated, expand_window)

if (value_filter) {
  bed_expanded_dt <- as.data.table(bed_expanded)[, .(seqnames, start, end, value)]
} else {
  bed_expanded_dt <- as.data.table(bed_expanded)[, .(seqnames, start, end)]
}

# Save results
timestamp_message("Saving output: ", output_bed_file)
fwrite(bed_expanded_dt, output_bed_file, sep = "\t", col.names = FALSE, quote = FALSE)
message("Isolated sites: ", nrow(bed_expanded_dt), " / ", nrow(bed_dt), " (", 
        round(nrow(bed_expanded_dt) / nrow(bed_dt) * 100, 2), "%)")

if (value_filter) {
  timestamp_message("Saving sorted output: ", output_sorted_bed_file)
  bed_expanded_sorted <- bed_expanded_dt[order(-value)]
  fwrite(bed_expanded_sorted, output_sorted_bed_file, sep = "\t", 
         col.names = FALSE, quote = FALSE)
}

