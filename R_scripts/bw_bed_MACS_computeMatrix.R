#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(tools)
})

# --------------------- Utils ---------------------
timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

# Read sequence lengths from FASTA .fai
read_length <- function(fasta_index) {
  if (!file.exists(fasta_index)) stop("FASTA index (.fai) not found: ", fasta_index)
  df <- read.delim(fasta_index, header = FALSE)
  setNames(as.integer(df[[2]]), df[[1]])
}

# ---------------- Core functions ----------------
read_bed3 <- function(file) {
  if (!file.exists(file)) stop("BED file not found: ", file)
  df <- read.table(file, header = FALSE, sep = "\t", comment.char = "", quote = "")
  if (ncol(df) < 3) stop("BED must have at least 3 columns: ", file)
  df <- df[, 1:3]
  colnames(df) <- c("chr", "start", "end")
  makeGRangesFromDataFrame(df)
}

# Summarize bigWig signal over a set of windows (extend by w on both sides)
compute1ValPerWin <- function(bed, wig, w = 0, seqlens, fun = "sum") {
  if (!inherits(wig, "SimpleRleList")) stop("wig must be SimpleRleList")
  n <- length(bed)
  out <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i %% 100 == 0) timestamp_message("Processed ", i, "/", n, " windows")
    bedi <- bed[i]
    chr  <- as.character(seqnames(bedi))
    if (!chr %in% names(wig)) { out[i] <- NA_real_; next }
    
    cov <- wig[[chr]]
    stW <- max(1L, start(bedi) - w)
    edW <- min(end(bedi) + w, seqlens[chr], length(cov))
    
    v <- Views(cov, start = stW, end = edW)
    out[i] <- switch(
      fun,
      sum  = sum(v),
      mean = mean(v),
      max  = max(v),
      sum_positive = {
        m <- as.matrix(v)
        m[m < 0] <- 0
        sum(m)
      },
      stop("Unknown fun: ", fun, " (use 'sum','mean','max','sum_positive')")
    )
  }
  out
}

overlap_bw <- function(bed, bw_file, window_size = 0, fun = "mean", seqlens) {
  if (!file.exists(bw_file)) stop("bigWig file not found: ", bw_file)
  bw <- import.bw(bw_file, as = "RleList")
  compute1ValPerWin(bed = bed, wig = bw, w = window_size, seqlens = seqlens, fun = fun)
}

# ---------------- Main Logic ----------------
# Compute a matrix of signals (rows = global sites, cols = bigWigs).
# If bed_files is provided, each bigWig uses its corresponding MACS window BED.
# Sites not present in a given window BED get NA for that column.
compute_matrix <- function(bed_file,
                           bw_files,
                           window_sizes,
                           funct = "mean",
                           seqlens,
                           bed_files = NULL) {
  if (length(bw_files) != length(window_sizes)) {
    stop("Length mismatch: length(bw_files) = ", length(bw_files),
         " but length(window_sizes) = ", length(window_sizes))
  }
  if (!is.null(bed_files) && length(bed_files) != length(bw_files)) {
    stop("If bed_files is provided, it must match bw_files in length. ",
         "length(bed_files) = ", length(bed_files),
         ", length(bw_files) = ", length(bw_files))
  }
  
  # Load global reference sites (these define the output row set)
  bed_sites <- read_bed3(bed_file)
  res <- vector("list", length(bw_files))
  
  for (i in seq_along(bw_files)) {
    bw <- bw_files[i]
    w  <- window_sizes[i]
    timestamp_message("Processing ", basename(bw), " with window +/-", w, " bp")
    
    # Pre-fill with NA to preserve row count
    signal_full <- rep(NA_real_, length(bed_sites))
    
    if (is.null(bed_files)) {
      # Global behavior: compute over the same site windows for all bws
      signal_full <- overlap_bw(bed_sites, bw, w, funct, seqlens)
    } else {
      # Per-bw MACS window behavior
      macs_bed <- read_bed3(bed_files[i])
      
      # Find which global sites fall in the MACS window(s)
      # "falls in a window" == any overlap of the site interval with the peak interval
      hits <- findOverlaps(bed_sites, macs_bed)
      if (length(hits) > 0) {
        idx_ref <- queryHits(hits)     # indices in global bed
        idx_macs <- subjectHits(hits)  # indices in macs_bed
        vals <- overlap_bw(macs_bed[idx_macs], bw, w, funct, seqlens)
        # Assign back to the matching site rows; others remain NA
        signal_full[idx_ref] <- vals
      }
    }
    
    colname <- paste0(file_path_sans_ext(basename(bw)), "_", w)
    res[[i]] <- signal_full
    names(res)[i] <- colname
  }
  
  out <- cbind(as.data.frame(bed_sites)[, 1:3, drop = FALSE], as.data.frame(res, check.names = FALSE))
  return(out)
}

# ---------------- Example usage ----------------
setwd("/Users/antienna/Documents/USC/ChioloLab/chiololab_DSBresponses_dm6")
seqlens <- seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)

bed_file <- "example_files/bed_files/AsiSIsite_dm6.bed"
bw_files <- c(
  "example_files/bw_files/END_seq_ExoVII_5hrVsUNT_A74A73_bwa_aligned_PE_SES_wdup.bw",
  "example_files/bw_files/END_seq_HO_rep2_B23B24_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
  "example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
  "example_files/bw_files/FHA_Ku80_rep2_5hrVsUNT_L2R_A30A29_PE_BWA_SES_nodup.bw",
  "example_files/bw_files/FHA_Ku80_rep3_5hrVsUNT_L2R_A44A43_PE_BWA_SES_nodup.bw",
  "example_files/bw_files/FHA_RPA70_rep2_5hrVsUNT_A12A11_L2R_PE_BWA_SES_nodup.bw",
  "example_files/bw_files/FHA_RPA70_rep3_5hrVsUNT_A28A27_L2R_PE_BWA_SES_nodup.bw",
  "example_files/bw_files/GFP_Rad51_rep2_A48A47_5hrVsUNT_bwa_aligned_PE_SES.bw",
  "example_files/bw_files/GFP_Rad51_rep3_A83A82_5hrVsUNT_bwa_aligned_PE_SES.bw"
)
win_files <- c(
  "example_files/bed_files/END_seq_ExoVII_broad_p0.05_peaks.broadPeak",
  "example_files/bed_files/END_seq_HO_rep2_broad_p0.05_peaks.broadPeak",
  "example_files/bed_files/END_seq_rep1_broad_p0.05_peaks.broadPeak",
  "example_files/bed_files/FHA_Ku80_rep2_broad_nolambda_p0.05_peaks.broadPeak",
  "example_files/bed_files/FHA_Ku80_rep3_broad_nolambda_p0.05_peaks.broadPeak",
  "example_files/bed_files/FHA_RPA70_rep2_broad_nolambda_p0.05_peaks.broadPeak",
  "example_files/bed_files/FHA_RPA70_rep3_broad_nolambda_p0.05_peaks.broadPeak",
  "example_files/bed_files/GFP_Rad51_rep2_broad_nolambda_p0.05_peaks.broadPeak",
  "example_files/bed_files/GFP_Rad51_rep3_broad_nolambda_p0.05_peaks.broadPeak"
)
window_sizes <- rep(0, length(bw_files))  # must match bw_files length

df <- compute_matrix(
  bed_file,
  bw_files,
  window_sizes,
  funct = "mean",
  seqlens,
  bed_files = win_files # or NULL
)
write.table(df, file = "example_files/output/signal_matrix_MACS.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
