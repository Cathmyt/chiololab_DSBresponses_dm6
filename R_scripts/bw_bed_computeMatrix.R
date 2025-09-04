#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(rtracklayer)
  library(GenomicRanges)
})

# --------------------- Utils ---------------------
timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

# ---------------- Core functions ----------------

compute1ValPerSite <- function(bed, wig, w = 20000, seqlens, fun = "sum") {
  if (!inherits(wig, "SimpleRleList")) stop("wig must be a SimpleRleList")
  if (!fun %in% c("sum", "mean", "max")) stop("fun must be 'sum','mean','max'")
  
  n <- length(bed)
  out <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i %% 100 == 0) timestamp_message(i, "/", n)
    bedi <- bed[i]
    chr  <- as.character(seqnames(bedi))
    if (!chr %in% names(wig)) {
      out[i] <- NA_real_
      next
    }
    cov <- wig[[chr]]
    center <- start(bedi) + 4L  # AsiSI-specific offset
    stW <- max(1L, center - w)
    edW <- min(center + w, seqlens[chr], length(cov))
    
    if (fun %in% c("sum", "mean")) {
      v <- Views(cov, start = stW, end = edW)
      out[i] <- if (fun == "sum") sum(v) else mean(v)
    } else {
      v1 <- Views(cov, start = stW, end = center)
      v2 <- Views(cov, start = center + 1, end = edW)
      out[i] <- mean(c(max(v1), max(v2)))
    }
  }
  out
}

overlap_bw <- function(bed, bw_file, window_size = 1000, fun = "mean", seqlens) {
  bw <- import.bw(bw_file, as = "RleList")
  compute1ValPerSite(bed = bed, wig = bw, w = window_size, seqlens = seqlens, fun = fun)
}

# ---------------- main logic ----------------

#' Compute coverage matrix for multiple bw files with different window sizes
#'
#' @param bed_file Path to BED file
#' @param bw_files Vector of bigWig file paths
#' @param window_sizes Vector of window sizes, same length as bw_files
#' @param funct Summarization function ("mean","sum","max")
#' @param seqlens Chromosome lengths
#' @return data.frame with one row per site, one column per bw file
compute_matrix <- function(bed_file, bw_files, window_sizes, funct = "mean", seqlens) {
  if (length(bw_files) != length(window_sizes)) {
    stop("bw_files and window_sizes must be the same length")
  }
  
  # load bed once
  bed_sites <- tryCatch(
    import.bed(bed_file),
    error = function(e) {
      bed_df <- read.table(bed_file, header = FALSE)[, 1:3]
      colnames(bed_df) <- c("chr", "start", "end")
      makeGRangesFromDataFrame(bed_df)
    }
  )
  
  res <- list()
  
  for (i in seq_along(bw_files)) {
    bw <- bw_files[i]
    w  <- window_sizes[i]
    timestamp_message("Processing ", basename(bw), " with window ", w)
    signal <- overlap_bw(bed_sites, bw, w, funct, seqlens)
    colname <- paste0(tools::file_path_sans_ext(basename(bw)), "_", w)
    res[[colname]] <- signal
  }
  
  df <- as.data.frame(res)
  # add site coordinates as reference
  df <- cbind(as.data.frame(bed_sites)[,1:3], df)
  return(df)
}

# ---------------- usage ----------------
setwd("/Users/antienna/Documents/USC/ChioloLab/chiololab_DSBresponses_dm6")
seqlens <- seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)

bed_file <- "example_files/bed_files/AsiSIsite_dm6.bed"
bw_files <- c("example_files/bw_files/END_seq_ExoVII_5hrVsUNT_A74A73_bwa_aligned_PE_SES_wdup.bw",
              "example_files/bw_files/END_seq_HO_rep2_B23B24_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
              "example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
              "example_files/bw_files/FHA_Ku80_rep2_5hrVsUNT_L2R_A30A29_PE_BWA_SES_nodup.bw",
              "example_files/bw_files/FHA_Ku80_rep3_5hrVsUNT_L2R_A44A43_PE_BWA_SES_nodup.bw",
              "example_files/bw_files/FHA_RPA70_rep2_5hrVsUNT_A12A11_L2R_PE_BWA_SES_nodup.bw",
              "example_files/bw_files/FHA_RPA70_rep3_5hrVsUNT_A28A27_L2R_PE_BWA_SES_nodup.bw",
              "example_files/bw_files/GFP_Rad51_rep2_A48A47_5hrVsUNT_bwa_aligned_PE_SES.bw",
              "example_files/bw_files/GFP_Rad51_rep3_A83A82_5hrVsUNT_bwa_aligned_PE_SES.bw")
window_sizes <- c(250, 250, 250, 250, 250, 2500, 2500, 2500, 2500)

df <- compute_matrix(bed_file, bw_files, window_sizes, funct = "mean", seqlens = seqlens)

write.table(df, file = "example_files/output/signal_matrix_allAsiSI.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



