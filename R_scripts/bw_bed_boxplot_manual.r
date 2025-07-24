#!/usr/bin/env Rscript

# ============================================================
# Boxplots of ChIP signal around BED sites from bigWig(s)
# ============================================================
# This script:
#   1) Reads a set of BED files (groups) and a list of bigWig files (replicates/conditions)
#   2) Computes signal (sum/mean/sum_positive/max) within each BED window +- extended window size defined
#   3) Produces boxplots (with jitter) and Wilcoxon significance bars

# Dependencies:
#   BSgenome.Dmelanogaster.UCSC.dm6
#   rtracklayer, GenomicRanges, ggplot2, ggsignif
# ============================================================

suppressPackageStartupMessages({
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(ggplot2)
  library(ggsignif)
  library(rtracklayer)
  library(GenomicRanges)
})

# --------------------- Utils ---------------------

timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

# ----------------- Core Functions ----------------

#' Compute one value per BED window
#'
#' @param bed GRanges object representing genomic windows
#' @param wig SimpleRleList from import.bw(as = "RleList")
#' @param w Extension around the window (bp)
#' @param seqlens Named integer vector of chromosome lengths
#' @param fun Summarization method: "sum", "mean", "max", or "sum_positive"
#'
#' @return Numeric vector of summarized signal per window
compute1ValPerWin <- function(bed, wig, w = 0, seqlens, fun = "sum") {
  if (!inherits(wig, "SimpleRleList")) stop("wig must be SimpleRleList")
  n <- length(bed)
  out <- numeric(n)
  
  for (i in seq_len(n)) {
    if (i %% 100 == 0) timestamp_message("Processed ", i, "/", n, " windows")
    bedi <- bed[i]
    chr <- as.character(seqnames(bedi))
    cov <- wig[[chr]]
    stW <- max(1, start(bedi) - w)
    edW <- min(end(bedi) + w, seqlens[chr], length(cov))
    
    v <- Views(cov, start = stW, end = edW)
    out[i] <- switch(fun,
                     sum = sum(v),
                     mean = mean(v),
                     max = max(v),
                     sum_positive = { # sums up only non-negative values within the specified window
                       m <- as.matrix(v)
                       m[m < 0] <- 0
                       sum(m)
                     },
                     stop("Unknown function: ", fun)
    )
  }
  out
}

#' Summarize bigWig signal over a BED file
#'
#' @param bed_file Path to BED file
#' @param bw_file Path to bigWig file
#' @param window_size Extension in bp
#' @param fun Summarization method
#'
#' @return Numeric vector of signals
overlap_bw <- function(bed_file, bw_file, window_size = 1000, fun = "mean") {
  if (!file.exists(bed_file)) stop("BED file not found: ", bed_file)
  if (!file.exists(bw_file)) stop("bigWig file not found: ", bw_file)
  
  timestamp_message("Reading BED: ", bed_file)
  bed_sites <- tryCatch(
    import.bed(bed_file),
    error = function(e) {
      bed_df <- read.table(bed_file)[, 1:3]
      colnames(bed_df) <- c("chr", "start", "end")
      makeGRangesFromDataFrame(bed_df)
    }
  )
  
  timestamp_message("Reading bigWig: ", bw_file)
  bw <- import.bw(bw_file, as = "RleList")
  
  compute1ValPerWin(bed = bed_sites, wig = bw, w = window_size, seqlens = seqlens, fun = fun)
}

#' Plot ChIP/ATAC boxplot
#'
#' @param results Data frame with Set, bwReads, and box_color columns
#' @param title Plot title
#' @param group_labels X-axis order
#' @param boxcenter_wMean Show mean line
#' @param boxcenter_noMedian Hide median
#' @param p_stat Group comparisons for Wilcoxon tests
#' @param show_pval_numeric Display numeric p-values
#'
#' @return ggplot2 object
plot_chipseq_boxplot <- function(results, title, group_labels = NULL,
                                 boxcenter_wMean = FALSE, boxcenter_noMedian = FALSE,
                                 p_stat = NULL, show_pval_numeric = TRUE) {
  timestamp_message("Generating plot: ", title)
  
  if (!is.null(group_labels)) results$Set <- factor(results$Set, levels = group_labels)
  
  p <- ggplot(results, aes(x = Set, y = bwReads, fill = box_color)) +
    geom_boxplot(outlier.shape = NA, width = 0.4) +
    geom_jitter(color = "black", size = 0.8, alpha = 0.3, width = 0.2) +
    scale_fill_identity() +
    theme_classic(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, y = "log2 ratio")
  
  if (boxcenter_wMean) {
    ltype <- if (boxcenter_noMedian) "solid" else "dashed"
    p <- p + stat_summary(geom = "errorbar", fun = mean, fun.min = mean, fun.max = mean,
                          width = 0.4, linetype = ltype)
  }
  
  # Add Wilcoxon tests
  if (!is.null(p_stat)) {
    comparisons <- lapply(p_stat, function(x) group_labels[x])
    ymax <- max(results$bwReads, na.rm = TRUE)
    y_positions <- ymax * 1.1 + seq(0, 0.05 * ymax, length.out = length(comparisons))
    
    pvals <- sapply(comparisons, function(grp) {
      wilcox.test(results$bwReads[results$Set == grp[1]],
                  results$bwReads[results$Set == grp[2]])$p.value
    })
    annotations <- sprintf("p = %.2g", pvals)
    
    for (i in seq_along(comparisons)) {
      p <- p + geom_signif(comparisons = list(comparisons[[i]]),
                           annotations = if (show_pval_numeric) annotations[i] else NULL,
                           y_position = y_positions[i],
                           tip_length = 0.03, textsize = 4,
                           map_signif_level = !show_pval_numeric)
    }
  }
  
  return(p)
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
# setwd("/Users/antienna/Documents/USC/ChioloLab/Chiolo_DSBresponses_dm6")

# genome lengths, either read from library or read from fasta index file
seqlens <- seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)
# seqlens <- read_length("~/wdir/Cathy/reference_genomes/hg19.fa.fai")

chipseq_file <- c(
  "../example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
  "../example_files/bw_files/END_seq_ExoVII_5hrVsUNT_A74A73_bwa_aligned_PE_SES_wdup.bw",
  "../example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
  "../example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
  "../example_files/bw_files/END_seq_ExoVII_5hrVsUNT_A74A73_bwa_aligned_PE_SES_wdup.bw"
)
bed_files <- c("../example_files/bed_files/EU_42DSBs.bed",
               "../example_files/bed_files/EU_42DSBs.bed",
               "../example_files/bed_files/dm6_random_200sites.bed",
               "../example_files/bed_files/HC_14DSBs.bed",
               "../example_files/bed_files/HC_14DSBs.bed")
group_labels <- c("rep1_EU", "rep2_EU", "reo1_random", "rep1_HC", "rep2_HC")
color_mapping <- c("darkgrey", "darkgrey", "orange", "#F8766D", "#F8766D")

p_stat <- list(c(1, 2), c(4,5))
window_sizes <- c(250, 2500)
funct <- "sum_positive"
boxcenter_wMean <- TRUE
boxcenter_noMedian <- TRUE
show_pval_numeric <- TRUE

# ----------------- Main -------------------
dir.create("plots", showWarnings = FALSE)

results <- data.frame()  # container for all computed signals

for (window_size in window_sizes) {
  timestamp_message("Processing window size: ", window_size)
  
  for (i in seq_along(chipseq_file)) {
    bw    <- chipseq_file[i]
    bed   <- bed_files[i]
    label <- group_labels[i]
    
    timestamp_message("Reading bigWig: ", basename(bw))
    signal <- overlap_bw(bed, bw, window_size, fun = funct)
    
    df <- data.frame(
      Set       = label,
      bwReads   = signal,
      box_color = color_mapping[i],
      stringsAsFactors = FALSE
    )
    results <- rbind(results, df)
  }
  
  # ----------------- Plot & Save -----------------
  plot_title <- paste0("boxplot_ENDseq_replicaes_acrossEUrandomHC_", funct, window_size, "bp")
  output_file <- file.path("plots", paste0(plot_title, "_pval.pdf"))
  
  timestamp_message("Saving plot to: ", output_file)
  pdf(output_file, width = 4.5, height = 6)
  p <- plot_chipseq_boxplot(
    results, plot_title, group_labels,
    boxcenter_wMean, boxcenter_noMedian, p_stat, show_pval_numeric
  )
  print(p)
  dev.off()
}

closeAllConnections()