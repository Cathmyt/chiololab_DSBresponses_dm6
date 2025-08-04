suppressPackageStartupMessages({
  library(dplyr)
  library(generics)
  library(plotly)
  library(tidymodels)
  library(GenomicRanges)
  library(rtracklayer)
})
# --------------------- Utils ---------------------

timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

# ----------------- Core Functions ----------------

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

#' Read a BED File with First Three Columns
#'
#' @param file Path to a BED file (must have at least 3 columns).
#'
#' @return A GRanges object with seqnames, start, and end positions.
#' @import GenomicRanges
read_bed3 <- function(file) {
  if (!file.exists(file)) stop("BED file not found: ", file)
  df <- read.table(file, header = FALSE, sep = "\t", comment.char = "", quote = "")
  df <- df[, 1:3]
  colnames(df) <- c("chr", "start", "end")
  makeGRangesFromDataFrame(df)
}



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
    read_bed3(bed_file),
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

#' Compute Averaged Signal Across Multiple bigWig Files
#'
#' @param bed_file Path to BED file.
#' @param bw_files Character vector of bigWig file paths.
#' @param window_size Window extension in base pairs.
#' @param funct Summarization function: "sum", "mean", "max", or "sum_positive".
#'
#' @return Numeric vector of averaged signal values across all bigWig files.
compute_bws <- function(bed_file, bw_files, window_size, funct = "sum") {
  sapply(bw_files, function(f) {
    overlap_bw(bed_file, f, window_size = window_size, fun = funct)
  }) %>% rowMeans(na.rm = TRUE)
}

#' Plot Signal Scatterplot Between Two Tracks over Intersecting Windows
#'
#' @param site_bed_file Path to BED file of sites.
#' @param win_bed1_file Path to first BED file of genomic windows (e.g., ChIP-seq peaks).
#' @param win_bed2_file Path to second BED file of genomic windows.
#' @param bw_files1 Vector of bigWig files corresponding to win_bed1_file.
#' @param bw_files2 Vector of bigWig files corresponding to win_bed2_file.
#' @param window_size1 Integer, extension size for bw_files1 signals.
#' @param window_size2 Integer, extension size for bw_files2 signals.
#' @param funct1 Summarization function for bw_files1.
#' @param funct2 Summarization function for bw_files2.
#' @param label1 X-axis label.
#' @param label2 Y-axis label.
#' @param plot_title Title of the scatter plot.
#'
#' @return A plotly scatter plot object with linear regression overlay.
#' @import plotly
#' @importFrom tidymodels linear_reg set_engine set_mode fit
#' @importFrom stats cor
plot_intersect_scatter <- function(site_bed_file, win_bed1_file, win_bed2_file,
                                   bw_files1, bw_files2,
                                   window_size1, window_size2,
                                   funct1 = "sum", funct2 = "sum",
                                   label1 = "X", label2 = "Y",
                                   plot_title = "Scatterplot") {

  site_gr <- read_bed3(site_bed_file)
  win1_gr <- read_bed3(win_bed1_file)
  win2_gr <- read_bed3(win_bed2_file)

  # Find intersection
  int1 <- subsetByOverlaps(site_gr, win1_gr)
  int2 <- subsetByOverlaps(site_gr, win2_gr)
  filtered_sites <- intersect(int1, int2)
  # Restrict windows to intersecting filtered_sites
  hits <- findOverlaps(filtered_sites, win1_gr)
  first_hits <- hits[!duplicated(queryHits(hits))]
  filtered_sites_matched <- filtered_sites[queryHits(first_hits)]
  filtered_win1  <- win1_gr[subjectHits(first_hits)]

  hits <- findOverlaps(filtered_sites, win2_gr)
  first_hits <- hits[!duplicated(queryHits(hits))]
  filtered_sites_matched <- filtered_sites[queryHits(first_hits)]
  filtered_win2  <- win2_gr[subjectHits(first_hits)]

  bed_sites_file <- tempfile(fileext = ".bed")
  bed_win1_file  <- tempfile(fileext = ".bed")
  bed_win2_file  <- tempfile(fileext = ".bed")
  export.bed(filtered_sites, bed_sites_file)
  export.bed(filtered_win1,  bed_win1_file)
  export.bed(filtered_win2,  bed_win2_file)

  # Compute enrichments
  X_avg <- compute_bws(bed_win1_file, bw_files1, window_size1, funct = funct1)
  Y_avg <- compute_bws(bed_win2_file, bw_files2, window_size2, funct = funct2)
  stopifnot(length(X_avg) == length(filtered_sites))
  stopifnot(length(Y_avg) == length(filtered_sites))
  
  df_all <- data.frame(chr = as.character(seqnames(filtered_sites)),
                       start = start(filtered_sites),
                       end = end(filtered_sites),
                       X = X_avg,
                       Y = Y_avg)

  # Regression
  lm_model <- linear_reg() %>% 
    set_engine('lm') %>% 
    set_mode('regression') %>%
    fit(Y ~ X, data = df_all)

  x_range <- seq(min(df_all$X), max(df_all$X), length.out = nrow(df_all))
  xdf <- data.frame(X = x_range)
  ydf <- predict(lm_model, xdf)
  xy <- cbind(xdf, ydf)
  colnames(xy)[2] <- "Y_e"

  # Plotting
  correlation <- cor(df_all$X, df_all$Y)
  str <- sprintf("correlation = %.4f", correlation)

  fig <- plot_ly()
  fig <- fig %>% add_trace(df_all , x = df_all$X, y = df_all$Y,
                           type = 'scatter', alpha = 1, mode = 'markers',
                           name = "Sites",
                           marker = list(size=4),
                           text = paste(df_all$chr, ":", df_all$start, "-", df_all$end, sep = ""),
                           hovertemplate = paste(
                             label2, ": %{y:e}<br>",
                             label1, ": %{x:e}<br>",
                             "site: %{text}",
                             "<extra></extra>", sep = ""))
  fig <- fig %>% add_trace(data = xy, x = ~X, y = ~Y_e,
                           name = 'Linear Fit', mode = 'lines', alpha = 1)
  fig <- fig %>% layout(xaxis = list(title = label1), 
                        yaxis = list(title = label2), 
                        title = plot_title,
                        legend = list(title = list(text = str)))
  return(fig)
}


# ----------------- User Config -------------------
setwd("/Users/antienna/Documents/USC/ChioloLab/chiololab_DSBresponses_dm6/R_scripts")
seqlens <- read_length("../example_files/ref_genome/dm6.fa.fai")
fig <- plot_intersect_scatter(
  site_bed_file = "../example_files/bed_files/AsiSIsite_dm6.bed",
  win_bed1_file = "../example_files/bed_files/PE_END_seq_rep1_broad_p0.05_dm6_larracuente_peaks.broadPeak",
  win_bed2_file = "../example_files/bed_files/H3K9me3_SNAP_ChIP_dm6_UNTvIN_broad_nolambda_bw300_p0.05_peaks.broadPeak",
  bw_files1 = c("../example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
   "../example_files/bw_files/END_seq_HO_rep2_B23B24_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
   "../example_files/bw_files/END_seq_ExoVII_5hrVsUNT_A74A73_bwa_aligned_PE_SES_wdup.bw"),
  bw_files2 = c("../example_files/bw_files/B26_H3K9me3_SNAP_ChIP_rep1_UNT_bwa_aligned_pe_dm6KMetStatPombe_CPM.bw"),
  window_size1 = 0,
  window_size2 = 0,
  funct1 = "sum", funct2 = "sum",
  label1 = "ENDseq", label2 = "H3K9me3",
  plot_title = "Scatterplot: ENDseq vs H3K9me3"
)
fig