#!/usr/bin/env Rscript

# ============================================================
# Boxplots of ChIP signal around BED sites from bigWig(s)
# ============================================================
# This script:
#   1) Reads a set of BED files (groups) and a list of bigWig files (replicates/conditions)
#   2) For each site, summarizes bigWig signal in a +/- window around the site center
#   3) Produces boxplots (with jitter) and Wilcoxon significance bars
#
# Dependencies: BSgenome.Dmelanogaster.UCSC.dm6, rtracklayer, GenomicRanges, ggplot2, ggsignif
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

# ----------------- Core functions ----------------

#' Compute one value per site per bigWig
#' specific to AsiSI, if using a different restriction enzyme, be sure to change the center
#'
#' @param bed GRanges of sites
#' @param wig SimpleRleList from import.bw(as = "RleList")
#' @param w half-window size in bp
#' @param seqlens named integer vector of chromosome lengths
#' @param fun summarization: "sum", "mean", "max"
#' @return numeric vector, length = length(bed)
compute1ValPerSite <- function(bed, wig, w = 20000, seqlens, fun = "sum") {
  if (!inherits(wig, "SimpleRleList")) {
    stop("wig must be a SimpleRleList (got: ", class(wig), ")")
  }
  if (!fun %in% c("sum", "mean", "max")) {
    stop("fun must be one of 'sum', 'mean', 'max'")
  }

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
    center <- start(bedi) + 4L
    stW <- max(1L, center - w)
    edW <- min(center + w, seqlens[chr], length(cov))

    if (fun %in% c("sum", "mean")) {
      v <- Views(cov, start = stW, end = edW)
      out[i] <- if (fun == "sum") sum(v) else mean(v)
    } else {
      v1 <- Views(cov, start = stW,       end = center)
      v2 <- Views(cov, start = center + 1, end = edW)
      out[i] <- mean(c(max(v1), max(v2)))
    }
  }
  out
}

#' Summarize bigWig signal over a BED for a single bigWig file
#'
#' @param bed_file path to BED
#' @param bw_file path to bigWig
#' @param window_size integer, half-window
#' @param fun "mean", "sum", "max"
#' @param seqlens chromosome lengths
#' @return numeric vector of signals
overlap_bw <- function(bed_file, bw_file, window_size = 1000, fun = "mean", seqlens) {
  bed_sites <- tryCatch(
    import.bed(bed_file),
    error = function(e) {
      bed_df <- read.table(bed_file, header = FALSE)[, 1:3]
      colnames(bed_df) <- c("chr", "start", "end")
      makeGRangesFromDataFrame(bed_df)
    }
  )

  bw <- import.bw(bw_file, as = "RleList")
  compute1ValPerSite(bed = bed_sites, wig = bw, w = window_size, seqlens = seqlens, fun = fun)
}

#' Compute coverage for multiple BED groups against one bigWig
#'
#' @param bed_files character vector of bed files (groups)
#' @param bw_file single bigWig file
#' @param window_size int
#' @param group_labels optional labels matching bed_files
#' @param funct summary function c("mean","sum","max")
#' @return data.frame with Set, label, bwReads
compute_cov_chipseq <- function(bed_files, bw_file, window_size = 1000,
                                group_labels = NULL, funct = "mean",
                                seqlens) {
  res <- vector("list", length(bed_files))
  for (i in seq_along(bed_files)) {
    bed_file <- bed_files[i]
    set_label <- if (is.null(group_labels)) {
      tools::file_path_sans_ext(basename(bed_file))
    } else {
      group_labels[i]
    }

    timestamp_message("Processing ", set_label, " with ", basename(bw_file))
    signal <- overlap_bw(bed_file, bw_file, window_size, funct, seqlens)
    res[[i]] <- data.frame(Set = set_label,
                           label = sub("^\\d+", "", set_label),
                           bwReads = signal,
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, res)
}

#' Plot boxplot with optional Wilcoxon tests
#'
#' @param results data.frame(Set, label, bwReads)
#' @param title plot title
#' @param group_labels ordering for x-axis
#' @param color_map_fill named vector for fill colors
#' @param boxcenter_wMean logical, add mean as dashed/solid line
#' @param boxcenter_noMedian logical, only show mean instead of median
#' @param show_pval_numeric logical, print numeric p-values
#' @return ggplot object
plot_chipseq_boxplot <- function(results, title = "Boxplot", group_labels = NULL,
                                 color_map_fill,
                                 boxcenter_wMean = FALSE,
                                 boxcenter_noMedian = FALSE,
                                 show_pval_numeric = FALSE) {

  if (!is.null(group_labels)) {
    results$Set <- factor(results$Set, levels = group_labels)
  } else {
    results$Set <- factor(results$Set, levels = unique(results$Set))
  }

  p <- ggplot(results, aes(x = Set, y = bwReads, fill = label))

  if (boxcenter_noMedian) {
    p <- p + geom_boxplot(fatten = NULL, outlier.shape = NA, width = 0.4)
  } else {
    p <- p + geom_boxplot(outlier.shape = NA, width = 0.4)
  }

  p <- p +
    scale_fill_manual(values = color_map_fill) +
    geom_jitter(color = "black", size = 0.8, alpha = 0.5, width = 0.2) +
    # geom_signif(comparisons = list(group_labels), map_signif_level=TRUE) + # default as wilcoxon test
    theme_classic() +
    # coord_cartesian(ylim=c(0, NA)) + # set x/y axis limits
    theme(axis.text = element_text(size = 15, color = "black"),
          text = element_text(size = 15, color = "black")) +
    ggtitle(title) +
    ylab("normalized read count")
  # ggtitle(tools::file_path_sans_ext(basename(bw_file)))
  # add p-values if multiple groups
  if (length(unique(results$Set)) > 1) {
    comparisons <- combn(levels(results$Set), 2, simplify = FALSE)
    # ####### Auto buffering from Leon #######
    ymax <- max(results$bwReads, na.rm = TRUE)
    ybase <- ymax * 1.1
    n_comp <- length(comparisons)
    offsets <- seq(0, 0.1 * ymax *(n_comp-1), length.out = n_comp)
    y_positions <- ybase + offsets
    # ########################################
    if (show_pval_numeric) {
      pvals <- sapply(comparisons, function(grp) {
        g1 <- results$bwReads[results$Set == grp[1]]
        g2 <- results$bwReads[results$Set == grp[2]]
        wilcox.test(g1, g2)$p.value
      })
      annotations <- sprintf("p = %.2g", pvals)
      p <- p + ggsignif::geom_signif(comparisons = comparisons,
                                     annotations = annotations,
                                     tip_length = 0.03,
                                     textsize = 5,
                                     y_position = y_positions,
                                     map_signif_level = FALSE)
    } else {
      p <- p + ggsignif::geom_signif(comparisons = comparisons,
                                     map_signif_level = TRUE,
                                     tip_length = 0.03,
                                     textsize = 5,
                                     y_position = y_positions)
    }
  }

  if (boxcenter_wMean) {
    ltype <- if (boxcenter_noMedian) "solid" else "dashed"
    p <- p + stat_summary(geom = "errorbar",
                          fun.min = mean, fun = mean, fun.max = mean,
                          width = 0.4, linetype = ltype)
  }

  p
}

#' Read reference sequence lengths from a FASTA index file
#'
#' @param fasta_index Path to the .fai file (FASTA index)
#' @return Named integer vector of contig lengths
#' @examples
#' lengths <- read_length("genome.fa.fai")
read_length <- function(fasta_index) {
  if (!file.exists(fasta_index)) {
    stop("FASTA index file not found: ", fasta_index)
  }
  # Read the .fai file (columns: contig, length, offset, linebases, linewidth)
  contig_df <- read.delim(
    fasta_index,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(contig_df)[1:2] <- c("contig", "length")

  contig_lengths <- setNames(as.integer(contig_df$length), contig_df$contig)
  return(contig_lengths)
}

# ----------------- User Config -------------------

# setwd("/Users/antienna/Documents/USC/ChioloLab/Chiolo_DSBresponses_dm6")

# genome lengths, either read from library or read from fasta index file
seqlens <- seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)
# seqlens <- read_length("~/wdir/Cathy/reference_genomes/hg19.fa.fai")

# plotting colors
color_map_fill <- c("EU" = "darkgrey", "HC" = "#F8766D", "FHC" = "#1f78b4")

# inputs
window_sizes   <- c(100, 250, 2500)
bed_files      <- c("../example_files/bed_files/EU_42DSBs.bed", 
                    "../example_files/bed_files/HC_14DSBs.bed", 
                    "../example_files/bed_files/FHC_10DSBs.bed")
group_labels   <- c("EU", "HC", "FHC")

bw_files  <- list(
  c("../example_files/bw_files/END_seq_ExoVII_5hrVsUNT_A74A73_bwa_aligned_PE_SES_wdup.bw",
    "../example_files/bw_files/END_seq_HO_rep2_B23B24_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
    "../example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw"),
  "../example_files/bw_files/B26_H3K9me3_SNAP_ChIP_rep1_UNT_bwa_aligned_pe_dm6KMetStatPombe_CPM.bw",
  "../example_files/bw_files/H3K4me3_IP_SRR392950_bwa_nodup_CPM.bw",
  "../example_files/bw_files/H3K27me3_SRR3452732_IP_bwa_realigned_pe_nodup_CPM.bw",
  "../example_files/bw_files/SRX4375842_IP2_RNAPOLII_trimmed_bwa.bw"
)
title_labels   <- c("ENDseq_log2r_A40A39rep1_A74A73rep3_B24B23rep4_wdup",
                    "H3K9me3_SNAP_ChIP_rep1_CPM", "H3K4me3_SRR392950_CPM", "H3K27me3_SRR3452731_CPM", "RNAPolII")

# settings
funct              <- "mean"
avg_btw_bw         <- TRUE
boxcenter_wMean    <- TRUE
boxcenter_noMedian <- TRUE
show_pval_numeric  <- TRUE

# ----------------- Main -------------------------

dir.create("../plots", showWarnings = FALSE)

all_plots <- list()

for (window_size in window_sizes) {
  for (i in seq_along(bw_files)) {

    bw_group    <- bw_files[[i]]
    title_label <- title_labels[i]

    # compute signal for each bw replicate/file in group
    res_list <- lapply(bw_group, function(bw) {
      compute_cov_chipseq(bed_files, bw, window_size, group_labels, funct, seqlens)
    })

    # bind and average (across bw files) per row
    res_df <- Reduce(function(x, y) {
      stopifnot(all(x$Set == y$Set), all(x$label == y$label))
      cbind(x, y$bwReads)
    }, res_list)

    # columns: Set, label, bwReads (first file), bwReads (2nd), ...
    if (ncol(res_df) == 3) {
      avg_reads <- as.numeric(res_df[, "bwReads"])
    } else {
      avg_reads <- rowMeans(res_df[, -(1:2)], na.rm = TRUE)
    }

    avgs <- data.frame(Set = res_df$Set,
                       label = res_df$label,
                       bwReads = avg_reads,
                       stringsAsFactors = FALSE)

    plot_title <- paste0(title_label, "_", window_size)
    file_stub  <- paste("boxplot",
                        title_label,
                        paste(group_labels, collapse = "_"),
                        paste0(window_size, "win"),
                        sep = "_")

    suffix <- paste0(
      if (avg_btw_bw) "_AVGbtwBW" else "",
      if (boxcenter_wMean && boxcenter_noMedian) "_asMean" else if (boxcenter_wMean) "_wMean" else if (boxcenter_noMedian) "_noMedian" else "",
      "_classic_scaled_",
      if (show_pval_numeric) "wPval" else "wPstar",
      ".pdf"
    )
    outfile <- file.path("plots", paste0(file_stub, suffix))

    timestamp_message("Writing ", outfile)
    pdf(outfile, width = 3.5, height = 6)
    p <- plot_chipseq_boxplot(avgs,
                              title = plot_title,
                              group_labels = group_labels,
                              color_map_fill = color_map_fill,
                              boxcenter_wMean = boxcenter_wMean,
                              boxcenter_noMedian = boxcenter_noMedian,
                              show_pval_numeric = show_pval_numeric)
    print(p)
    dev.off()

    all_plots[[paste0(plot_title, "_", window_size)]] <- p
  }
}

closeAllConnections()

invisible(all_plots)
