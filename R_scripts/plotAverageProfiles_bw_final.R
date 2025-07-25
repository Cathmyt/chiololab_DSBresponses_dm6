#!/usr/bin/env Rscript

# ============================================================
# Average Profile Plotting of ChIP signal from bigWig(s)
# ============================================================
# This script:
#   1) Extends and bins BED file regions (e.g., ±5kb around cut sites).
#   2) Extracts bigWig signal over these binned regions.
#   3) Computes average signal profiles with optional smoothing.
#   4) Plots profiles with optional error bars (SD/SEM).
#
# Dependencies: BSgenome.Dmelanogaster.UCSC.dm6, rtracklayer,
#               GenomicRanges, ggplot2, zoo, parallel
# ============================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(zoo)
  library(parallel)
})

# --------------------- Utils ---------------------

timestamp_message <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

# ----------------- Core Functions ----------------

#' Read reference genome lengths
#'
#' @param genome Either a BSgenome object or a FASTA index (.fai).
#' @return Named integer vector of chromosome sizes.
read_genome_lengths <- function(genome) {
  if (is.character(genome) && file.exists(genome)) {
    df <- read.delim(genome, header = FALSE)
    return(setNames(as.integer(df[[2]]), df[[1]]))
  } else if (inherits(genome, "BSgenome")) {
    return(seqlengths(genome))
  } else {
    stop("Invalid genome input. Must be BSgenome object or .fai file.")
  }
}

#' Extend BED regions by upstream/downstream.
#'
#' @param bed_file Path to BED file.
#' @param upstream Bases upstream to add.
#' @param downstream Bases downstream to add.
#' @param seqlens Chromosome lengths for clipping.
#' @return GRanges object.
tidy_import_bed <- function(bed_file, upstream, downstream, seqlens) {
  timestamp_message("Importing BED: ", basename(bed_file))
  bed <- import.bed(bed_file)
  bed_extended <- resize(resize(bed, width(bed) + upstream, fix = "end"),
                         width(bed) + upstream + downstream, fix = "start")
  start(bed_extended) <- pmax(start(bed_extended), 1)
  end(bed_extended) <- pmin(end(bed_extended), seqlens[as.character(seqnames(bed_extended))])
  bed_extended
}


#' Bin BED regions into fixed number of bins
#'
#' @param bed GRanges object.
#' @param nbins Number of bins per region.
#' @return GRanges object with tiled regions.
bin_bed <- function(bed, nbins) {
  return(unlist(tile(bed, nbins)))
}

#' Compute coverage for binned regions from bigWig
#'
#' @param regions GRanges of binned regions.
#' @param bw_file Path to bigWig file.
#' @param nbins Number of bins.
#' @param plotEBwithinBed Optional: "SD" or "SEM" for error bars.
#' @return Vector or list of coverage values (and SD if requested).
compute_bw_coverage <- function(regions, bw_file, nbins, plotEBwithinBed = NULL) {
  bw <- import.bw(bw_file, which = regions)
  coverage_vals <- sapply(seq_along(regions), function(i) {
    bin_region <- regions[i]
    cov_values <- bw[overlapsAny(bw, bin_region)]
    mean(score(cov_values), na.rm = TRUE)
  })
  
  cov_matrix <- matrix(coverage_vals, ncol = nbins, byrow = TRUE)
  avg_coverage <- colMeans(cov_matrix, na.rm = TRUE)
  
  if (!is.null(plotEBwithinBed)) {
    if (plotEBwithinBed == "SD") cov_sd <- apply(cov_matrix, 2, sd, na.rm = TRUE)
    else if (plotEBwithinBed == "SEM") {cov_sd <- apply(cov_matrix, 2, sd, na.rm = TRUE); cov_sd <- cov_sd / sqrt(nrow(cov_matrix))}
    else {message(plotEBwithinBed, "cannot be plotted as errorbars.")}
    return(list(avg_coverage, cov_sd))
  }
  return(avg_coverage)
}


#' Plot average coverage profile
#'
#' @param df Data frame with Position, Coverage, SD, and Group columns.
#' @param binsize Bin size in bp.
#' @param upstream Upstream extension.
#' @param downstream Downstream extension.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param title_label Plot title.
#' @param zoom_above_zero Logical: Zoom y-axis to >= 0.
#' @param y_upper_lim Optional y-axis upper limit.
#' @param filling Logical: Use ribbons under the line plot.
#' @return ggplot object.
plot_profile <- function(df, binsize, upstream, downstream, x_label, y_label, title_label, 
                         zoom_above_zero, plotSDbtwBed, y_upper_lim, filling) {
  color_map <- c("EU" = "dimgrey", "HC" = "#F8766D", "FHC" = "#1f78b4", "AsiSI" = "#fcba03")
  df$Group <- factor(df$Group, levels = c("EU", "FHC", "HC", "AsiSI"))
  
  p <- ggplot()
  
  for (group in levels(df$Group)) {
    df_subset <- df[df$Group == group, ]
    if (filling) p <- p + geom_ribbon(data = df_subset, aes(x = Position, ymin = 0, ymax = Coverage, fill = Group), alpha = 0.3)
    p <- p + geom_errorbar(data = df_subset, aes(x = Position, ymin = Coverage - SD, ymax = Coverage + SD, color = Group), width = 0, alpha = 0.6) +
      geom_line(data = df_subset, aes(x = Position, y = Coverage, color = Group), linewidth = 1.5)
  }
  
  p <- p +
    scale_colour_manual(values = color_map) +
    scale_fill_manual(values = color_map) +
    theme_classic() +
    xlab(x_label) +
    ylab(y_label) +
    ggtitle(title_label) +
    theme(axis.text = element_text(size = 12, color = "black"), text = element_text(size = 14, color = "black"))
  if (zoom_above_zero) {
    if (is.na(y_upper_lim)) {
      p <- p + coord_cartesian(ylim = c(0, NA)) + 
        scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) 
    } else {
      p <- p + coord_cartesian(ylim = c(0, y_upper_lim)) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), breaks = seq(0.0, y_upper_lim, by = 0.4))
        # scale_y_continuous(expand = c(0,0))
    } } else { p <- p + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) }
  return(p)
}



# --------------------- Main Function ---------------------
plot_average_profile_bw <- function(bed_files, bw_files, binsize, smoothsize, upstream, downstream, 
                                    zoom_above_zero, plotSDbtwBed, plotEBwithinBed, smooth_cut_for_strandsep, y_upper_lim, filling,
                                    x_label, y_label, title_label, seqlens) {
  nbins <- ceiling((upstream + downstream) / binsize)
  positions <- seq(-upstream + binsize/2, downstream - binsize/2, binsize)
  df <- data.frame(Position = numeric(), Coverage = numeric(), SD = numeric(), Group = character())
  
  for (bed_file in bed_files) {
    wins <- tidy_import_bed(bed_file, upstream, downstream, seqlens)
    binned_wins <- bin_bed(wins, nbins)

    quantifications <- parallel::mclapply(bw_files, function(bw) {
      compute_bw_coverage(binned_wins, bw, nbins, plotEBwithinBed)
    })
    
    if (!is.null(plotEBwithinBed)) {
      sd_quants <- quantifications[[1]][[2]]
      quantifications <- quantifications[[1]][1]
    }
    
    combined <- do.call(cbind, quantifications)
    avg_profile <- rowMeans(combined, na.rm = TRUE)
    sd_profile <- apply(combined, 1, sd, na.rm = TRUE)
    
    if (!is.null(plotEBwithinBed)) sd_profile <- sd_quants
    if (smoothsize > binsize) {
      smooth_nbins <- ceiling(smoothsize / binsize)
    } else {message("Smoothing window is smaller than binsize, smoothing is not applied.")}
    if (!smooth_cut_for_strandsep) {
      avg_profile <- rollapply(avg_profile, smooth_nbins, mean, partial = TRUE, align = "center")
      sd_profile <- rollapply(sd_profile, smooth_nbins, mean, partial = TRUE, align = "center")
    } else {
      message("Reads close to cut sites are not smoothed.")
      mid_point <- ceiling(length(avg_profile)/2)
      avg_profile_p1 <- rollapply(avg_profile[1:mid_point], smooth_nbins, mean, partial = T, align = "center")
      avg_profile_p2 <- rollapply(avg_profile[(mid_point+1):length(avg_profile)], smooth_nbins, mean, partial = T, align = "center")
      avg_profile <- c(avg_profile_p1, avg_profile_p2)
      sd_profile_p1 <- rollapply(sd_profile[1:mid_point], smooth_nbins, mean, partial = T, align = "center")
      sd_profile_p2 <- rollapply(sd_profile[(mid_point+1):length(sd_profile)], smooth_nbins, mean, partial = T, align = "center")
      sd_profile <- c(sd_profile_p1, sd_profile_p2)
    }
    
    bed_label <- regmatches(tools::file_path_sans_ext(basename(bed_file)), 
                            regexpr("(FHC|HC|EU|AsiSI)", tools::file_path_sans_ext(basename(bed_file))))
    df <- rbind(df, data.frame(Position = positions, Coverage = avg_profile, SD = sd_profile, Group = bed_label))
  }
  
  p <- plot_profile(df[-1,], binsize, upstream, downstream, x_label, y_label, title_label, zoom_above_zero, plotSDbtwBed, y_upper_lim, filling)
  return(p)
}


# ----------------- User Config -------------------
setwd("/Users/antienna/Documents/USC/ChioloLab/chiololab_DSBresponses_dm6/R_scripts")

# Reference genome: can be BSgenome object or FASTA .fai index
# seqlens <- read_genome_lengths("~/ref/hg19.fa.fai")
seqlens <- read_genome_lengths(BSgenome.Dmelanogaster.UCSC.dm6)

# Parameters
upstream <- 5000; downstream <- 5000; binsize <- 50; smoothsize <- 250
y_label <- "normalized read counts"; x_label <- "Positions relative to cut sites (bp)"
zoom_above_zero <- TRUE; smooth_cut_for_strandsep <- FALSE; y_upper_lim <- NA; filling <- TRUE
plotSDbtwBed <- FALSE; plotEBwithinBed <- NULL # NULL/"SD"/"SEM"

bed_files      <- c("../example_files/bed_files/EU_42DSBs.bed",
                    "../example_files/bed_files/HC_14DSBs.bed")
group_labels   <- c("EU", "HC")

bw_files  <- list(
  c("../example_files/bw_files/END_seq_ExoVII_5hrVsUNT_A74A73_bwa_aligned_PE_SES_wdup.bw",
    "../example_files/bw_files/END_seq_HO_rep2_B23B24_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw",
    "../example_files/bw_files/END_seq_rep1_A40A39_5hrVsUNT_bwa_aligned_PE_SES_wdup.bw")
)
title_labels   <- c("ENDseq_log2r_A40A39rep1_A74A73rep3_B24B23rep4_wdup")
  

dir.create("../plots", showWarnings = FALSE)
# ----------------- Main --------------------------
timestamp_message("Starting average profile computation...")

for (i in seq_along(bw_files)) {
  bw_files_i <- bw_files[[i]]
  title_label <- title_labels[i]
  pdf_title <- file.path("../plots", sprintf(
    "AverageProfile_%s_%dbin_%dsmooth_%dwin_%s.pdf",
    title_label, binsize, smoothsize, upstream + downstream, paste(group_labels, collapse = "_")
  ))
  if (!is.null(plotEBwithinBed)) pdf_title <- paste(strsplit(pdf_title, ".pdf")[[1]], "_", plotEBwithinBed, "asEB", ".pdf", sep = "")
  pdf(pdf_title, width = 5, height = 4)
  print(plot_average_profile_bw(bed_files, bw_files_i, 
                                binsize, smoothsize, upstream, downstream, 
                                zoom_above_zero, plotSDbtwBed, plotEBwithinBed, 
                                smooth_cut_for_strandsep, y_upper_lim, filling, 
                                x_label, y_label, title_label, seqlens))
  dev.off()
}
closeAllConnections()
