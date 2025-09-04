#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(dplyr)
})

# ---------------- Functions ----------------

#' Barplot of site signals
#'
#' @param df dataframe (rows=sites, cols=signals)
#' @param order_col index of column used to order sites (default: 4)
#' @param value_cols vector of column indices for plotting
#' @param normalize_percent logical, if TRUE scale to % of max in each column
#' @param site_label_cols optional indices for chr/start/end to form site labels
#' @return ggplot object
plot_site_barplot <- function(df,
                              order_col = 4,
                              value_cols = c(4,5),
                              normalize_percent = FALSE,
                              site_label_cols = 1:3) {
  df[is.na(df)] <- 0
  
  # site labels
  if (!is.null(site_label_cols)) {
    df$site_id <- apply(df[, site_label_cols], 1, paste, collapse = "_")
  } else {
    df$site_id <- seq_len(nrow(df))
  }
  
  # order sites by chosen column
  df <- df[order(df[[order_col]], decreasing = FALSE), ]
  df$site_id <- factor(df$site_id, levels = df$site_id)
  
  # subset and reshape
  df_sub <- df[, c(ncol(df), value_cols)]
  df_long <- melt(df_sub, id.vars = "site_id", variable.name = "bw_file", value.name = "signal")
  
  # normalize if needed
  if (normalize_percent) {
    df_long <- df_long %>%
      group_by(bw_file) %>%
      mutate(signal = 100 * signal / max(signal, na.rm = TRUE)) %>%
      ungroup()
    ylab <- "Signal (% of max)"
  } else {
    ylab <- "Signal (absolute reads)"
  }
  
  # sort bars so tallest is drawn last (back)
  # df_long <- df_long[order(df_long$signal), ]
  df_long <- df_long %>%
    group_by(site_id) %>%
    arrange(desc(signal), .by_group = TRUE) %>%
    mutate(bw_file = factor(bw_file, levels = bw_file)) %>%
    ungroup()
  
  # plot
  p <- ggplot(df_long, aes(x = site_id, y = signal, fill = bw_file)) +
    geom_bar(stat = "identity", position = position_identity(), alpha = 1, width = 0.8) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 2)) +
    ylab(ylab) +
    xlab("Sites") +
    ggtitle("Signal per site")
  
  return(p)
}

# ---------------- usage ----------------
setwd("/Users/antienna/Documents/USC/ChioloLab/chiololab_DSBresponses_dm6")

df <- read.delim("example_files/output/signal_matrix_MACS.tsv", stringsAsFactors = FALSE)

p <- plot_site_barplot(df, order_col = 4, value_cols = c(4,5,6,7), normalize_percent = FALSE)
ggsave("example_files/output/barplot_absolute_allAsiSI_MACS.pdf", p, width = 45, height = 10)

p <- plot_site_barplot(df, order_col = 4, value_cols = c(4,5,6,7), normalize_percent = TRUE)
ggsave("example_files/output/barplot_percentage_allAsiSI_MACS.pdf", p, width = 45, height = 10)

