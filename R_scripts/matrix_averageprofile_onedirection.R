library(ggplot2)
library(zoo)

# ---- Parse header line to get metadata ----
parse_matrix_metadata <- function(file_path) {
  lines <- readLines(file_path, n = 3)
  meta_line <- gsub("^#+", "", lines[2])
  parts <- unlist(strsplit(meta_line, "\\t"))
  info <- list()
  for (p in parts) {
    if (grepl(":", p)) {
      kv <- unlist(strsplit(p, ":"))
      key <- trimws(kv[1])
      value <- as.numeric(trimws(kv[2]))
      info[[key]] <- value
    }
  }
  return(list(upstream = info[["upstream"]],
              downstream = info[["downstream"]],
              binsize = info[["bin size"]]))
}

# ---- Load matrix excluding header ----
load_matrix_body <- function(file_path) {
  data <- read.table(file_path, skip = 3, header = FALSE)
  return(as.matrix(data))
}

# ---- Process full matrix into (n_replicates x bins) ----
process_matrix <- function(mat, upstream, downstream, binsize, drop_idx = NULL) {
  n_bins <- (upstream + downstream) / binsize
  total_cols <- ncol(mat)
  n_reps <- total_cols / n_bins
  stopifnot(total_cols %% n_bins == 0)
  
  # Reshape into replicates
  reshaped <- matrix(NA, nrow = n_reps, ncol = n_bins)
  for (i in 1:n_reps) {
    reshaped[i, ] <- colMeans(mat[ , ((i-1)*n_bins + 1):(i*n_bins)], na.rm = TRUE)
  }
  
  # Drop specified replicates
  if (!is.null(drop_idx)) {
    reshaped <- reshaped[-drop_idx, , drop = FALSE]
  }
  
  return(reshaped)
}


# ---- Flip left matrix and average with right ----
merge_left_right <- function(left_mat, right_mat) {
  flipped_left <- left_mat[, ncol(left_mat):1]
  merged <- (flipped_left + right_mat) / 2
  avg <- colMeans(merged)
  sd  <- apply(merged, 2, sd)
  return(list(avg = avg, sd = sd))
}

# ---- Plot all merged profiles ----
plot_multiple_profiles <- function(results_list, binsize, upstream, group_labels, color_mapping, 
                                   title = "Merged One-Directional Profiles", zoom_above_zero = FALSE) {
  df <- data.frame(Position = numeric(), Average = numeric(), SD = numeric(), Group = character())
  
  for (i in seq_along(results_list)) {
    res <- results_list[[i]]
    pos <- seq(0, by = binsize, length.out = length(res$avg))
    df <- rbind(df, data.frame(Position = pos, Average = res$avg, SD = res$sd, Group = group_labels[i]))
  }
  
  p <- ggplot(df, aes(x = Position, y = Average, color = Group, fill = Group)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = Average - SD, ymax = Average + SD), alpha = 0.3, linetype = 0) +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme_classic() +
    labs(x = "Distance from site (bp)", y = "Signal", title = title)
  
  if (zoom_above_zero) {
    p <- p +
      coord_cartesian(ylim = c(0, NA)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
  } else {
    p <- p +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
  }
  
  
  return(p)
}

# ---- Main Wrapper ----
run_profile_plot <- function(matrix_files, group_labels, color_mapping, smooth_bp = 0,
                               zoom_above_zero = TRUE, title = "Merged One-Directional Profiles", drop_reps = list()) {
  results_list <- list()
  
  for (i in seq_along(matrix_files)) {
    left_file <- matrix_files[[i]][1]
    right_file <- matrix_files[[i]][2]
    
    left_meta <- parse_matrix_metadata(left_file)
    right_meta <- parse_matrix_metadata(right_file)
    
    if (left_meta$upstream != right_meta$downstream) {
      stop(sprintf("Mismatch: left upstream (%d) != right downstream (%d)", left_meta$upstream, right_meta$downstream))
    }
    if (left_meta$binsize != right_meta$binsize) {
      stop(sprintf("Mismatch: bin size (%d vs %d)", left_meta$binsize, right_meta$binsize))
    }
    
    message(sprintf("Processing: %s vs %s", basename(left_file), basename(right_file)))
    message(sprintf("  Binsize = %d bp | Upstream = %d bp | Downstream = %d bp", 
                    left_meta$binsize, left_meta$upstream, left_meta$downstream))
    if (smooth_bp > 0) message(sprintf("  Smoothing window = %d bp", smooth_bp))
    
    left_mat <- load_matrix_body(left_file)
    right_mat <- load_matrix_body(right_file)
    
    # left_proc <- process_matrix(left_mat, left_meta$upstream, left_meta$downstream, left_meta$binsize)
    # right_proc <- process_matrix(right_mat, right_meta$upstream, right_meta$downstream, right_meta$binsize)
    left_proc <- process_matrix(left_mat, 
                                left_meta$upstream, left_meta$downstream, left_meta$binsize, 
                                drop_idx = drop_reps[[i]])
    right_proc <- process_matrix(right_mat, 
                                 right_meta$upstream, right_meta$downstream, right_meta$binsize, 
                                 drop_idx = drop_reps[[i]])
    
    
    merged <- merge_left_right(left_proc, right_proc)
    
    if (smooth_bp > 0 && smooth_bp > left_meta$binsize) {
      smooth_n <- round(smooth_bp / left_meta$binsize)
      merged$avg <- rollapply(merged$avg, width = smooth_n, FUN = mean, fill = "extend", align = "center")
    }
    
    results_list[[i]] <- merged
  }
  
  plot <- plot_multiple_profiles(results_list, left_meta$binsize, left_meta$upstream, group_labels, color_mapping, zoom_above_zero = zoom_above_zero, title = title)
  return(plot)
}




# ---- Usage ----
setwd("/Users/antienna/Documents/USC/ChioloLab/scripts/oneDirProfile")
matrix_files <- list(
  c("pH2Av_TRTvsUNT_SES_3reps_across55EU_-500kb.tab", "pH2Av_TRTvsUNT_SES_3reps_across55EU_+500kb.tab"),
  c("pH2Av_TRTvsUNT_SES_3reps_across15HC_-500kb.tab", "pH2Av_TRTvsUNT_SES_3reps_across15HC_+500kb.tab")
)
drop_reps <- list(c(3), c(3)) # Drop 3rd replicate from both groups
group_labels <- c("55EU", "15HC")
color_mapping <- c("55EU" = "dimgrey", "15HC" = "#F8766D")
smooth_bp <- 5000
zoom_above_zero <- TRUE
title <- "pH2Av avg profile, two sides merged"

# Generate plot
p <- run_profile_plot(matrix_files, group_labels, color_mapping, 
                      smooth_bp, zoom_above_zero, title, drop_reps)
p
ggsave("one_directional_zoomAboveZero.pdf", p, width = 6, height = 4)

