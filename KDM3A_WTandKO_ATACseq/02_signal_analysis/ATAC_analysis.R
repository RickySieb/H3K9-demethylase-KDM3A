##========================================================================================
## loading required packages
##========================================================================================
# Time tracking
start_time <- Sys.time()
cat("Starting to load packages", "\n")
# Source the packages for the setup
library(parallel)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(tidyr)

##========================================================================================
## Enable parallel processing
##========================================================================================
# Time tracking
start_time <- Sys.time()
cat("Enabling parallel processing", "\n")
# Source and set the packages for the processing
library(future)
plan("multicore", workers = 15)
# Time tracking
end_time <- Sys.time()
cat("Finished enabling paralell processing, it took", (end_time - start_time), "seconds", "\n")
cat("\n")

##========================================================================================
## Wrapper function for Rscript compatability
##========================================================================================
# Define a wrapper function for Diffbind object handling
manage_Robjects <- function(object_name, output_dir, create_fn) {
  # log the current time
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S ", "\n"))
  start_time <- Sys.time()
  # Construct the output path
  output_path <- file.path(output_dir, paste0(object_name, ".rds"))
  
  # Check if the object already exists in the environment
  if (!exists(object_name, envir = .GlobalEnv)) {
    
    # Check if the file exists
    if (file.exists(output_path)) {
      # Load the object from file
      assign(object_name, readRDS(output_path), envir = .GlobalEnv)
      cat(object_name, "loaded from file:", output_path, "\n")
    } else {
      # File does not exist, create the object using the provided function
      new_object <- create_fn()
      assign(object_name, new_object, envir = .GlobalEnv)
      
      # Save the object to the file
      saveRDS(new_object, output_path)
      cat(object_name, "created and saved to:", output_path, "\n")
    }
    
  } else {
    cat(object_name, "already exists in the environment. Skipping unnecessary processing.\n")
  }
  end_time <- Sys.time()
  cat("Finished finding/loading/processing ->", object_name,", it took", (end_time - start_time), "seconds", "\n")
  cat("\n")
}

##========================================================================================
## Accept command-line arguments
##========================================================================================
# Time tracking
start_time <- Sys.time()
cat("Integrating command line arguments", "\n")
# command line argument integration
args <- commandArgs(trailingOnly = TRUE) 
if (length(args) < 2) { 
  stop("Usage: Rscript ...script.R <data_dir> <output_dir>")
}

data_dir    <- args[1]
output_dir  <- args[2]

## Ensure the output directory exists
if (!dir.exists(output_dir)) { 
  cat("Creating output directory:", output_dir, "\n") 
  dir.create(output_dir, recursive = TRUE)
}
# Time tracking
end_time <- Sys.time()
cat("Finished integrating command line arguments, it took", (end_time - start_time), "seconds", "\n")
cat("\n")

##========================================================================================
## Interactive session data paths — edit before running interactively
##========================================================================================
if (interactive()) {
  data_dir   <- "/path/to/KDM3A/mm10/"
  output_dir <- "/path/to/KDM3A/mm10/R_analyses/"
}

##========================================================================================
## analysis - load KDM3A peaksets
##========================================================================================
# BED files produced by the bedtools intersect/unique step (see chipseq/03_bedtools_peaks/)
BED_CTRL_SPECIFIC  <- file.path(data_dir, "KDM3A/bedtools/KDM3A_CTRL_specific.bed")
BED_LPS_SPECIFIC   <- file.path(data_dir, "KDM3A/bedtools/KDM3A_LPS_specific.bed")
BED_CONSTITUTIVE   <- file.path(data_dir, "KDM3A/bedtools/intersect_KDM3A_CTRL_LPS.bed")

create_combined_peaks <- function() {
  df1 <- read.delim(BED_CTRL_SPECIFIC, header = FALSE, stringsAsFactors = FALSE)
  df2 <- read.delim(BED_LPS_SPECIFIC,  header = FALSE, stringsAsFactors = FALSE)
  df3 <- read.delim(BED_CONSTITUTIVE,  header = FALSE, stringsAsFactors = FALSE)
  
  df1$source <- "CTRL_specific"
  df2$source <- "LPS_specific"
  df3$source <- "constitutive"
  
  df_combined <- rbind(df1, df2, df3)
  
  peaks <- GRanges(
    seqnames = df_combined[[1]],
    ranges = IRanges(
      start = as.integer(df_combined[[2]]),
      end   = as.integer(df_combined[[3]])
    ),
    source = df_combined$source
  )
}

cat("Finding/loading/processing combined peaks into GRanges object\n")
manage_Robjects("peaks", output_dir, create_combined_peaks)

##========================================================================================
## analysis from BW
##========================================================================================
# ATAC BigWig files — edit paths and sample names to match your data directory
ATAC_BW_DIR <- file.path(data_dir, "ATAC/BigWigs")

create_ATAC_signal_tibble <- function() {
  # Vector of bigWig files
  bw_files <- c(
    file.path(ATAC_BW_DIR, "ATAC_WT_S4-trimmed.sorted.dups_marked.filtered.bw"),
    file.path(ATAC_BW_DIR, "ATAC_WT_S80-trimmed.sorted.dups_marked.filtered.bw"),
    file.path(ATAC_BW_DIR, "ATAC_KO_S1-trimmed.sorted.dups_marked.filtered.bw"),
    file.path(ATAC_BW_DIR, "ATAC_KO_S77-trimmed.sorted.dups_marked.filtered.bw"),
    file.path(ATAC_BW_DIR, "ATAC_LPSWT_S3-trimmed.sorted.dups_marked.filtered.bw"),
    file.path(ATAC_BW_DIR, "ATAC_LPSWT_S79-trimmed.sorted.dups_marked.filtered.bw"),
    file.path(ATAC_BW_DIR, "ATAC_LPSKO_S2-trimmed.sorted.dups_marked.filtered.bw"),
    file.path(ATAC_BW_DIR, "ATAC_LPSKO_S78-trimmed.sorted.dups_marked.filtered.bw")
  )
  
  sample_names <- c("CTRL_1", "CTRL_2", "KO_1", "KO_2", "CTRLLPS_1", "CTRLLPS_2", "KOLPS_1", "KOLPS_2")
  
  # Inner function for one file
  get_coverage_per_peak <- function(bw_file, peaks) {
    coverage <- import(bw_file, which = peaks, as = "RleList")
    mean_cov <- numeric(length(peaks))
    for (i in seq_along(peaks)) {
      chr <- as.character(seqnames(peaks)[i])
      start <- start(peaks)[i]
      end <- end(peaks)[i]
      if (chr %in% names(coverage)) {
        region_cov <- coverage[[chr]][start:end]
        mean_cov[i] <- mean(region_cov, na.rm = TRUE)
      } else {
        mean_cov[i] <- NA
      }
    }
    return(mean_cov)
  }
  
  # Apply to all files
  coverage_matrix <- sapply(bw_files, get_coverage_per_peak, peaks = peaks)
  colnames(coverage_matrix) <- sample_names
  
  peak_info <- as.data.frame(peaks)
  
  combined <- cbind(peak_info, coverage_matrix)
  
  ATAC_signal_tibble <- combined %>%
    pivot_longer(
      cols = all_of(sample_names),
      names_to = "sample",
      values_to = "signal"
    ) %>%
    mutate(
      group = case_when(
        grepl("^CTRL_", sample) ~ "CTRL",
        grepl("^KO_", sample) ~ "KO",
        grepl("^CTRLLPS_", sample) ~ "CTRLLPS",
        grepl("^KOLPS_", sample) ~ "KOLPS",
        TRUE ~ NA_character_
      ),
      group = factor(group, levels = c("CTRL", "KO", "CTRLLPS", "KOLPS")),
      source = factor(source, levels = c("CTRL_specific", "LPS_specific", "constitutive"))
    )
  
  return(ATAC_signal_tibble)
}

cat("Finding/loading/processing long_df from bigWig coverage\n")
manage_Robjects("ATAC_signal_tibble", output_dir, create_ATAC_signal_tibble)

##========================================================================================
## analysis from BW - figure ATAC
##========================================================================================
if (interactive()) {
  long_df_noNA <- na.omit(ATAC_signal_tibble)
  ggplot(long_df_noNA, aes(x = group, y = log2(signal + 1), fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ source, nrow = 1) +
    ylim(0, 6) +
    labs(
      title = "ATAC at KDM3A",
      y = "log2(Tag density + 1)", x = "Sample Group"
    ) +
    # Define colors here
    scale_fill_manual(
      values = c(
        "CTRL" = "white",     # blue
        "KO" = "#1f77b4",       # orange
        "CTRLLPS" = "white",  # green
        "KOLPS" = "#1f77b4"     # red
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 14),
      panel.border = element_blank(),  # No border box
      axis.line = element_line(color = "black", linewidth = 0.5),  # Add axis lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white")
    )
}

if (interactive()) {
  # Test all pairwise comparisons within each peak set
  stat_test <- long_df_noNA %>%
    group_by(source) %>%
    wilcox_test(signal ~ group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
}
