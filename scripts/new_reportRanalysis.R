#!/usr/bin/env Rscript

#### Header ####
# Report Generation Script
# Author: Denis Rivard, modified by GitHub Copilot
# Description: 

library(here)
source(here("scripts", "utils.R"))
library(purrr)
library(tidyr)
library(dplyr)

# Load parallel libraries conditionally
if ("--parallel" %in% commandArgs(trailingOnly = TRUE)) {
  library(foreach)
  library(doParallel)
}

# Function to display help
show_help <- function() {
  cat("Usage: Rscript reportRanalysis.R [OPTIONS]\n\n")
  cat("Inputs:\n")
  cat("  --unaligned-kreports <DIR> Path to unaligned kreports folder (default: NULL)\n")
  cat("  --unaligned-bracken <DIR>  Path to unaligned bracken folder (default: NULL)\n")
  cat("  --nonhuman-kreports <DIR>  Path to nonhuman kreports folder (default: NULL)\n")
  cat("  --nonhuman-bracken <DIR>   Path to nonhuman bracken folder (default: NULL)\n")
  cat("  --runtime-dir <DIR>       Path to runtime folder for parsing runtime files (default: NULL)\n")
  cat("  --rrstats-dir <DIR>       Path to read statistics folder for parsing rrstats files (default: NULL)\n")
  cat("  --metadata-dir <DIR>      Path to directory containing sample metadata CSV files (default: NULL)\n")
  cat("Options:\n")
  cat("  --no-bracken              Skip bracken file processing\n")
  cat("  --subspecies              Include subspecies (S1, S2, S3) in addition to species (default: FALSE)\n")
  cat("  --params                  Add confidence_levels, minimum_hit_groups, and human_reads columns to long data (default: FALSE)\n")
  cat("  --parallel                Use parallel processing with cluster from global environment 'cl' (default: FALSE)\n")
  cat("  --output                  Save results to CSV files (default: FALSE)\n")
  cat("  --output-dir <DIR>         Output directory (default: outputs/full_run)\n")
  cat("Filtering:\n")
  cat("  --exclude-taxid <ID>       Exclude species with this taxonomy ID (e.g., 9606 for Homo sapiens)\n")
  cat("  --min-reads <N>           Minimum clade reads threshold (default: 0)\n")
  cat("  --top-n <N>               The n number of top species to count in batch frequency analysis (default: 3)\n")
  cat("  --minimizer-ratio <R>     Filter by minimum ratio of distinct_minimizers/cladeReads (default: NULL)\n")
  cat("  --minimizer-threshold <N> Filter by minimum distinct_minimizers threshold (default: NULL)\n")
}
# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Check for help request
if ("--help" %in% args || "-h" %in% args) {
  show_help()
  stop("Help requested. Execution stopped.")
}

options(scipen = 999)

# Configuration variables with default values
MIN_CLADE_READS <- 0
EXCLUDE_TAXID <- NULL  # Set to NULL by default (no exclusion)
TOP_SPECIES_N <- 3
PROCESS_BRACKEN <- TRUE   # Set to TRUE to enable bracken processing by default
SAVE_OUTPUT <- FALSE  # Set to FALSE to skip saving CSV files
ADD_PARAMS <- FALSE  # Set to FALSE to skip adding parameter columns
INCLUDE_SUBSPECIES <- FALSE  # Set to FALSE to exclude subspecies
MINIMIZER_RATIO <- NULL  # Set to NULL to skip minimizer ratio filtering
MINIMIZER_THRESHOLD <- NULL  # Set to NULL to skip minimizer threshold filtering
USE_PARALLEL <- FALSE  # Set to FALSE to skip parallel processing

# Default input directory paths
UNALIGNED_KREPORTS_DIR <- NULL
UNALIGNED_BRACKEN_DIR <- NULL
NONHUMAN_KREPORTS_DIR <- NULL
NONHUMAN_BRACKEN_DIR <- NULL
RUNTIME_DIR <- NULL
RRSTATS_DIR <- NULL
METADATA_DIR <- NULL
OUTPUT_DIR <- here("outputs", "new_samples")


run_as_script(here("scripts", "output_processing.R"), 
              "--unaligned-kreports", "/Volumes/Promise/Denis_Summer_2025/k2RNAsq_pipeline/RNASeqPipeline/kraken2/c15m3stdk2dbdl_unaligned",
              "--unaligned-bracken", "/Volumes/Promise/Denis_Summer_2025/k2RNAsq_pipeline/RNASeqPipeline/bracken/c15m3stdk2dbdl_unaligned",
              "--nonhuman-kreports", "/Volumes/Promise/Denis_Summer_2025/k2RNAsq_pipeline/RNASeqPipeline/kraken2_nonhuman/c15m3stdk2dbdl_nonhuman",
              "--nonhuman-bracken", "/Volumes/Promise/Denis_Summer_2025/k2RNAsq_pipeline/RNASeqPipeline/bracken/c15m3stdk2dbdl_nonhuman",
              "--runtime-dir", "/Volumes/Promise/Denis_Summer_2025/k2RNAsq_pipeline/RNASeqPipeline/",
              "--rrstats-dir", "/Volumes/Promise/Denis_Summer_2025/k2RNAsq_pipeline/RNASeqPipeline/QC",
              "--parallel", "--subspecies", "--params", 
              "--output", "--output-dir", "/Volumes/Promise/Denis_Summer_2025/k2RNAsq_pipeline/RNASeqPipeline/Ranalysis/outputs")

species_df <- unaligned_results$species_list
run_as_script(here("scripts", "annotate_species.R"), "--df", "species_df", "--output", "annotated_species")

# data 1: metadata


# data 2: kraken2 results

create_correlation_analysis(
  data1 = unaligned_results$merged,
  data2 = nonhuman_results$merged,
  data1_name = "Unaligned",
  data2_name = "Non-Human",
  data1_col_suffix = "cladeReads",
  data2_col_suffix = "cladeReads",
  output_dir = here("reports")
)

host_contamination <- setdiff(unaligned_results$merged$name, nonhuman_results$merged$name)
host_contamination
host_contaminated_results <- unaligned_results$merged %>%
  filter(name %in% host_contamination)
host_contaminated_results[, c("name", "cladeReads_max")]

# Function to create correlation analysis
create_correlation_analysis <- function(data1, data2, data1_name, data2_name, data1_col_suffix, data2_col_suffix, output_dir) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get sample columns
  cols1 <- grep(paste0("_", data1_col_suffix, "$"), colnames(data1), value = TRUE)
  cols2 <- grep(paste0("_", data2_col_suffix, "$"), colnames(data2), value = TRUE)
  
  samples1 <- str_remove(cols1, paste0("_", data1_col_suffix, "$"))
  samples2 <- str_remove(cols2, paste0("_", data2_col_suffix, "$"))
  
  common_samples <- intersect(samples1, samples2)
  
  if (length(common_samples) == 0) {
    warning("No common samples found between datasets")
    return(NULL)
  }
  
  cat("Found", length(common_samples), "common samples for correlation analysis\n")
  
  start_time <- Sys.time()
  
  # Process samples sequentially (simplified for reportRanalysis)
  cat("Processing correlation analysis sequentially...\n")
  sample_results <- list()
  for (sample_id in common_samples) {
    col1 <- paste0(sample_id, "_", data1_col_suffix)
    col2 <- paste0(sample_id, "_", data2_col_suffix)
    
    if (col1 %in% colnames(data1) && col2 %in% colnames(data2)) {
      # Select data from both datasets
      data1_subset <- data1 %>% select(taxID, name, all_of(col1))
      data2_subset <- data2 %>% select(taxID, name, all_of(col2))
      
      # Rename the columns to generic names before joining
      colnames(data1_subset)[3] <- "reads1"
      colnames(data2_subset)[3] <- "reads2"
      
      sample_comparison <- inner_join(
        data1_subset,
        data2_subset,
        by = c("taxID", "name")
      ) %>%
        filter(reads1 > 0 & reads2 > 0)
      
      if (nrow(sample_comparison) > 10) {
        cor_coef <- cor(log10(sample_comparison$reads1 + 1), 
                        log10(sample_comparison$reads2 + 1), 
                        use = "complete.obs")
        
        sample_comparison <- sample_comparison %>%
          mutate(
            log_reads1 = log10(reads1 + 1),
            log_reads2 = log10(reads2 + 1),
            label = ifelse(abs(log_reads1 - log_reads2) > 2 | 
                             log_reads1 > 4 | log_reads2 > 4, name, NA)
          )
        # Create the correlation plot
        library(ggplot2)
        library(ggrepel)
        gg_correlation <- ggplot(sample_comparison, aes(x = log_reads1, y = log_reads2)) +
          geom_point(alpha = 0.6, size = 1.5) + 
          geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
          geom_text_repel(aes(label = label), na.rm = TRUE, size = 2.5, max.overlaps = 10) +
          labs(
            x = paste0("Log10(", data1_name, data1_col_suffix, " + 1)"),
            y = paste0("Log10(", data2_name, data1_col_suffix, " + 1)"),
            title = paste0("Correlation: ", data1_name, " vs ", data2_name, " (", sample_id, ")"),
            subtitle = paste0("Pearson r = ", round(cor_coef, 5), " (n = ", nrow(sample_comparison), " species)")
          ) + my_ggplot_theme + coord_fixed()
        
        sample_results[[length(sample_results) + 1]] <- list(
          sample = sample_id,
          pearson_r = cor_coef,
          n_species = nrow(sample_comparison),
          plot = gg_correlation
        )
      }
    }
  }
  
  end_time <- Sys.time()
  cat("Finished correlation analysis in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
  
  # Remove NULL entries
  sample_results <- sample_results[!sapply(sample_results, is.null)]
  
  # Extract correlation results and plots
  correlation_results <- tibble(
    sample = character(),
    pearson_r = numeric(),
    n_species = integer()
  )
  
  plots_list <- list()
  last_plot <- NULL
  
  for (result in sample_results) {
    if (!is.null(result)) {
      correlation_results <- bind_rows(
        correlation_results,
        tibble(
          sample = result$sample,
          pearson_r = result$pearson_r,
          n_species = result$n_species
        )
      )
      plots_list[[result$sample]] <- result$plot
      last_plot <- result$plot
    }
  }
  
  cat("Saving plots to PDF\n")
  # Save plots to PDF
  pdf(file.path(output_dir, paste0("correlation_", data1_name, "_vs_", data2_name, "_by_sample.pdf")), 
      width = 12, height = 8)
  for (plot in plots_list) {
    print(plot)
  }
  dev.off()
  cat("Correlation plots saved to PDF\n")
  
  # Print correlation results table
  cat("\n=== Correlation Analysis Results ===\n")
  correlation_results <- correlation_results %>%
    mutate(pearson_r = round(pearson_r, 5))
  print(correlation_results)
  cat("=====================================\n")
  
  return(list(common_samples = common_samples, correlation_results = correlation_results))
}


