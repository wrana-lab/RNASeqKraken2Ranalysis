#!/usr/bin/env Rscript

#### Header ####
# Kraken/Bracken Analysis Orchestration Script
# Author: Denis Rivard, modified by GitHub Copilot
# Date: August 19, 2025
# Description: Orchestrates downstream analysis of kraken2/bracken results including
#              processing kreports, species annotation with risk groups and HOMD data,
#              and correlation analysis between datasets
# Usage: Rscript new_reportRanalysis.R --input-dir <path> --proj-name <name> [options]

library(here)
source(here("scripts", "utils.R"))
library(parallel)
library(doParallel)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display help
show_help <- function() {
  cat("Usage: Rscript reportRanalysis.R [OPTIONS]\n\n")
  cat("Description:\n")
  cat("  Orchestrates comprehensive analysis of kraken2/bracken outputs including:\n")
  cat("  - Processing kreport and bracken files from multiple datasets\n")
  cat("  - Species annotation with risk groups and HOMD classifications\n")
  cat("  - Correlation analysis between unaligned and non-human datasets\n")
  cat("  - Generation of summary reports and visualizations\n\n")
  cat("Required Options:\n")
  cat("  --input-dir <DIR>         Path to input directory containing kreport/bracken folders\n")
  cat("  --proj-name <NAME>        Project name for output directory naming\n\n")
  cat("Report Options:\n")
  cat("  --all                     Include all report Analyses\n")
  cat("  --rt-stats                Include Runtime Analysis\n")
  cat("  --rr-stats                Include Read/Alignment Report Statistics\n")
  cat("  --kr-stats                Include Kraken Report Statistics\n")
  cat("  --sa-stats                Include Species Annotation Statistics\n")
  cat("  --pd-report               Pathogen Detection Analysis (RG3&4 calls) \n")
  cat("  --split-report            Split each per sample report into its own pdf document\n")
  cat("  --no-per-sample           Skip analysis and report generation for each sample\n")
  cat("  --no-corr                 Skip correlation analysis and report generation\n")
  cat("  --no-params               Skip adding parameter columns\n")
  cat("  --subspecies              Include subspecies in analysis\n")
  cat("  --top-n-plot <N>          Number of top species to show in plots (default: 50)\n\n")
  cat("Optional Options:\n")
  cat("  --output-base-dir <DIR>   Base output directory (default: outputs)\n")
  cat("  --databases-dir <DIR>     Path to databases directory (default: Ranalysis/databases)\n")
  cat("  --top-n-freq <N>          Number of top species for frequency analysis (default: 25)\n")
  cat("  --min-reads <N>           Minimum clade reads threshold (default: 0)\n")
  cat("  --exclude-taxid <ID>      Exclude species with this taxonomy ID\n")
  cat("  --minimizer-ratio <R>     Filter by minimum ratio of distinct_minimizers/cladeReads\n")
  cat("  --minimizer-threshold <N> Filter by minimum distinct_minimizers threshold\n")
  cat("  --cores <N>               Number of cores to use (default: detectCores())\n")
  cat("  --help, -h                Show this help message\n\n")
  cat("Optional Manual Output Processing Inputs:\n")
  cat("  --unaligned-kreports <DIR> Path to unaligned kreports folder (default: NULL)\n")
  cat("  --unaligned-bracken <DIR>  Path to unaligned bracken folder (default: NULL)\n")
  cat("  --nonhuman-kreports <DIR>  Path to nonhuman kreports folder (default: NULL)\n")
  cat("  --nonhuman-bracken <DIR>   Path to nonhuman bracken folder (default: NULL)\n")
  cat("  --runtime-dir <DIR>       Path to runtime folder for parsing runtime files (default: NULL)\n")
  cat("  --rrstats-dir <DIR>       Path to read statistics folder for parsing rrstats files (default: NULL)\n")
  cat("  --metadata-dir <DIR>      Path to directory containing sample metadata CSV files (default: NULL)\n")
  cat("Examples:\n")
  cat("  Rscript new_reportRanalysis.R --input-dir /path/to/data --proj-name my_project\n")
  cat("  Rscript new_reportRanalysis.R --input-dir /path/to/data --proj-name my_project --top-n 50 --cores 8\n")
  cat("  Rscript new_reportRanalysis.R --input-dir /path/to/data --proj-name my_project --no-corr --subspecies\n\n")
}

# Check for help request
if ("--help" %in% args || "-h" %in% args) {
  show_help()
  stop("Help requested. Execution stopped.")
} else {rm(show_help)}

# Configuration structure
create_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- list(
    INPUT_DIR = NULL,
    PROJ_NAME = NULL,
    OUTPUT_BASE_DIR = "outputs",
    DATABASES_DIR = here("databases"),
    TOP_N_FREQ = 25,
    TOP_N_PLOT = 50,
    MIN_READS = 0,
    EXCLUDE_TAXID = NULL,
    PERFORM_CORRELATION = TRUE,
    ADD_PARAMS = TRUE,
    INCLUDE_SUBSPECIES = FALSE,
    MINIMIZER_RATIO = NULL,
    MINIMIZER_THRESHOLD = NULL,
    CORES = detectCores(),
    # Report options
    ALL_REPORTS = FALSE,
    RT_STATS = FALSE,
    RR_STATS = FALSE,
    KR_STATS = FALSE,
    SA_STATS = FALSE,
    PD_STATS = FALSE,
    SPLIT_REPORT = FALSE,
    NO_PER_SAMPLE = FALSE,
    # Manual input paths
    UNALIGNED_KREPORTS_DIR = NULL,
    UNALIGNED_BRACKEN_DIR = NULL,
    NONHUMAN_KREPORTS_DIR = NULL,
    NONHUMAN_BRACKEN_DIR = NULL,
    RUNTIME_DIR = NULL,
    RRSTATS_DIR = NULL,
    METADATA_DIR = NULL
  )
  
  # Parse command line arguments
  if (length(args) > 0) {
    i <- 1
    while (i <= length(args)) {
      if (args[i] == "--input-dir" && i < length(args)) {
        config$INPUT_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--proj-name" && i < length(args)) {
        config$PROJ_NAME <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--output-base-dir" && i < length(args)) {
        config$OUTPUT_BASE_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--databases-dir" && i < length(args)) {
        config$DATABASES_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--top-n-freq" && i < length(args)) {
        config$TOP_N_FREQ <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--top-n-plot" && i < length(args)) {
        config$TOP_N_PLOT <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--min-reads" && i < length(args)) {
        config$MIN_READS <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--exclude-taxid" && i < length(args)) {
        config$EXCLUDE_TAXID <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--no-corr") {
        config$PERFORM_CORRELATION <- FALSE
        i <- i + 1
      } else if (args[i] == "--no-params") {
        config$ADD_PARAMS <- FALSE
        i <- i + 1
      } else if (args[i] == "--subspecies") {
        config$INCLUDE_SUBSPECIES <- TRUE
        i <- i + 1
      } else if (args[i] == "--minimizer-ratio" && i < length(args)) {
        config$MINIMIZER_RATIO <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--minimizer-threshold" && i < length(args)) {
        config$MINIMIZER_THRESHOLD <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--cores" && i < length(args)) {
        config$CORES <- as.numeric(args[i + 1])
        i <- i + 2
      # Report options
      } else if (args[i] == "--all") {
        config$ALL_REPORTS <- TRUE
        i <- i + 1
      } else if (args[i] == "--rt-stats") {
        config$RT_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--rr-stats") {
        config$RR_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--kr-stats") {
        config$KR_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--sa-stats") {
        config$SA_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--pd-stats") {
        config$PD_STATS <- TRUE
        i <- i + 1
      } else if (args[i] == "--split-report") {
        config$SPLIT_REPORT <- TRUE
        i <- i + 1
      } else if (args[i] == "--no-per-sample") {
        config$NO_PER_SAMPLE <- TRUE
        i <- i + 1
      # Manual input paths
      } else if (args[i] == "--unaligned-kreports" && i < length(args)) {
        config$UNALIGNED_KREPORTS_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--unaligned-bracken" && i < length(args)) {
        config$UNALIGNED_BRACKEN_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--nonhuman-kreports" && i < length(args)) {
        config$NONHUMAN_KREPORTS_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--nonhuman-bracken" && i < length(args)) {
        config$NONHUMAN_BRACKEN_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--runtime-dir" && i < length(args)) {
        config$RUNTIME_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--rrstats-dir" && i < length(args)) {
        config$RRSTATS_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--metadata-dir" && i < length(args)) {
        config$METADATA_DIR <- args[i + 1]
        i <- i + 2
      } else if (args[i] == "--help" || args[i] == "-h") {
        show_help()
        stop("Help requested. Execution stopped.")
      } else {
        i <- i + 1
      }
    }
  }
  
  # Validate required arguments
  if (is.null(config$INPUT_DIR)) {
    stop("ERROR: --input-dir is required.")
  }
  if (is.null(config$PROJ_NAME)) {
    stop("ERROR: --proj-name is required.")
  }
  if (!dir.exists(config$INPUT_DIR)) {
    stop("ERROR: Input directory does not exist: ", config$INPUT_DIR)
  }
  if (!dir.exists(config$DATABASES_DIR)) {
    stop("ERROR: Databases directory does not exist: ", config$DATABASES_DIR)
  }
  
  # Create output paths
  config$output_dir <- file.path(config$OUTPUT_BASE_DIR, config$PROJ_NAME, "outputs")
  config$reports_dir <- file.path(config$OUTPUT_BASE_DIR, config$PROJ_NAME, "reports")
  
  return(config)
}

# Setup parallel processing
setup_parallel <- function(cores) {
  cat("Setting up parallel processing with", cores, "cores\n")
  start_time <- Sys.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  assign("cl", cl, envir = .GlobalEnv)
  end_time <- Sys.time()
  cat("Parallel cluster initialized with", getDoParWorkers(), "workers, in ", round(difftime(end_time, start_time, units = "secs")), "seconds\n")
  return(cl)
}

# Build arguments for output_processing.R (updated to work with new version)
build_output_processing_args <- function(config) {
  # Use manual paths if provided, otherwise auto-detect
  unaligned_kreports_dir <- config$UNALIGNED_KREPORTS_DIR
  nonhuman_kreports_dir <- config$NONHUMAN_KREPORTS_DIR
  unaligned_bracken_dir <- config$UNALIGNED_BRACKEN_DIR
  nonhuman_bracken_dir <- config$NONHUMAN_BRACKEN_DIR
  
  # Auto-detect directories only if manual paths are not provided
  if (is.null(unaligned_kreports_dir) || is.null(nonhuman_kreports_dir)) {
    # Try different expected directory structures
    if (dir.exists(file.path(config$INPUT_DIR, "kraken2"))) {
      # Standard pipeline structure - look for subdirectories
      kraken2_dirs <- list.dirs(file.path(config$INPUT_DIR, "kraken2"), full.names = TRUE, recursive = FALSE)
      unaligned_dirs <- grep("unaligned", kraken2_dirs, value = TRUE)
      nonhuman_dirs <- grep("nonhuman", kraken2_dirs, value = TRUE)
      
      if (is.null(unaligned_kreports_dir) && length(unaligned_dirs) > 0) unaligned_kreports_dir <- unaligned_dirs[1]
      if (is.null(nonhuman_kreports_dir) && length(nonhuman_dirs) > 0) nonhuman_kreports_dir <- nonhuman_dirs[1]
    } else if (dir.exists(file.path(config$INPUT_DIR, "unaligned_kreports"))) {
      if (is.null(unaligned_kreports_dir)) unaligned_kreports_dir <- file.path(config$INPUT_DIR, "unaligned_kreports")
    } else if (dir.exists(file.path(config$INPUT_DIR, "kreports"))) {
      # Use general kreports directory for both
      if (is.null(unaligned_kreports_dir)) unaligned_kreports_dir <- file.path(config$INPUT_DIR, "kreports")
    }
    
    if (is.null(nonhuman_kreports_dir)) {
      if (dir.exists(file.path(config$INPUT_DIR, "kraken2"))) {
        # Look for nonhuman subdirectories within the main kraken2 folder
        kraken2_dirs <- list.dirs(file.path(config$INPUT_DIR, "kraken2"), full.names = TRUE, recursive = FALSE)
        nonhuman_dirs <- grep("nonhuman", kraken2_dirs, value = TRUE)
        if (length(nonhuman_dirs) > 0) nonhuman_kreports_dir <- nonhuman_dirs[1]
      } else if (dir.exists(file.path(config$INPUT_DIR, "kraken2_nonhuman"))) {
        # Fallback: check for legacy structure (kraken2_nonhuman folder) - deprecated
        cat("Warning: Using legacy kraken2_nonhuman directory structure. Consider updating to new kraken2/{param_label}_nonhuman structure.\n")
        kraken2_nonhuman_dirs <- list.dirs(file.path(config$INPUT_DIR, "kraken2_nonhuman"), full.names = TRUE, recursive = FALSE)
        nonhuman_dirs <- grep("nonhuman", kraken2_nonhuman_dirs, value = TRUE)
        if (length(nonhuman_dirs) > 0) nonhuman_kreports_dir <- nonhuman_dirs[1]
      } else if (dir.exists(file.path(config$INPUT_DIR, "nonhuman_kreports"))) {
        nonhuman_kreports_dir <- file.path(config$INPUT_DIR, "nonhuman_kreports")
      } else if (dir.exists(file.path(config$INPUT_DIR, "kreports"))) {
        nonhuman_kreports_dir <- file.path(config$INPUT_DIR, "kreports")
      }
    }
  }
  
  # Check if we found valid kreport directories
  if (is.null(unaligned_kreports_dir)) {
    stop("ERROR: Could not find unaligned kreports directory in: ", config$INPUT_DIR)
  }
  if (is.null(nonhuman_kreports_dir)) {
    stop("ERROR: Could not find nonhuman kreports directory in: ", config$INPUT_DIR)
  }
  
  args <- c(
    "--unaligned-kreports", unaligned_kreports_dir,
    "--nonhuman-kreports", nonhuman_kreports_dir,
    "--output-dir", config$output_dir,
    "--top-n", as.character(config$TOP_N_FREQ),
    "--min-reads", as.character(config$MIN_READS),
    "--output",
    "--parallel"
  )
  
  # Add bracken directories (use manual paths if provided, otherwise auto-detect)
  if (!is.null(unaligned_bracken_dir)) {
    args <- c(args, "--unaligned-bracken", unaligned_bracken_dir)
  } else if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
    # Standard pipeline structure - look for subdirectories
    bracken_dirs <- list.dirs(file.path(config$INPUT_DIR, "bracken"), full.names = TRUE, recursive = FALSE)
    unaligned_bracken_dirs <- grep("unaligned", bracken_dirs, value = TRUE)
    
    if (length(unaligned_bracken_dirs) > 0) {
      args <- c(args, "--unaligned-bracken", unaligned_bracken_dirs[1])
    }
  } else {
    # Try alternative directory structures
    if (dir.exists(file.path(config$INPUT_DIR, "unaligned_bracken"))) {
      args <- c(args, "--unaligned-bracken", file.path(config$INPUT_DIR, "unaligned_bracken"))
    }
  }
  
  if (!is.null(nonhuman_bracken_dir)) {
    args <- c(args, "--nonhuman-bracken", nonhuman_bracken_dir)
  } else if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
    # Standard pipeline structure - look for subdirectories
    bracken_dirs <- list.dirs(file.path(config$INPUT_DIR, "bracken"), full.names = TRUE, recursive = FALSE)
    nonhuman_bracken_dirs <- grep("nonhuman", bracken_dirs, value = TRUE)
    
    if (length(nonhuman_bracken_dirs) > 0) {
      args <- c(args, "--nonhuman-bracken", nonhuman_bracken_dirs[1])
    }
  } else {
    # Try alternative directory structures
    if (dir.exists(file.path(config$INPUT_DIR, "nonhuman_bracken"))) {
      args <- c(args, "--nonhuman-bracken", file.path(config$INPUT_DIR, "nonhuman_bracken"))
    }
    
    # Also try legacy bracken directory structure 
    if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
      bracken_dir <- file.path(config$INPUT_DIR, "bracken")
      # Note: param_label detection would need to be implemented for this fallback
      unaligned_bracken_files <- list.files(bracken_dir, pattern = "_unaligned_.*\\.bracken$", full.names = TRUE)
      nonhuman_bracken_files <- list.files(bracken_dir, pattern = "_nonhuman_.*\\.bracken$", full.names = TRUE)
      
      if (length(unaligned_bracken_files) > 0) {
        args <- c(args, "--unaligned-bracken", bracken_dir)
      }
      if (length(nonhuman_bracken_files) > 0) {
        args <- c(args, "--nonhuman-bracken", bracken_dir)
      }
    }
  }
  
  # Add metadata arguments (use manual paths if provided, otherwise auto-detect)
  if (!is.null(config$RUNTIME_DIR)) {
    args <- c(args, "--runtime-dir", config$RUNTIME_DIR)
  } else if (dir.exists(config$INPUT_DIR) && file.exists(file.path(config$INPUT_DIR, "runtime"))) {
    args <- c(args, "--runtime-dir", config$INPUT_DIR)
  }
  
  if (!is.null(config$RRSTATS_DIR)) {
    args <- c(args, "--rrstats-dir", config$RRSTATS_DIR)
  } else if (dir.exists(config$INPUT_DIR) && dir.exists(file.path(config$INPUT_DIR, "QC"))) {
    args <- c(args, "--rrstats-dir", file.path(config$INPUT_DIR, "QC"))
  }
  
  if (!is.null(config$METADATA_DIR)) {
    args <- c(args, "--metadata-dir", config$METADATA_DIR)
  } else if (dir.exists(config$INPUT_DIR)) {
    metadata_files <- list.files(config$INPUT_DIR, pattern = "sample_metadata.*\\.csv$", full.names = TRUE, recursive = FALSE)
    if (length(metadata_files) > 0) {
      args <- c(args, "--metadata-dir", config$INPUT_DIR)
    }
  }
  
  # Add optional flags
  if (config$ADD_PARAMS) {
    args <- c(args, "--params")
  }
  if (config$INCLUDE_SUBSPECIES) {
    args <- c(args, "--subspecies")
  }
  if (!is.null(config$EXCLUDE_TAXID)) {
    args <- c(args, "--exclude-taxid", as.character(config$EXCLUDE_TAXID))
  }
  if (!is.null(config$MINIMIZER_RATIO)) {
    args <- c(args, "--minimizer-ratio", as.character(config$MINIMIZER_RATIO))
  }
  if (!is.null(config$MINIMIZER_THRESHOLD)) {
    args <- c(args, "--minimizer-threshold", as.character(config$MINIMIZER_THRESHOLD))
  }
  
  return(args)
}

# ============================================================================
# Analysis Reports
# ============================================================================

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
  
  end_time <- Sys.time()
  cat("Finished correlation analysis in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
  
  # Extract correlation results and plots
  correlation_results <- tibble(
    sample = character(),
    pearson_r = numeric(),
    n_species = integer()
  )
  
  plots_list <- list()
  last_plot <- NULL
  
  # TODO: Add Parallel Processing
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

# Function to create runtime analysis
runtime_analysis_report <- function(runtime_data, report_dir) {
  pdf(file.path(report_dir, "runtime_analysis.pdf"), width = 8, height = 6)
  plot.new()
  
  title_text <- "Part 1: Runtime Information"
  text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
  
  # Create reference text instead of displaying the table
  ref_text <- paste0(
    "Runtime information has been compiled and saved as CSV files.\n\n",
    "The runtime data includes:\n",
    "- Sample names (* indicates DB loading samples)\n",
    "- Cutadapt/STAR/Samtools runtime (seconds)\n",
    "- Kraken2 unaligned runtime (seconds)\n",
    "- Kraken2 non-human runtime (seconds)\n",
    "- Merged input read counts from RRstats\n",
    "- Merged percent mapped reads from RRstats\n",
    "- Database loading flags (first samples)\n\n",
    "Data source: Parsed from runtime*.txt files in runtimes/ folder\n",
    "Merging: Combined with read statistics by sample name\n",
    "Output location: Check outputs/ directory for CSV files containing this data\n\n",
    "Key findings:\n",
    "- Total samples with runtime data: ", nrow(runtime_data), "\n",
    "- Average Cutadapt/STAR/Samtools runtime: ", round(mean(runtime_data$cutadapt_star_samtools_rt, na.rm = TRUE)), " seconds\n",
    "- Average Kraken2 unaligned runtime: ", round(mean(runtime_data$kraken2_unaligned_rt, na.rm = TRUE)), " seconds\n",
    "- Average Kraken2 non-human runtime: ", round(mean(runtime_data$kraken2_nonhuman_rt, na.rm = TRUE)), " seconds\n",
    "- Samples involved in DB loading: ", sum(runtime_data$is_first_k2_unaligned | runtime_data$is_first_k2_nonhuman, na.rm = TRUE)
  )
  
  text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
  
  # Prepare data for plotting
  runtime_plot_data <- runtime_data
  
  # Add asterisks to sample names for first samples
  # TODO: Fix this to ensure it occurs for both unaligned and nonhuman samples
  runtime_plot_data$sample_display <- ifelse(
    runtime_plot_data$is_first_k2_unaligned | runtime_plot_data$is_first_k2_nonhuman,
    paste0(runtime_plot_data$sample, "*"),
    runtime_plot_data$sample
  )
  
  # Add read statistics if available
  if ("num_input_reads" %in% colnames(runtime_data)) {
      merge_cols <- c("sample", "num_input_reads")
      if ("percent_mapped_reads" %in% colnames(runtime_data)) {
          merge_cols <- c(merge_cols, "percent_mapped_reads")
      }
      runtime_plot_data <- merge(runtime_plot_data, 
                                 runtime_data[, merge_cols, drop = FALSE], 
                                 by = "sample", all.x = TRUE)
  }
  
  # Calculate total runtime
  runtime_plot_data$total_runtime <- rowSums(runtime_plot_data[, c("cutadapt_star_samtools_rt", "kraken2_unaligned_rt", "kraken2_nonhuman_rt")], na.rm = TRUE)
  
  # Calculate averages for non-database loading samples
  non_db_samples_unaligned <- !runtime_plot_data$is_first_k2_unaligned
  non_db_samples_nonhuman <- !runtime_plot_data$is_first_k2_nonhuman
  
  avg_k2_unaligned_non_db <- mean(runtime_plot_data$kraken2_unaligned_rt[non_db_samples_unaligned], na.rm = TRUE)
  avg_k2_nonhuman_non_db <- mean(runtime_plot_data$kraken2_nonhuman_rt[non_db_samples_nonhuman], na.rm = TRUE)
  
  # Plot 1: Total Runtime
  p1 <- ggplot(runtime_plot_data, aes(x = reorder(sample_display, -total_runtime), y = total_runtime)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = mean(runtime_plot_data$total_runtime, na.rm = TRUE),
               color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = nrow(runtime_plot_data), y = mean(runtime_plot_data$total_runtime, na.rm = TRUE), 
             label = paste("Avg:", round(mean(runtime_plot_data$total_runtime, na.rm = TRUE)), "s"), 
             hjust = 1, vjust = -0.5, color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    labs(title = "Total Runtime per Sample", x = "Sample", y = "Total Runtime (seconds)")
  
  print(p1)
  
  # Plot 2: Cutadapt/STAR/Samtools Runtime
  p2 <- ggplot(runtime_plot_data, aes(x = reorder(sample_display, -cutadapt_star_samtools_rt), y = cutadapt_star_samtools_rt)) +
    geom_bar(stat = "identity", fill = "forestgreen", alpha = 0.7) +
    geom_hline(yintercept = mean(runtime_plot_data$cutadapt_star_samtools_rt, na.rm = TRUE), 
               color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = nrow(runtime_plot_data), y = mean(runtime_plot_data$cutadapt_star_samtools_rt, na.rm = TRUE), 
             label = paste("Avg:", round(mean(runtime_plot_data$cutadapt_star_samtools_rt, na.rm = TRUE)), "s"), 
             hjust = 1, vjust = -0.5, color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    labs(title = "Cutadapt/STAR/Samtools Runtime", x = "Sample", y = "Runtime (seconds)")
  
  print(p2)
  
  # Plot 3: Kraken2 Unaligned Runtime (with DB loading stacked)
  # Prepare stacked data for K2 unaligned
  # Prepare K2 unaligned stacked components correctly
  k2_unaligned_data <- runtime_plot_data
  k2_nonhuman_data <- runtime_plot_data
  
  k2_unaligned_data$db_loading_time <- ifelse(
      k2_unaligned_data$is_first_k2_unaligned, pmax(0, k2_nonhuman_data$kraken2_nonhuman_rt - avg_k2_nonhuman_non_db),
      0
  )
  k2_unaligned_data$classification_time <- k2_unaligned_data$kraken2_unaligned_rt - k2_unaligned_data$db_loading_time
  
  # Prepare K2 non-human stacked components
  k2_nonhuman_data$db_loading_time <- ifelse(
      k2_nonhuman_data$is_first_k2_nonhuman,
      pmax(0, k2_nonhuman_data$kraken2_nonhuman_rt - avg_k2_nonhuman_non_db),
      0
  )
  k2_nonhuman_data$classification_time <- k2_nonhuman_data$kraken2_nonhuman_rt - k2_nonhuman_data$db_loading_time
  
  # Reshape for stacked bar plot - unaligned
  k2_unaligned_long <- bind_rows(
      k2_unaligned_data %>% select(sample_display, classification_time, db_loading_time, kraken2_unaligned_rt) %>%
          rename(sample = sample_display, Classification = classification_time, `DB Loading Estimate` = db_loading_time, total_runtime = kraken2_unaligned_rt) %>%
          pivot_longer(cols = c("Classification", "DB Loading Estimate"), names_to = "runtime_type", values_to = "runtime")
  )
  k2_unaligned_long$runtime_type <- factor(k2_unaligned_long$runtime_type, levels = c("Classification", "DB Loading Estimate"))
  
  p3 <- ggplot(k2_unaligned_long, aes(x = reorder(sample, -total_runtime), y = runtime, fill = runtime_type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Classification" = "orange", "DB Loading Estimate" = "darkred")) +
      geom_hline(yintercept = mean(runtime_plot_data$kraken2_unaligned_rt, na.rm = TRUE),
                           color = "blue", linetype = "dashed", linewidth = 1) +
      annotate("text", x = nrow(runtime_plot_data), y = mean(runtime_plot_data$kraken2_unaligned_rt, na.rm = TRUE),
                       label = paste("Avg:", round(mean(runtime_plot_data$kraken2_unaligned_rt, na.rm = TRUE)), "s"),
                       hjust = 1, vjust = -0.5, color = "blue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                  legend.position = "bottom") +
      labs(title = "Kraken2 Unaligned Runtime", x = "Sample", y = "Runtime (seconds)", fill = "Runtime Type")
  
  print(p3)
  
  # Reshape for stacked bar plot - nonhuman
  k2_nonhuman_long <- bind_rows(
      k2_nonhuman_data %>% select(sample_display, classification_time, db_loading_time, kraken2_nonhuman_rt) %>%
          rename(sample = sample_display, Classification = classification_time, `DB Loading Estimate` = db_loading_time, total_runtime = kraken2_nonhuman_rt) %>%
          pivot_longer(cols = c("Classification", "DB Loading Estimate"), names_to = "runtime_type", values_to = "runtime")
  )
  k2_nonhuman_long$runtime_type <- factor(k2_nonhuman_long$runtime_type, levels = c("Classification", "DB Loading Estimate"))
  
  p4 <- ggplot(k2_nonhuman_long, aes(x = reorder(sample, -total_runtime), y = runtime, fill = runtime_type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Classification" = "purple", "DB Loading Estimate" = "darkred")) +
      geom_hline(yintercept = mean(runtime_plot_data$kraken2_nonhuman_rt, na.rm = TRUE),
                           color = "blue", linetype = "dashed", linewidth = 1) +
      annotate("text", x = nrow(runtime_plot_data), y = mean(runtime_plot_data$kraken2_nonhuman_rt, na.rm = TRUE),
                       label = paste("Avg:", round(mean(runtime_plot_data$kraken2_nonhuman_rt, na.rm = TRUE)), "s"),
                       hjust = 1, vjust = -0.5, color = "blue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                  legend.position = "bottom") +
      labs(title = "Kraken2 Non-human Runtime", x = "Sample", y = "Runtime (seconds)", fill = "Runtime Type")
  
  print(p4)
  dev.off()
}

runread_stats_report <- function(runread_data, report_dir) {
  pdf(file.path(report_dir, "run_read_stats_report.pdf"), width = 8, height = 6)
  plot.new()
  
  title_text <- "Part 2: Batch RNAseq Statistics"
  text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
  
  # Create reference text instead of displaying the table
  ref_text <- paste0(
    "RNAseq statistics have been compiled and saved as CSV files.\n\n",
    "The RNAseq statistics data includes:\n",
    "- Sample names\n",
    "- Number of input reads\n",
    "- Average input read length\n",
    "- Uniquely mapped reads\n",
    "- Average mapped length\n", 
    "- Reads assigned to genes\n",
    "- Percent mapped reads\n",
    "- Uniquely mapped percent\n",
    "- Flag for samples with <10M unique reads (highlighted)\n\n",
    "Data source: Parsed from run_read_stats*.txt files in rrstats/ folder\n",
    "Output location: Check outputs/ directory for CSV files containing this data\n\n",
    "Key findings:\n",
    "- Total samples with read statistics: ", nrow(rrstats_data), "\n",
    "- Samples with <10M unique reads: ", sum(rrstats_data$less_than_10M_unique, na.rm = TRUE), "\n",
    "- Average input reads: ", format(round(mean(rrstats_data$num_input_reads, na.rm = TRUE)), big.mark = ","), "\n",
    "- Average mapping percentage: ", round(mean(rrstats_data$percent_mapped_reads, na.rm = TRUE), 1), "%\n",
    "- Average unique mapping percentage: ", round(mean(rrstats_data$uniquely_mapped_percent, na.rm = TRUE), 1), "%"
  )
  
  text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
  dev.off()
}

kraken_stats_report <- function(unaligned_results, nonhuman_results, combined_report_data, config, report_dir, batch) {
  if (batch) {
    pdf(file.path(report_dir, "kraken_stats_report.pdf"), width = 8, height = 6)
    plot.new()
    
    title_text <- "Part 3: Kraken2 Statistics"
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
    
    # Create reference text instead of displaying the table
    ref_text <- paste0(
      "Kraken2 classification results have been compiled and saved as CSV files.\n\n",
      "The Kraken2 results data includes:\n",
      "- Sample names\n",
      "- Species classifications\n",
      "- Unaligned clade reads\n",
      "- Non-human clade reads\n",
      "- Dataset information (unaligned/nonhuman)\n",
      "- Merged input read counts from RRstats\n",
      "- Additional optional columns (confidence, minimum hit groups, etc.)\n\n",
      "Data source: Combined from unaligned_kreports/ and nonhuman_kreports/\n",
      "Processing: output_processing.R with species filtering and merging\n",
      "Merging: Combined with runtime and read statistics by sample name\n",
      "Output location: Check outputs/ directory for:\n",
      "  - sample_report_data_with_metadata*op.csv (combined data)\n",
      "  - sample_report_data_[dataset]_*op.csv (dataset-specific files)\n\n",
      "Key findings:\n",
      "- Unique species across all unaligned samples: ", nrow(unaligned_results$species_list), "\n",
      "- Average number of species found per sample before secondary filtering: ", mean(combined_report_data$species_count), "\n",
      "- Average number of species found per sample after secondary filtering: ", mean(combined_report_data$filtered_species_count), "\n"
    )
    # TODO: add config checks to add info such as subspecies info etc. 
    
    top_clades <- unaligned_results$merged %>% 
      arrange(desc(cladeReads_mean)) %>%
      slice_head(n = config$TOP_N_FREQ)
    
    ref_text <- paste(ref_text, paste(top_clades$name, top_clades$cladeReads_mean, sep = ": ", collapse = "\n"), sep = "\n")
    
    text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    
    dev.off()
  } else {
    
  }
}

species_annotation_report <- function(annotated_species, report_dir, batch) {
  if (batch) {
    pdf(file.path(report_dir, "species_annotation_report.pdf"), width = 8, height = 6)
    plot.new()
    
    title_text <- "Part 4: Species Annotation"
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
    
    # Create reference text instead of displaying the table
    ref_text <- paste0(
      "Species annotation results have been compiled and saved as CSV files.\n\n",
      "The species annotation data includes:\n",
      "- Sample names\n",
      "- Annotated species names\n",
      "- TaxIDs and clade information\n",
      "- Additional metadata columns (if available)\n\n",
      "Data source: Parsed from annotated_species object\n",
      "Output location: Check outputs/ directory for:\n",
      "  - annotated_species.csv (annotated species list)\n\n",
      "Key findings:\n",
      "- Total annotated species: ", nrow(annotated_species), "\n"
    )
    
    text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    
    dev.off()
    
  } else {
    
  }
}

pathogen_detection_report <- function(annotated_species, report_dir, batch) {
  if (batch) {
    
  } else {
    
  }
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Create configuration from command line arguments
config <- create_config(args)

# Setup parallel processing
cl <- setup_parallel(config$CORES)

# Build arguments for output_processing.R
output_args <- build_output_processing_args(config)

cat("Running output_processing.R with arguments:\n")

# Run output processing
run_as_script(here("scripts", "output_processing.R"), output_args)

# Run correlation analysis if requested
if (config$PERFORM_CORRELATION && !is.null(unaligned_results$merged) && !is.null(nonhuman_results$merged)) {
  cat("Running correlation analysis...\n")
  
  # Ensure reports directory exists
  if (!dir.exists(config$reports_dir)) {
    dir.create(config$reports_dir, recursive = TRUE)
  }
  
  correlation_results <- create_correlation_analysis(
    data1 = unaligned_results$merged,
    data2 = nonhuman_results$merged,
    data1_name = "Unaligned",
    data2_name = "Non-Human",
    data1_col_suffix = "cladeReads",
    data2_col_suffix = "cladeReads",
    output_dir = config$reports_dir
  )
  
  # Analysis of host contamination
  cat("Analyzing host contamination...\n")
  host_contamination <- setdiff(unaligned_results$merged$name, nonhuman_results$merged$name)
  cat("Species found only in unaligned data (potential host contamination):", length(host_contamination), "\n")
  
  if (length(host_contamination) > 0) {
    host_contaminated_results <- unaligned_results$merged %>%
      filter(name %in% host_contamination) %>%
      select(name, contains("cladeReads")) %>%
      arrange(desc(rowSums(select(., contains("cladeReads")), na.rm = TRUE)))
    
    cat("Top host-contaminated species:\n")
    print(head(host_contaminated_results, 10))
    
    # Save host contamination results
    write.csv(host_contaminated_results, 
              file.path(config$reports_dir, "host_contamination_species.csv"), 
              row.names = FALSE)
  }
}

#### Batch Report Exporting ####

# Run runtime statistics report
if (config$RT_STATS & exists("runtime_data", envir = .GlobalEnv)) {
  runtime_data <- get("runtime_data", envir = .GlobalEnv)
  cat("Runtime results loaded successfully\n")
  runtime_analysis_report(runtime_data, config$reports_dir)
} else {
  read.csv(get_latest_timestamped_file(input_dir = config$output_dir, pattern = "runtime_info"), 
           stringsAsFactors = FALSE) -> runtime_data
  runtime_analysis_report(runtime_data, config$reports_dir)
}

# Run read statistics report
if (config$RR_STATS & exists("rrstats_data", envir = .GlobalEnv)) {
  rrstats_data <- get("rrstats_data", envir = .GlobalEnv)
  cat("Read statistics data loaded successfully\n")
  runread_stats_report(rrstats_data, config$reports_dir)
} else {
  read.csv(get_latest_timestamped_file(input_dir = config$output_dir, pattern = "read_statistics"), 
           stringsAsFactors = FALSE) -> rrstats_data
}

if (exists("unaligned_results", envir = .GlobalEnv)) {
  unaligned_results <- get("unaligned_results", envir = .GlobalEnv)
  cat("Unaligned results loaded successfully\n")
  
  if (!is.null(config$NONHUMAN_KREPORTS_DIR) & exists("nonhuman_results", envir = .GlobalEnv)) {
    nonhuman_results <- get("nonhuman_results", envir = .GlobalEnv)
    cat("Non-human results loaded successfully\n")
    if (config$KR_STATS & exists(combined_report_data, envir = .GlobalEnv)) {
      combined_report_data <- get("combined_report_data", envir = .GlobalEnv)
      cat("Combined report data loaded successfully\n")
      kraken_stats_report(unaligned_results, nonhuman_results, combined_report_data, config, config$reports_dir)
    } else {
      combined_report_data <- read.csv(get_latest_timestamped_file(input_dir = config$output_dir, pattern = "sample_report_data"), 
                                       stringsAsFactors = FALSE)
    }
  } else {
    stop("ERROR: nonhuman_results not found after running output_processing.R")
  }
  
  # Run species annotation if species_list is available
  if (config$SA_STATS & !is.null(unaligned_results$species_list)) {
    cat("Running species annotation...\n")
    species_df <- unaligned_results$species_list
    assign("species_df", species_df, envir = .GlobalEnv)
    
    # Run annotation with custom databases directory
    run_as_script(here("scripts", "annotate_species.R"), 
                  "--df", "species_df", 
                  "--output", "annotated_species",
                  "--databases-dir", config$DATABASES_DIR)
    
    if (exists("annotated_species", envir = .GlobalEnv)) {
      annotated_species <- get("annotated_species", envir = .GlobalEnv)
      cat("Species annotation completed successfully\n")
      species_annotation_report(annotated_species, config$reports_dir)
    } else {
      stop("ERROR: annotated_species not found after running annotate_species.R")
    }
  }
} else {
stop("ERROR: unaligned_results not found after running output_processing.R")
}

#### Per Sample Report Generation ####

if (!config$NO_PER_SAMPLE) {
  # split merged_long by sample
  trimmed_long <- nonhuman_results$merged_long %>% filter(!is.na(cladeReads)) # remove rows with na values in cladeReads of unaligned_results$merged_lon
  trimmed_by_sample <- split(trimmed_long, trimmed_long$sample) # split trimmed_long into a named list of data.frames by sample
  trimmed_by_sample <- lapply(trimmed_by_sample, as.data.frame) # ensure each element is a plain data.frame (optional)
  
  
  # Create per-sample reports
  cat("Generating per-sample reports...\n")
  sample_dir <- file.path(config$reports_dir, "per_sample_reports")
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }
  if (config$SPLIT_REPORT) {
    for (sample_name in names(trimmed_by_sample)) {
      sample_data <- trimmed_by_sample[[sample_name]]
      pdf(file.path(sample_dir, paste0(sample_name, "_report.pdf")), width = 12, height = 8)
      
      
      dev.off()
    }
  } else {
    # Create a single report for all samples
    pdf(file.path(sample_dir, "all_samples_report.pdf"), width = 12, height = 8)
    for (sample_name in names(trimmed_by_sample)) {
      sample_data <- trimmed_by_sample[[sample_name]]
      
    }
    dev.off()
  }
}

# Clean up parallel cluster
if (exists("cl")) {
  stopCluster(cl)
  cat("Parallel cluster stopped\n")
}

cat("Analysis completed successfully!\n")
cat("Output directory:", config$output_dir, "\n")
cat("Reports directory:", config$reports_dir, "\n")


