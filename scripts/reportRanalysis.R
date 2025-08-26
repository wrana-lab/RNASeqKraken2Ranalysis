#!/usr/bin/env Rscript

#### Header ####
# Kraken/Bracken Analysis Orchestration Script
# Author: Denis Rivard, modified by GitHub Copilot
# Last updated: August 26th 2025
# Description: Orchestrates downstream analysis of kraken2/bracken results including
#              processing kreports, species annotation with risk groups and HOMD data,
#              and correlation analysis between datasets
# Usage: Rscript reportRanalysis.R --input-dir <path> --proj-name <name> [options]

library(here)
library(parallel)
library(doParallel)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

# Utility functions from utils.R
run_as_script <- function() {
  return(!interactive() && length(sys.calls()) == 0)
}

get_latest_timestamped_file <- function(pattern, path = ".") {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No files found matching pattern:", pattern, "in path:", path, "\n")
    return(NULL)
  }
  
  # Extract timestamps from the files
  # Assuming format: prefix_YYMMDD_suffix
  timestamp_pattern <- "([0-9]{6})"
  timestamps <- str_extract(files, timestamp_pattern)
  
  if (all(is.na(timestamps))) {
    cat("Warning: No timestamps found in file names. Returning first file.\n")
    return(files[1])
  }
  
  # Handle files without timestamps
  valid_files <- files[!is.na(timestamps)]
  valid_timestamps <- timestamps[!is.na(timestamps)]
  
  if (length(valid_files) == 0) {
    cat("Warning: No files with valid timestamps found. Returning first file.\n")
    return(files[1])
  }
  
  # Find the most recent file
  latest_index <- which.max(valid_timestamps)
  latest_file <- valid_files[latest_index]
  
  cat("Latest file found:", basename(latest_file), "with timestamp:", valid_timestamps[latest_index], "\n")
  
  return(latest_file)
}

# Function to run scripts with command-line arguments
run_script_with_args <- function(script_path, args_vector) {
  # Handle vector of arguments
  assign("commandArgs", function(trailingOnly = TRUE) args_vector, envir = .GlobalEnv)
  # Run the script
  source(script_path)
}

# Try to load openxlsx for Excel export functionality
excel_available <- FALSE
tryCatch({
  library(openxlsx)
  excel_available <- TRUE
}, error = function(e) {
  cat("Warning: openxlsx not available, Excel export will be skipped\n")
})

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
  cat("  --all-reports             Include all report Analyses\n")
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
    RANALYSIS_DIR = here(),  # Default to current here() location
    OUTPUT_BASE_DIR = "outputs",
    DATABASES_DIR = here("databases"),
    TOP_N_FREQ = 25,
    TOP_N_PLOT = 50,
    MIN_READS = 10,
    EXCLUDE_TAXID = NULL,
    PERFORM_CORRELATION = TRUE,
    ADD_PARAMS = TRUE,
    INCLUDE_SUBSPECIES = FALSE,
    MINIMIZER_RATIO = NULL,
    MINIMIZER_THRESHOLD = 500,
    CORES = detectCores(),
    # Report options
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
      } else if (args[i] == "--ranalysis-dir" && i < length(args)) {
        config$RANALYSIS_DIR <- args[i + 1]
        config$DATABASES_DIR <- file.path(args[i + 1], "databases")
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
      } else if (args[i] == "--all-reports") {
        config$RT_STATS <- TRUE
        config$RR_STATS <- TRUE
        config$KR_STATS <- TRUE
        config$SA_STATS <- TRUE
        config$PD_STATS <- TRUE
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
      } else if (args[i] == "--pd-report") {
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
  if (config$SA_STATS && !dir.exists(config$DATABASES_DIR)) {
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
  # Reuse existing cluster if present and valid
  if (exists("cl", envir = .GlobalEnv)) {
    cl_existing <- tryCatch(get("cl", envir = .GlobalEnv), error = function(e) NULL)
    if (!is.null(cl_existing) && inherits(cl_existing, "cluster")) {
      workers <- tryCatch(getDoParWorkers(), error = function(e) NA)
      if (!is.na(workers) && workers > 0) {
        cat("Existing cluster detected with", workers, "workers. Reusing it.\n")
        registerDoParallel(cl_existing)
        return(cl_existing)
      } else {
        cat("Existing cluster object found but not registered. Registering it now.\n")
        registerDoParallel(cl_existing)
        return(cl_existing)
      }
    }
  }
  # No valid existing cluster found -> create a new one
  start_time <- Sys.time()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  assign("cl", cl, envir = .GlobalEnv)
  end_time <- Sys.time()
  workers <- tryCatch(getDoParWorkers(), error = function(e) NA)
  cat("Parallel cluster initialized with", workers, "workers, in", 
      round(difftime(end_time, start_time, units = "secs")), "seconds\n")
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
    } else {
      # Check if INPUT_DIR itself contains kreport files (direct directory case)
      kreport_files <- list.files(config$INPUT_DIR, pattern = "\\.(kreport|k2report)$", recursive = TRUE)
      if (length(kreport_files) > 0) {
        # Check subdirectories for unaligned and nonhuman patterns
        subdirs <- list.dirs(config$INPUT_DIR, full.names = TRUE, recursive = FALSE)
        unaligned_subdirs <- grep("unaligned", subdirs, value = TRUE)
        nonhuman_subdirs <- grep("nonhuman", subdirs, value = TRUE)
        
        if (is.null(unaligned_kreports_dir) && length(unaligned_subdirs) > 0) {
          unaligned_kreports_dir <- unaligned_subdirs[1]
        }
        if (is.null(nonhuman_kreports_dir) && length(nonhuman_subdirs) > 0) {
          nonhuman_kreports_dir <- nonhuman_subdirs[1]
        }
      }
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
  
  # Count samples in unaligned directory
  unaligned_samples <- 0
  if (dir.exists(unaligned_kreports_dir)) {
    unaligned_files <- list.files(unaligned_kreports_dir, pattern = "\\.(kreport|k2report)$")
    unaligned_samples <- length(unaligned_files)
  }
  
  # Check nonhuman directory and validate sample count match
  nonhuman_valid <- FALSE
  if (!is.null(nonhuman_kreports_dir) && dir.exists(nonhuman_kreports_dir)) {
    nonhuman_files <- list.files(nonhuman_kreports_dir, pattern = "\\.(kreport|k2report)$")
    nonhuman_samples <- length(nonhuman_files)
    
    if (nonhuman_samples == unaligned_samples && nonhuman_samples > 0) {
      nonhuman_valid <- TRUE
      cat("Found matching nonhuman kreports directory with", nonhuman_samples, "samples\n")
    } else {
      cat("Warning: Nonhuman kreports directory found but sample count mismatch (unaligned:", unaligned_samples, ", nonhuman:", nonhuman_samples, ")\n")
      nonhuman_kreports_dir <- NULL
    }
  } else {
    cat("Warning: No valid nonhuman kreports directory found. Analysis will use only unaligned data.\n")
    nonhuman_kreports_dir <- NULL
  }
  
  args <- c(
    "--unaligned-kreports", unaligned_kreports_dir,
    "--output-dir", config$output_dir,
    "--top-n", as.character(config$TOP_N_FREQ),
    "--min-reads", as.character(config$MIN_READS),
    "--output"
  )
  
  if (config$CORES > 1) {
    args <- c(args, "--parallel")
  } else {
    cat("Note: Running in single-core mode. For faster processing, consider using multiple cores with --cores <N> > 1.\n")
  }
  
  # Only add nonhuman kreports if valid
  if (nonhuman_valid && !is.null(nonhuman_kreports_dir)) {
    args <- c(args, "--nonhuman-kreports", nonhuman_kreports_dir)
  }
  
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
  
  # Only add nonhuman bracken if nonhuman kreports are valid
  if (nonhuman_valid) {
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
create_correlation_analysis <- function(data1, data2, data1_name, data2_name, data1_col_suffix, data2_col_suffix, output_dir, top_n_label = config$TOP_N_PLOT) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  top_n_label <- max(1, floor(config$TOP_N_PLOT / 2))
  
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
  
  # Process samples with parallel processing when available
  if (exists("cl", envir = .GlobalEnv) && !is.null(get("cl", envir = .GlobalEnv))) {
    cat("Processing correlation analysis in parallel using", length(get("cl", envir = .GlobalEnv)), "cores...\n")
    
    # Set up parallel backend
    cluster <- get("cl", envir = .GlobalEnv)
    
    # Export required variables to cluster nodes
    clusterExport(cluster, c("data1", "data2", "data1_col_suffix", "data2_col_suffix", 
                             "data1_name", "data2_name"), envir = environment())
    
    # Parallel processing using foreach
    sample_results <- foreach(sample_id = common_samples, 
                              .packages = c("dplyr", "ggplot2", "ggrepel"),
                              .combine = 'c') %dopar% {
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
        
        # Create labels for only the top_n_label most different points
        sample_comparison <- sample_comparison %>%
          mutate(
            log_reads1 = log10(reads1 + 1),
            log_reads2 = log10(reads2 + 1),
            diff = abs(log_reads1 - log_reads2)
          ) %>%
          arrange(desc(diff)) %>%
          mutate(
            rank = row_number(),
            label = ifelse(rank <= top_n_label, name, NA_character_)
          ) %>%
          select(-diff, -rank)
        
        # Create the correlation plot
        gg_correlation <- ggplot(sample_comparison, aes(x = log_reads1, y = log_reads2)) +
          geom_point(alpha = 0.6, size = 1.5) + 
          geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
          geom_text_repel(aes(label = label), na.rm = TRUE, size = 2.5, max.overlaps = top_n_label) +
          labs(
            x = paste0("Log10(", data1_name, " ", data1_col_suffix, " + 1)"),
            y = paste0("Log10(", data2_name, " ", data2_col_suffix, " + 1)"),
            title = paste0("Correlation: ", data1_name, " vs ", data2_name, " (", sample_id, ")"),
            subtitle = paste0("Pearson r = ", round(cor_coef, 5), " (n = ", nrow(sample_comparison), " species)")
          ) + coord_fixed()
        
        list(list(
          sample = sample_id,
          pearson_r = cor_coef,
          n_species = nrow(sample_comparison),
          plot = gg_correlation
        ))
      } else {
        list(NULL)
      }
    }
    
    # Remove NULL entries
    sample_results <- sample_results[!sapply(sample_results, is.null)]
    
  } else {
    # Fallback to sequential processing
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
        
        # Create labels for only the top_n_label most different points
        sample_comparison <- sample_comparison %>%
          mutate(
            log_reads1 = log10(reads1 + 1),
            log_reads2 = log10(reads2 + 1),
            diff = abs(log_reads1 - log_reads2)
          ) %>%
          arrange(desc(diff)) %>%
          mutate(
            rank = row_number(),
            label = ifelse(rank <= top_n_label, name, NA_character_)
          ) %>%
          select(-diff, -rank)
        
        # Create the correlation plot
        gg_correlation <- ggplot(sample_comparison, aes(x = log_reads1, y = log_reads2)) +
            geom_point(alpha = 0.6, size = 1.5) + 
            geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
            geom_text_repel(aes(label = label), na.rm = TRUE, size = 2.5, max.overlaps = top_n_label) +
            labs(
              x = paste0("Log10(", data1_name, " ", data1_col_suffix, " + 1)"),
              y = paste0("Log10(", data2_name, " ", data2_col_suffix, " + 1)"),
              title = paste0("Correlation: ", data1_name, " vs ", data2_name, " (", sample_id, ")"),
              subtitle = paste0("Pearson r = ", round(cor_coef, 5), " (n = ", nrow(sample_comparison), " species)")
            ) + coord_fixed()
        
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
  
  # Extract correlation results and plots
  correlation_results <- tibble(
    sample = character(),
    pearson_r = numeric(),
    n_species = integer()
  )
  
  plots_list <- list()
  last_plot <- NULL
  
  # Process results (compatible with both parallel and sequential)
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
  
  # Create reference text based on available data
  has_nonhuman_runtime <- "kraken2_nonhuman_rt" %in% colnames(runtime_data) && 
                          any(!is.na(runtime_data$kraken2_nonhuman_rt))
  
  ref_text <- paste0(
    "Runtime information has been compiled and saved as CSV files.\n\n",
    "The runtime data includes:\n",
    "- Sample names (* indicates DB loading samples)\n",
    "- Cutadapt/STAR/Samtools runtime (seconds)\n",
    "- Kraken2 unaligned runtime (seconds)\n"
  )
  
  if (has_nonhuman_runtime) {
    ref_text <- paste0(ref_text, "- Kraken2 non-human runtime (seconds)\n")
  } else {
    ref_text <- paste0(ref_text, "- Kraken2 non-human runtime (not available - unaligned only)\n")
  }
  
  ref_text <- paste0(ref_text,
    "- Merged input read counts from RRstats\n",
    "- Merged percent mapped reads from RRstats\n",
    "- Database loading flags (first samples)\n\n",
    "Data source: Parsed from runtime*.txt files in runtimes/ folder\n",
    "Merging: Combined with read statistics by sample name\n",
    "Output location: Check outputs/ directory for CSV files containing this data\n\n",
    "Key findings:\n",
    "- Total samples with runtime data: ", nrow(runtime_data), "\n",
    "- Average Cutadapt/STAR/Samtools runtime: ", round(mean(runtime_data$cutadapt_star_samtools_rt, na.rm = TRUE)), " seconds\n",
    "- Average Kraken2 unaligned runtime: ", round(mean(runtime_data$kraken2_unaligned_rt, na.rm = TRUE)), " seconds\n"
  )
  
  if (has_nonhuman_runtime) {
    ref_text <- paste0(ref_text, 
      "- Average Kraken2 non-human runtime: ", round(mean(runtime_data$kraken2_nonhuman_rt, na.rm = TRUE)), " seconds\n",
      "- Samples involved in DB loading: ", 
      sum(runtime_data$is_first_k2_unaligned | 
          (if("is_first_k2_nonhuman" %in% colnames(runtime_data)) runtime_data$is_first_k2_nonhuman else FALSE), 
          na.rm = TRUE), "\n"
    )
  } else {
    ref_text <- paste0(ref_text, 
      "- Average Kraken2 non-human runtime: N/A (unaligned analysis only)\n",
      "- Samples involved in DB loading: ", sum(runtime_data$is_first_k2_unaligned, na.rm = TRUE), "\n"
    )
  }
  
  # Add runtime per input reads statistics if data is available
  if ("num_input_reads" %in% colnames(runtime_data)) {
    unaligned_rt_per_input <- runtime_data$kraken2_unaligned_rt / runtime_data$num_input_reads
    ref_text <- paste0(ref_text, 
      "- Average Kraken2 unaligned runtime per input read: ", 
      round(mean(unaligned_rt_per_input, na.rm = TRUE), 6), " seconds/read\n"
    )
    
    if (has_nonhuman_runtime && "percent_mapped_reads" %in% colnames(runtime_data)) {
      # For non-human: runtime per (input_reads * percent_mapped_reads)
      nonhuman_effective_reads <- runtime_data$num_input_reads * (runtime_data$percent_mapped_reads / 100)
      nonhuman_rt_per_effective <- runtime_data$kraken2_nonhuman_rt / nonhuman_effective_reads
      ref_text <- paste0(ref_text, 
        "- Average Kraken2 non-human runtime per effective read: ", 
        round(mean(nonhuman_rt_per_effective, na.rm = TRUE), 6), " seconds/read\n"
      )
    }
  }
  
  text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
  
  # Prepare data for plotting
  runtime_plot_data <- runtime_data
  
  # Check if we have non-human data
  has_nonhuman_data <- "kraken2_nonhuman_rt" %in% colnames(runtime_plot_data) && 
    any(!is.na(runtime_plot_data$kraken2_nonhuman_rt))
  
  # Add asterisks to sample names for first samples
  # Handle cases where non-human data might not be available
  if (has_nonhuman_data && "is_first_k2_nonhuman" %in% colnames(runtime_plot_data)) {
    runtime_plot_data$sample_display <- ifelse(
      runtime_plot_data$is_first_k2_unaligned | runtime_plot_data$is_first_k2_nonhuman,
      paste0(runtime_plot_data$sample, "*"),
      runtime_plot_data$sample
    )
  } else {
    runtime_plot_data$sample_display <- ifelse(
      runtime_plot_data$is_first_k2_unaligned,
      paste0(runtime_plot_data$sample, "*"),
      runtime_plot_data$sample
    )
  }
  
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
  
  # Calculate total runtime (handle missing non-human data)
  if (has_nonhuman_data) {
    runtime_plot_data$total_runtime <- rowSums(runtime_plot_data[, c("cutadapt_star_samtools_rt", "kraken2_unaligned_rt", "kraken2_nonhuman_rt")], na.rm = TRUE)
  } else {
    runtime_plot_data$total_runtime <- rowSums(runtime_plot_data[, c("cutadapt_star_samtools_rt", "kraken2_unaligned_rt")], na.rm = TRUE)
  }
  
  # Calculate averages for non-database loading samples
  non_db_samples_unaligned <- !runtime_plot_data$is_first_k2_unaligned
  avg_k2_unaligned_non_db <- mean(runtime_plot_data$kraken2_unaligned_rt[non_db_samples_unaligned], na.rm = TRUE)
  
  if (has_nonhuman_data) {
    non_db_samples_nonhuman <- !runtime_plot_data$is_first_k2_nonhuman
    avg_k2_nonhuman_non_db <- mean(runtime_plot_data$kraken2_nonhuman_rt[non_db_samples_nonhuman], na.rm = TRUE)
  }
  
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
  
  k2_unaligned_data$db_loading_time <- ifelse(
      k2_unaligned_data$is_first_k2_unaligned, 
      pmax(0, k2_unaligned_data$kraken2_unaligned_rt - avg_k2_nonhuman_non_db),
      0
  )
  k2_unaligned_data$classification_time <- k2_unaligned_data$kraken2_unaligned_rt - k2_unaligned_data$db_loading_time
  
  # Only prepare K2 non-human stacked components if data is available
  if (has_nonhuman_data) {
    k2_nonhuman_data <- runtime_plot_data
    k2_nonhuman_data$db_loading_time <- ifelse(
        k2_nonhuman_data$is_first_k2_nonhuman,
        pmax(0, k2_nonhuman_data$kraken2_nonhuman_rt - avg_k2_nonhuman_non_db),
        0
    )
    k2_nonhuman_data$classification_time <- k2_nonhuman_data$kraken2_nonhuman_rt - k2_nonhuman_data$db_loading_time
    # Use more accurate K2 non-human runtime data for DB loading estimate
    k2_unaligned_data$db_loading_time <- ifelse(
      k2_unaligned_data$is_first_k2_unaligned, 
      max(k2_nonhuman_data$db_loading_time),
      0
    )
    k2_unaligned_data$classification_time <- k2_unaligned_data$kraken2_unaligned_rt - k2_unaligned_data$db_loading_time
  }
  
  # Reshape for stacked bar plot - unaligned
  k2_unaligned_long <- k2_unaligned_data %>% 
    select(sample_display, classification_time, db_loading_time, kraken2_unaligned_rt) %>%
    rename(sample = sample_display, total_runtime = kraken2_unaligned_rt) %>%
    pivot_longer(cols = c("classification_time", "db_loading_time"), 
                 names_to = "runtime_type", values_to = "runtime") %>%
    mutate(runtime_type = case_when(
      runtime_type == "classification_time" ~ "Classification",
      runtime_type == "db_loading_time" ~ "DB Loading Estimate"
    )) %>%
    mutate(runtime_type = factor(runtime_type, levels = c("Classification", "DB Loading Estimate")))
  
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
  
  # Only create non-human plot if data is available
  if (has_nonhuman_data) {
    # Reshape for stacked bar plot - nonhuman
    k2_nonhuman_long <- k2_nonhuman_data %>% 
      select(sample_display, classification_time, db_loading_time, kraken2_nonhuman_rt) %>%
      rename(sample = sample_display, total_runtime = kraken2_nonhuman_rt) %>%
      pivot_longer(cols = c("classification_time", "db_loading_time"), 
                   names_to = "runtime_type", values_to = "runtime") %>%
      mutate(runtime_type = case_when(
        runtime_type == "classification_time" ~ "Classification",
        runtime_type == "db_loading_time" ~ "DB Loading Estimate"
      )) %>%
      mutate(runtime_type = factor(runtime_type, levels = c("Classification", "DB Loading Estimate")))
    
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
  } else {
    # Create a placeholder plot indicating no non-human data
    plot.new()
    title_text <- "Kraken2 Non-human Runtime"
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
    text(0.5, 0.5, "No non-human Kraken2 data available\n(Only unaligned analysis was performed)", 
         adj = c(0.5, 0.5), cex = 1, col = "gray50")
  }
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
    "- Flag for samples with <10M unique reads\n\n",
    "Data source: Parsed from run_read_stats*.txt files in rrstats/ folder\n",
    "Output location: Check outputs/ directory for CSV files containing this data\n\n",
    "Key findings:\n",
    "- Total samples with read statistics: ", nrow(runread_data), "\n",
    "- Samples with <10M unique reads: ", sum(runread_data$less_than_10M_unique, na.rm = TRUE), "\n",
    "- Average input reads: ", format(round(mean(runread_data$num_input_reads, na.rm = TRUE)), big.mark = ","), "\n",
    "- Average mapping percentage: ", round(mean(runread_data$percent_mapped_reads, na.rm = TRUE), 1), "%\n",
    "- Average unique mapping percentage: ", round(mean(runread_data$uniquely_mapped_percent, na.rm = TRUE), 1), "%"
  )
  
  text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
  dev.off()
  
  # Export to Excel with conditional formatting if openxlsx is available
  tryCatch({
    if (excel_available) {
      excel_file <- file.path(report_dir, "run_read_stats_report.xlsx")
      
      # Create workbook
      wb <- createWorkbook()
      addWorksheet(wb, "Read_Stats")

      # write.csv(runread_data, file = file.path(report_dir, "run_read_stats_report.csv"), row.names = FALSE)
      # runread_data <- read.table(file.path(report_dir, "run_read_stats_report.csv"), header = TRUE, sep = ",")
      
      # Create headers without underscores
      headers <- colnames(runread_data)
      clean_headers <- gsub("_", " ", headers)
      clean_headers <- tools::toTitleCase(clean_headers)
      
      # Prepare data for Excel export - ensure proper data types
      excel_data <- runread_data
      
      # Convert numeric columns to ensure they're properly formatted
      numeric_cols <- c("num_input_reads", "avg_input_read_length", "uniquely_mapped_reads", 
                       "avg_mapped_length", "reads_assigned_to_genes", "percent_mapped_reads", 
                       "uniquely_mapped_percent")
      
      for (col in numeric_cols) {
        if (col %in% colnames(excel_data)) {
          excel_data[[col]] <- as.numeric(excel_data[[col]])
        }
      }
      
      # Ensure logical column is properly formatted
      if ("less_than_10M_unique" %in% colnames(excel_data)) {
        excel_data[["less_than_10M_unique"]] <- as.logical(excel_data[["less_than_10M_unique"]])
      }
      
      # Ensure sample column is character
      if ("sample" %in% colnames(excel_data)) {
        excel_data[["sample"]] <- as.character(excel_data[["sample"]])
      }
      
      # Write headers and data separately
      writeData(wb, "Read_Stats", t(clean_headers), startRow = 1, colNames = FALSE)
      writeData(wb, "Read_Stats", excel_data, startRow = 2, colNames = FALSE)
      
      # Add conditional formatting
      data_rows <- 2:(nrow(excel_data) + 1)
      # Data bars for num_input_reads
      if ("num_input_reads" %in% headers) {
        num_reads_col <- which(headers == "num_input_reads")
        conditionalFormatting(wb, "Read_Stats",
                            cols = num_reads_col,
                            rows = data_rows,
                            type = "dataBar",
                            style = "#0070C0")
        cat("Applied data bars to num_input_reads column\n")
      }
      
      # THESE BREAK THE EXCEL FILE - DISABLED FOR NOW
      # # Green data bars for uniquely_mapped_reads
      # if ("uniquely_mapped_reads" %in% headers) {
      #   unique_col <- which(headers == "uniquely_mapped_reads")
      #   conditionalFormatting(wb, "Read_Stats", 
      #                       cols = unique_col, 
      #                       rows = data_rows,
      #                       type = "dataBar",
      #                       style = "#00B050")
      #   cat("Applied green data bars to uniquely_mapped_reads column\n")
      # }
      # # Orange data bars for reads_assigned_to_genes
      # if ("reads_assigned_to_genes" %in% headers) {
      #   assigned_col <- which(headers == "reads_assigned_to_genes")
      #   conditionalFormatting(wb, "Read_Stats",
      #                       cols = assigned_col,
      #                       rows = data_rows,
      #                       type = "dataBar",
      #                       style = "#FF8C00")
      #   cat("Applied orange data bars to reads_assigned_to_genes column\n")
      # }
    
      # Red-White-Blue color scale for avg_mapped_length
      if ("avg_mapped_length" %in% headers) {
        avg_length_col <- which(headers == "avg_mapped_length")
        conditionalFormatting(wb, "Read_Stats",
                            cols = avg_length_col,
                            rows = data_rows,
                            type = "colorScale",
                            style = c("#FF0000", "#FFFFFF", "#0000FF"))
        # Add borders to color scale cells
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        addStyle(wb, "Read_Stats", border_style, rows = data_rows, cols = avg_length_col)
        cat("Applied red-white-blue color scale to avg_mapped_length column\n")
      }
      
      # Green-Yellow color scale for percent_mapped_reads
      if ("percent_mapped_reads" %in% headers) {
        percent_mapped_col <- which(headers == "percent_mapped_reads")
        conditionalFormatting(wb, "Read_Stats",
                            cols = percent_mapped_col,
                            rows = data_rows,
                            type = "colorScale",
                            style = c("#00FF00", "#FFFF00"))
        # Add borders to color scale cells
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        addStyle(wb, "Read_Stats", border_style, rows = data_rows, cols = percent_mapped_col)
        cat("Applied green-yellow color scale to percent_mapped_reads column\n")
      }

      # Green-Yellow color scale for uniquely_mapped_percent
      if ("uniquely_mapped_percent" %in% headers) {
        unique_percent_col <- which(headers == "uniquely_mapped_percent")
        conditionalFormatting(wb, "Read_Stats",
                            cols = unique_percent_col,
                            rows = data_rows,
                            type = "colorScale",
                            style = c("#00FF00", "#FFFF00"))
        # Add borders to color scale cells
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        addStyle(wb, "Read_Stats", border_style, rows = data_rows, cols = unique_percent_col)
        cat("Applied green-yellow color scale to uniquely_mapped_percent column\n")
      }
      
      # Auto-size columns
      setColWidths(wb, "Read_Stats", cols = 1:length(headers), widths = "auto")
      
      # Add header formatting
      header_style <- createStyle(
        textDecoration = "bold",
        fgFill = "#E6E6FA",
        border = "TopBottomLeftRight",
        borderStyle = "thin"
      )
      addStyle(wb, "Read_Stats", header_style, rows = 1, cols = 1:length(headers))
      
      # Save workbook
      saveWorkbook(wb, excel_file, overwrite = TRUE)
      cat("Read statistics exported to Excel:", excel_file, "\n")
    } else {
      cat("Excel export skipped - openxlsx not available\n")
    }
  }, error = function(e) {
    cat("Warning: Failed to export read statistics to Excel:", e$message, "\n")
    print(e)
  })
}

kraken_stats_report <- function(unaligned_results, nonhuman_results = NULL, combined_report_data, config, report_dir, batch) {
  if (batch) {
    pdf(file.path(report_dir, "kraken_stats_report.pdf"), width = 8, height = 6)
    plot.new()
    
    title_text <- "Part 3: Kraken2 Statistics"
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
    
    # Check if nonhuman data is available
    has_nonhuman <- !is.null(nonhuman_results) && !is.null(nonhuman_results$merged)
    
    # Create reference text instead of displaying the table
    ref_text <- paste0(
      "Kraken2 classification results have been compiled and saved as CSV files.\n\n",
      "The Kraken2 results data includes:\n",
      "- Sample names\n",
      "- Species classifications\n",
      "- Unaligned clade reads\n"
    )
    
    if (has_nonhuman) {
      ref_text <- paste0(ref_text, "- Non-human clade reads\n")
      ref_text <- paste0(ref_text, "- Dataset information (unaligned/nonhuman)\n")
    } else {
      ref_text <- paste0(ref_text, "- Dataset information (unaligned only - no nonhuman data available)\n")
    }
    
    ref_text <- paste0(ref_text,
      "- Merged input read counts from RRstats\n",
      "- Additional optional columns (confidence, minimum hit groups, etc.)\n\n",
      "Data source: Combined from unaligned_kreports/", if (has_nonhuman) " and nonhuman_kreports/" else "", "\n",
      "Processing: output_processing.R with species filtering and merging\n",
      "Merging: Combined with runtime and read statistics by sample name\n",
      "Output location: Check outputs/ directory for:\n",
      "  - sample_report_data_with_metadata*op.csv (combined data)\n",
      "  - sample_report_data_[dataset]_*op.csv (dataset-specific files)\n\n",
      "Key findings:\n",
      "- Average number of species found per sample before secondary filtering: ", mean(combined_report_data$species_count, na.rm = TRUE), "\n",
      "- Average number of species found per sample after secondary filtering: ", mean(combined_report_data$filtered_species_count, na.rm = TRUE), "\n",
      "- Number of unique species across all unaligned samples: ", nrow(unaligned_results$species_list), "\n"
    )
    
    if (config$INCLUDE_SUBSPECIES) {
      ref_text <- paste(ref_text, "- Average subspecies across all samples before secondary filtering: ", mean(combined_report_data$subspecies_count), "\n", sep = "")
      ref_text <- paste(ref_text, "- Average subspecies across all samples after secondary filtering: ", mean(combined_report_data$subspecies_count), "\n", sep = "")
    }
    
    top_clades <- unaligned_results$merged %>% 
      arrange(desc(cladeReads_mean)) %>%
      slice_head(n = config$TOP_N_FREQ)
    ref_text <- paste(ref_text, paste(top_clades$name, top_clades$cladeReads_mean, sep = ": ", collapse = "\n"), sep = "\n")
    
    text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    
    if (has_nonhuman) {
      cat("Creating Kraken2 statistics report using non-human sample results...\n")
      species_list <- nonhuman_results$species_list
    } else {
      species_list <- unaligned_results$species_list
    }
    
    # Assign to global environment for the plot scripts to access
    assign("species_list", species_list, envir = .GlobalEnv)
    
    # Source the function-based plotting (disable main execution)
    assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
    source(file.path(config$RANALYSIS_DIR, "scripts", "annotate_species_plots.R"))
    
    # Kingdom basic plot
    kingdom_basic <- generate_annotation_plot(
      df = species_list,
      plot_type = "kingdom"
    )
    print(kingdom_basic)
    
    # Kingdom detailed plots
    detailed_kingdom_plot <- generate_annotation_plot(
      df = species_list,
      plot_type = "kingdom",
      detailed = "Freq",
      log = TRUE
    )
    print(detailed_kingdom_plot)
    
    # Check if bracken data is available for advanced plotting
    if (has_nonhuman) {
      bracken_cols <- grep("bracken_reads_mean", colnames(nonhuman_results$species_list), value = TRUE)
    } else {
      # For unaligned results, check the same way
      bracken_cols <- grep("bracken_reads_mean", colnames(unaligned_results$species_list), value = TRUE)
    }
    
    if (length(bracken_cols) > 0) {
      # First: Show log mean bracken reads kingdom plot
      bracken_kingdom_plot <- generate_annotation_plot(
        df = species_list,
        plot_type = "kingdom",
        detailed = "bracken_reads_mean",
        log = TRUE
      )
      print(bracken_kingdom_plot)
      
      # Second: Create side-by-side comparison plot (clade reads vs bracken reads)
      cat("Creating side-by-side clade vs bracken reads comparison plot...\n")
      
      # Calculate dynamic N based on config$TOP_N_PLOT
      dynamic_n <- max(1, floor(config$TOP_N_PLOT / 2))
      
      # Prepare data for side-by-side comparison (ranked left-to-right by species)
      top_species_comparison <- unaligned_results$species_list %>%
        filter(!is.na(cladeReads_mean) & !is.na(bracken_reads_mean)) %>%
        arrange(desc(cladeReads_mean)) %>%
        slice_head(n = dynamic_n) %>%
        mutate(species_rank = row_number()) %>%
        select(name, cladeReads_mean, bracken_reads_mean, species_rank)
      
      if (nrow(top_species_comparison) > 0) {
        # Transform data for side-by-side plotting
        comparison_long <- top_species_comparison %>%
          pivot_longer(cols = c("cladeReads_mean", "bracken_reads_mean"), 
                       names_to = "read_type", values_to = "reads") %>%
          mutate(read_type = case_when(
            read_type == "cladeReads_mean" ~ "Kraken2 Clade",
            read_type == "bracken_reads_mean" ~ "Bracken"
          )) %>%
          mutate(read_type = factor(read_type, levels = c("Kraken2 Clade", "Bracken")))
        
        # Create side-by-side bar plot with improved ranking and log scale
        comparison_plot <- ggplot(comparison_long, 
                                 aes(x = reorder(name, species_rank), y = log10(reads + 1), fill = read_type)) +
          geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
          scale_fill_manual(values = c("Kraken2 Clade" = "steelblue", "Bracken" = "orange")) +
          scale_y_continuous(trans = "log10", 
                             labels = function(x) format(10^x, scientific = FALSE)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                legend.position = "bottom") +
          labs(
            title = paste0("Kraken2 Clade Reads vs Bracken Reads (Top ", dynamic_n, " Species)"),
            subtitle = paste("Ranked by Kraken2 clade reads (log10 scale)"),
            x = "Species (ranked by clade reads)",
            y = "Log10(Reads + 1)",
            fill = "Read Type"
          )
        
        print(comparison_plot)
      } else {
        plot.new()
        text(0.5, 0.5, "Insufficient data for Bracken comparison plot", 
             adj = c(0.5, 0.5), cex = 1, col = "gray50")
      }
    } else {
      # No bracken data - show the original clade reads mean plot
      mean_kingdom_plot <- generate_annotation_plot(
        df = species_list,
        plot_type = "kingdom",
        detailed = "cladeReads_mean",
        log = TRUE
      )
      print(mean_kingdom_plot)
      
      plot.new()
      text(0.5, 0.5, "No Bracken data available for comparison", 
           adj = c(0.5, 0.5), cex = 1, col = "gray50")
    }
    
    dev.off()
    
    # Add new columns to combined_report_data and export to Excel
    tryCatch({
      # Export to Excel with conditional formatting if openxlsx is available
      if (excel_available) {
        excel_file <- file.path(report_dir, "kraken_stats_report.xlsx")
        
        # Create workbook
        wb <- createWorkbook()
        addWorksheet(wb, "Kraken_Stats")
        
        # Create headers without underscores
        headers <- colnames(combined_report_data)
        clean_headers <- gsub("_", " ", headers)
        clean_headers <- tools::toTitleCase(clean_headers)
        
        # Write headers and data separately
        writeData(wb, "Kraken_Stats", t(clean_headers), startRow = 1, colNames = FALSE)
        writeData(wb, "Kraken_Stats", combined_report_data, startRow = 2, colNames = FALSE)
        
        # Add conditional formatting
        data_rows <- 2:(nrow(combined_report_data) + 1)
        
        # Red-white color scale for unclassified counts
        unclassified_patterns <- c("unclassified_counts", "unclassified.count", "unclassified")
        unclassified_col <- NA
        for (pattern in unclassified_patterns) {
          if (pattern %in% headers) {
            unclassified_col <- which(headers == pattern)
            break
          }
        }
        
        if (!is.na(unclassified_col)) {
          conditionalFormatting(wb, "Kraken_Stats", 
                              cols = unclassified_col, 
                              rows = data_rows,
                              type = "colorScale", 
                              style = c("#FFFFFF", "#FF0000"))
          cat("Applied red-white color scale to unclassified column\n")
        }
        
        # Data bars for species_count
        species_patterns <- c("species_count", "species.count", "species")
        species_col <- NA
        for (pattern in species_patterns) {
          if (pattern %in% headers) {
            species_col <- which(headers == pattern)
            break
          }
        }
        
        if (!is.na(species_col)) {
          conditionalFormatting(wb, "Kraken_Stats", 
                              cols = species_col, 
                              rows = data_rows,
                              type = "dataBar",
                              style = "#0070C0")
          cat("Applied data bars to species_count column\n")
        }
        
        # Auto-size columns
        setColWidths(wb, "Kraken_Stats", cols = 1:length(headers), widths = "auto")
        
        # Add some basic formatting to headers
        header_style <- createStyle(
          textDecoration = "bold",
          fgFill = "#E6E6FA",
          border = "TopBottomLeftRight",
          borderStyle = "thin"
        )
        addStyle(wb, "Kraken_Stats", header_style, rows = 1, cols = 1:length(headers))
        border_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "thin")
        # Add borders to all data cells
        for (col in 1:length(headers)) {
          addStyle(wb, "Kraken_Stats", border_style, rows = data_rows, cols = col)
        }
        
        # Save workbook
        saveWorkbook(wb, excel_file, overwrite = TRUE)
        cat("Kraken statistics exported to Excel:", excel_file, "\n")
      } else {
        cat("Excel export skipped - openxlsx not available\n")
      }
    }, error = function(e) {
      cat("Warning: Failed to modify data or export to Excel:", e$message, "\n")
      print(e)
    })
    
  } else {
    # Per-sample mode - generate plots for individual sample
    sample_merged_long <- unaligned_results  # This is the per-sample data
    sample_name <- sample_merged_long$sample[1]
    kraken_plots <- list()
    has_nonhuman <- nonhuman_results

    # Summary plot
    plot.new()
    title_text <- paste("Part 3: Kraken2 Statistics -", sample_name)
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)

    # Create reference text for individual sample
    ref_text <- paste0(
      "The Kraken2 results for sample ", sample_name, ":\n"
    )

    if (has_nonhuman) {
      ref_text <- paste0(ref_text, "- Using Non-human sample results\n")
    } else {
      ref_text <- paste0(ref_text, "- Using Unaligned sample results only (no nonhuman data available)\n")
    }

    ref_text <- paste0(ref_text,
      "Key findings:\n",
      "- Number of species found before secondary filtering: ", combined_report_data$species_count[1], "\n",
      "- Number of species found after secondary filtering: ", combined_report_data$filtered_species_count[1], "\n"
    )

    if (config$INCLUDE_SUBSPECIES) {
      ref_text <- paste0(ref_text, 
        "- Subspecies before secondary filtering: ", combined_report_data$subspecies_count[1], "\n",
        "- Subspecies after secondary filtering: ", combined_report_data$filtered_subspecies_count[1], "\n"
      )
    }

    top_clades <- sample_merged_long %>% 
      filter(!is.na(cladeReads)) %>%
      arrange(desc(cladeReads)) %>%
      slice_head(n = min(config$TOP_N_FREQ, nrow(.)))

    ref_text <- paste0(ref_text, 
      "\nTop species in this sample:\n",
      paste(top_clades$name, top_clades$cladeReads, sep = ": ", collapse = "\n")
    )

    text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    kraken_plots[[paste0(sample_name, "_summary")]] <- recordPlot()
    
    # Filter out TaxID 9606 (human) and calculate percentages for kingdom plot
    sample_filtered <- sample_merged_long %>%
      filter(taxID != 9606 & !is.na(cladeReads))
    
    if (nrow(sample_filtered) > 0) {
      # Use bracken reads if available, otherwise cladeReads
      read_col <- if ("bracken_reads" %in% colnames(sample_filtered) && 
                     any(!is.na(sample_filtered$bracken_reads))) {
        "bracken_reads"
      } else {
        "cladeReads"
      }
      
      # Calculate percentages
      total_reads <- sum(sample_filtered[[read_col]], na.rm = TRUE)
      sample_filtered$percentage <- (sample_filtered[[read_col]] / total_reads) * 100
      
      # Source the function-based plotting instead of run_as_script (disable main execution)
      assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
      source(file.path(config$RANALYSIS_DIR, "scripts", "annotate_species_plots.R"))
      
      # Kingdom basic plot using function interface
      kingdom_basic <- generate_annotation_plot(
        df = sample_filtered,
        plot_type = "kingdom",
        sample_name = sample_name
      )
      kraken_plots[[paste0(sample_name, "_kingdom_basic")]] <- kingdom_basic
      
      # Kingdom detailed plot using percentage
      detailed_kingdom_plot <- generate_annotation_plot(
        df = sample_filtered,
        plot_type = "kingdom",
        detailed = "percentage",
        sample_name = sample_name,
        top_n = config$TOP_N_PLOT
      )
      kraken_plots[[paste0(sample_name, "_kingdom_detailed")]] <- detailed_kingdom_plot
      
    } else {
      # No data available - create placeholder plots
      empty_plot <- ggplot() + 
        theme_void() + 
        labs(title = paste("No data available for", sample_name))
      
      kraken_plots[[paste0(sample_name, "_kingdom_basic")]] <- empty_plot
      kraken_plots[[paste0(sample_name, "_kingdom_detailed")]] <- empty_plot
    }
    
    return(kraken_plots)
  }
}

species_annotation_report <- function(annotated_species, report_dir, batch = TRUE, config) {
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
    
    # Generate Summary Plots
    tryCatch({
      # Source the function-based plotting (disable main execution)
      assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
      source(file.path(config$RANALYSIS_DIR, "scripts", "annotate_species_plots.R"))
      
      # HOMD basic plot
      homd_basic <- generate_annotation_plot(
        df = annotated_species,
        plot_type = "homd"
      )
      print(homd_basic)
      
      # Risk group basic plot
      risk_group_basic <- generate_annotation_plot(
        df = annotated_species,
        plot_type = "risk_group"
      )
      print(risk_group_basic)
      
      # Summary basic plot
      summary_basic <- generate_annotation_plot(
        df = annotated_species,
        plot_type = "summary"
      )
      print(summary_basic)
      
    }, error = function(e) {
      cat("Warning: Error generating annotation plots:", e$message, "\n")
      plot.new()
      text(0.5, 0.5, "Annotation plots could not be generated", cex = 1.2, adj = 0.5)
    })
    
    # Generate Detailed Plots
    tryCatch({
      # Detailed HOMD plot
      detailed_homd_plot <- generate_annotation_plot(
        df = annotated_species,
        plot_type = "homd",
        detailed = "Freq"
      )
      print(detailed_homd_plot)
      
      # Detailed risk group plot
      detailed_risk_group_plot <- generate_annotation_plot(
        df = annotated_species,
        plot_type = "risk_group",
        detailed = "Freq"
      )
      print(detailed_risk_group_plot)
      
    }, error = function(e) {
      cat("Warning: Error generating detailed annotation plots:", e$message, "\n")
      plot.new()
      text(0.5, 0.5, "Detailed annotation plots could not be generated", cex = 1.2, adj = 0.5)
    })
    
    dev.off()
    
    # Export annotated species to Excel with formatting
    tryCatch({
      if (excel_available) {
        excel_file <- file.path(report_dir, "species_annotation_report.xlsx")
        
        # Create workbook
        wb <- createWorkbook()
        addWorksheet(wb, "Species_Annotation")
        
        # Create headers without underscores and periods
        headers <- colnames(annotated_species)
        clean_headers <- gsub("[._]", " ", headers)
        clean_headers <- tools::toTitleCase(clean_headers)
        
        # Prepare data for Excel export - ensure proper data types
        excel_data <- annotated_species
        
        # Convert numeric columns to ensure they're properly formatted
        numeric_cols <- c("Freq", "cladeReads_mean", "cladeReads_max", "bracken_reads_mean", "bracken_reads_max")
        for (col in numeric_cols) {
          if (col %in% colnames(excel_data)) {
            excel_data[[col]] <- as.numeric(excel_data[[col]])
          }
        }
        
        # Ensure character columns are properly formatted
        character_cols <- c("name", "taxID", "taxRank", "taxLineage", "RiskGroup", "HOMD", "HOMD.Category")
        for (col in character_cols) {
          if (col %in% colnames(excel_data)) {
            excel_data[[col]] <- as.character(excel_data[[col]])
          }
        }
        
        # Write headers and data separately
        writeData(wb, "Species_Annotation", t(clean_headers), startRow = 1, colNames = FALSE)
        writeData(wb, "Species_Annotation", excel_data, startRow = 2, colNames = FALSE)
        
        # Auto-size columns
        setColWidths(wb, "Species_Annotation", cols = 1:length(headers), widths = "auto")
        
        # Add header formatting
        header_style <- createStyle(
          textDecoration = "bold",
          fgFill = "#E6E6FA",
          border = "TopBottomLeftRight",
          borderStyle = "thin"
        )
        addStyle(wb, "Species_Annotation", header_style, rows = 1, cols = 1:length(headers))
        
        # Add borders to data cells (apply to each row individually)
        data_rows <- 2:(nrow(excel_data) + 1)
        data_style <- createStyle(
          border = "TopBottomLeftRight",
          borderStyle = "thin"
        )
        for (row in data_rows) {
          addStyle(wb, "Species_Annotation", data_style, rows = row, cols = 1:length(headers))
        }
        
        # Save workbook
        saveWorkbook(wb, excel_file, overwrite = TRUE)
        cat("Species annotation exported to Excel:", excel_file, "\n")
      } else {
        cat("Excel export skipped - openxlsx not available\n")
      }
    }, error = function(e) {
      cat("Warning: Failed to export species annotation to Excel:", e$message, "\n")
      print(e)
    })
    
  } else {
    # Per-sample mode - generate plots for individual sample
    sample_merged_long <- annotated_species  # This is the per-sample data
    sample_name <- sample_merged_long$sample[1]
    annotation_plots <- list()
    
    # Get annotation data for this sample by merging with the global annotated_species
    if (exists("annotated_species", envir = .GlobalEnv)) {
      global_annotated_species <- get("annotated_species", envir = .GlobalEnv)
      
      # Merge annotation columns from global annotated_species by taxID
      sample_annotated <- sample_merged_long %>%
        left_join(
          global_annotated_species %>% 
            select(taxID, RiskGroup, HOMD, HOMD.Category) %>%
            distinct(),
          by = "taxID"
        ) %>%
        # Fill missing annotations
        mutate(
          RiskGroup = ifelse(is.na(RiskGroup), "NotAnnotated", RiskGroup),
          HOMD.Category = ifelse(is.na(HOMD.Category), "NotAnnotated", HOMD.Category)
        )
    } else {
      # Fallback if no global annotations available
      sample_annotated <- sample_merged_long %>%
        mutate(
          RiskGroup = "NotAnnotated",
          HOMD.Category = "NotAnnotated"
        )
    }
    
    # Summary plot
    plot.new()
    title_text <- paste("Part 4: Species Annotation -", sample_name)
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
    
    ref_text <- paste0(
      "Species annotation results for sample ", sample_name, ":\n\n",
      "Key findings:\n",
      "- Total species in sample: ", nrow(sample_annotated), "\n",
      "- Species with risk annotations: ", sum(sample_annotated$RiskGroup != "NotAnnotated"), "\n",
      "- Species with HOMD annotations: ", sum(sample_annotated$HOMD.Category != "NotAnnotated"), "\n",
      "- High-risk species (RG3/4): ", sum(sample_annotated$RiskGroup %in% c("RG3", "RG4")), "\n",
      "- Pathogenic species (HOMD): ", sum(sample_annotated$HOMD.Category %in% c("Pathogen", "Opportunist")), "\n"
    )
    
    text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    annotation_plots[[paste0(sample_name, "_summary")]] <- recordPlot()
    
    # Source the function-based plotting (disable main execution)
    assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
    source(file.path(config$RANALYSIS_DIR, "scripts", "annotate_species_plots.R"))
    
    # Generate annotation plots if we have annotated data
    if (nrow(sample_annotated) > 0) {
      tryCatch({
        # HOMD basic plot
        homd_basic <- generate_annotation_plot(
          df = sample_annotated,
          plot_type = "homd",
          sample_name = sample_name
        )
        annotation_plots[[paste0(sample_name, "_homd_basic")]] <- homd_basic
        
        # Risk group basic plot
        risk_group_basic <- generate_annotation_plot(
          df = sample_annotated,
          plot_type = "risk_group", 
          sample_name = sample_name
        )
        annotation_plots[[paste0(sample_name, "_risk_group_basic")]] <- risk_group_basic
        
        # Summary basic plot
        summary_basic <- generate_annotation_plot(
          df = sample_annotated,
          plot_type = "summary",
          sample_name = sample_name
        )
        annotation_plots[[paste0(sample_name, "_summary_basic")]] <- summary_basic
        
        # Check for read data for detailed plots
        read_col <- if ("brackenReads" %in% colnames(sample_annotated) && 
                       any(!is.na(sample_annotated$brackenReads))) {
          "brackenReads"
        } else if ("cladeReads" %in% colnames(sample_annotated) && 
                  any(!is.na(sample_annotated$cladeReads))) {
          "cladeReads"
        } else {
          NULL
        }
        
        if (!is.null(read_col)) {
          # Detailed HOMD plot
          detailed_homd_plot <- generate_annotation_plot(
            df = sample_annotated,
            plot_type = "homd",
            detailed = read_col,
            sample_name = sample_name,
            top_n = config$TOP_N_PLOT,
            log = TRUE
          )
          annotation_plots[[paste0(sample_name, "_homd_detailed")]] <- detailed_homd_plot
          
          # Detailed risk group plot
          detailed_risk_group_plot <- generate_annotation_plot(
            df = sample_annotated,
            plot_type = "risk_group",
            detailed = read_col,
            sample_name = sample_name,
            top_n = config$TOP_N_PLOT,
            log = TRUE
          )
          annotation_plots[[paste0(sample_name, "_risk_group_detailed")]] <- detailed_risk_group_plot
        }
        
      }, error = function(e) {
        cat("Warning: Error generating annotation plots for", sample_name, ":", e$message, "\n")
        empty_plot <- ggplot() + 
          theme_void() + 
          labs(title = paste("Annotation plots not available for", sample_name))
        annotation_plots[[paste0(sample_name, "_error")]] <- empty_plot
      })
    } else {
      # No annotated data available
      empty_plot <- ggplot() + 
        theme_void() + 
        labs(title = paste("No annotation data available for", sample_name))
      annotation_plots[[paste0(sample_name, "_no_data")]] <- empty_plot
    }
    
    return(annotation_plots)
  }
}

pathogen_detection_report <- function(annotated_species, merged_long_data, report_dir, batch, config) {
  if (batch) {
    pdf(file.path(report_dir, "pathogen_detection_report.pdf"), width = 8, height = 6)
    plot.new()
    
    title_text <- "Part 5: Pathogen Detection (RG3/4 Species)"
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
    
    # Filter for pathogenic species (RG3/4 and non-microbiome)
    pathogen_species <- annotated_species %>%
      filter(
        RiskGroup %in% c("RG3", "RG4") | 
        (is.na(RiskGroup) & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME")) |
        (!RiskGroup %in% c("RG1", "RG2") & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME"))
      ) %>%
      filter(!is.na(name))
    
    # Create reference text
    ref_text <- paste0(
      "Pathogen detection analysis focuses on RG3/4 species and non-microbiome organisms.\n\n",
      "Filtering criteria:\n",
      "- Include: RG3, RG4, risk groups and NotAnnotated species\n",
      "- Exclude: RG1, RG2, risk groups\n",
      "- Exclude: Species annotated as 'Microbiome' in HOMD.Category\n\n",
      "Key findings:\n",
      "- Total species after annotation: ", nrow(annotated_species), "\n",
      "- Potential pathogenic species identified: ", nrow(pathogen_species), "\n"
    )
    
    if (nrow(pathogen_species) > 0) {
      ref_text <- paste0(ref_text,
        "- Most abundant pathogen: ", pathogen_species$name[which.max(pathogen_species$cladeReads_mean)], "\n",
        "- Average pathogen clade reads: ", round(mean(pathogen_species$cladeReads_mean, na.rm = TRUE), 1), "\n"
      )
    }
    
    ref_text <- paste0(ref_text,
      "\nData source: Filtered from annotated_species object\n",
      "Output location: Check outputs/ directory for:\n",
      "  - pathogen_detection_report.xlsx (pathogen summary by sample)\n"
    )
    
    text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    
    # Generate the same annotation plots but filtered for pathogens only
    if (nrow(pathogen_species) > 0) {
      tryCatch({
        # Source the function-based plotting (disable main execution)
        assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
        source(file.path(config$RANALYSIS_DIR, "scripts", "annotate_species_plots.R"))
        
        # Risk group plot for pathogens
        pathogen_risk_plot <- generate_annotation_plot(
          df = pathogen_species,
          plot_type = "risk_group"
        )
        print(pathogen_risk_plot)
        
        # Detailed risk group plot
        detailed_pathogen_risk_plot <- generate_annotation_plot(
          df = pathogen_species,
          plot_type = "risk_group",
          detailed = "cladeReads_mean",
          log = TRUE
        )
        print(detailed_pathogen_risk_plot)
        
      }, error = function(e) {
        cat("Warning: Error generating pathogen plots:", e$message, "\n")
        plot.new()
        text(0.5, 0.5, "Pathogen plots could not be generated", cex = 1.2, adj = 0.5)
      })
    } else {
      plot.new()
      text(0.5, 0.5, "No pathogenic species identified", cex = 1.5, adj = 0.5, col = "darkgreen")
    }
    
    dev.off()
    
    # Generate Excel report with sample-level pathogen percentages
    tryCatch({
      if (excel_available && nrow(pathogen_species) > 0) {
        # Calculate sample-level percentages
        pathogen_sample_data <- calculate_pathogen_percentages(pathogen_species, merged_long_data)
        
        if (nrow(pathogen_sample_data) > 0) {
          excel_file <- file.path(report_dir, "pathogen_detection_report.xlsx")
          
          # Create workbook
          wb <- createWorkbook()
          addWorksheet(wb, "Pathogen_Detection")
          
          # Create clean headers
          headers <- c("Sample Name", "RG3/4 and NotAnnotated Species", "Percentage of Sample Reads")
          
          # Write headers and data
          writeData(wb, "Pathogen_Detection", t(headers), startRow = 1, colNames = FALSE)
          writeData(wb, "Pathogen_Detection", pathogen_sample_data, startRow = 2, colNames = FALSE)
          
          # Auto-size columns
          setColWidths(wb, "Pathogen_Detection", cols = 1:length(headers), widths = "auto")
          
          # Add header formatting
          header_style <- createStyle(
            textDecoration = "bold",
            fgFill = "#FFE6E6",  # Light red for pathogen report
            border = "TopBottomLeftRight",
            borderStyle = "thin"
          )
          addStyle(wb, "Pathogen_Detection", header_style, rows = 1, cols = 1:length(headers))
          
          # Add borders to data cells
          data_rows <- 2:(nrow(pathogen_sample_data) + 1)
          data_style <- createStyle(
            border = "TopBottomLeftRight",
            borderStyle = "thin"
          )
          for (row in data_rows) {
            addStyle(wb, "Pathogen_Detection", data_style, rows = row, cols = 1:length(headers))
          }
          
          # Add conditional formatting for high percentages
          if (nrow(pathogen_sample_data) > 1) {
            conditionalFormatting(wb, "Pathogen_Detection", 
                                cols = 3, 
                                rows = data_rows,
                                type = "colorScale", 
                                style = c("#FFFFFF", "#FF6B6B"))
          }
          
          # Save workbook
          saveWorkbook(wb, excel_file, overwrite = TRUE)
          cat("Pathogen detection report exported to Excel:", excel_file, "\n")
        } else {
          cat("No pathogen detections found in samples\n")
        }
      } else if (!excel_available) {
        cat("Excel export skipped - openxlsx not available\n")
      } else {
        cat("No pathogenic species to report\n")
      }
    }, error = function(e) {
      cat("Warning: Failed to export pathogen detection to Excel:", e$message, "\n")
      print(e)
    })
    
  } else {
    # Per-sample mode - generate plots for individual sample
    sample_merged_long <- merged_long_data  # This is the per-sample data
    sample_name <- sample_merged_long$sample[1]
    pathogen_plots <- list()
    
    # Get annotation data for this sample by merging with the global annotated_species
    if (exists("annotated_species", envir = .GlobalEnv)) {
      global_annotated_species <- get("annotated_species", envir = .GlobalEnv)
      
      # Merge annotation columns from global annotated_species by taxID
      sample_annotated <- sample_merged_long %>%
        left_join(
          global_annotated_species %>% 
            select(taxID, RiskGroup, HOMD, HOMD.Category) %>%
            distinct(),
          by = "taxID"
        ) %>%
        # Fill missing annotations
        mutate(
          RiskGroup = ifelse(is.na(RiskGroup), "NotAnnotated", RiskGroup),
          HOMD.Category = ifelse(is.na(HOMD.Category), "NotAnnotated", HOMD.Category)
        )
    } else {
      # Fallback if no global annotations available
      sample_annotated <- sample_merged_long %>%
        mutate(
          RiskGroup = "NotAnnotated",
          HOMD.Category = "NotAnnotated"
        )
    }
    
    # Filter for pathogenic species (same criteria as batch mode)
    sample_pathogen_species <- sample_annotated %>%
      filter(
        RiskGroup %in% c("RG3", "RG4") | 
        (is.na(RiskGroup) & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME")) |
        (!RiskGroup %in% c("RG1", "RG2") & !HOMD.Category %in% c("Microbiome", "microbiome", "MICROBIOME"))
      ) %>%
      filter(!is.na(name))
    
    # Summary plot
    plot.new()
    title_text <- paste("Part 5: Pathogen Detection (RG3/4 Species) -", sample_name)
    text(0.5, 0.95, title_text, adj = c(0.5, 1), cex = 1.2, font = 2)
    
    ref_text <- paste0(
      "Pathogen detection analysis for sample ", sample_name, ":\n\n",
      "Filtering criteria:\n",
      "- Include: RG3, RG4 risk groups\n",
      "- Exclude: RG1, RG2 risk groups\n",
      "- Exclude: Species annotated as 'Microbiome' in HOMD.Category\n\n",
      "Key findings:\n",
      "- Total species in sample: ", nrow(sample_annotated), "\n",
      "- Potential pathogenic species identified: ", nrow(sample_pathogen_species), "\n"
    )
    
    if (nrow(sample_pathogen_species) > 0) {
      # Find the most abundant pathogen based on available read data
      read_col <- if ("brackenReads" %in% colnames(sample_pathogen_species) && 
                     any(!is.na(sample_pathogen_species$brackenReads))) {
        "brackenReads"
      } else if ("cladeReads" %in% colnames(sample_pathogen_species) && 
                any(!is.na(sample_pathogen_species$cladeReads))) {
        "cladeReads"
      } else {
        NULL
      }
      
      if (!is.null(read_col)) {
        most_abundant_idx <- which.max(sample_pathogen_species[[read_col]])
        ref_text <- paste0(ref_text,
          "- Most abundant pathogen: ", sample_pathogen_species$name[most_abundant_idx], "\n",
          "- Average pathogen reads: ", round(mean(sample_pathogen_species[[read_col]], na.rm = TRUE), 1), "\n"
        )
      }
    }
    
    text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    pathogen_plots[[paste0(sample_name, "_summary")]] <- recordPlot()
    
    # Generate pathogen plots if we have pathogenic species
    if (nrow(sample_pathogen_species) > 0) {
      tryCatch({
        # Source the function-based plotting (disable main execution)
        assign("..main_execution_disabled..", TRUE, envir = .GlobalEnv)
        source(file.path(config$RANALYSIS_DIR, "scripts", "annotate_species_plots.R"))
        
        # Risk group plot for pathogens
        pathogen_risk_plot <- generate_annotation_plot(
          df = sample_pathogen_species,
          plot_type = "risk_group",
          sample_name = paste(sample_name, "- Pathogens")
        )
        pathogen_plots[[paste0(sample_name, "_pathogen_risk_basic")]] <- pathogen_risk_plot
        
        # Check for read data for detailed plots
        read_col <- if ("brackenReads" %in% colnames(sample_pathogen_species) && 
                       any(!is.na(sample_pathogen_species$brackenReads))) {
          "brackenReads"
        } else if ("cladeReads" %in% colnames(sample_pathogen_species) && 
                  any(!is.na(sample_pathogen_species$cladeReads))) {
          "cladeReads"
        } else {
          NULL
        }
        
        if (!is.null(read_col)) {
          # Detailed risk group plot for pathogens
          detailed_pathogen_risk_plot <- generate_annotation_plot(
            df = sample_pathogen_species,
            plot_type = "risk_group",
            detailed = read_col,
            sample_name = paste(sample_name, "- Pathogens"),
            top_n = config$TOP_N_PLOT,
            log = TRUE
          )
          pathogen_plots[[paste0(sample_name, "_pathogen_risk_detailed")]] <- detailed_pathogen_risk_plot
        }
        
      }, error = function(e) {
        cat("Warning: Error generating pathogen plots for", sample_name, ":", e$message, "\n")
        empty_plot <- ggplot() + 
          theme_void() + 
          labs(title = paste("Pathogen plots not available for", sample_name))
        pathogen_plots[[paste0(sample_name, "_error")]] <- empty_plot
      })
    } else {
      # No pathogenic species identified
      no_pathogen_plot <- ggplot() + 
        theme_void() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("No pathogenic species identified in", sample_name), 
                size = 6, color = "darkgreen", hjust = 0.5, vjust = 0.5)
      pathogen_plots[[paste0(sample_name, "_no_pathogens")]] <- no_pathogen_plot
    }
    
    return(pathogen_plots)
  }
}

# Helper function to calculate species percentages by sample
calculate_species_percentages <- function(species_data, merged_long_data) {
  tryCatch({
    # Get list of species taxIDs from the input species data
    species_taxids <- unique(species_data$taxID)
    
    # Determine which read column to use from merged_long_data
    if ("bracken_reads" %in% colnames(merged_long_data) && 
        any(!is.na(merged_long_data$bracken_reads))) {
      read_col <- "bracken_reads"
    } else {
      read_col <- "cladeReads"
    }
    
    # Filter merged_long_data for species and exclude human (taxID 9606) 
    species_reads <- merged_long_data %>%
      filter(taxID %in% species_taxids, taxID != "9606") %>%
      select(sample, name, taxID, all_of(read_col)) %>%
      filter(!is.na(.data[[read_col]]), .data[[read_col]] > 0)
    
    # Calculate total non-human reads per sample
    total_reads_per_sample <- merged_long_data %>%
      filter(taxID != "9606", !is.na(.data[[read_col]])) %>%
      group_by(sample) %>%
      summarise(total_reads = sum(.data[[read_col]], na.rm = TRUE), .groups = 'drop')
    
    # Calculate species percentages
    species_summary <- species_reads %>%
      left_join(total_reads_per_sample, by = "sample") %>%
      mutate(
        percentage = round((.data[[read_col]] / total_reads) * 100, 2)
      ) %>%
      filter(percentage > 0) %>%
      select(sample, name, reads = all_of(read_col), percentage) %>%
      arrange(sample, desc(percentage))
    
    return(species_summary)
    
  }, error = function(e) {
    cat("Error calculating species percentages:", e$message, "\n")
    return(data.frame())
  })
}

# Helper function to calculate pathogen percentages by sample (wrapper for backward compatibility)
calculate_pathogen_percentages <- function(pathogen_species, merged_long_data) {
  pathogen_summary <- calculate_species_percentages(pathogen_species, merged_long_data)
  
  if (nrow(pathogen_summary) > 0) {
    # Format for pathogen detection report
    pathogen_summary <- pathogen_summary %>%
      mutate(percentage_formatted = paste0(percentage, "%")) %>%
      select(`Sample Name` = sample, `RG3/4 and NotAnnotated Species` = name, `Percentage of Sample Reads` = percentage_formatted) %>%
      arrange(desc(as.numeric(gsub("%", "", `Percentage of Sample Reads`))))
  }
  
  return(pathogen_summary)
}

#### MAIN EXECUTION #### 

# Create configuration from command line arguments
config <- create_config(args)

# Setup parallel processing
if (config$CORES > 1) {
  cl <- setup_parallel(config$CORES)
} else {
  cat("Running in single-core mode\n")
}

# Build arguments for output_processing.R
output_args <- build_output_processing_args(config)

cat("Running output_processing.R with arguments:\n")

# Run output processing
run_script_with_args(file.path(config$RANALYSIS_DIR, "scripts", "output_processing.R"), output_args)

if (exists("unaligned_results", envir = .GlobalEnv)) {
  unaligned_results <- get("unaligned_results", envir = .GlobalEnv)
  cat("Unaligned results loaded successfully\n")
} else {
  # Attempt to load in from files
  stop("ERROR: unaligned_results not found after running output_processing.R")
}

if (exists("nonhuman_results", envir = .GlobalEnv)) {
  nonhuman_results <- get("nonhuman_results", envir = .GlobalEnv)
  cat("Non-human results loaded successfully\n")
} else {
  cat("WARNING: nonhuman_results not found after running output_processing.R\n")
  nonhuman_results <- NULL
}

# Run correlation analysis if requested
if (config$PERFORM_CORRELATION) {
  cat("Running correlation analysis...\n")
  
  # Ensure reports directory exists
  if (!dir.exists(config$reports_dir)) {
    dir.create(config$reports_dir, recursive = TRUE)
  }
  
  if (!is.null(nonhuman_results) && exists("nonhuman_results", envir = .GlobalEnv)) {
    correlation_results <- create_correlation_analysis(
      data1 = unaligned_results$merged,
      data2 = nonhuman_results$merged,
      data1_name = "Unaligned",
      data2_name = "Non-Human",
      data1_col_suffix = "cladeReads",
      data2_col_suffix = "cladeReads",
      output_dir = config$reports_dir
    )
    
    # Check if bracken data is available and perform clade vs bracken correlation analysis
    bracken_cols_unaligned <- grep("_bracken_reads$", colnames(unaligned_results$merged), value = TRUE)
    bracken_cols_nonhuman <- if(!is.null(nonhuman_results)) grep("_bracken_reads$", colnames(nonhuman_results$merged), value = TRUE) else character(0)
    
    if (length(bracken_cols_unaligned) > 0) {
      cat("Running clade reads vs bracken reads correlation analysis for unaligned data...\n")
      bracken_correlation_unaligned <- create_correlation_analysis(
        data1 = unaligned_results$merged,
        data2 = unaligned_results$merged,
        data1_name = "Unaligned_clade_reads",
        data2_name = "Unaligned_bracken_reads",
        data1_col_suffix = "cladeReads",
        data2_col_suffix = "bracken_reads",
        output_dir = config$reports_dir
      )
    }
    
    if (length(bracken_cols_nonhuman) > 0) {
      cat("Running clade reads vs bracken reads correlation analysis for non-human data...\n")
      bracken_correlation_nonhuman <- create_correlation_analysis(
        data1 = nonhuman_results$merged,
        data2 = nonhuman_results$merged,
        data1_name = "Non-Human_clade_reads",
        data2_name = "Non-Human_bracken_reads",
        data1_col_suffix = "cladeReads",
        data2_col_suffix = "bracken_reads",
        output_dir = config$reports_dir
      )
    }
    
    # Analysis of host contamination
    cat("Analyzing host contamination...\n")
    host_contamination <- setdiff(unaligned_results$merged$name, nonhuman_results$merged$name)
    cat("Species found only in unaligned data (potential host contamination):", length(host_contamination), "\n")
    
    if (length(host_contamination) > 0) {
      host_contaminated_results <- unaligned_results$merged %>%
        filter(name %in% host_contamination) %>%
        select(name, contains("cladeReads")) %>%
        arrange(desc(rowSums(select(., contains("cladeReads")), na.rm = TRUE)))
      
      cat("Host-contaminated species:\n")
      print(host_contaminated_results$name)
      
      # Save host contamination results
      write.csv(host_contaminated_results, 
                file.path(config$reports_dir, "host_contamination_species.csv"), 
                row.names = FALSE)
    }
  } else {
    cat("Skipping correlation analysis - nonhuman results not available\n")
  }
}

#### Batch Report Exporting ####

# Run runtime statistics report 
if (config$RT_STATS) {
  tryCatch({
    if (exists("runtime_data", envir = .GlobalEnv)) {
      runtime_data <- get("runtime_data", envir = .GlobalEnv)
      cat("Runtime results loaded successfully\n")
    } else {
      runtime_data <- read.csv(get_latest_timestamped_file(input_dir = config$output_dir, pattern = "runtime_info"), 
                               stringsAsFactors = FALSE)
    }
    tryCatch({
      # Add read statistics if available
      if (exists("rrstats_data", envir = .GlobalEnv)) {
        # Add num_input_reads from rrstats_data matching by sample
        rrstats_data_for_rt <- get("rrstats_data", envir = .GlobalEnv)
        runtime_data <- merge(runtime_data, rrstats_data_for_rt[, c("sample", "num_input_reads", "percent_mapped_reads")], 
                              by = "sample", all.x = TRUE)
      } else {
        # Load read statistics data from the latest file
        rrstats_data_for_rt <- read.csv(get_latest_timestamped_file(input_dir = config$output_dir, pattern = "read_statistics"), 
                                 stringsAsFactors = FALSE)
        runtime_data <- merge(runtime_data, rrstats_data_for_rt[, c("sample", "num_input_reads", "percent_mapped_reads")], 
                              by = "sample", all.x = TRUE)
      }
      rm(rrstats_data_for_rt)
    }, error = function(e) {
      cat("Warning: Could not merge read statistics with runtime data:", e$message, "\n")
    })
    
    runtime_analysis_report(runtime_data, config$reports_dir)
  }, error = function(e) {
    cat("Warning: Runtime statistics report skipped - no runtime data available:", e$message, "\n")
  })
}

# Run read statistics report
if (config$RR_STATS) {
  tryCatch({
    if (exists("rrstats_data", envir = .GlobalEnv)) {
      rrstats_data <- get("rrstats_data", envir = .GlobalEnv)
      cat("Read statistics data loaded successfully\n")
    } else {
      # Load read statistics data from the latest file
      cat("Loading read statistics data from the latest generated file...\n")
      rrstats_data <- read.csv(get_latest_timestamped_file(input_dir = config$output_dir, pattern = "read_statistics"), 
                               stringsAsFactors = FALSE)
    }
    runread_stats_report(rrstats_data, config$reports_dir)
  }, error = function(e) {
    cat("Warning: Read statistics report skipped - no read statistics data available:", e$message, "\n")
  })
}

# Run kraken statistics report
if (config$KR_STATS) {
  tryCatch({
    if (exists("combined_report_data", envir = .GlobalEnv)) {
      combined_report_data <- get("combined_report_data", envir = .GlobalEnv)
      cat("Combined report data loaded successfully\n")
    } else {
      combined_report_data <- read.csv(get_latest_timestamped_file(input_dir = config$output_dir, pattern = "sample_report_data"), 
                                       stringsAsFactors = FALSE)
    }
    
    kraken_stats_report(unaligned_results, nonhuman_results, combined_report_data, config, config$reports_dir, TRUE)
  }, error = function(e) {
    cat("Warning: Kraken statistics report skipped - no kraken data available:", e$message, "\n")
  })
}

species_df <- unaligned_results$species_list
assign("species_df", species_df, envir = .GlobalEnv)
# Run annotation with custom databases directory
run_script_with_args(file.path(config$RANALYSIS_DIR, "scripts", "annotate_species.R"), 
                     c("--df", "species_df", "--output", "annotated_species", 
                       "--databases-dir", config$DATABASES_DIR))

# Run species annotation if annotated species_list is available
if (config$SA_STATS) {
  cat("Running species annotation...\n")

  if (exists("annotated_species", envir = .GlobalEnv)) {
    annotated_species <- get("annotated_species", envir = .GlobalEnv)
    species_annotation_report(annotated_species, config$reports_dir, TRUE, config)
    cat("Species annotation completed successfully\n")

  } else {
    stop("ERROR: annotated_species not found after running annotate_species.R")
  }
}

# Run pathogen detection if requested & annotated species_list is available
if (config$PD_STATS) {
  if (exists("annotated_species", envir = .GlobalEnv)) {
    annotated_species <- get("annotated_species", envir = .GlobalEnv)
    cat("Running pathogen detection analysis...\n")
    # Get merged_long data for sample-level calculations
    merged_long_data <- if (!is.null(nonhuman_results) && !is.null(nonhuman_results$merged_long)) {
      nonhuman_results$merged_long
    } else {
      unaligned_results$merged_long
    }
    pathogen_detection_report(annotated_species, merged_long_data, config$reports_dir, TRUE, config)
  } else {
    stop("ERROR: annotated_species not found after running annotate_species.R")
  }
}

#### Per Sample Report Generation ####
if (!config$NO_PER_SAMPLE) {
  cat("Generating per-sample reports...\n")
  
  # Use nonhuman results if available, otherwise unaligned
  if (!is.null(nonhuman_results) && !is.null(nonhuman_results$merged_long)) {
    merged_long_data <- nonhuman_results$merged_long
    using_nonhuman_per_sample <- TRUE
  } else {
    merged_long_data <- unaligned_results$merged_long
    using_nonhuman_per_sample <- FALSE
  }
  
  # Filter out NA values in cladeReads
  trimmed_long <- merged_long_data %>% filter(!is.na(cladeReads))
  # filter out human reads (taxID 9606)
  trimmed_long <- trimmed_long %>% filter(taxID != "9606")
  
  # Split by sample
  trimmed_by_sample <- split(trimmed_long, trimmed_long$sample)
  trimmed_by_sample <- lapply(trimmed_by_sample, as.data.frame)

  # Create per-sample report directories
  per_sample_dir <- file.path(config$reports_dir, "per_sample_reports")
  if (!dir.exists(per_sample_dir)) {
    dir.create(per_sample_dir, recursive = TRUE)
  }

  # Helper functions to generate per-sample reports
  generate_per_sample_kraken_report_plots <- function(combined_report_data, sample_name, sample_data, config) {
    sample_report_data <- as.data.frame(combined_report_data) %>%
      filter(sample == sample_name)
    return(kraken_stats_report(
      unaligned_results = sample_data,
      nonhuman_results = using_nonhuman_per_sample,
      combined_report_data = sample_report_data,
      config = config,
      report_dir = NULL,
      batch = FALSE)
    )
  }
  
  generate_per_sample_species_annotation_plots <- function(sample_name, sample_data, config) {
    return(species_annotation_report(
      annotated_species = sample_data,
      report_dir = NULL,
      batch = FALSE,
      config = config)
    )
  }
  
  generate_per_sample_pathogen_detection_plots <- function(sample_name, sample_data, config) {
    return(pathogen_detection_report(
      annotated_species = NULL,  # Will be pulled from global
      merged_long_data = sample_data,
      report_dir = NULL,
      batch = FALSE,
      config = config)
    )
  }

  if (config$SPLIT_REPORT) {
    # Generate individual PDF reports per sample
    
    if (exists("cl")) {
      # Generating reports in Parallel
      cat("Generating per-sample reports in parallel using", length(get("cl", envir = .GlobalEnv)), "cores...\n")
      
      cluster <- get("cl", envir = .GlobalEnv)
      
      # Export required variables to cluster nodes
      clusterExport(cluster, c("trimmed_by_sample", "config", "combined_report_data", 
                   "per_sample_dir", "generate_per_sample_kraken_report_plots",
                   "generate_per_sample_species_annotation_plots", 
                   "generate_per_sample_pathogen_detection_plots",
                   "using_nonhuman_per_sample", "annotated_species"), 
            envir = environment())
      
      # Process samples in parallel using foreach
      foreach(sample_name = names(trimmed_by_sample), 
          .packages = c("dplyr", "ggplot2", "here"),
          .combine = 'c') %dopar% {
      
      cat("Generating report for sample:", sample_name, "\n")
      sample_data <- trimmed_by_sample[[sample_name]]
      sample_report_dir <- file.path(per_sample_dir, sample_name)
      if (!dir.exists(sample_report_dir)) {
        dir.create(sample_report_dir, recursive = TRUE)
      }
      
      # Generate Kraken report
      if (config$KR_STATS) {
        tryCatch({
        kraken_plots <- generate_per_sample_kraken_report_plots(combined_report_data, sample_name, sample_data, config)
        
        pdf(file = file.path(sample_report_dir, paste0(sample_name, "_kraken_report.pdf")), width = 8, height = 6)
        if (!is.null(kraken_plots) && length(kraken_plots) > 0) {
          for (plot in kraken_plots) {
          print(plot)
          }
        } else {
          plot.new()
          text(0.5, 0.5, "Kraken plots not available for this sample", cex = 1.2, adj = 0.5)
        }
        dev.off()
        }, error = function(e) {
        cat("Error generating Kraken report for", sample_name, ":", e$message, "\n")
        })
      }
      
      # Generate Species Annotation report
      if (config$SA_STATS) {
        tryCatch({
        annotation_plots <- generate_per_sample_species_annotation_plots(sample_name, sample_data, config)
        
        pdf(file = file.path(sample_report_dir, paste0(sample_name, "_species_annotation_report.pdf")), width = 8, height = 6)
        if (!is.null(annotation_plots) && length(annotation_plots) > 0) {
          for (plot in annotation_plots) {
          print(plot)
          }
        } else {
          plot.new()
          text(0.5, 0.5, "Species annotation plots not available for this sample", cex = 1.2, adj = 0.5)
        }
        dev.off()
        }, error = function(e) {
        cat("Error generating species annotation report for", sample_name, ":", e$message, "\n")
        })
      }
      
      # Generate Pathogen Detection report
      if (config$PD_STATS) {
        tryCatch({
        pathogen_plots <- generate_per_sample_pathogen_detection_plots(sample_name, sample_data, config)
        
        pdf(file = file.path(sample_report_dir, paste0(sample_name, "_pathogen_detection_report.pdf")), width = 8, height = 6)
        if (!is.null(pathogen_plots) && length(pathogen_plots) > 0) {
          for (plot in pathogen_plots) {
          print(plot)
          }
        } else {
          plot.new()
          text(0.5, 0.5, "Pathogen detection plots not available for this sample", cex = 1.2, adj = 0.5)
        }
        dev.off()
        }, error = function(e) {
        cat("Error generating pathogen detection report for", sample_name, ":", e$message, "\n")
        })
      }
      
      # Return sample name for tracking completion
      sample_name
      }
      
      cat("Parallel per-sample report generation completed\n")
      
    } else {
      for (sample_name in names(trimmed_by_sample)) {
      cat("Generating report for sample:", sample_name, "\n")
      sample_data <- trimmed_by_sample[[sample_name]]
      sample_report_dir <- file.path(per_sample_dir, sample_name)
      if (!dir.exists(sample_report_dir)) {
        dir.create(sample_report_dir, recursive = TRUE)
      }
      
      # Generate Kraken report
      if (config$KR_STATS) {
        kraken_plots <- generate_per_sample_kraken_report_plots(combined_report_data, sample_name, sample_data, config)
        
        pdf(file = file.path(sample_report_dir, paste0(sample_name, "_kraken_report.pdf")), width = 8, height = 6)
        if (!is.null(kraken_plots) && length(kraken_plots) > 0) {
          for (plot in kraken_plots) {
            print(plot)
          }
        } else {
        plot.new()
        text(0.5, 0.5, "Kraken plots not available for this sample", cex = 1.2, adj = 0.5)
        }
        dev.off()
      }
      
      # Generate Species Annotation report
      if (config$SA_STATS) {
        annotation_plots <- generate_per_sample_species_annotation_plots(sample_name, sample_data, config)
        
        pdf(file = file.path(sample_report_dir, paste0(sample_name, "_species_annotation_report.pdf")), width = 8, height = 6)
        if (!is.null(annotation_plots) && length(annotation_plots) > 0) {
        for (plot in annotation_plots) {
          print(plot)
        }
        } else {
        plot.new()
        text(0.5, 0.5, "Species annotation plots not available for this sample", cex = 1.2, adj = 0.5)
        }
        dev.off()
      }
      
      # Generate Pathogen Detection report
      if (config$PD_STATS) {
        pathogen_plots <- generate_per_sample_pathogen_detection_plots(sample_name, sample_data, config)
        
        pdf(file = file.path(sample_report_dir, paste0(sample_name, "_pathogen_detection_report.pdf")), width = 8, height = 6)
        if (!is.null(pathogen_plots) && length(pathogen_plots) > 0) {
        for (plot in pathogen_plots) {
          print(plot)
        }
        } else {
        plot.new()
        text(0.5, 0.5, "Pathogen detection plots not available for this sample", cex = 1.2, adj = 0.5)
        }
        dev.off()
      }
      }
    }
  } else { 
    # Combine all reports into one PDF per report type
    cat("Generating combined per-sample reports...\n")
    
    # Generate all sample plots for each report type
    if (config$KR_STATS) {
      if (exists("cl") && !is.null(get("cl", envir = .GlobalEnv))) {
        # Parallel processing for Kraken plots
        cat("Generating Kraken plots in parallel using", length(get("cl", envir = .GlobalEnv)), "cores...\n")
        cluster <- get("cl", envir = .GlobalEnv)
        # Export required variables to cluster nodes
        clusterExport(cluster, c("trimmed_by_sample", "config", "combined_report_data", 
                                 "generate_per_sample_kraken_report_plots"), 
                      envir = environment())
        # Process samples in parallel using foreach
        all_samples_kraken_plots <- foreach(sample_name = names(trimmed_by_sample), 
                                           .packages = c("dplyr", "ggplot2", "here"),
                                           .combine = 'c') %dopar% {
          sample_data <- trimmed_by_sample[[sample_name]]
          sample_plots <- generate_per_sample_kraken_report_plots(combined_report_data, sample_name, sample_data, config)
          list(sample_plots)
        }
        names(all_samples_kraken_plots) <- names(trimmed_by_sample)
      } else {
        # Sequential processing fallback
        all_samples_kraken_plots <- list()
        for (sample_name in names(trimmed_by_sample)) {
          sample_data <- trimmed_by_sample[[sample_name]]
          sample_plots <- generate_per_sample_kraken_report_plots(combined_report_data, sample_name, sample_data, config)
          all_samples_kraken_plots[[sample_name]] <- sample_plots
        }
      }
      
      # Save all kraken plots to one PDF
      cat("Saving combined Kraken report to PDF...\n")
      pdf(file = file.path(per_sample_dir, "all_samples_kraken_report.pdf"), width = 8, height = 6)
      for (sample_name in names(all_samples_kraken_plots)) {
        sample_plots <- all_samples_kraken_plots[[sample_name]]
        if (!is.null(sample_plots) && length(sample_plots) > 0) {
          for (plot in sample_plots) {
            print(plot)
          }
        } else {
          plot.new()
          text(0.5, 0.5, paste("Kraken plots not available for", sample_name), cex = 1.2, adj = 0.5)
        }
      }
      dev.off()
      cat("Combined Kraken report saved to PDF\n")
    }
    
    if (config$SA_STATS) {
      if (exists("cl") && !is.null(get("cl", envir = .GlobalEnv))) {
        # Parallel processing for Species Annotation plots
        cat("Generating Species Annotation plots in parallel using", length(get("cl", envir = .GlobalEnv)), "cores...\n")
        cluster <- get("cl", envir = .GlobalEnv)
        # Export required variables to cluster nodes
        clusterExport(cluster, c("trimmed_by_sample", "config", 
                                 "generate_per_sample_species_annotation_plots", "annotated_species"), 
                      envir = environment())
        # Process samples in parallel using foreach
        all_samples_species_annotation_plots <- foreach(sample_name = names(trimmed_by_sample), 
                                           .packages = c("dplyr", "ggplot2", "here"),
                                           .combine = 'c') %dopar% {
          sample_data <- trimmed_by_sample[[sample_name]]
          sample_plots <- generate_per_sample_species_annotation_plots(sample_name, sample_data, config)
          list(sample_plots)
        }
        names(all_samples_species_annotation_plots) <- names(trimmed_by_sample)
      } else {
        all_samples_species_annotation_plots <- list()
        for (sample_name in names(trimmed_by_sample)) {
          sample_data <- trimmed_by_sample[[sample_name]]
          sample_plots <- generate_per_sample_species_annotation_plots(sample_name, sample_data, config)
          all_samples_species_annotation_plots[[sample_name]] <- sample_plots
        }
      }
      
      # Save all annotation plots to one PDF
      cat("Saving combined Species Annotation report to PDF...\n")
      pdf(file = file.path(per_sample_dir, "all_samples_species_annotation_report.pdf"), width = 8, height = 6)
      for (sample_name in names(all_samples_species_annotation_plots)) {
        sample_plots <- all_samples_species_annotation_plots[[sample_name]]
        if (!is.null(sample_plots) && length(sample_plots) > 0) {
          for (plot in sample_plots) {
            print(plot)
          }
        } else {
          plot.new()
          text(0.5, 0.5, paste("Species annotation plots not available for", sample_name), cex = 1.2, adj = 0.5)
        }
      }
      dev.off()
      cat("Species Annotation report saved to PDF\n")
    }
    
    if (config$PD_STATS) {
      if (exists("cl") && !is.null(get("cl", envir = .GlobalEnv))) {
        # Parallel processing for Pathogen Detection plots
        cat("Generating Pathogen Detection plots in parallel using", length(get("cl", envir = .GlobalEnv)), "cores...\n")
        cluster <- get("cl", envir = .GlobalEnv)
        # Export required variables to cluster nodes
        clusterExport(cluster, c("trimmed_by_sample", "config", 
                                 "generate_per_sample_pathogen_detection_plots", "annotated_species"), 
                      envir = environment())
        # Process samples in parallel using foreach
        all_samples_pathogen_plots <- foreach(sample_name = names(trimmed_by_sample), 
                                               .packages = c("dplyr", "ggplot2", "here"),
                                               .combine = 'c') %dopar% {
          sample_data <- trimmed_by_sample[[sample_name]]
          sample_plots <- generate_per_sample_pathogen_detection_plots(sample_name, sample_data, config)
          list(sample_plots)
        }
        names(all_samples_pathogen_plots) <- names(trimmed_by_sample)
      } else {
        all_samples_pathogen_plots <- list()
        for (sample_name in names(trimmed_by_sample)) {
          sample_data <- trimmed_by_sample[[sample_name]]
          sample_plots <- generate_per_sample_pathogen_detection_plots(sample_name, sample_data, config)
          all_samples_pathogen_plots[[sample_name]] <- sample_plots
        }
      }
      # Save all pathogen plots to one PDF
      cat("Saving combined Pathogen Detection report to PDF...\n")
      pdf(file = file.path(per_sample_dir, "all_samples_pathogen_detection_report.pdf"), width = 8, height = 6)
      for (sample_name in names(all_samples_pathogen_plots)) {
        sample_plots <- all_samples_pathogen_plots[[sample_name]]
        if (!is.null(sample_plots) && length(sample_plots) > 0) {
          for (plot in sample_plots) {
            print(plot)
          }
        } else {
          plot.new()
          text(0.5, 0.5, paste("Pathogen detection plots not available for", sample_name), cex = 1.2, adj = 0.5)
        }
      }
      dev.off()
      cat("Pathogen Detection report saved to PDF\n")
    }
  }
}

# Clean up parallel cluster
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
  cat("Parallel cluster stopped\n")
}

cat("Analysis completed successfully!\n")
cat("Output directory:", config$output_dir, "\n")
cat("Reports directory:", config$reports_dir, "\n")


