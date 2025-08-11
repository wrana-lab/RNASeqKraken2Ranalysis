#!/usr/bin/env Rscript

#### Header ####
# Output Processing Script
# Author: Denis Rivard, modified by GitHub Copilot
# Description: Improved and simplified kraken/bracken output processing script
# Compatible with various kreport file naming patterns

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
  cat("Usage: Rscript output_processing_v5.R [OPTIONS]\n\n")
  cat("Options:\n")
  cat("  --exclude-taxid <ID>       Exclude species with this taxonomy ID (e.g., 9606 for Homo sapiens)\n")
  cat("  --no-bracken              Skip bracken file processing\n")
  cat("  --min-reads <N>           Minimum clade reads threshold (default: 0)\n")
  cat("  --top-n <N>               Number of top species to analyze (default: 3)\n")
  cat("  --unaligned-kreports <DIR> Path to unaligned kreports folder (default: NULL)\n")
  cat("  --unaligned-bracken <DIR>  Path to unaligned bracken folder (default: NULL)\n")
  cat("  --nonhuman-kreports <DIR>  Path to nonhuman kreports folder (default: NULL)\n")
  cat("  --nonhuman-bracken <DIR>   Path to nonhuman bracken folder (default: NULL)\n")
  cat("  --runtime-dir <DIR>       Path to runtime folder for parsing runtime files (default: NULL)\n")
  cat("  --rrstats-dir <DIR>       Path to read statistics folder for parsing rrstats files (default: NULL)\n")
  cat("  --metadata-dir <DIR>      Path to directory containing sample metadata CSV files (default: NULL)\n")
  cat("  --output-dir <DIR>         Output directory (default: outputs/full_run)\n")
  cat("  --output                  Save results to CSV files (default: FALSE)\n")
  cat("  --params                  Add confidence_levels, minimum_hit_groups, and human_reads columns to long data (default: FALSE)\n")
  cat("  --subspecies              Include subspecies (S1, S2, S3) in addition to species (default: FALSE)\n")
  cat("  --minimizer-ratio <R>     Filter by minimum ratio of distinct_minimizers/cladeReads (default: NULL)\n")
  cat("  --minimizer-threshold <N> Filter by minimum distinct_minimizers threshold (default: NULL)\n")
  cat("  --parallel                Use parallel processing with cluster from global environment 'cl' (default: FALSE)\n")
  cat("  --help, -h                Show this help message\n\n")
  cat("Examples:\n")
  cat("  Rscript output_processing_v5.R                          # Default settings\n")
  cat("  Rscript output_processing_v5.R --exclude-taxid 9606     # Exclude Homo sapiens\n")
  cat("  Rscript output_processing_v5.R --no-bracken             # Skip bracken processing\n")
  cat("  Rscript output_processing_v5.R --unaligned-kreports /path/to/kreports # Custom kreports path\n")
  cat("  Rscript output_processing_v5.R --output-dir /path/to/output # Custom output directory\n")
  cat("  Rscript output_processing_v5.R --output                 # Save results to CSV files\n")
  cat("  Rscript output_processing_v5.R --params                 # Add parameter columns to long data\n")
  cat("  Rscript output_processing_v5.Rs --subspecies             # Include subspecies in analysis\n")
  cat("  Rscript output_processing_v5.R --spec-list              # Output complete species list\n")
  cat("  Rscript output_processing_v5.R --minimizer-ratio 0.1    # Filter by minimizer ratio >= 0.1\n")
  cat("  Rscript output_processing_v5.R --minimizer-threshold 5  # Filter by distinct_minimizers >= 5\n")
  cat("  Rscript output_processing_v5.R --parallel               # Use parallel processing\n")
  cat("  Rscript output_processing_v5.R --runtime-dir /path/to/runtimes # Parse runtime files\n")
  cat("  Rscript output_processing_v5.R --rrstats-dir /path/to/rrstats # Parse read statistics files\n")
  cat("  Rscript output_processing_v5.R --metadata-dir /path/to/metadata # Parse sample metadata files\n\n")
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

# Parse command line arguments for configuration
if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--exclude-taxid" && i < length(args)) {
      EXCLUDE_TAXID <- as.numeric(args[i + 1])
      cat("Excluding taxID:", EXCLUDE_TAXID, "\n")
      i <- i + 2
    } else if (args[i] == "--no-bracken") {
      PROCESS_BRACKEN <- FALSE
      cat("Skipping bracken processing\n")
      i <- i + 1
    } else if (args[i] == "--min-reads" && i < length(args)) {
      MIN_CLADE_READS <- as.numeric(args[i + 1])
      cat("Minimum clade reads set to:", MIN_CLADE_READS, "\n")
      i <- i + 2
    } else if (args[i] == "--top-n" && i < length(args)) {
      TOP_SPECIES_N <- as.numeric(args[i + 1])
      cat("Top N species set to:", TOP_SPECIES_N, "\n")
      i <- i + 2
    } else if (args[i] == "--unaligned-kreports" && i < length(args)) {
      UNALIGNED_KREPORTS_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Unaligned kreports directory set to:", UNALIGNED_KREPORTS_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--unaligned-bracken" && i < length(args)) {
      UNALIGNED_BRACKEN_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Unaligned bracken directory set to:", UNALIGNED_BRACKEN_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--nonhuman-kreports" && i < length(args)) {
      NONHUMAN_KREPORTS_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Nonhuman kreports directory set to:", NONHUMAN_KREPORTS_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--nonhuman-bracken" && i < length(args)) {
      NONHUMAN_BRACKEN_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Nonhuman bracken directory set to:", NONHUMAN_BRACKEN_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--runtime-dir" && i < length(args)) {
      RUNTIME_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Runtime directory set to:", RUNTIME_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--rrstats-dir" && i < length(args)) {
      RRSTATS_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Read statistics directory set to:", RRSTATS_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--metadata-dir" && i < length(args)) {
      METADATA_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Sample metadata directory set to:", METADATA_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--output-dir" && i < length(args)) {
      OUTPUT_DIR <- normalizePath(args[i + 1], mustWork = FALSE)
      cat("Output directory set to:", OUTPUT_DIR, "\n")
      i <- i + 2
    } else if (args[i] == "--output") {
      SAVE_OUTPUT <- TRUE
      cat("CSV output enabled\n")
      i <- i + 1
    } else if (args[i] == "--params") {
      ADD_PARAMS <- TRUE
      cat("Parameter columns will be added to long data\n")
      i <- i + 1
    } else if (args[i] == "--subspecies") {
      INCLUDE_SUBSPECIES <- TRUE
      cat("Subspecies will be included in analysis\n")
      i <- i + 1
    } else if (args[i] == "--minimizer-ratio" && i < length(args)) {
      MINIMIZER_RATIO <- as.numeric(args[i + 1])
      cat("Minimizer ratio filtering set to:", MINIMIZER_RATIO, "\n")
      i <- i + 2
    } else if (args[i] == "--minimizer-threshold" && i < length(args)) {
      MINIMIZER_THRESHOLD <- as.numeric(args[i + 1])
      cat("Minimizer threshold filtering set to:", MINIMIZER_THRESHOLD, "\n")
      i <- i + 2
    } else if (args[i] == "--parallel") {
      USE_PARALLEL <- TRUE
      cat("Parallel processing enabled\n")
      i <- i + 1
    } else {
      cat("Unknown argument:", args[i], "\n")
      show_help()
      stop("Help requested. Execution stopped.")
    }
  }
}
rm(show_help)

# Display final configuration
cat("\n=== Configuration ===\n")
cat("Minimum clade reads:", MIN_CLADE_READS, "\n")
cat("Top N species:", TOP_SPECIES_N, "\n")
cat("Process bracken:", PROCESS_BRACKEN, "\n")
cat("Save CSV output:", SAVE_OUTPUT, "\n")
cat("Add parameter columns:", ADD_PARAMS, "\n")
cat("Include subspecies:", INCLUDE_SUBSPECIES, "\n")
cat("Minimizer ratio threshold:", if(is.null(MINIMIZER_RATIO)) "None" else MINIMIZER_RATIO, "\n")
cat("Minimizer absolute threshold:", if(is.null(MINIMIZER_THRESHOLD)) "None" else MINIMIZER_THRESHOLD, "\n")
cat("Use parallel processing:", USE_PARALLEL, "\n")
if (is.null(EXCLUDE_TAXID)) {
  cat("Exclude taxID: None (including all species)\n")
} else {
  cat("Exclude taxID:", EXCLUDE_TAXID, "\n")
}
cat("\nInput directories:\n")
cat("  Unaligned kreports:", UNALIGNED_KREPORTS_DIR, "\n")
cat("  Unaligned bracken:", UNALIGNED_BRACKEN_DIR, "\n")
cat("  Nonhuman kreports:", NONHUMAN_KREPORTS_DIR, "\n")
cat("  Nonhuman bracken:", NONHUMAN_BRACKEN_DIR, "\n")
cat("  Runtime directory:", RUNTIME_DIR, "\n")
cat("  Read statistics directory:", RRSTATS_DIR, "\n")
cat("  Sample metadata directory:", METADATA_DIR, "\n")
cat("  Output directory:", OUTPUT_DIR, "\n")
cat("=====================\n\n")

# Setup parallel processing if requested
if (USE_PARALLEL) {
  if (exists("cl", envir = .GlobalEnv)) {
    cat("Setting up parallel processing with cluster from global environment\n")
    registerDoParallel(cl)
    cat("Number of workers:", getDoParWorkers(), "\n")
  } else {
    warning("Parallel processing requested but no cluster 'cl' found in global environment")
    cat("To create a cluster, run: cl <- parallel::makeCluster(8) # or desired number of cores\n")
    USE_PARALLEL <- FALSE
    cat("Falling back to sequential processing\n")
  }
}

#### Logic ####
# Enhanced file parsing function that handles multiple naming patterns
parse_sample_info <- function(filename) {
  # Remove file extension
  sample_name <- tools::file_path_sans_ext(basename(filename))
  
  # Initialize default values
  result <- list(
    original_name = sample_name,
    trimmed_name = sample_name,
    confidence = NA,
    minimum_hits = NA,
    database_used = "unknown"
  )
  
  # Pattern 1: Standard pattern like "W50504632_S109_c45m30stdk2dbdl_report"
  if (grepl("_S\\d+_c\\d+m\\d+", sample_name)) {
    result$trimmed_name <- sub("(_S\\d+).*", "", sample_name)
    conf_min_db <- regmatches(sample_name, regexpr("c\\d+m\\d+[A-Za-z0-9_]*", sample_name))
    if (length(conf_min_db) > 0) {
      result$confidence <- as.integer(sub("c(\\d+).*", "\\1", conf_min_db))
      result$minimum_hits <- as.integer(sub(".*m(\\d+).*", "\\1", conf_min_db))
      db_part <- sub(".*m\\d+([A-Za-z0-9_]+).*", "\\1", conf_min_db)
      result$database_used <- sub("_report$", "", db_part)  # Remove trailing _report
    }
  }
  # Pattern 2: Pattern like "temp_ZQ3002179_B_c0m45stdk2dbdl_report" or "ZQ1400570_A_Pos_c45m30stdk2dbdl_report"
  else if (grepl("_[A-Z](_[A-Za-z]*)?_c\\d+m\\d+", sample_name)) {
    result$trimmed_name <- sub("_c\\d+.*", "", sample_name)
    conf_min_db <- regmatches(sample_name, regexpr("c\\d+m\\d+[A-Za-z0-9_]*", sample_name))
    if (length(conf_min_db) > 0) {
      result$confidence <- as.integer(sub("c(\\d+).*", "\\1", conf_min_db))
      result$minimum_hits <- as.integer(sub(".*m(\\d+).*", "\\1", conf_min_db))
      db_part <- sub(".*m\\d+([A-Za-z0-9_]+).*", "\\1", conf_min_db)
      result$database_used <- sub("_report$", "", db_part)  # Remove trailing _report
    }
  }
  # Pattern 3: General pattern that captures any filename with c#m# format
  else if (grepl("c\\d+m\\d+", sample_name)) {
    result$trimmed_name <- sub("_c\\d+.*", "", sample_name)
    conf_min_db <- regmatches(sample_name, regexpr("c\\d+m\\d+[A-Za-z0-9_]*", sample_name))
    if (length(conf_min_db) > 0) {
      result$confidence <- as.integer(sub("c(\\d+).*", "\\1", conf_min_db))
      result$minimum_hits <- as.integer(sub(".*m(\\d+).*", "\\1", conf_min_db))
      db_part <- sub(".*m\\d+([A-Za-z0-9_]+).*", "\\1", conf_min_db)
      result$database_used <- sub("_report$", "", db_part)  # Remove trailing _report
    }
  }
  # Pattern 4: Specific pattern like "W50504632_S109_bracken_S"
  else if (grepl("_S\\d+_bracken_S$", sample_name)) {
    result$trimmed_name <- sub("_S\\d+_bracken_S$", "", sample_name)
    result$confidence <- NA
    result$minimum_hits <- NA
    result$database_used <- "unknown"
  }
  # Pattern 5: Specific pattern like "Neg_Control_W_S90_bracken_S"
  else if (grepl("_S\\d+_bracken_S$", sample_name)) {
    result$trimmed_name <- sub("_S\\d+_bracken_S$", "", sample_name)
    result$confidence <- NA
    result$minimum_hits <- NA
    result$database_used <- "unknown"
  }
  # Pattern 6: Specific pattern like "W50504632_bracken_S"
  else if (grepl("_bracken_S$", sample_name)) {
    result$trimmed_name <- sub("_bracken_S$", "", sample_name)
    result$confidence <- NA
    result$minimum_hits <- NA
    result$database_used <- "unknown"
  }
  # Pattern 7: Specific pattern like "Neg_Control_W_bracken_S"
  else if (grepl("_bracken_S$", sample_name)) {
    result$trimmed_name <- sub("_bracken_S$", "", sample_name)
    result$confidence <- NA
    result$minimum_hits <- NA
    result$database_used <- "unknown"
  }
  # Pattern 8: Simple pattern with just confidence (fallback)
  else if (grepl("c\\d+", sample_name)) {
    result$trimmed_name <- sub("_c\\d+.*", "", sample_name)
    conf_match <- regmatches(sample_name, regexpr("c\\d+", sample_name))
    if (length(conf_match) > 0) {
      result$confidence <- as.integer(sub("c(\\d+)", "\\1", conf_match))
    }
  }
  # Pattern 7: No specific pattern - use full name
  else {
    result$trimmed_name <- sample_name
  }
  return(result)
} # Denis Rivard: LGTM 30/7/2025

# Function to filter kreport data based on clade reads and minimizer criteria
filter_kreport_data <- function(kreport_df, sample_name, min_clade_reads = MIN_CLADE_READS, 
                                minimizer_ratio = MINIMIZER_RATIO, minimizer_threshold = MINIMIZER_THRESHOLD,
                                exclude_taxid = EXCLUDE_TAXID) {
  initial_count <- nrow(kreport_df)
  
  # Apply taxID exclusion filter first if specified
  if (!is.null(exclude_taxid)) {
    excluded_count <- kreport_df %>% filter(taxID == exclude_taxid) %>% nrow()
    kreport_df <- kreport_df %>%
      filter(taxID != exclude_taxid)
    if (excluded_count > 0) {
      cat("Sample", sample_name, ": excluded", excluded_count, "entries with taxID", exclude_taxid, "\n")
    }
  }
  
  # Apply minimum clade reads filter
  if (min_clade_reads > 0) {
    pre_reads_count <- nrow(kreport_df)
    kreport_df <- kreport_df %>%
      filter(cladeReads >= min_clade_reads)
    cat("Sample", sample_name, ": filtered from", pre_reads_count, "to", nrow(kreport_df), 
        "species with >=", min_clade_reads, "clade reads\n")
  }
  
  # Apply minimizer filtering if either parameter is specified
  if (!is.null(minimizer_ratio) || !is.null(minimizer_threshold)) {
    pre_minimizer_count <- nrow(kreport_df)
    
    # Apply ratio filter if specified
    if (!is.null(minimizer_ratio)) {
      kreport_df <- kreport_df %>%
        filter(is.na(distinct_minimizers) | distinct_minimizers == 0 | 
               (cladeReads / pmax(distinct_minimizers, 1)) >= (1/minimizer_ratio))
      cat("Sample", sample_name, ": applied minimizer ratio filter (>=", minimizer_ratio, 
          "): from", pre_minimizer_count, "to", nrow(kreport_df), "species\n")
    }
    
    # Apply threshold filter if specified
    if (!is.null(minimizer_threshold)) {
      threshold_pre_count <- nrow(kreport_df)
      kreport_df <- kreport_df %>%
        filter(is.na(distinct_minimizers) | distinct_minimizers >= minimizer_threshold)
      cat("Sample", sample_name, ": applied minimizer threshold filter (>=", minimizer_threshold, 
          "): from", threshold_pre_count, "to", nrow(kreport_df), "species\n")
    }
  }
  
  return(kreport_df)
}

# Simplified function to process kreport files
process_kreport_folder <- function(kreports_folder, include_subspecies = INCLUDE_SUBSPECIES) {
  if (!dir.exists(kreports_folder)) {
    stop("Directory does not exist: ", kreports_folder)
  }
  
  # Initialize fresh report_data for this folder
  if (include_subspecies) {
    report_data_df <- tibble(
      sample = character(), 
      unclassified_counts = integer(), 
      confidence = integer(), 
      minimum_hits = integer(), 
      database_used = character(),
      species_count = integer(),
      subspecies_count = integer(),
      total_species_count = integer(),
      filtered_species_count = integer(),
      filtered_subspecies_count = integer(),
      filtered_total_species_count = integer()
    )
  } else {
    report_data_df <- tibble(
      sample = character(), 
      unclassified_counts = integer(), 
      confidence = integer(), 
      minimum_hits = integer(), 
      database_used = character(),
      species_count = integer(),
      filtered_species_count = integer()
    )
  }
  
  kreports_files <- list.files(kreports_folder, pattern = "\\.(kreport|k2report)$", full.names = TRUE)
  
  if (length(kreports_files) == 0) {
    warning("No .kreport or .k2report files found in: ", kreports_folder)
    return(list(report_data = report_data_df, kreports_combined = tibble()))
  }
  
  cat("Processing", length(kreports_files), "kreport files from", kreports_folder, "\n")
  
  start_time <- Sys.time()
  # Initialize list to store raw kreports
  kreports_raw_list <- list()
  
  # Read all kreport files (sequential or parallel)
  if (USE_PARALLEL && exists("cl", envir = .GlobalEnv)) {
    cat("Reading kreport files in parallel...\n")
    kreports_raw_list <- foreach(i = seq_along(kreports_files), 
                                .export = c("read_report2", "build_kraken_tree", "collapse.taxRanks", "delete_taxRanks_below", "parse_sample_info")) %dopar% {
      file <- kreports_files[i]
      
      # Read kreport file
      kreport_df <- read_report2(myfile = file, has_header = FALSE, add_taxRank_columns = FALSE, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1", "S2", "S3"))
      
      if (is.null(kreport_df) || nrow(kreport_df) == 0) {
        return(NULL)
      }
      
      # Store raw kreport with file info
      list(
        file = file,
        data = kreport_df,
        sample_info = parse_sample_info(file)
      )
    }
  } else {
    cat("Reading kreport files sequentially...\n")
    for (i in seq_along(kreports_files)) {
      file <- kreports_files[i]
      cat("Processing kreport:", file, "\n")
      
      # Read kreport file
      kreport_df <- read_report2(myfile = file, has_header = FALSE, add_taxRank_columns = FALSE, keep_taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S", "S1", "S2", "S3"))

      if (is.null(kreport_df) || nrow(kreport_df) == 0) {
        warning("Failed to read or empty file: ", file)
        next
      }
      
      # Store raw kreport with file info
      kreports_raw_list[[i]] <- list(
        file = file,
        data = kreport_df,
        sample_info = parse_sample_info(file)
      )
    }
  }

  # Remove NULL entries
  kreports_raw_list <- kreports_raw_list[!sapply(kreports_raw_list, is.null)]
  
  if (length(kreports_raw_list) == 0) {
    warning("No valid kreport files processed")
    return(list(report_data = report_data_df, kreports_combined = tibble()))
  }
  end_time <- Sys.time()
  cat("Finished reading kreport files in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")


  
  # Now process each kreport for report data and filtering
  kreports_list <- list()
  
  for (i in seq_along(kreports_raw_list)) {
    kreport_info <- kreports_raw_list[[i]]
    file <- kreport_info$file
    kreport_df <- kreport_info$data
    sample_info <- kreport_info$sample_info
    
    # Get unclassified reads count
    unclassified_reads <- kreport_df$cladeReads[kreport_df$taxID == 0]
    if (length(unclassified_reads) == 0) unclassified_reads <- 0
    
    # Add to report data
    if (include_subspecies) {
      report_data_df <- bind_rows(
        report_data_df,
        tibble(
          sample = sample_info$trimmed_name,
          unclassified_counts = unclassified_reads,
          confidence = sample_info$confidence,
          minimum_hits = sample_info$minimum_hits,
          database_used = sample_info$database_used,
          species_count = NA_integer_,
          subspecies_count = NA_integer_,
          total_species_count = NA_integer_,
          filtered_species_count = NA_integer_,
          filtered_subspecies_count = NA_integer_,
          filtered_total_species_count = NA_integer_
        )
      )
    } else {
      report_data_df <- bind_rows(
        report_data_df,
        tibble(
          sample = sample_info$trimmed_name,
          unclassified_counts = unclassified_reads,
          confidence = sample_info$confidence,
          minimum_hits = sample_info$minimum_hits,
          database_used = sample_info$database_used,
          species_count = NA_integer_,
          filtered_species_count = NA_integer_
        )
      )
    }
    
    # Process kreport data - always read all levels, then filter based on include_subspecies
    # First, get all species and subspecies data
    all_species_data <- kreport_df %>%
      filter(taxRank %in% c("S", "S1", "S2", "S3")) %>%
      mutate(
        name = str_remove(name, "^s\\d*_"),  # remove leading s_, s1_, s2_, s3_, etc.
        name = str_remove_all(name, "'") # remove apostrophes
      )
    
    # Filter based on include_subspecies parameter
    if (include_subspecies) {
      # Include both species and subspecies levels
      species_level_data <- all_species_data
    } else {
      # Include only species level
      species_level_data <- all_species_data %>%
        filter(taxRank == "S")
    }
    
    # Calculate initial species counts (before filtering)
    if (include_subspecies) {
      species_count <- all_species_data %>% filter(taxRank == "S") %>% nrow()
      subspecies_count <- all_species_data %>% filter(taxRank %in% c("S1", "S2", "S3")) %>% nrow()
      total_species_count <- nrow(all_species_data)
    } else {
      species_count <- all_species_data %>% filter(taxRank == "S") %>% nrow()
    }
    
    # Apply filtering to species-level data
    filtered_data <- filter_kreport_data(species_level_data, sample_info$trimmed_name)
    
    # Calculate filtered species counts
    if (include_subspecies) {
      filtered_species_count <- filtered_data %>% filter(taxRank == "S") %>% nrow()
      filtered_subspecies_count <- filtered_data %>% filter(taxRank %in% c("S1", "S2", "S3")) %>% nrow()
      filtered_total_species_count <- nrow(filtered_data)
      
      # Update the last row of report_data with species counts
      report_data_df[nrow(report_data_df), "species_count"] <- species_count
      report_data_df[nrow(report_data_df), "subspecies_count"] <- subspecies_count
      report_data_df[nrow(report_data_df), "total_species_count"] <- total_species_count
      report_data_df[nrow(report_data_df), "filtered_species_count"] <- filtered_species_count
      report_data_df[nrow(report_data_df), "filtered_subspecies_count"] <- filtered_subspecies_count
      report_data_df[nrow(report_data_df), "filtered_total_species_count"] <- filtered_total_species_count
    } else {
      filtered_species_count <- nrow(filtered_data)
      
      # Update the last row of report_data with species counts
      report_data_df[nrow(report_data_df), "species_count"] <- species_count
      report_data_df[nrow(report_data_df), "filtered_species_count"] <- filtered_species_count
    }
    
    # Prepare final processed data (common for both subspecies and species-only)
    processed_kreport <- filtered_data %>%
      select(cladeReads, minimizers, distinct_minimizers, taxID, name, taxRank, taxLineage) %>%
      rename_with(~ paste0(sample_info$trimmed_name, "_", .), 
                  c(cladeReads, minimizers, distinct_minimizers))
    
    # Check for duplicate taxIDs
    duplicate_taxids <- processed_kreport %>% 
      count(taxID) %>% 
      filter(n > 1)
    
    if (nrow(duplicate_taxids) > 0) {
      warning("Sample ", sample_info$trimmed_name, " has ", nrow(duplicate_taxids), 
             " duplicate taxIDs: ", paste(duplicate_taxids$taxID, collapse = ", "))
    }

    # Add processed kreport to list
    kreports_list[[i]] <- processed_kreport
  }
  
  # Remove NULL entries
  kreports_list <- kreports_list[!sapply(kreports_list, is.null)]
  
  if (length(kreports_list) == 0) {
    warning("No valid kreport files processed")
    # Return properly structured empty tibble
    empty_kreports <- tibble(
      taxID = integer(0),
      name = character(0),
      taxRank = character(0),
      taxLineage = character(0)
    )
    return(list(report_data = report_data_df, kreports_combined = empty_kreports))
  }
  
  # Combine all kreports
  kreports_combined <- reduce(kreports_list, full_join, by = c("taxID", "name", "taxRank", "taxLineage")) %>%
    select(taxID, name, taxRank, taxLineage, everything())
  
  return(list(report_data = report_data_df, kreports_combined = kreports_combined))
}

# Simplified function to process bracken files
process_bracken_folder <- function(bracken_folder) {
  if (!dir.exists(bracken_folder)) {
    warning("Directory does not exist: ", bracken_folder)
    return(tibble())
  }
  
  bracken_files <- list.files(bracken_folder, pattern = ".*bracken.*\\.txt$", full.names = TRUE)
  
  if (length(bracken_files) == 0) {
    warning("No bracken .txt files found in: ", bracken_folder)
    return(tibble())
  }
  
  cat("Processing", length(bracken_files), "bracken files from", bracken_folder, "\n")
  
  bracken_list <- map(bracken_files, function(file) {
    tryCatch({
      # Ensure libraries are available for parallel processing
      
      sample_info <- parse_sample_info(file)
      
      bracken_df <- read.table(file, sep = "\t", header = T, comment.char = "", quote = "") %>%
        select(name, taxonomy_id, new_est_reads) %>%
        rename(taxID = taxonomy_id) %>%
        rename_with(~ paste0(sample_info$trimmed_name, "_bracken_reads"), 
                    .cols = "new_est_reads")
      
      return(bracken_df)
    }, error = function(e) {
      cat("ERROR processing bracken file:", basename(file), "-", e$message, "\n")
      return(NULL)
    })
  })
  
  # Combine all bracken reports
  # Remove NULL entries first
  bracken_list <- bracken_list[!sapply(bracken_list, is.null)]
  
  if (length(bracken_list) > 0) {
    bracken_combined <- reduce(bracken_list, full_join, by = c("name", "taxID"))
  } else {
    bracken_combined <- tibble()
  }
  
  return(bracken_combined)
}

# Function to merge kreport and bracken data
merge_data <- function(kreports_data, bracken_data = NULL, add_params = ADD_PARAMS, report_data = NULL) {
  # Merge with bracken data if available
  if (!is.null(bracken_data) && nrow(bracken_data) > 0) {
    merged_data <- kreports_data %>%
      full_join(bracken_data, by = c("name", "taxID"))
  } else {
    merged_data <- kreports_data
  }
  
  # Reorder columns
  merged_data <- merged_data %>%
    select(taxID, name, taxRank, taxLineage, everything())
  
  # Create long format
  merged_long <- merged_data %>%
    pivot_longer(
      cols = -c(taxID, name, taxRank, taxLineage),
      names_to = c("sample", "measurement"),
      names_pattern = "^(.*?)_(cladeReads|minimizers|distinct_minimizers|bracken_reads)$",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = measurement,
      values_from = value
    )
  
  # Add parameter columns if requested
  if (add_params) {
    # Join with report_data to get parameter information
    if (!is.null(report_data) && nrow(report_data) > 0) {
      # Create a lookup table from report_data for parameter information
      param_lookup <- report_data %>%
        select(sample, confidence, minimum_hits, database_used) %>%
        distinct() %>%
        rename(
          confidence_levels = confidence,
          minimum_hit_groups = minimum_hits
        ) %>%
        mutate(
          # Determine human_reads based on dataset or database information
          # This is a placeholder - adjust logic based on your specific criteria
          human_reads = grepl("human|nonhuman", database_used, ignore.case = TRUE)
        )
      
      # Join the parameter information with merged_long
      merged_long <- merged_long %>%
        left_join(param_lookup, by = "sample")
      
      cat("Added parameter columns: confidence_levels, minimum_hit_groups, human_reads\n")
    } else {
      warning("report_data is NULL or empty - cannot add parameter columns")
    }
  }
  
  return(list(merged = merged_data, merged_long = merged_long))
}

# Function to add summary statistics
add_summary_stats <- function(merged_data) {
  cladeReads_cols <- grep("_cladeReads$", colnames(merged_data), value = TRUE)
  bracken_reads_cols <- grep("_bracken_reads$", colnames(merged_data), value = TRUE)
  
  result <- merged_data
  
  if (length(cladeReads_cols) > 0) {
    result <- result %>%
      mutate(
        cladeReads_mean = rowMeans(select(., all_of(cladeReads_cols)), na.rm = TRUE),
        cladeReads_max = apply(select(., all_of(cladeReads_cols)), 1, max, na.rm = TRUE)
      )
  }
  
  if (length(bracken_reads_cols) > 0) {
    result <- result %>%
      mutate(
        bracken_reads_mean = rowMeans(select(., all_of(bracken_reads_cols)), na.rm = TRUE),
        bracken_reads_max = apply(select(., all_of(bracken_reads_cols)), 1, max, na.rm = TRUE)
      )
  }
  
  return(result)
}

# Function to get complete species list with topN frequency included
get_species_list_with_frequency <- function(merged_data, count_column_pattern = "_cladeReads$", 
                                           top_n = 10, exclude_taxid = NULL) {
  count_cols <- grep(count_column_pattern, colnames(merged_data), value = TRUE)
  
  if (length(count_cols) == 0) {
    warning("No columns matching pattern found: ", count_column_pattern)
    # Return properly structured empty tibble with required columns
    return(tibble(
      name = character(0),
      taxID = integer(0),
      taxRank = character(0),
      taxLineage = character(0),
      Freq = integer(0)
    ))
  }
  
  # Apply taxID exclusion only if exclude_taxid is specified
  if (!is.null(exclude_taxid)) {
    filtered_data <- merged_data %>%
      filter(taxID != exclude_taxid)
    cat("Excluding taxID", exclude_taxid, "from species list\n")
  } else {
    filtered_data <- merged_data
    cat("No taxID exclusion applied in species list\n")
  }
  
  # Get top species from each sample
  top_species_list <- map_dfr(count_cols, function(sample_col) {
    temp_df <- filtered_data %>%
      select(name, taxID, taxRank, taxLineage, count = all_of(sample_col)) %>%
      arrange(desc(count)) %>%
      slice_head(n = top_n)
    
    return(temp_df)
  })
  
  # Create frequency table for top species
  top_species_freq <- top_species_list %>%
    count(taxID, name = "topN_freq") %>%
    filter(topN_freq >= 1)
  
  # Create complete species list with same format but include topN frequency
  complete_list_df <- filtered_data %>%
    select(name, taxID, taxRank, taxLineage) %>%
    distinct() %>%
    left_join(top_species_freq, by = "taxID") %>%
    mutate(Freq = ifelse(is.na(topN_freq), 0, topN_freq)) %>%
    select(-topN_freq) %>%
    arrange(desc(Freq), name)
  
  # Add cladeReads_mean and cladeReads_max if they exist in the data
  if ("cladeReads_mean" %in% colnames(filtered_data) && "cladeReads_max" %in% colnames(filtered_data)) {
    # Get the mean and max values for each species
    summary_stats <- filtered_data %>%
      select(taxID, cladeReads_mean, cladeReads_max) %>%
      distinct()
    
    # Join with the species list
    complete_list_df <- complete_list_df %>%
      left_join(summary_stats, by = "taxID")
  }
  
  return(complete_list_df)
}

# Function to parse sample metadata files
parse_sample_metadata <- function(input_dir) {
  # Look for sample metadata CSV files
  metadata_files <- list.files(input_dir, pattern = "sample_metadata.*\\.csv$", full.names = TRUE, recursive = FALSE)
  
  if (length(metadata_files) == 0) {
    cat("No sample metadata files found in input directory\n")
    return(NULL)
  }
  
  cat("Found", length(metadata_files), "sample metadata files\n")
  
  all_metadata <- list()
  
  for (metadata_file in metadata_files) {
    cat("Processing metadata file:", basename(metadata_file), "\n")
    
    tryCatch({
      metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
      
      # Check for required columns (case-insensitive)
      required_cols <- c("sample", "condition")
      col_names_lower <- tolower(colnames(metadata))
      
      # Find matching columns
      sample_col <- which(col_names_lower %in% c("sample", "samples"))
      condition_col <- which(col_names_lower %in% c("condition", "conditions"))
      
      if (length(sample_col) == 0 || length(condition_col) == 0) {
        cat("Warning: Metadata file missing required columns. Expected 'Sample' and 'Condition'\n")
        cat("Available columns:", paste(colnames(metadata), collapse = ", "), "\n")
        next
      }
      
      # Standardize column names
      colnames(metadata)[sample_col[1]] <- "sample"
      colnames(metadata)[condition_col[1]] <- "condition"
      
      # Remove rows with missing sample names
      metadata <- metadata[!is.na(metadata$sample) & metadata$sample != "", ]
      
      if (nrow(metadata) > 0) {
        all_metadata[[basename(metadata_file)]] <- metadata
      }
      
    }, error = function(e) {
      cat("Error reading", basename(metadata_file), ":", e$message, "\n")
    })
  }
  
  # Combine all metadata
  if (length(all_metadata) > 0) {
    combined_metadata <- do.call(rbind, all_metadata)
    
    # Remove duplicates, keeping latest values
    final_metadata <- combined_metadata %>%
      group_by(sample) %>%
      slice_tail(n = 1) %>%
      ungroup()
    
    return(final_metadata)
  } else {
    return(NULL)
  }
}

# Function to parse runtime files
parse_runtime_files <- function(runtime_folder) {
  if (!dir.exists(runtime_folder)) {
    cat("Runtime folder not found:", runtime_folder, "\n")
    return(NULL)
  }
  
  runtime_files <- list.files(runtime_folder, pattern = "^runtime", full.names = TRUE)
  
  if (length(runtime_files) == 0) {
    cat("No runtime files found in:", runtime_folder, "\n")
    return(NULL)
  }
  
  cat("Found", length(runtime_files), "runtime files\n")
  
  all_runtime_data <- list()
  
  for (runtime_file in runtime_files) {
    cat("Processing runtime file:", basename(runtime_file), "\n")
    
    # Read the file
    lines <- readLines(runtime_file, warn = FALSE)
    
    # Initialize data structures for this file
    runtime_data <- data.frame(
      sample = character(),
      cutadapt_star_samtools_rt = numeric(),
      kraken2_unaligned_rt = numeric(),
      kraken2_nonhuman_rt = numeric(),
      is_first_k2_unaligned = logical(),
      is_first_k2_nonhuman = logical(),
      stringsAsFactors = FALSE
    )
    
    # Track first samples for database loading
    first_k2_unaligned <- TRUE
    first_k2_nonhuman <- TRUE
    
    # Parse each line
    for (line in lines) {
      # Parse Cutadapt/STAR/Samtools runtime
      if (grepl("^Cutadapt/STAR/Samtools \\(.*\\); RT:", line)) {
        sample_match <- regmatches(line, regexpr("\\(([^)]+)\\)", line))
        sample_name <- gsub("[()]", "", sample_match)
        rt_match <- regmatches(line, regexpr("RT: (\\d+) seconds", line))
        rt_seconds <- as.numeric(gsub("RT: | seconds", "", rt_match))
        
        # Find existing row or create new one
        existing_row <- which(runtime_data$sample == sample_name)
        if (length(existing_row) > 0) {
          # Update existing row (take latest runtime)
          runtime_data$cutadapt_star_samtools_rt[existing_row] <- rt_seconds
        } else {
          # Add new row
          new_row <- data.frame(
            sample = sample_name,
            cutadapt_star_samtools_rt = rt_seconds,
            kraken2_unaligned_rt = NA,
            kraken2_nonhuman_rt = NA,
            is_first_k2_unaligned = FALSE,
            is_first_k2_nonhuman = FALSE,
            stringsAsFactors = FALSE
          )
          runtime_data <- rbind(runtime_data, new_row)
        }
      }
      
      # Parse Kraken2 unaligned runtime
      if (grepl("^Kraken2 unaligned \\(.*\\); RT:", line)) {
        sample_match <- regmatches(line, regexpr("\\(([^)]+)\\)", line))
        sample_name <- gsub("[()]", "", sample_match)
        # Remove suffix pattern like "_S109" to match sample names
        sample_name <- sub("_S\\d+$", "", sample_name)
        rt_match <- regmatches(line, regexpr("RT: (\\d+) seconds", line))
        rt_seconds <- as.numeric(gsub("RT: | seconds", "", rt_match))
        
        # Find existing row or create new one
        existing_row <- which(runtime_data$sample == sample_name)
        if (length(existing_row) > 0) {
          # Update existing row (take latest runtime)
          runtime_data$kraken2_unaligned_rt[existing_row] <- rt_seconds
          runtime_data$is_first_k2_unaligned[existing_row] <- first_k2_unaligned
        } else {
          # Add new row
          new_row <- data.frame(
            sample = sample_name,
            cutadapt_star_samtools_rt = NA,
            kraken2_unaligned_rt = rt_seconds,
            kraken2_nonhuman_rt = NA,
            is_first_k2_unaligned = first_k2_unaligned,
            is_first_k2_nonhuman = FALSE,
            stringsAsFactors = FALSE
          )
          runtime_data <- rbind(runtime_data, new_row)
        }
        
        # Mark that we've seen the first K2 unaligned sample
        if (first_k2_unaligned) {
          cat("  First K2 unaligned sample in", basename(runtime_file), ":", sample_name, "\n")
          first_k2_unaligned <- FALSE
        }
      }
      
      # Parse Kraken2 non-human runtime
      if (grepl("^Kraken2 non-human \\(.*\\); RT:", line)) {
        sample_match <- regmatches(line, regexpr("\\(([^)]+)\\)", line))
        sample_name <- gsub("[()]", "", sample_match)
        rt_match <- regmatches(line, regexpr("RT: (\\d+) seconds", line))
        rt_seconds <- as.numeric(gsub("RT: | seconds", "", rt_match))
        
        # Find existing row or create new one
        existing_row <- which(runtime_data$sample == sample_name)
        if (length(existing_row) > 0) {
          # Update existing row (take latest runtime)
          runtime_data$kraken2_nonhuman_rt[existing_row] <- rt_seconds
          runtime_data$is_first_k2_nonhuman[existing_row] <- first_k2_nonhuman
        } else {
          # Add new row
          new_row <- data.frame(
            sample = sample_name,
            cutadapt_star_samtools_rt = NA,
            kraken2_unaligned_rt = NA,
            kraken2_nonhuman_rt = rt_seconds,
            is_first_k2_unaligned = FALSE,
            is_first_k2_nonhuman = first_k2_nonhuman,
            stringsAsFactors = FALSE
          )
          runtime_data <- rbind(runtime_data, new_row)
        }
        
        # Mark that we've seen the first K2 non-human sample
        if (first_k2_nonhuman) {
          cat("  First K2 non-human sample in", basename(runtime_file), ":", sample_name, "\n")
          first_k2_nonhuman <- FALSE
        }
      }
    }
    
    all_runtime_data[[basename(runtime_file)]] <- runtime_data
  }
  
  # Combine all runtime data from different files
  if (length(all_runtime_data) > 0) {
    combined_runtime <- do.call(rbind, all_runtime_data)
    
    # Remove any rows with empty sample names
    combined_runtime <- combined_runtime[!is.na(combined_runtime$sample) & combined_runtime$sample != "", ]
    
    # If there are duplicates (same sample from multiple files), keep the latest values
    # Group by sample and take the non-NA values, preferring later files
    # For boolean flags, use OR logic (TRUE if any file marks it as first)
    suppressWarnings({
      final_runtime <- combined_runtime %>%
        group_by(sample) %>%
        summarise(
          cutadapt_star_samtools_rt = last(cutadapt_star_samtools_rt[!is.na(cutadapt_star_samtools_rt)]),
          kraken2_unaligned_rt = last(kraken2_unaligned_rt[!is.na(kraken2_unaligned_rt)]),
          kraken2_nonhuman_rt = last(kraken2_nonhuman_rt[!is.na(kraken2_nonhuman_rt)]),
          is_first_k2_unaligned = any(is_first_k2_unaligned, na.rm = TRUE),
          is_first_k2_nonhuman = any(is_first_k2_nonhuman, na.rm = TRUE),
          .groups = 'drop'
        )
    })
    
    return(final_runtime)
  } else {
    return(NULL)
  }
}

# Function to parse read statistics files
parse_rrstats_files <- function(rrstats_folder) {
  if (!dir.exists(rrstats_folder)) {
    cat("Read stats folder not found:", rrstats_folder, "\n")
    return(NULL)
  }
  
  rrstats_files <- list.files(rrstats_folder, pattern = "run_read_stats.*\\.txt$", full.names = TRUE)
  
  if (length(rrstats_files) == 0) {
    cat("No read stats files found in:", rrstats_folder, "\n")
    return(NULL)
  }
  
  cat("Found", length(rrstats_files), "read stats files\n")
  
  all_rrstats_data <- list()
  
  for (rrstats_file in rrstats_files) {
    cat("Processing read stats file:", basename(rrstats_file), "\n")
    
    # Read the file as tab-separated
    tryCatch({
      rrstats_data <- read.table(rrstats_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                                skip = 0, fill = TRUE, comment.char = "")
      
      # Remove any completely empty rows
      rrstats_data <- rrstats_data[!apply(rrstats_data == "" | is.na(rrstats_data), 1, all), ]
      
      # Skip if no data after cleaning
      if (nrow(rrstats_data) == 0) {
        cat("No valid data found in", basename(rrstats_file), "\n")
        next
      }
      
      # Set column names based on the header structure
      colnames(rrstats_data) <- c("sample", "num_input_reads", "avg_input_read_length", 
                                 "uniquely_mapped_reads", "avg_mapped_length", 
                                 "reads_assigned_to_genes", "percent_mapped_reads", 
                                 "uniquely_mapped_percent", "less_than_10M_unique")
      
      # Convert numeric columns with better error handling
      numeric_cols <- c("num_input_reads", "avg_input_read_length", "uniquely_mapped_reads", 
                       "avg_mapped_length", "reads_assigned_to_genes", "percent_mapped_reads", 
                       "uniquely_mapped_percent")
      
      for (col in numeric_cols) {
        if (col %in% colnames(rrstats_data)) {
          # Suppress coercion warnings and handle them properly
          suppressWarnings({
            rrstats_data[[col]] <- as.numeric(as.character(rrstats_data[[col]]))
          })
        }
      }
      
      # Convert boolean column with better handling
      if ("less_than_10M_unique" %in% colnames(rrstats_data)) {
        rrstats_data$less_than_10M_unique <- as.logical(rrstats_data$less_than_10M_unique)
      }
      
      # Remove rows where sample name is NA or empty
      rrstats_data <- rrstats_data[!is.na(rrstats_data$sample) & rrstats_data$sample != "", ]
      
      if (nrow(rrstats_data) > 0) {
        all_rrstats_data[[basename(rrstats_file)]] <- rrstats_data
      }
      
    }, error = function(e) {
      cat("Error reading", basename(rrstats_file), ":", e$message, "\n")
    })
  }
  
  # Combine all read stats data from different files
  if (length(all_rrstats_data) > 0) {
    combined_rrstats <- do.call(rbind, all_rrstats_data)
    
    # Remove any remaining NA sample names
    combined_rrstats <- combined_rrstats[!is.na(combined_rrstats$sample) & combined_rrstats$sample != "", ]
    
    # If there are duplicates (same sample from multiple files), keep the latest values
    suppressWarnings({
      final_rrstats <- combined_rrstats %>%
        group_by(sample) %>%
        summarise(
          num_input_reads = last(num_input_reads[!is.na(num_input_reads)]),
          avg_input_read_length = last(avg_input_read_length[!is.na(avg_input_read_length)]),
          uniquely_mapped_reads = last(uniquely_mapped_reads[!is.na(uniquely_mapped_reads)]),
          avg_mapped_length = last(avg_mapped_length[!is.na(avg_mapped_length)]),
          reads_assigned_to_genes = last(reads_assigned_to_genes[!is.na(reads_assigned_to_genes)]),
          percent_mapped_reads = last(percent_mapped_reads[!is.na(percent_mapped_reads)]),
          uniquely_mapped_percent = last(uniquely_mapped_percent[!is.na(uniquely_mapped_percent)]),
          less_than_10M_unique = last(less_than_10M_unique[!is.na(less_than_10M_unique)]),
          .groups = 'drop'
        )
    })
    
    return(final_rrstats)
  } else {
    return(NULL)
  }
}

# Main processing function
process_dataset <- function(dataset_name, kreports_folder, bracken_folder = NULL, 
                           process_bracken = PROCESS_BRACKEN, add_params = ADD_PARAMS, 
                           include_subspecies = INCLUDE_SUBSPECIES) {
  cat("\n=== Processing", dataset_name, "dataset ===\n")
  
  # Process kreports
  kreport_results <- process_kreport_folder(kreports_folder, include_subspecies = include_subspecies)
  
  # Process bracken if folder provided and processing is enabled
  bracken_data <- NULL
  if (process_bracken && !is.null(bracken_folder) && dir.exists(bracken_folder)) {
    bracken_data <- process_bracken_folder(bracken_folder)
    cat("Bracken processing enabled for", dataset_name, "\n")
  } else {
    cat("Bracken processing skipped for", dataset_name, "\n")
  }
  
  # Merge data
  merge_results <- merge_data(kreport_results$kreports_combined, bracken_data, 
                             add_params = add_params, report_data = kreport_results$report_data)
  
  # Add summary statistics
  merged_with_stats <- add_summary_stats(merge_results$merged)
  
  # Get species list with topN frequency
  species_list <- get_species_list_with_frequency(merged_with_stats, top_n = TOP_SPECIES_N, 
                                                 exclude_taxid = EXCLUDE_TAXID)
  cat("Generated species list with", nrow(species_list), "species\n")
  
  result <- list(
    report_data = kreport_results$report_data,
    merged = merged_with_stats,
    merged_long = merge_results$merged_long,
    species_list = species_list,
    bracken_data = bracken_data
  )
  return(result)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cat("Starting kraken/bracken output processing...\n")

# Check if unaligned kreports directory is specified
if (is.null(UNALIGNED_KREPORTS_DIR)) {
  cat("ERROR: No unaligned kreports directory specified. Use --unaligned-kreports to specify the path.\n")
  stop("Required unaligned kreports directory not specified.")
}

# Process unaligned dataset
unaligned_bracken_folder <- if (PROCESS_BRACKEN) UNALIGNED_BRACKEN_DIR else NULL
unaligned_results <- process_dataset(
  dataset_name = "Unaligned",
  kreports_folder = UNALIGNED_KREPORTS_DIR,
  bracken_folder = unaligned_bracken_folder,
  process_bracken = PROCESS_BRACKEN,
  add_params = ADD_PARAMS,
  include_subspecies = INCLUDE_SUBSPECIES
)

# Skip processing nonhuman dataset if NONHUMAN_KREPORTS_DIR is NULL
if (!is.null(NONHUMAN_KREPORTS_DIR)) {
  # Process nonhuman dataset  
  nonhuman_bracken_folder <- if (PROCESS_BRACKEN) NONHUMAN_BRACKEN_DIR else NULL
  nonhuman_results <- process_dataset(
    dataset_name = "Nonhuman",
    kreports_folder = NONHUMAN_KREPORTS_DIR, 
    bracken_folder = nonhuman_bracken_folder,
    process_bracken = PROCESS_BRACKEN,
    add_params = ADD_PARAMS,
    include_subspecies = INCLUDE_SUBSPECIES
  )
} else {
  cat("\nSkipping nonhuman dataset processing as NONHUMAN_KREPORTS_DIR is NULL\n")
}

# Only save if SAVE_OUTPUT is TRUE
if (SAVE_OUTPUT) {
  cat("\n=== Saving outputs ===\n")
  output_dir <- OUTPUT_DIR

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Save unaligned results
  output_csv_file(unaligned_results$merged, "unaligned_merged", output_dir, "op")
  output_csv_file(unaligned_results$merged_long, "unaligned_merged_long", output_dir, "op")
  output_csv_file(unaligned_results$species_list, "species_list_unaligned", output_dir, "op")
  
  # Save unaligned bracken data if it exists
  if (!is.null(unaligned_results$bracken_data) && nrow(unaligned_results$bracken_data) > 0) {
    output_csv_file(unaligned_results$bracken_data, "unaligned_bracken", output_dir, "op")
  }

  # Save nonhuman results only if they exist
  if (!is.null(NONHUMAN_KREPORTS_DIR) && exists("nonhuman_results")) {
    output_csv_file(nonhuman_results$merged, "nonhuman_merged", output_dir, "op")
    output_csv_file(nonhuman_results$merged_long, "nonhuman_merged_long", output_dir, "op")
    output_csv_file(nonhuman_results$species_list, "species_list_nonhuman", output_dir, "op")
    
    # Save nonhuman bracken data if it exists
    if (!is.null(nonhuman_results$bracken_data) && nrow(nonhuman_results$bracken_data) > 0) {
      output_csv_file(nonhuman_results$bracken_data, "nonhuman_bracken", output_dir, "op")
    }
  }

  # Save combined report data
  if (!is.null(NONHUMAN_KREPORTS_DIR) && exists("nonhuman_results")) {
    combined_report_data <- bind_rows(
      unaligned_results$report_data %>% mutate(dataset = "unaligned"),
      nonhuman_results$report_data %>% mutate(dataset = "nonhuman")
    )
  } else {
    combined_report_data <- unaligned_results$report_data %>% mutate(dataset = "unaligned")
  }
  output_csv_file(combined_report_data, "sample_report_data", output_dir, "op")
  
  # Parse and save metadata files if requested
  cat("\n=== Parsing additional metadata ===\n")
  
  # Parse runtime files if directory is provided
  if (!is.null(RUNTIME_DIR)) {
    cat("Parsing runtime files...\n")
    runtime_data <- parse_runtime_files(RUNTIME_DIR)
    if (!is.null(runtime_data)) {
      output_csv_file(runtime_data, "runtime_info", output_dir, "op")
      cat("Runtime data saved successfully\n")
    }
  }
  
  # Parse read statistics files if directory is provided
  if (!is.null(RRSTATS_DIR)) {
    cat("Parsing read statistics files...\n")
    rrstats_data <- parse_rrstats_files(RRSTATS_DIR)
    if (!is.null(rrstats_data)) {
      output_csv_file(rrstats_data, "read_statistics", output_dir, "op")
      cat("Read statistics data saved successfully\n")
    }
  }
  
  # Parse sample metadata files if directory is provided
  if (!is.null(METADATA_DIR)) {
    cat("Parsing sample metadata files...\n")
    metadata_data <- parse_sample_metadata(METADATA_DIR)
    if (!is.null(metadata_data)) {
      output_csv_file(metadata_data, "sample_metadata", output_dir, "op")
      cat("Sample metadata saved successfully\n")
    }
  }
  
  # If multiple metadata types were parsed, create a combined metadata file
  parsed_metadata <- list()
  if (!is.null(RUNTIME_DIR) && exists("runtime_data") && !is.null(runtime_data)) {
    parsed_metadata$runtime <- runtime_data
  }
  if (!is.null(RRSTATS_DIR) && exists("rrstats_data") && !is.null(rrstats_data)) {
    parsed_metadata$rrstats <- rrstats_data
  }
  if (!is.null(METADATA_DIR) && exists("metadata_data") && !is.null(metadata_data)) {
    parsed_metadata$metadata <- metadata_data
  }
  
  if (length(parsed_metadata) > 1) {
    cat("Creating combined metadata file...\n")
    # Start with the first dataset as base
    combined_metadata <- parsed_metadata[[1]]
    
    # Merge with other datasets
    for (i in 2:length(parsed_metadata)) {
      combined_metadata <- merge(combined_metadata, parsed_metadata[[i]], by = "sample", all = TRUE)
    }
    
    output_csv_file(combined_metadata, "combined_metadata", output_dir, "op")
    cat("Combined metadata saved successfully\n")
  }
  
  cat("Processing completed successfully!\n")
  cat("Output files saved to:", output_dir, "\n")
} else {
  cat("CSV output skipped (--output not specified)\n")
  cat("Processing completed successfully!\n")
}

# Clean up configuration variables
rm(args, PROCESS_BRACKEN, SAVE_OUTPUT, ADD_PARAMS, INCLUDE_SUBSPECIES, USE_PARALLEL, OUTPUT_DIR, MIN_CLADE_READS,
   RUNTIME_DIR, RRSTATS_DIR, METADATA_DIR, parsed_metadata)

