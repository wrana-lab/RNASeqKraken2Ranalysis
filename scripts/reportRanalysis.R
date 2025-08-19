#!/usr/bin/env Rscript

#### Header ####
# Report Analysis Orchestration Script
# Author: Denis Rivard
# Date: Aug 5th, 2025
# Description: Modular orchestration of downstream analysis of kreports
# Usage: Rscript reportRanalysis.R --input-dir <path> --proj-name <name> [options]

library(here)
source(here("Ranalysis", "scripts", "utils.R"))
library(parallel)
library(doParallel)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(purrr)
library(tidyr)
library(stringr)

# Function to display help
show_help <- function() {
  cat("Usage: Rscript reportRanalysis.R [OPTIONS]\n\n")
  cat("Required Options:\n")
  cat("  --input-dir <DIR>          Path to input directory containing kreport folders\n")
  cat("  --proj-name <NAME>         Project name for output directory naming\n\n")
  cat("Optional Options:\n")
  cat("  --output-base-dir <DIR>    Base output directory (default: outputs)\n")
  cat("  --top-n <N>               Number of top species to analyze (default: 25)\n")
  cat("  --min-reads <N>           Minimum clade reads threshold (default: 0)\n")
  cat("  --exclude-taxid <ID>      Exclude species with this taxonomy ID\n")
  cat("  --no-bracken              Skip bracken file processing\n")
  cat("  --no-corr                 Skip correlation analysis\n")
  cat("  --no-params               Skip adding parameter columns\n")
  cat("  --subspecies              Include subspecies in analysis\n")
  cat("  --minimizer-ratio <R>     Filter by minimum ratio of distinct_minimizers/cladeReads\n")
  cat("  --minimizer-threshold <N> Filter by minimum distinct_minimizers threshold\n")
  cat("  --cores <N>               Number of cores to use (default: detectCores())\n")
  cat("  --help, -h                Show this help message\n\n")
}

# Configuration structure
create_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- list(
    INPUT_DIR = NULL,
    PROJ_NAME = NULL,
    OUTPUT_BASE_DIR = "outputs",
    TOP_N = 25,
    MIN_READS = 0,
    EXCLUDE_TAXID = NULL,
    PROCESS_BRACKEN = TRUE,
    PERFORM_CORRELATION = TRUE,
    ADD_PARAMS = TRUE,
    INCLUDE_SUBSPECIES = FALSE,
    MINIMIZER_RATIO = NULL,
    MINIMIZER_THRESHOLD = NULL,
    CORES = detectCores()
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
      } else if (args[i] == "--top-n" && i < length(args)) {
        config$TOP_N <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--min-reads" && i < length(args)) {
        config$MIN_READS <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--exclude-taxid" && i < length(args)) {
        config$EXCLUDE_TAXID <- as.numeric(args[i + 1])
        i <- i + 2
      } else if (args[i] == "--no-bracken") {
        config$PROCESS_BRACKEN <- FALSE
        i <- i + 1
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
  
  # Create output paths - separate outputs and reports directories
  config$output_dir <- file.path(config$OUTPUT_BASE_DIR, paste0(config$PROJ_NAME, "_outputs"))
  config$reports_dir <- file.path(config$output_dir, "reports")
  
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

#### Logic ####
# Build arguments for output_processing.R
build_output_processing_args <- function(config) {
  # Check for different possible directory structures
  unaligned_kreports_dir <- NULL
  nonhuman_kreports_dir <- NULL
  unaligned_bracken_dir <- NULL
  nonhuman_bracken_dir <- NULL
  
  # Try the expected structure first
  if (dir.exists(file.path(config$INPUT_DIR, "unaligned_kreports"))) {
    unaligned_kreports_dir <- file.path(config$INPUT_DIR, "unaligned_kreports")
  } else if (dir.exists(file.path(config$INPUT_DIR, "kreports"))) {
    # Use general kreports directory for both unaligned and nonhuman
    unaligned_kreports_dir <- file.path(config$INPUT_DIR, "kreports")
  }
  
  if (dir.exists(file.path(config$INPUT_DIR, "nonhuman_kreports"))) {
    nonhuman_kreports_dir <- file.path(config$INPUT_DIR, "nonhuman_kreports")
  } else if (dir.exists(file.path(config$INPUT_DIR, "kreports"))) {
    # Use general kreports directory for both unaligned and nonhuman
    nonhuman_kreports_dir <- file.path(config$INPUT_DIR, "kreports")
  }
  
  # Check if we found valid kreport directories
  if (is.null(unaligned_kreports_dir) || is.null(nonhuman_kreports_dir)) {
    stop("ERROR: Could not find kreports directory in: ", config$INPUT_DIR)
  }
  
  args <- c(
    "--unaligned-kreports", unaligned_kreports_dir,
    "--nonhuman-kreports", nonhuman_kreports_dir,
    "--output-dir", config$output_dir,
    "--top-n", as.character(config$TOP_N),
    "--min-reads", as.character(config$MIN_READS),
    "--output",
    "--parallel"
  )
  
  # Add bracken directories if processing is enabled
  if (config$PROCESS_BRACKEN) {
    if (dir.exists(file.path(config$INPUT_DIR, "unaligned_bracken"))) {
      unaligned_bracken_dir <- file.path(config$INPUT_DIR, "unaligned_bracken")
    } else if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
      unaligned_bracken_dir <- file.path(config$INPUT_DIR, "bracken")
    }
    
    if (dir.exists(file.path(config$INPUT_DIR, "nonhuman_bracken"))) {
      nonhuman_bracken_dir <- file.path(config$INPUT_DIR, "nonhuman_bracken")
    } else if (dir.exists(file.path(config$INPUT_DIR, "bracken"))) {
      nonhuman_bracken_dir <- file.path(config$INPUT_DIR, "bracken")
    }
    
    if (!is.null(unaligned_bracken_dir) && !is.null(nonhuman_bracken_dir)) {
      args <- c(args, "--unaligned-bracken", unaligned_bracken_dir)
      args <- c(args, "--nonhuman-bracken", nonhuman_bracken_dir)
    } else {
      cat("Warning: Bracken directories not found, skipping bracken processing\n")
      args <- c(args, "--no-bracken")
    }
  } else {
    args <- c(args, "--no-bracken")
  }
  
  # Add metadata arguments if directories exist
  runtime_folder <- file.path(config$INPUT_DIR, "runtimes")
  rrstats_folder <- file.path(config$INPUT_DIR, "rrstats")
  metadata_files <- list.files(config$INPUT_DIR, pattern = "sample_metadata.*\\.csv$", full.names = TRUE, recursive = FALSE)
  
  if (dir.exists(runtime_folder)) {
    args <- c(args, "--runtime-dir", runtime_folder)
  }
  if (dir.exists(rrstats_folder)) {
    args <- c(args, "--rrstats-dir", rrstats_folder)
  }
  if (length(metadata_files) > 0) {
    args <- c(args, "--metadata-dir", config$INPUT_DIR)
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

# Run output processing
run_output_processing <- function(config, args = NULL) {
  if (is.null(args)) {
    args <- build_output_processing_args(config)
  }
  
  cat("--- Running output_processing.R ---\n")
  cat("Command arguments:\n")
  cat(paste(args, collapse = " "), "\n\n")
  
  start_time <- Sys.time()
  
  run_as_script(here("Ranalysis", "scripts", "output_processing.R"), args)
  
  end_time <- Sys.time()
  processing_time <- difftime(end_time, start_time, units = "mins")
  cat("\n--- output_processing.R completed successfully ---\n")
  cat("Processing time:", round(processing_time, 2), "minutes\n")
  
  return(TRUE)
}

# Load metadata files
load_metadata_files <- function(config) {
  cat("\n--- Loading metadata files ---\n")
  cat("Looking for metadata files in:", config$output_dir, "\n")
  
  # List all CSV files in output directory for debugging
  all_csv_files <- list.files(config$output_dir, pattern = "\\.csv$", full.names = FALSE)
  cat("All CSV files found:", paste(all_csv_files, collapse = ", "), "\n")
  
  metadata <- list(
    runtime_data = NULL,
    rrstats_data = NULL,
    sample_metadata = NULL,
    combined_metadata = NULL
  )
  
  tryCatch({
    # Try to load each type of metadata file - note the pattern includes timestamp
    file_patterns <- list(
      runtime_data = "runtime_info_.*op\\.csv$",
      rrstats_data = "read_statistics_.*op\\.csv$", 
      sample_metadata = "sample_metadata_.*op\\.csv$",
      combined_metadata = "combined_metadata_.*op\\.csv$"
    )
    
    for (data_type in names(file_patterns)) {
      pattern <- file_patterns[[data_type]]
      files <- list.files(config$output_dir, pattern = pattern, full.names = TRUE)
      
      if (length(files) > 0) {
        # Use the most recent file if multiple exist
        file_path <- files[which.max(file.info(files)$mtime)]
        metadata[[data_type]] <- read.csv(file_path, stringsAsFactors = FALSE)
        cat("Loaded", data_type, "from", basename(file_path), "for", nrow(metadata[[data_type]]), "samples\n")
      } else {
        cat("No files found for", data_type, "with pattern:", pattern, "\n")
      }
    }
    
  }, error = function(e) {
    cat("Warning: Could not load metadata files:", e$message, "\n")
  })
  
  return(metadata)
}

# Run annotation and generate plots
run_annotation_analysis <- function(config) {
  cat("\n--- Running annotation of species ---\n")
  
  # Look for species list files created by output_processing.R
  species_list_files <- list.files(config$output_dir, pattern = "species_list_unaligned.*op\\.csv$", full.names = TRUE)
  
  if (length(species_list_files) > 0) {
    # Load the species list for annotation
    species_list_file <- species_list_files[1]
    cat("Loading species list from:", basename(species_list_file), "\n")
    species_list_unaligned <- read.csv(species_list_file, stringsAsFactors = FALSE)
    
    cat("Loaded", nrow(species_list_unaligned), "species for annotation\n")
    
    # Put the data in global environment with expected name
    assign("species_list_unaligned", species_list_unaligned, envir = .GlobalEnv)
    
    # Run annotation script
    tryCatch({
      run_as_script(here("Ranalysis", "scripts", "annotate_species.R"), 
                    "--df", "species_list_unaligned",
                    "--output", "annotated_species_list")
      
      cat("Species annotation completed\n")
      
      # Generate annotation plots if annotated species data exists
      if (exists("annotated_species_list", envir = .GlobalEnv)) {
        cat("\n--- Generating annotation plots ---\n")
        # Get the annotated data
        annotated_data <- get("annotated_species_list", envir = .GlobalEnv)
        
        # Check if required annotation columns exist
        required_annotation_cols <- c("RiskGroup", "HOMD.Category")
        missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotated_data))
        
        if (length(missing_annotation_cols) == 0) {
          cat("Generating comprehensive annotation plots...\n")
          
          # Generate basic annotation plots
          tryCatch({
            # HOMD basic plot
            run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                          "--df", "annotated_species_list",
                          "--plot-type", "homd")
            
            # Risk group basic plot
            run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                          "--df", "annotated_species_list",
                          "--plot-type", "risk_group")
            
            # Kingdom basic plot
            run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                          "--df", "annotated_species_list",
                          "--plot-type", "kingdom")
            
            # Summary basic plot
            run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                          "--df", "annotated_species_list",
                          "--plot-type", "summary")
            
            cat("Basic annotation plots generated successfully\n")
            
            # Generate detailed plots if cladeReads_mean exists
            if ("cladeReads_mean" %in% colnames(annotated_data)) {
              cat("Generating detailed annotation plots with cladeReads data...\n")
              
              # Add log transformed cladeReads_mean and cladeReads_max
              annotated_data$cladeReads_mean_log <- log10(annotated_data$cladeReads_mean + 1)
              if ("cladeReads_max" %in% colnames(annotated_data)) {
                annotated_data$cladeReads_max_log <- log10(annotated_data$cladeReads_max + 1)
              }
              
              # Update the global environment with log-transformed data
              assign("annotated_species_list", annotated_data, envir = .GlobalEnv)
              
              # Generate detailed plots with log cladeReads_mean
              run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                            "--df", "annotated_species_list",
                            "--plot-type", "homd",
                            "--detailed", "cladeReads_mean_log")
              
              run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                            "--df", "annotated_species_list",
                            "--plot-type", "risk_group",
                            "--detailed", "cladeReads_mean_log")
              
              # Generate detailed plots with log cladeReads_max if available
              if ("cladeReads_max_log" %in% colnames(annotated_data)) {
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "homd",
                              "--detailed", "cladeReads_max_log")
                
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "risk_group",
                              "--detailed", "cladeReads_max_log")
              }
              
              # Generate detailed plots with Frequency if available
              if ("Freq" %in% colnames(annotated_data)) {
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "homd",
                              "--detailed", "Freq")
                
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "risk_group",
                              "--detailed", "Freq")
              }
              
              # Generate detailed plots for Bracken data if available
              bracken_cols <- grep("_bracken_reads$", colnames(annotated_data), value = TRUE)
              if (length(bracken_cols) > 0) {
                cat("Generating detailed annotation plots with Bracken data...\n")
                
                # Calculate mean and max Bracken reads
                bracken_matrix <- annotated_data[, bracken_cols, drop = FALSE]
                annotated_data$bracken_reads_mean <- rowMeans(bracken_matrix, na.rm = TRUE)
                annotated_data$bracken_reads_max <- apply(bracken_matrix, 1, max, na.rm = TRUE)
                
                # Add log transformed Bracken reads  
                annotated_data$bracken_reads_mean_log <- log10(annotated_data$bracken_reads_mean + 1)
                annotated_data$bracken_reads_max_log <- log10(annotated_data$bracken_reads_max + 1)
                
                # Update the global environment
                assign("annotated_species_list", annotated_data, envir = .GlobalEnv)
                
                # Generate detailed plots with log Bracken reads mean
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "homd",
                              "--detailed", "bracken_reads_mean_log")
                
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "risk_group",
                              "--detailed", "bracken_reads_mean_log")
                
                # Generate detailed plots with log Bracken reads max
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "homd",
                              "--detailed", "bracken_reads_max_log")
                
                run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                              "--df", "annotated_species_list",
                              "--plot-type", "risk_group",
                              "--detailed", "bracken_reads_max_log")
              } else {
                cat("Note: No Bracken data columns found - skipping Bracken detailed plots\n")
              }
            }
            
          }, error = function(e) {
            cat("Warning: Error generating annotation plots:", e$message, "\n")
          })
          
        } else {
          cat("Warning: Missing required annotation columns:", paste(missing_annotation_cols, collapse = ", "), "\n")
          cat("Skipping annotation plots generation\n")
        }
        
        # Save annotated data to output directory
        annotated_file <- file.path(config$output_dir, "annotated_species_list.csv")
        write.csv(annotated_data, annotated_file, row.names = FALSE)
        cat("Annotated species data saved to", basename(annotated_file), "\n")
        
      } else {
        cat("Warning: annotated_species_list not found - skipping annotation plots\n")
      }
      
    }, error = function(e) {
      cat("Warning: Could not run annotation analysis:", e$message, "\n")
    })
    
  } else {
    cat("Warning: No species list files found for annotation\n")
  }
}

# Function to create correlation analysis
create_correlation_analysis <- function(data1, data2, data1_name, data2_name, output_dir) {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get sample columns
  cols1 <- grep("_cladeReads$", colnames(data1), value = TRUE)
  cols2 <- grep("_cladeReads$", colnames(data2), value = TRUE)

  samples1 <- str_remove(cols1, "_cladeReads$")
  samples2 <- str_remove(cols2, "_cladeReads$")

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
    col1 <- paste0(sample_id, "_cladeReads")
    col2 <- paste0(sample_id, "_cladeReads")

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
        gg_correlation <- ggplot(sample_comparison, aes(x = log_reads1, y = log_reads2)) +
          geom_point(alpha = 0.6, size = 1.5) + 
          geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
          geom_text_repel(aes(label = label), na.rm = TRUE, size = 2.5, max.overlaps = 10) +
          labs(
            x = paste0("Log10(", data1_name, " cladeReads + 1)"),
            y = paste0("Log10(", data2_name, " cladeReads + 1)"),
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

# Function to create per-sample Kraken vs Bracken correlation analysis
create_kraken_bracken_correlation <- function(config) {
  cat("\n--- Running per-sample Kraken vs Bracken correlation analysis ---\n")
  
  # Look for merged files that contain both Kraken and Bracken data (exclude _long files)
  merged_files <- list.files(config$output_dir, pattern = "merged.*op\\.csv$", full.names = TRUE)
  merged_files <- merged_files[!grepl("_long", merged_files)]  # Exclude long format files
  
  if (length(merged_files) == 0) {
    cat("Warning: No merged files found for Kraken vs Bracken correlation analysis\n")
    return(NULL)
  }
  
  # Use the first merged file and check if it has both Kraken and Bracken columns
  merged_file <- merged_files[1]
  cat("Loading merged data from:", basename(merged_file), "\n")
  
  merged_data <- read.csv(merged_file, stringsAsFactors = FALSE)
  
  # Check for Kraken columns (_cladeReads) and Bracken columns (_bracken_reads)
  kraken_cols <- grep("_cladeReads$", colnames(merged_data), value = TRUE)
  bracken_cols <- grep("_bracken_.*_bracken_reads$", colnames(merged_data), value = TRUE)
  
  if (length(kraken_cols) == 0) {
    cat("Warning: No Kraken cladeReads columns found in merged file\n")
    return(NULL)
  }
  
  if (length(bracken_cols) == 0) {
    cat("Warning: No Bracken read columns found in merged file\n")
    return(NULL)
  }
  
  cat("Found", length(kraken_cols), "Kraken columns and", length(bracken_cols), "Bracken columns\n")
  
  # Extract sample names from both column types
  kraken_samples <- str_remove(kraken_cols, "_cladeReads$")
  bracken_samples <- str_extract(bracken_cols, "^[^_]+")  # Extract sample ID before first underscore
  
  # Find common samples
  common_samples <- intersect(kraken_samples, bracken_samples)
  
  if (length(common_samples) == 0) {
    cat("Warning: No common samples found between Kraken and Bracken data\n")
    cat("Kraken samples:", length(kraken_samples), "found\n")
    cat("Bracken samples:", length(bracken_samples), "found\n")
    return(NULL)
  }
  
  cat("Found", length(common_samples), "common samples for Kraken vs Bracken correlation\n")
  
  # Create pseudo-datasets for the correlation function
  # Prepare Kraken data (keep existing structure)
  kraken_data <- merged_data[, c("taxID", "name", kraken_cols)]
  
  # Prepare Bracken data (rename columns to match Kraken structure)
  bracken_data <- merged_data[, c("taxID", "name", bracken_cols)]
  
  # Rename bracken columns to match kraken pattern for common samples
  for (sample in common_samples) {
    kraken_col <- paste0(sample, "_cladeReads")
    bracken_col <- bracken_cols[grepl(paste0("^", sample, "_"), bracken_cols)]
    
    if (length(bracken_col) > 0) {
      # Rename the bracken column to match kraken pattern
      colnames(bracken_data)[colnames(bracken_data) == bracken_col[1]] <- kraken_col
    }
  }
  
  # Determine dataset type from filename
  dataset_type <- if (grepl("unaligned", merged_file)) "Unaligned" else "Nonhuman"
  
  # Create analysis names
  kraken_name <- paste0("Kraken_", dataset_type)
  bracken_name <- paste0("Bracken_", dataset_type)
  
  cat("Comparing", kraken_name, "vs", bracken_name, "using", length(common_samples), "samples\n")
  
  result <- create_correlation_analysis(kraken_data, bracken_data, kraken_name, bracken_name, config$reports_dir)
  
  if (!is.null(result)) {
    cat("Kraken vs Bracken correlation analysis completed successfully\n")
  }
  
  return(result)
}

# Load correlation data and run analysis
run_correlation_analysis <- function(config) {
  if (!config$PERFORM_CORRELATION) {
    return(NULL)
  }
  
  cat("\n--- Running correlation analysis ---\n")
  
  # Set USE_PARALLEL if it doesn't exist (for compatibility with output_processing.R functions)
  if (!exists("USE_PARALLEL", envir = .GlobalEnv)) {
    assign("USE_PARALLEL", FALSE, envir = .GlobalEnv)
  }
  
  # 1. Run traditional unaligned vs nonhuman correlation
  unaligned_files <- list.files(config$output_dir, pattern = "unaligned_merged.*op\\.csv$", full.names = TRUE)
  nonhuman_files <- list.files(config$output_dir, pattern = "nonhuman_merged.*op\\.csv$", full.names = TRUE)
  
  # Exclude _long files
  unaligned_files <- unaligned_files[!grepl("_long", unaligned_files)]
  nonhuman_files <- nonhuman_files[!grepl("_long", nonhuman_files)]
  
  correlation_result <- NULL
  if (length(unaligned_files) > 0 && length(nonhuman_files) > 0) {
    tryCatch({
      unaligned_data <- read.csv(unaligned_files[1], stringsAsFactors = FALSE)
      nonhuman_data <- read.csv(nonhuman_files[1], stringsAsFactors = FALSE)
      
      correlation_result <- create_correlation_analysis(unaligned_data, nonhuman_data, "Unaligned", "Nonhuman", config$reports_dir)
    }, error = function(e) {
      cat("Warning: Unaligned vs Nonhuman correlation analysis failed:", e$message, "\n")
    })
  } else {
    cat("Warning: Could not find merged data files for unaligned vs nonhuman correlation analysis\n")
  }
  
  # 2. Run Kraken vs Bracken correlation if Bracken files are available
  if (config$PROCESS_BRACKEN) {
    tryCatch({
      kraken_bracken_result <- create_kraken_bracken_correlation(config)
    }, error = function(e) {
      cat("Warning: Kraken vs Bracken correlation analysis failed:", e$message, "\n")
    })
  }
  
  return(correlation_result)
}

# Generate batch report
generate_batch_report <- function(config, metadata) {
  
  # Extract combined_metadata from metadata parameter (if available)
  combined_metadata <- metadata$combined_metadata
  
  cat("\n--- Generating batch report ---\n")
  
  # Load the combined_report_data
  combined_report_files <- list.files(config$output_dir, pattern = "sample_report_data.*op\\.csv$", full.names = TRUE)
  
  if (length(combined_report_files) == 0) {
    cat("No sample report data files found for batch report generation\n")
    return(NULL)
  }
  
  combined_report_file <- combined_report_files[1]  # Use the first (most recent) file
  cat("Loading combined report data from:", basename(combined_report_file), "\n")
  combined_report_data <- read.csv(combined_report_file, stringsAsFactors = FALSE)
  
  # Display info about the loaded data
  cat("Loaded report data with", nrow(combined_report_data), "rows and", ncol(combined_report_data), "columns\n")
  
  # Show which optional columns are present
  optional_cols <- c("confidence_levels", "minimum_hit_groups", "human_reads", "database_used")
  present_optional <- optional_cols[optional_cols %in% colnames(combined_report_data)]
  if (length(present_optional) > 0) {
    cat("Optional columns detected:", paste(present_optional, collapse = ", "), "\n")
  }
  
  # Check for subspecies data
  subspecies_cols <- grep("^S[123]_", colnames(combined_report_data), value = TRUE)
  if (length(subspecies_cols) > 0) {
    cat("Subspecies columns detected:", length(subspecies_cols), "columns\n")
  }
  
  # Initialize combined_report_with_metadata as the basic report data
  combined_report_with_metadata <- combined_report_data
  
  # Merge metadata if available
  if (!is.null(combined_metadata)) {
    cat("Successfully parsed metadata for", nrow(combined_metadata), "samples\n")
    cat("Merging metadata with combined report data...\n")
    
    combined_report_with_metadata <- merge(combined_report_data, combined_metadata, by = "sample", all.x = TRUE)
    cat("Merged data has", nrow(combined_report_with_metadata), "rows\n")
  } else {
    cat("No metadata available - generating basic batch report\n")
  }
  
  # Save the updated combined report data
  updated_file <- file.path(config$output_dir, paste0("sample_report_data_with_metadata", format(Sys.time(), "%y%m%d"), "op.csv"))
  write.csv(combined_report_with_metadata, updated_file, row.names = FALSE)
  cat("Saved combined report data to:", basename(updated_file), "\n")
  
  # Export uncombined report dataframes for each dataset
  cat("\nExporting dataset-specific report dataframes...\n")
  
  # Check if dataset column exists
  if ("dataset" %in% colnames(combined_report_with_metadata)) {
        unique_datasets <- unique(combined_report_with_metadata$dataset)
        cat("Found datasets:", paste(unique_datasets, collapse = ", "), "\n")
        
        for (dataset in unique_datasets) {
          if (!is.na(dataset) && dataset != "") {
            dataset_data <- combined_report_with_metadata[combined_report_with_metadata$dataset == dataset, ]
            
            # Remove completely empty columns
            dataset_data <- dataset_data[, !apply(is.na(dataset_data) | dataset_data == "", 2, all)]
            
            # Save dataset-specific file
            dataset_file <- file.path(config$output_dir, paste0("sample_report_data_", dataset, "_", format(Sys.time(), "%y%m%d"), "op.csv"))
            write.csv(dataset_data, dataset_file, row.names = FALSE)
            cat("Saved", dataset, "dataset report to:", basename(dataset_file), "\n")
            cat("  - Contains", nrow(dataset_data), "samples and", ncol(dataset_data), "columns\n")
            
            # Show sample of columns for verification
            runtime_cols <- grep("_rt$", colnames(dataset_data), value = TRUE)
            if (length(runtime_cols) > 0) {
              cat("  - Runtime columns:", paste(runtime_cols, collapse = ", "), "\n")
            }
            
            if ("condition" %in% colnames(dataset_data)) {
              conditions <- unique(dataset_data$condition[!is.na(dataset_data$condition)])
              if (length(conditions) > 0) {
                cat("  - Conditions:", paste(conditions, collapse = ", "), "\n")
              }
            }
          }
        }
      } else {
        cat("No dataset column found in combined data. Saving single combined file only.\n")
      }
      
      # Display summary
      if (!is.null(combined_metadata) && !is.null(metadata$runtime_data)) {
        cat("\nRuntime data summary:\n")
        runtime_summary_cols <- c("cutadapt_star_samtools_rt", "kraken2_unaligned_rt", "kraken2_nonhuman_rt")
        print(summary(combined_metadata[, runtime_summary_cols[runtime_summary_cols %in% colnames(combined_metadata)]]))
        
        # Create Kraken database loading samples dataframe
        kraken_db_load_samples <- combined_metadata[combined_metadata$is_first_k2_unaligned | combined_metadata$is_first_k2_nonhuman, ]
        
        if (nrow(kraken_db_load_samples) > 0) {
          cat("\nKraken database loading samples:\n")
          
          # Create clean display with only relevant columns
          kraken_display <- data.frame(
            sample = kraken_db_load_samples$sample,
            kraken2_unaligned_rt = ifelse(kraken_db_load_samples$is_first_k2_unaligned, 
                                          kraken_db_load_samples$kraken2_unaligned_rt, NA),
            kraken2_nonhuman_rt = ifelse(kraken_db_load_samples$is_first_k2_nonhuman, 
                                         kraken_db_load_samples$kraken2_nonhuman_rt, NA),
            stringsAsFactors = FALSE
          )
          
          # Remove NA columns for cleaner display
          kraken_display <- kraken_display[, !apply(is.na(kraken_display), 2, all)]
          
          print(kraken_display)
          
          # Also show comparison with other samples
          if (any(kraken_db_load_samples$is_first_k2_unaligned)) {
            first_unaligned_rt <- kraken_db_load_samples$kraken2_unaligned_rt[kraken_db_load_samples$is_first_k2_unaligned][1]
            other_unaligned_rt <- median(combined_metadata$kraken2_unaligned_rt[!combined_metadata$is_first_k2_unaligned], na.rm = TRUE)
            if (!is.na(first_unaligned_rt) && !is.na(other_unaligned_rt)) {
              cat("K2 unaligned: First sample RT =", first_unaligned_rt, "s, Median other samples RT =", round(other_unaligned_rt), "s\n")
            }
          }
          
          if (any(kraken_db_load_samples$is_first_k2_nonhuman)) {
            first_nonhuman_rt <- kraken_db_load_samples$kraken2_nonhuman_rt[kraken_db_load_samples$is_first_k2_nonhuman][1]
            other_nonhuman_rt <- median(combined_metadata$kraken2_nonhuman_rt[!combined_metadata$is_first_k2_nonhuman], na.rm = TRUE)
            if (!is.na(first_nonhuman_rt) && !is.na(other_nonhuman_rt)) {
              cat("K2 non-human: First sample RT =", first_nonhuman_rt, "s, Median other samples RT =", round(other_nonhuman_rt), "s\n")
            }
          }
        }
      }
      
      if (!is.null(combined_metadata) && !is.null(metadata$rrstats_data)) {
        cat("\nRead statistics summary:\n")
        stats_summary_cols <- c("num_input_reads", "uniquely_mapped_reads", "percent_mapped_reads", "uniquely_mapped_percent")
        print(summary(combined_metadata[, stats_summary_cols[stats_summary_cols %in% colnames(combined_metadata)]]))
      }
      
      # Show first few rows
      cat("\nFirst 5 rows of combined data with metadata:\n")
      print(head(combined_report_with_metadata, 5))
  
  
  # Generate batch report PDF if we have data (metadata is optional)
  if (!is.null(combined_report_with_metadata)) {
    
    # Load required libraries for PDF generation
    suppressPackageStartupMessages({
      library(ggplot2)
      library(gridExtra)
      library(grid)
      library(knitr)
      library(kableExtra)
    })
    
    runtime_folder <- file.path(config$INPUT_DIR, "runtimes")
    rrstats_folder <- file.path(config$INPUT_DIR, "rrstats")
    metadata_files <- list.files(config$INPUT_DIR, pattern = "sample_metadata.*\\.csv$", full.names = TRUE, recursive = FALSE)
    
    # Create PDF file path
    pdf_file <- file.path(config$reports_dir, paste0("batch_report_", config$PROJ_NAME, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
    
    # Start PDF device with larger pages to accommodate plots
    pdf(pdf_file, width = 11, height = 8.5, onefile = TRUE)
    
    # Calculate metadata for Part 1
    earliest_timestamp <- NULL
    latest_timestamp <- Sys.time()
    
    if (!is.null(combined_metadata) && !is.null(metadata$runtime_data)) {
      # Try to extract timestamps from runtime files to get processing start time
      runtime_files <- list.files(runtime_folder, pattern = "^runtime", full.names = TRUE)
      if (length(runtime_files) > 0) {
        file_times <- file.info(runtime_files)$mtime
        if (length(file_times) > 0) {
          earliest_timestamp <- min(file_times, na.rm = TRUE)
        }
      }
    }
    
    sample_metadata <- metadata$sample_metadata
    
    num_samples <- if (!is.null(combined_metadata)) nrow(combined_metadata) else nrow(combined_report_with_metadata)
    num_positive <- if (!is.null(sample_metadata) && "condition" %in% colnames(sample_metadata)) {
      sum(grepl("positive|pos", sample_metadata$condition, ignore.case = TRUE), na.rm = TRUE)
    } else {
      "N/A"
    }
    
    overall_process_time <- if (!is.null(earliest_timestamp)) {
      round(difftime(latest_timestamp, earliest_timestamp, units = "hours"), 2)
    } else {
      "N/A"
    }
    
    # PAGE 1: Metadata and Summary
    plot.new()
    
    # Create text plots for metadata
    metadata_text <- paste0(
      "SPARKEN Batch Report - Sample Processing Summary\n\n",
      "Part 1: Metadata\n",
      "Samples collection from: June 5 to June 9 2025\n",
      "Sequencing date: \n",
      "Processed: ", if (!is.null(earliest_timestamp)) format(earliest_timestamp, "%Y-%m-%d %H:%M:%S") else "N/A", "\n",
      "Overall Process time: ", overall_process_time, " hours\n",
      "by: Denis Rivard\n",
      "Number of samples: ", num_samples, "\n",
      "Number of positive samples: ", num_positive, "\n",
      "RNA kit used: COVID-19 RNA kit\n\n"
    )
    
    # Add RNAseq statistics
    if (!is.null(combined_metadata) && !is.null(metadata$rrstats_data)) {
      rnaseq_text <- paste0(
        "Part 2: Batch RNAseq Statistics\n",
        "Total samples with read statistics: ", nrow(metadata$rrstats_data), "\n",
        "Samples with <10M unique reads: ", sum(metadata$rrstats_data$less_than_10M_unique, na.rm = TRUE), "\n",
        "Average input reads: ", format(round(mean(metadata$rrstats_data$num_input_reads, na.rm = TRUE)), big.mark = ","), "\n",
        "Average mapping percentage: ", round(mean(metadata$rrstats_data$percent_mapped_reads, na.rm = TRUE), 1), "%\n",
        "Average unique mapping percentage: ", round(mean(metadata$rrstats_data$uniquely_mapped_percent, na.rm = TRUE), 1), "%\n\n"
      )
    } else {
      rnaseq_text <- "Part 2: Batch RNAseq Statistics\nNo read statistics data available\n\n"
    }
    
    # Add runtime information
    if (!is.null(metadata$runtime_data)) {
      runtime_text <- paste0(
        "Part 3: Runtime Information\n",
        "Total samples with runtime data: ", nrow(metadata$runtime_data), "\n",
        "Average Cutadapt/STAR/Samtools runtime: ", round(mean(metadata$runtime_data$cutadapt_star_samtools_rt, na.rm = TRUE)), " seconds\n",
        "Average Kraken2 unaligned runtime: ", round(mean(metadata$runtime_data$kraken2_unaligned_rt, na.rm = TRUE)), " seconds\n",
        "Average Kraken2 non-human runtime: ", round(mean(metadata$runtime_data$kraken2_nonhuman_rt, na.rm = TRUE)), " seconds\n",
        "Samples involved in DB loading: ", sum(metadata$runtime_data$is_first_k2_unaligned | metadata$runtime_data$is_first_k2_nonhuman, na.rm = TRUE), "\n\n"
      )
    } else {
      runtime_text <- "Part 3: Runtime Information\nNo runtime data available\n\n"
    }
    
    full_text <- paste0(metadata_text, rnaseq_text, runtime_text, 
                        "Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    
    text(0.05, 0.95, full_text, adj = c(0, 1), cex = 0.8, family = "mono")
    
    # PAGE 2: RNAseq Statistics Table Reference
    if (!is.null(metadata$rrstats_data)) {
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
        "- Total samples with read statistics: ", nrow(metadata$rrstats_data), "\n",
        "- Samples with <10M unique reads: ", sum(metadata$rrstats_data$less_than_10M_unique, na.rm = TRUE), "\n",
        "- Average input reads: ", format(round(mean(metadata$rrstats_data$num_input_reads, na.rm = TRUE)), big.mark = ","), "\n",
        "- Average mapping percentage: ", round(mean(metadata$rrstats_data$percent_mapped_reads, na.rm = TRUE), 1), "%\n",
        "- Average unique mapping percentage: ", round(mean(metadata$rrstats_data$uniquely_mapped_percent, na.rm = TRUE), 1), "%"
      )
      
      text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    }
    
    # PAGE 3: Runtime Information Table Reference 
    if (!is.null(metadata$runtime_data)) {
      plot.new()
      
      title_text <- "Part 3: Runtime Information"
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
        "- Total samples with runtime data: ", nrow(metadata$runtime_data), "\n",
        "- Average Cutadapt/STAR/Samtools runtime: ", round(mean(metadata$runtime_data$cutadapt_star_samtools_rt, na.rm = TRUE)), " seconds\n",
        "- Average Kraken2 unaligned runtime: ", round(mean(metadata$runtime_data$kraken2_unaligned_rt, na.rm = TRUE)), " seconds\n",
        "- Average Kraken2 non-human runtime: ", round(mean(metadata$runtime_data$kraken2_nonhuman_rt, na.rm = TRUE)), " seconds\n",
        "- Samples involved in DB loading: ", sum(metadata$runtime_data$is_first_k2_unaligned | metadata$runtime_data$is_first_k2_nonhuman, na.rm = TRUE)
      )
      
      text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    }
    
    # PAGE 4: Kraken2 Results Table Reference
    if (!is.null(combined_report_data)) {
      plot.new()
      
      title_text <- "Part 4: Kraken2 Results"
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
        "- Total records: ", nrow(combined_report_data), "\n",
        "- Unique samples: ", length(unique(combined_report_data$sample)), "\n",
        "- Unique species: ", length(unique(combined_report_data$species)), "\n"
      )
      
      # Add dataset info if available
      if ("dataset" %in% colnames(combined_report_data)) {
        datasets <- unique(combined_report_data$dataset)
        dataset_text <- paste0("- Datasets: ", paste(datasets, collapse = ", "))
        ref_text <- paste0(ref_text, dataset_text)
      }
      
      text(0.05, 0.85, ref_text, adj = c(0, 1), cex = 0.7, family = "mono")
    }
    
    # Runtime Plots - each on its own page
    if (!is.null(metadata$runtime_data)) {
      
      # Prepare data for plotting
      runtime_plot_data <- metadata$runtime_data
      
      # Add asterisks to sample names for first samples
      runtime_plot_data$sample_display <- ifelse(
        runtime_plot_data$is_first_k2_unaligned | runtime_plot_data$is_first_k2_nonhuman,
        paste0(runtime_plot_data$sample, "*"),
        runtime_plot_data$sample
      )
      
      # Add read statistics if available
      if (!is.null(metadata$rrstats_data)) {
        runtime_plot_data <- merge(runtime_plot_data, 
                                   metadata$rrstats_data[, c("sample", "num_input_reads", "percent_mapped_reads")], 
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
      k2_unaligned_data <- runtime_plot_data
      k2_unaligned_data$db_loading_time <- ifelse(k2_unaligned_data$is_first_k2_unaligned, 
                                                  pmax(0, k2_unaligned_data$kraken2_unaligned_rt - avg_k2_unaligned_non_db), 0)
      k2_unaligned_data$classification_time <- ifelse(k2_unaligned_data$is_first_k2_unaligned,
                                                      avg_k2_unaligned_non_db,
                                                      k2_unaligned_data$kraken2_unaligned_rt)
      
      # Reshape for stacked bar plot
      k2_unaligned_long <- data.frame(
        sample = rep(k2_unaligned_data$sample_display, 2),
        runtime_type = rep(c("Classification", "DB Loading"), each = nrow(k2_unaligned_data)),
        runtime = c(k2_unaligned_data$classification_time, k2_unaligned_data$db_loading_time),
        total_runtime = rep(k2_unaligned_data$kraken2_unaligned_rt, 2)
      )
      
      p3 <- ggplot(k2_unaligned_long, aes(x = reorder(sample, -total_runtime), y = runtime, fill = runtime_type)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("Classification" = "orange", "DB Loading" = "darkred")) +
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
      
      # Plot 4: Kraken2 Non-human Runtime (with DB loading stacked)
      # Prepare stacked data for K2 non-human
      k2_nonhuman_data <- runtime_plot_data
      k2_nonhuman_data$db_loading_time <- ifelse(k2_nonhuman_data$is_first_k2_nonhuman, 
                                                 pmax(0, k2_nonhuman_data$kraken2_nonhuman_rt - avg_k2_nonhuman_non_db), 0)
      k2_nonhuman_data$classification_time <- ifelse(k2_nonhuman_data$is_first_k2_nonhuman,
                                                     avg_k2_nonhuman_non_db,
                                                     k2_nonhuman_data$kraken2_nonhuman_rt)
      
      # Reshape for stacked bar plot
      k2_nonhuman_long <- data.frame(
        sample = rep(k2_nonhuman_data$sample_display, 2),
        runtime_type = rep(c("Classification", "DB Loading"), each = nrow(k2_nonhuman_data)),
        runtime = c(k2_nonhuman_data$classification_time, k2_nonhuman_data$db_loading_time),
        total_runtime = rep(k2_nonhuman_data$kraken2_nonhuman_rt, 2)
      )
      
      p4 <- ggplot(k2_nonhuman_long, aes(x = reorder(sample, -total_runtime), y = runtime, fill = runtime_type)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("Classification" = "purple", "DB Loading" = "darkred")) +
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
    }
    
    # Annotation Plots Section - Add comprehensive annotation plots like experiment_comparison_v2.R
    if (exists("homd_plot", envir = .GlobalEnv) && exists("risk_group_plot", envir = .GlobalEnv)) {
      cat("Adding annotation plots to PDF...\n")
      
      # PAGE: Basic Annotation Plots (2x2 grid like experiment_comparison_v2.R)
      tryCatch({
        # Get the basic plots from global environment
        homd_basic <- get("homd_plot", envir = .GlobalEnv)
        risk_basic <- get("risk_group_plot", envir = .GlobalEnv)
        
        # Check for kingdom and summary plots
        kingdom_basic <- if (exists("kingdom_plot", envir = .GlobalEnv)) {
          get("kingdom_plot", envir = .GlobalEnv)
        } else {
          # Create empty plot if kingdom plot doesn't exist
          ggplot() + theme_void() + labs(title = "Kingdom Plot Not Available")
        }
        
        summary_basic <- if (exists("summary_plot", envir = .GlobalEnv)) {
          get("summary_plot", envir = .GlobalEnv)
        } else {
          # Create empty plot if summary plot doesn't exist
          ggplot() + theme_void() + labs(title = "Summary Plot Not Available")
        }
        
        # Create 2x2 grid of basic annotation plots with smaller size to prevent cutoff
        suppressPackageStartupMessages(library(gridExtra))
        basic_grid <- grid.arrange(
          homd_basic + labs(title = paste("HOMD Categories -", config$PROJ_NAME)) + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
          risk_basic + labs(title = paste("Risk Groups -", config$PROJ_NAME)) + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
          kingdom_basic + labs(title = paste("Kingdom Distribution -", config$PROJ_NAME)) + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
          summary_basic + labs(title = paste("Annotation Summary -", config$PROJ_NAME)) + theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")),
          nrow = 2, ncol = 2,
          top = paste("Species Annotation Overview -", config$PROJ_NAME),
          heights = c(0.45, 0.45),
          widths = c(0.48, 0.48)
        )
        
        cat("Basic annotation plots page added\n")
        
      }, error = function(e) {
        cat("Warning: Error creating basic annotation plots page:", e$message, "\n")
        # Create placeholder page
        plot.new()
        text(0.5, 0.5, "Annotation plots could not be generated", cex = 1.2, adj = 0.5)
      })
      
      # PAGE: Detailed Annotation Plots with log cladeReads_mean (if available)
      if (exists("detailed_homd_plot", envir = .GlobalEnv) && exists("detailed_risk_group_plot", envir = .GlobalEnv)) {
        homd_detailed_mean <- get("detailed_homd_plot", envir = .GlobalEnv) + 
          labs(title = paste(config$PROJ_NAME, "- HOMD Categories (log cladeReads_mean)"))
        
        risk_detailed_mean <- get("detailed_risk_group_plot", envir = .GlobalEnv) + 
          labs(title = paste(config$PROJ_NAME, "- Risk Groups (log cladeReads_mean)"))
        
        # Create 1x2 grid for detailed mean plots
        detailed_mean_grid <- grid.arrange(
          homd_detailed_mean,
          risk_detailed_mean,
          nrow = 2, ncol = 1,
          top = paste("Detailed Species Analysis (log cladeReads_mean) -", config$PROJ_NAME)
        )
        
        cat("Detailed annotation plots (cladeReads_mean) page added\n")
      }
      
      # PAGE: Detailed Annotation Plots with log cladeReads_max (if available)
      # Check if we have plots with cladeReads_max data by checking if the plot objects were updated
      
      # The detailed plots might have been overwritten, so check if we have max data
      if (exists("annotated_species_list", envir = .GlobalEnv)) {
        annotated_data <- get("annotated_species_list", envir = .GlobalEnv)
        
        if ("cladeReads_max_log" %in% colnames(annotated_data) && 
            exists("detailed_homd_plot", envir = .GlobalEnv) && 
            exists("detailed_risk_group_plot", envir = .GlobalEnv)) {
          
          # Re-generate plots specifically for cladeReads_max to ensure we have the right plots
          run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                        "--df", "annotated_species_list",
                        "--plot-type", "homd",
                        "--detailed", "cladeReads_max_log")
          homd_detailed_max <- get("detailed_homd_plot", envir = .GlobalEnv) + 
            labs(title = paste(config$PROJ_NAME, "- HOMD Categories (log cladeReads_max)"))
          
          run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                        "--df", "annotated_species_list",
                        "--plot-type", "risk_group",
                        "--detailed", "cladeReads_max_log")
          risk_detailed_max <- get("detailed_risk_group_plot", envir = .GlobalEnv) + 
            labs(title = paste(config$PROJ_NAME, "- Risk Groups (log cladeReads_max)"))
          
          # Create 1x2 grid for detailed max plots
          detailed_max_grid <- grid.arrange(
            homd_detailed_max,
            risk_detailed_max,
            nrow = 2, ncol = 1,
            top = paste("Detailed Species Analysis (log cladeReads_max) -", config$PROJ_NAME)
          )
          
          cat("Detailed annotation plots (cladeReads_max) page added\n")
        }
      }
      
      # PAGE: Detailed Annotation Plots with Bracken data (if available)
      
      if (exists("annotated_species_list", envir = .GlobalEnv)) {
        annotated_data <- get("annotated_species_list", envir = .GlobalEnv)
        
        # Check if we have Bracken data columns
        bracken_cols <- grep("_bracken_reads$", colnames(annotated_data), value = TRUE)
        
        if (length(bracken_cols) > 0 && "bracken_reads_mean_log" %in% colnames(annotated_data)) {
          cat("Adding Bracken detailed annotation plots to PDF...\n")
          
          # Generate plots for Bracken mean reads
          run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                        "--df", "annotated_species_list",
                        "--plot-type", "homd",
                        "--detailed", "bracken_reads_mean_log")
          homd_bracken_mean <- get("detailed_homd_plot", envir = .GlobalEnv) + 
            labs(title = paste(config$PROJ_NAME, "- HOMD Categories (log Bracken reads mean)"))
          
          run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                        "--df", "annotated_species_list",
                        "--plot-type", "risk_group",
                        "--detailed", "bracken_reads_mean_log")
          risk_bracken_mean <- get("detailed_risk_group_plot", envir = .GlobalEnv) + 
            labs(title = paste(config$PROJ_NAME, "- Risk Groups (log Bracken reads mean)"))
          
          # Create 1x2 grid for Bracken mean plots
          bracken_mean_grid <- grid.arrange(
            homd_bracken_mean,
            risk_bracken_mean,
            nrow = 2, ncol = 1,
            top = paste("Detailed Species Analysis (log Bracken reads mean) -", config$PROJ_NAME)
          )
          
          cat("Detailed annotation plots (Bracken reads mean) page added\n")
          
          # Generate plots for Bracken max reads if available
          if ("bracken_reads_max_log" %in% colnames(annotated_data)) {
            run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                          "--df", "annotated_species_list",
                          "--plot-type", "homd",
                          "--detailed", "bracken_reads_max_log")
            homd_bracken_max <- get("detailed_homd_plot", envir = .GlobalEnv) + 
              labs(title = paste(config$PROJ_NAME, "- HOMD Categories (log Bracken reads max)"))
            
            run_as_script(here("Ranalysis", "scripts", "annotate_species_plots.R"), 
                          "--df", "annotated_species_list",
                          "--plot-type", "risk_group",
                          "--detailed", "bracken_reads_max_log")
            risk_bracken_max <- get("detailed_risk_group_plot", envir = .GlobalEnv) + 
              labs(title = paste(config$PROJ_NAME, "- Risk Groups (log Bracken reads max)"))
            
            # Create 1x2 grid for Bracken max plots
            bracken_max_grid <- grid.arrange(
              homd_bracken_max,
              risk_bracken_max,
              nrow = 2, ncol = 1,
              top = paste("Detailed Species Analysis (log Bracken reads max) -", config$PROJ_NAME)
            )
            
            cat("Detailed annotation plots (Bracken reads max) page added\n")
          }
        } else {
          cat("Note: No Bracken data available for detailed plots\n")
        }
      }
    } else {
      cat("Note: Annotation plots not available - basic plots not found\n")
    }
    
    # Close PDF device
    tryCatch({
      dev.off()
    }, error = function(e) {
      cat("Warning: Error closing PDF device:", e$message, "\n")
      # Try to clean up any remaining graphics devices
      while (dev.cur() > 1) {
        dev.off()
      }
    })
    
    cat("Multi-page batch report saved to:", pdf_file, "\n")
    
  } else {
    cat("Insufficient data to generate batch report PDF\n")
  }
}

# Main execution function
main <- function() {
  # Check for help request first
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args || "-h" %in% args || length(args) == 0) {
    show_help()
    if (length(args) == 0) {
      stop("No arguments provided. Use --help for usage information.")
    } else {
      stop("Help requested. Execution stopped.")
    }
  }
  
  # Create configuration
  config <- create_config(args)
  
  # Create output directories
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
  }
  if (!dir.exists(config$reports_dir)) {
    dir.create(config$reports_dir, recursive = TRUE)
  }
  
  # Display configuration
  cat("=== Report Analysis Configuration ===\n")
  cat("Input directory:", config$INPUT_DIR, "\n")
  cat("Project name:", config$PROJ_NAME, "\n")
  cat("Output directory:", config$output_dir, "\n")
  cat("Reports directory:", config$reports_dir, "\n")
  cat("Top N species:", config$TOP_N, "\n")
  cat("Process bracken:", config$PROCESS_BRACKEN, "\n")
  cat("Perform correlation:", config$PERFORM_CORRELATION, "\n")
  cat("Cores to use:", config$CORES, "\n")
  cat("=====================================\n\n")
  
  # Setup parallel processing
  if (!exists("cl")) {
    cl <- setup_parallel(config$CORES)
  }
  
  # Start timing
  start_time <- Sys.time()
  
  cat("=== Starting Report Analysis Pipeline ===\n")
  cat("Project:", config$PROJ_NAME, "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Run main pipeline
  # 1. Run output processing
  run_output_processing(config)
  
  # 2. Load metadata
  metadata <- load_metadata_files(config)
  
  # 3. Run correlation analysis
  run_correlation_analysis(config)
  
  # 4. Run annotation analysis
  run_annotation_analysis(config)
  
  # 5. Generate batch report
  generate_batch_report(config, metadata)
  
  # Final summary
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  
  cat("\n=== Pipeline Completed Successfully ===\n")
  cat("Project:", config$PROJ_NAME, "\n")
  cat("Total runtime:", round(total_time, 2), "minutes\n")
  cat("Output directory:", config$output_dir, "\n")
  cat("Reports directory:", config$reports_dir, "\n")
  cat("==========================================\n")
  
  # Clean up parallel processing
  if (exists("cl") && !is.null(cl)) {
      # stopCluster(cl)
  }
}

# Execute main function if script is run directly
if (!interactive()) {
  main()
}
