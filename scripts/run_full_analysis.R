#!/usr/bin/env Rscript

# ==============================================================================
#
#          SINGLE-SCRIPT K2RNASEQ PIPELINE ANALYSIS & REPORTING
#
# ==============================================================================
#
# Author: Gemini 2.5 Pro
# Date: August 15, 2025
#
# Description:
# This script consolidates the entire R-based analysis pipeline for the
# K2RNASeq project into a single, executable file. It handles everything
# from initial data ingestion of Kraken2/Bracken reports to final PDF report
# generation, including data filtering, correlation analysis, species
# annotation, metadata processing, and comprehensive plotting.
#
# By integrating all steps, this script aims to:
#   - Eliminate redundant file I/O between steps.
#   - Improve performance and maintainability.
#   - Provide a single point of execution for the entire analysis.
#
# Usage:
#   Rscript run_full_analysis.R --input-dir <path> --run-name <name> [OPTIONS]
#
# ==============================================================================
#
# SCRIPT STRUCTURE:
#
#   1. HEADER & SETUP
#      - Package Loading
#      - Argument Parsing & Help Message
#      - Global Configuration
#
#   2. HELPER FUNCTIONS
#      - Core utilities for data reading, processing, and plotting.
#
#   3. DATA INGESTION
#      - Reading Kraken2 reports, Bracken abundance files, and metadata.
#
#   4. DATA FILTERING & PRE-PROCESSING
#      - Applying thresholds (clade reads, minimizers).
#      - Merging datasets and creating long-format data.
#
#   5. DATA ANALYSIS
#      - Correlation analysis (Unaligned vs. Non-human, Kraken vs. Bracken).
#      - Species frequency analysis.
#
#   6. SPECIES ANNOTATION
#      - Annotating unique species with risk group and HOMD information.
#
#   7. METADATA ANALYSIS
#      - Parsing and summarizing runtime and read statistics.
#
#   8. PLOTTING
#      - Generating all analytical and summary plots.
#
#   9. DATA & REPORT SAVING
#      - Saving processed dataframes (CSV).
#      - Generating all PDF reports (Batch Report, Correlations, etc.).
#
#   10. BATCH REPORT GENERATION
#      - Compiling all results and plots into a comprehensive multi-page PDF report.
#
# ==============================================================================


# ==============================================================================
# SECTION 1: HEADER & SETUP
# ==============================================================================

# --- 1.1: Package Loading ---
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(gridExtra)
  library(knitr)
  library(kableExtra)
  library(doParallel)
  library(foreach)
  library(iterators)
  library(RColorBrewer)
  library(tibble)
})

# --- 1.2: Argument Parsing & Help Message ---
option_list <- list(
  # Required
  make_option(c("--input_dir"), type = "character", default = NULL,
              help = "Path to the main input directory containing kreport and bracken subdirectories."),
  make_option(c("--run_name"), type = "character", default = NULL,
              help = "Project name used for naming output files and directories."),

  # Optional Paths
  make_option(c("--output_base_dir"), type = "character", default = "Ranalysis",
              help = "Base directory for outputs [default: %default]"),
  make_option(c("--output_prefix"), type = "character", default = NULL,
              help = "Optional prefix for all output files."),
  make_option(c("--kreport_glob"), type = "character", default = NULL,
              help = "Glob pattern to find k-report files (e.g., 'path/to/*-kreport.txt'). Overrides --input_dir for k-reports."),
  make_option(c("--db_dir"), type = "character", default = "Ranalysis/databases",
              help = "Directory containing annotation databases (epathogen, HOMD) [default: %default]"),
  make_option(c("--metadata_dir"), type = "character", default = NULL,
              help = "Optional directory containing sample metadata files (CSV)."),
  make_option(c("--runtime_dir"), type = "character", default = "runtime",
              help = "Directory containing pipeline runtime files [default: %default]"),
  make_option(c("--rrstats_dir"), type = "character", default = "QC",
              help = "Directory containing read statistics files (run_read_stats) [default: %default]"),

  # Filtering & Analysis Parameters
  make_option(c("--top_n"), type = "integer", default = 25,
              help = "Number of top species to include in analysis [default: %default]"),
  make_option(c("--min_reads"), type = "integer", default = 0,
              help = "Minimum clade reads threshold for a species to be included [default: %default]"),
  make_option(c("--exclude_taxids"), type = "character", default = NULL,
              help = "Comma-separated list of taxonomy IDs to exclude."),
  make_option(c("--minimizer_ratio"), type = "double", default = NULL,
              help = "Filter by minimum ratio of distinct_minimizers/cladeReads."),
  make_option(c("--minimizer_threshold"), type = "integer", default = NULL,
              help = "Filter by minimum absolute distinct_minimizers threshold."),
  make_option(c("--no_bracken"), action = "store_true", default = FALSE,
              help = "Skip all Bracken file processing."),
  make_option(c("--no_corr"), action = "store_true", default = FALSE,
              help = "Skip correlation analysis."),
  make_option(c("--subspecies"), action = "store_true", default = FALSE,
              help = "Include subspecies-level (S1, S2, etc.) data in the analysis."),
  make_option(c("--generate_report"), action = "store_true", default = FALSE,
              help = "Generate the final multi-page PDF batch report."),

  # Execution
  make_option(c("--cores"), type = "integer", default = 0,
              help = "Number of cores for parallel processing. 0 detects automatically [default: %default]")
)

parser <- OptionParser(
  usage = "%prog --input_dir <DIR> --run_name <NAME> [OPTIONS]",
  option_list = option_list,
  description = "A single-script pipeline for K2RNASeq analysis and reporting."
)

args <- parse_args(parser)

# --- 1.3: Global Configuration & Validation ---
main_start_time <- Sys.time()

if (is.null(args$input_dir) || is.null(args$run_name)) {
  print_help(parser)
  stop("Missing required arguments: --input_dir and --run_name", call. = FALSE)
}

# Set up parallel processing
if (args$cores == 0) {
  cores <- detectCores()
  cat("Auto-detected", cores, "cores.\n")
} else {
  cores <- args$cores
  cat("Using", cores, "cores as specified.\n")
}
if (cores > 1) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  cat("Parallel cluster initialized with", length(cl), "workers.\n")
} else {
  registerDoParallel() # Use sequential backend
  cat("Running in sequential mode.\n")
}

# Create output directories
output_dir <- file.path(args$output_base_dir, paste0(args$run_name, "_outputs"))
reports_dir <- file.path(output_dir, "reports")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(reports_dir)) dir.create(reports_dir, recursive = TRUE)

# Display final configuration
cat("\n=== K2RNASeq Analysis Configuration ===\n")
cat("Run Name:           ", args$run_name, "\n")
cat("Input Directory:    ", args$input_dir, "\n")
cat("Output Directory:   ", output_dir, "\n")
cat("Reports Directory:  ", reports_dir, "\n")
cat("Annotation DBs:     ", args$db_dir, "\n")
cat("Top N Species:      ", args$top_n, "\n")
cat("Min Clade Reads:    ", args$min_reads, "\n")
cat("Process Bracken:    ", !args$no_bracken, "\n")
cat("Run Correlations:   ", !args$no_corr, "\n")
cat("Generate Report:    ", args$generate_report, "\n")
cat("Cores for Parallel: ", cores, "\n")
cat("=======================================\n\n")


# ==============================================================================
# SECTION 2: HELPER FUNCTIONS
# ==============================================================================
cat("Initializing helper functions...\n")

#' Read Kraken-style Report File
#'
#' Reads a standard Kraken2 report file, focusing on species-level data.
#'
#' @param myfile Path to the report file.
#' @return A tibble with columns: `name`, `taxID`, `cladeReads`, `taxonReads`.
read_report <- function(myfile) {
  tryCatch({
    report <- utils::read.table(myfile, sep = "\t", header = FALSE,
                      col.names = c("percentage", "cladeReads", "taxonReads", "taxRank", "taxID", "name"),
                      quote = "", stringsAsFactors = FALSE, comment.char = "#")
    
    if (nrow(report) == 0) return(tibble())

    report$name <- gsub("^ *", "", report$name)
    
    report %>% 
      filter(taxRank == "S") %>%
      select(name, taxID, cladeReads, taxonReads)
  }, error = function(e) {
     warning("Failed to read file ", myfile, ": ", e$message)
     return(tibble())
  })
}

#' Get Parameters from Filename
#'
#' Extracts key-value parameters from a Kraken2 report filename.
#'
#' @param filename The name of the file.
#' @return A tibble with columns: `confidence_levels`, `minimum_hit_groups`, `human_reads`.
get_params_from_filename <- function(filename) {
  confidence <- as.numeric(str_extract(filename, "c(?<val>[0-9.]+)", group = "val"))
  min_hits <- as.numeric(str_extract(filename, "m(?<val>[0-9.]+)", group = "val"))
  human_reads <- str_extract(filename, "(?<val>human|nonhuman)", group = "val")
  
  tibble(
    confidence_levels = confidence,
    minimum_hit_groups = min_hits,
    human_reads = human_reads
  )
}


#' Process a Folder of Kraken2 Reports
#'
#' Reads all k-reports in a folder and combines them into a single dataframe.
#'
#' @param folder_path Path to the directory containing report files.
#' @return A tibble containing the combined data from all reports.
process_kreport_folder <- function(folder_path) {
  files <- list.files(folder_path, pattern = "\\.k(2)?report$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No .kreport or .k2report files found in: ", folder_path)
    return(tibble())
  }
  
  cat("Processing", length(files), "k-report files from", basename(folder_path), "...\n")
  
  report_data <- map_dfr(files, ~{
    sample_name <- str_remove(basename(.x), "-.*$")
    report <- read_report(.x) 
    if (nrow(report) > 0) {
      report %>% mutate(sample = sample_name, .before = 1)
    } else {
      NULL
    }
  })
  
  return(report_data)
}

#' Process a set of Kraken2 reports based on a glob pattern.
#'
#' Reads all k-reports matching a glob pattern and combines them.
#' It determines if a report is 'unaligned' or 'nonhuman' from its filename.
#'
#' @param glob_pattern A glob pattern to find k-report files.
#' @return A list of two tibbles: `unaligned` and `nonhuman`.
process_kreport_glob <- function(glob_pattern) {
  files <- Sys.glob(glob_pattern)
  if (length(files) == 0) {
    warning("No files found matching glob pattern: ", glob_pattern)
    return(list(unaligned = tibble(), nonhuman = tibble()))
  }

  cat("Processing", length(files), "k-report files from glob pattern...\n")

  all_data <- map_dfr(files, ~{
    file_path <- .x
    sample_name <- str_remove(basename(file_path), "-.*$")
    report_type <- if (grepl("unaligned", basename(file_path), ignore.case = TRUE)) {
      "unaligned"
    } else if (grepl("nonhuman", basename(file_path), ignore.case = TRUE)) {
      "nonhuman"
    } else {
      "unknown"
    }

    report <- read_report(file_path)
    if (nrow(report) > 0) {
      if (report_type == "unknown") {
        warning("Skipping file with unknown type (not 'unaligned' or 'nonhuman' in filename): ", basename(file_path))
        return(NULL)
      }
      report %>% mutate(sample = sample_name, type = report_type, .before = 1)
    } else {
      NULL
    }
  })

  if (nrow(all_data) > 0) {
    unaligned_data <- all_data %>% filter(type == "unaligned") %>% select(-type)
    nonhuman_data <- all_data %>% filter(type == "nonhuman") %>% select(-type)
    return(list(unaligned = unaligned_data, nonhuman = nonhuman_data))
  } else {
    return(list(unaligned = tibble(), nonhuman = tibble()))
  }
}


#' Process a Folder of Bracken Abundance Files
#'
#' Reads all Bracken output files (.txt) in a folder and combines them.
#'
#' @param folder_path Path to the directory containing Bracken files.
#' @return A tibble with combined Bracken data.
process_bracken_folder <- function(folder_path) {
  files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)
   if (length(files) == 0) {
    warning("No .txt files found in Bracken directory: ", folder_path)
    return(tibble())
  }
  
  cat("Processing", length(files), "Bracken files from", basename(folder_path), "...\n")
  
  bracken_data <- map_dfr(files, ~{
    sample_name <- str_remove(basename(.x), "-.*$")
    tryCatch({
      readr::read_tsv(.x, col_types = readr::cols(), progress = FALSE) %>%
        mutate(sample = sample_name, .before = 1)
    }, error = function(e) {
      warning("Could not read Bracken file: ", .x)
      return(NULL)
    })
  })
  
  if (nrow(bracken_data) > 0) {
    bracken_data <- bracken_data %>%
      rename(
        bracken_reads = new_est_reads,
        bracken_fraction_total = fraction_total_reads
      ) %>%
      select(sample, name, taxID = taxonomy_id, bracken_reads, bracken_fraction_total)
  }
  
  return(bracken_data)
}

#' Parse and combine sample metadata files from a directory.
#'
#' @param input_dir Directory containing sample metadata CSV files.
#' @return A single tibble with combined and deduplicated metadata, or NULL.
parse_sample_metadata <- function(input_dir) {
  if (is.null(input_dir) || !dir.exists(input_dir)) {
    return(NULL)
  }
  
  metadata_files <- list.files(input_dir, pattern = "sample_metadata.*\\.csv$", full.names = TRUE, recursive = FALSE)
  if (length(metadata_files) == 0) {
    return(NULL)
  }
  
  cat("Processing", length(metadata_files), "sample metadata files...\n")
  
  all_metadata <- map_dfr(metadata_files, ~{
    tryCatch({
      df <- read.csv(.x, stringsAsFactors = FALSE)
      sample_col <- names(df)[tolower(names(df)) %in% c("sample", "samples")]
      condition_col <- names(df)[tolower(names(df)) %in% c("condition", "conditions")]
      
      if (length(sample_col) == 1 && length(condition_col) == 1) {
        df %>%
          rename(sample = all_of(sample_col), condition = all_of(condition_col)) %>%
          select(sample, condition, everything()) %>%
          filter(!is.na(sample) & sample != "")
      } else {
        warning("Skipping metadata file ", basename(.x), " due to missing 'sample' or 'condition' columns.")
        NULL
      }
    }, error = function(e) {
      warning("Error reading metadata file ", basename(.x), ": ", e$message)
      NULL
    })
  })
  
  if (nrow(all_metadata) > 0) {
    final_metadata <- all_metadata %>%
      group_by(sample) %>%
      slice_tail(n = 1) %>%
      ungroup()
    cat("Successfully parsed and combined metadata for", nrow(final_metadata), "unique samples.\n")
    return(final_metadata)
  }
  
  return(NULL)
}

#' Parse pipeline runtime files.
#'
#' @param runtime_folder Path to the directory with runtime files.
#' @return A tibble with parsed runtime data, or NULL.
parse_runtime_files <- function(runtime_folder) {
  if (is.null(runtime_folder) || !dir.exists(runtime_folder)) {
    return(NULL)
  }
  
  runtime_files <- list.files(runtime_folder, pattern = "^runtime", full.names = TRUE)
  if (length(runtime_files) == 0) {
    return(NULL)
  }
  
  cat("Processing", length(runtime_files), "runtime files...\n")
  
  all_runtime_data <- map_dfr(runtime_files, ~{
    lines <- readLines(.x, warn = FALSE)
    map_dfr(lines, function(line) {
      sample_name <- str_extract(line, "(?<=\\()([^)]+)(?=\\))")
      rt_seconds <- as.numeric(str_extract(line, "(?<=RT: )\\d+"))
      if (is.na(sample_name) || is.na(rt_seconds)) return(NULL)
      process_type <- case_when(
        grepl("^Cutadapt/STAR/Samtools", line) ~ "cutadapt_star_samtools_rt",
        grepl("^Kraken2 unaligned", line) ~ "kraken2_unaligned_rt",
        grepl("^Kraken2 non-human", line) ~ "kraken2_nonhuman_rt",
        TRUE ~ "unknown_rt"
      )
      sample_name_clean <- str_remove(sample_name, "_S\\d+$")
      tibble(sample = sample_name_clean, process = process_type, runtime = rt_seconds)
    })
  })
  
  if (nrow(all_runtime_data) > 0) {
    final_runtime <- all_runtime_data %>%
      group_by(sample, process) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      pivot_wider(names_from = process, values_from = runtime)
    cat("Successfully parsed runtime data for", nrow(final_runtime), "samples.\n")
    return(final_runtime)
  }
  
  return(NULL)
}

#' Parse read statistics files (run_read_stats.txt).
#'
#' @param rrstats_folder Path to the directory with read stat files.
#' @return A tibble with parsed read statistics, or NULL.
parse_rrstats_files <- function(rrstats_folder) {
  if (is.null(rrstats_folder) || !dir.exists(rrstats_folder)) {
    return(NULL)
  }
  
  rrstats_files <- list.files(rrstats_folder, pattern = "run_read_stats.*\\.txt$", full.names = TRUE)
  if (length(rrstats_files) == 0) {
    return(NULL)
  }
  
  cat("Processing", length(rrstats_files), "read stats files...\n")
  
  col_names <- c("sample", "num_input_reads", "avg_input_read_length", 
                 "uniquely_mapped_reads", "avg_mapped_length", 
                 "reads_assigned_to_genes", "percent_mapped_reads", 
                 "uniquely_mapped_percent", "less_than_10M_unique")
                 
  all_rrstats_data <- map_dfr(rrstats_files, ~{
    tryCatch({
      read.table(.x, header = FALSE, sep = "\t", col.names = col_names,
                 fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
    }, error = function(e) {
      warning("Error reading read stats file ", basename(.x), ": ", e$message)
      NULL
    })
  })
  
  if (nrow(all_rrstats_data) > 0) {
    final_rrstats <- all_rrstats_data %>%
      filter(!is.na(sample) & sample != "") %>%
      group_by(sample) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      mutate(across(where(is.character) & !matches("sample"), as.numeric))
    cat("Successfully parsed read statistics for", nrow(final_rrstats), "unique samples.\n")
    return(final_rrstats)
  }
  
  return(NULL)
}

#' Apply user-defined filters to the processed data.
#'
#' @param data A dataframe of combined k-report data.
#' @param min_reads Minimum clade reads threshold.
#' @param exclude_taxids A vector of taxIDs to exclude.
#' @param minimizer_ratio Minimum distinct_minimizers/cladeReads ratio.
#' @param minimizer_threshold Minimum distinct_minimizers count.
#' @return A filtered tibble.
apply_filters <- function(data, min_reads, exclude_taxids, minimizer_ratio, minimizer_threshold) {
  cat("Applying filters...\n")
  original_rows <- nrow(data)
  
  # 1. Minimum clade reads filter
  if (min_reads > 0) {
    data <- data %>% filter(cladeReads >= min_reads)
    cat("  - Applied min-reads filter (>= ", min_reads, "): ", nrow(data), " rows remaining.\n", sep = "")
  }
  
  # 2. Exclude taxID filter
  if (!is.null(exclude_taxids)) {
    exclude_ids <- as.numeric(strsplit(exclude_taxids, ",")[[1]])
    data <- data %>% filter(!taxID %in% exclude_ids)
    cat("  - Applied exclude-taxid filter: ", nrow(data), " rows remaining.\n", sep = "")
  }
  
  # 3. Minimizer ratio filter
  if (!is.null(minimizer_ratio) && "distinct_minimizers" %in% names(data)) {
    data <- data %>%
      filter(is.na(distinct_minimizers) | distinct_minimizers == 0 | 
             (cladeReads / pmax(distinct_minimizers, 1)) >= (1 / minimizer_ratio))
    cat("  - Applied minimizer-ratio filter (>= ", minimizer_ratio, "): ", nrow(data), " rows remaining.\n", sep = "")
  }
  
  # 4. Minimizer threshold filter
  if (!is.null(minimizer_threshold) && "distinct_minimizers" %in% names(data)) {
    data <- data %>%
      filter(is.na(distinct_minimizers) | distinct_minimizers >= minimizer_threshold)
    cat("  - Applied minimizer-threshold filter (>= ", minimizer_threshold, "): ", nrow(data), " rows remaining.\n", sep = "")
  }
  
  cat("Filtering complete. Total rows removed:", original_rows - nrow(data), "\n")
  return(data)
}

#' Save a Dataframe to CSV
#'
#' @param data The dataframe to save.
#' @param filename_prefix A descriptive prefix for the filename.
#' @param run_name The project name.
#' @param output_dir The directory to save the file in.
#' @param output_prefix An optional prefix for the filename.
save_data <- function(data, filename_prefix, run_name, output_dir, output_prefix = NULL) {
  if (nrow(data) > 0) {
    timestamp <- format(Sys.time(), "%Y%m%d")
    final_prefix <- if (!is.null(output_prefix)) paste0(output_prefix, "_") else ""
    file_path <- file.path(output_dir, paste0(final_prefix, filename_prefix, "_", run_name, "_", timestamp, ".csv"))
    write.csv(data, file_path, row.names = FALSE)
    cat("Successfully saved", nrow(data), "rows to", basename(file_path), "\n")
  } else {
    cat("Skipped saving", filename_prefix, "because it was empty.\n")
  }
}

#' Standard ggplot2 Theme
my_ggplot_theme <- theme_bw() + 
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

#' Create a bar chart for annotation categories.
#'
#' @param data The annotated dataframe.
#' @param category_col The column name of the category to plot.
#' @param title The title for the plot.
#' @return A ggplot object.
create_annotation_bar_chart <- function(data, category_col, title) {
  if (!category_col %in% names(data)) return(NULL)
  
  plot_data <- data %>%
    count(!!sym(category_col), name = "count") %>%
    mutate(percentage = count / sum(count) * 100) %>%
    arrange(desc(count))
    
  ggplot(plot_data, aes(x = reorder(!!sym(category_col), -count), y = count, fill = !!sym(category_col))) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = -0.5, size = 3) +
    labs(
      title = title,
      x = category_col,
      y = "Number of Species"
    ) +
    my_ggplot_theme +
    theme(legend.position = "none")
}

#' Create a plot of top N species based on a metric.
#'
#' @param data The long-format annotated dataframe.
#' @param metric_col The metric to rank species by (e.g., "cladeReads").
#' @param top_n The number of species to show.
#' @param color_by The annotation column to color bars by (e.g., "RiskGroup").
#' @param title The plot title.
#' @return A ggplot object.
create_top_species_plot <- function(data, metric_col, top_n, color_by, title) {
  if (!metric_col %in% names(data)) return(NULL)
  
  plot_data <- data %>%
    group_by(name, !!sym(color_by)) %>%
    summarise(mean_metric = mean(!!sym(metric_col), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_metric)) %>%
    slice_head(n = top_n)
    
  ggplot(plot_data, aes(x = reorder(name, mean_metric), y = mean_metric, fill = !!sym(color_by))) +
    geom_col() +
    coord_flip() +
    labs(
      title = title,
      subtitle = paste("Top", top_n, "species, colored by", color_by),
      x = "Species",
      y = paste("Mean", metric_col)
    ) +
    my_ggplot_theme +
    theme(legend.position = "right")
}

#' Save a list of ggplot objects to a multi-page PDF.
#'
#' @param plots_list A named list of ggplot objects.
#' @param filename The name of the output PDF file.
#' @param output_dir The directory to save the PDF in.
#' @param output_prefix An optional prefix for the filename.
save_plots_to_pdf <- function(plots_list, filename, output_dir, output_prefix = NULL) {
  final_prefix <- if (!is.null(output_prefix)) paste0(output_prefix, "_") else ""
  pdf_path <- file.path(output_dir, paste0(final_prefix, filename))
  
  pdf(pdf_path, width = 11, height = 8.5)
  for (plot_name in names(plots_list)) {
    if (!is.null(plots_list[[plot_name]])) {
      print(plots_list[[plot_name]])
    }
  }
  dev.off()
  
  cat("Saved plots to", basename(pdf_path), "\n")
}

cat("Helper functions initialized.\n\n")


# ==============================================================================
# SECTION 3: DATA INGESTION
# ==============================================================================
cat("\n\n=======================================\n")
cat("  SECTION 3: DATA INGESTION\n")
cat("=======================================\n\n")

# --- 3.1: Define input paths ---
unaligned_bracken_dir  <- file.path(args$input_dir, "bracken")
nonhuman_bracken_dir   <- file.path(args$input_dir, "bracken")

# --- 3.2: Process Kraken2 Reports ---
if (!is.null(args$kreport_glob)) {
  cat("Processing k-reports using glob pattern:", args$kreport_glob, "\n")
  glob_results <- process_kreport_glob(args$kreport_glob)
  unaligned_kreport_data <- glob_results$unaligned
  nonhuman_kreport_data <- glob_results$nonhuman
} else {
  unaligned_kreports_dir <- file.path(args$input_dir, "kraken2")
  nonhuman_kreports_dir  <- file.path(args$input_dir, "kraken2_nonhuman")
  unaligned_kreport_data <- process_kreport_folder(unaligned_kreports_dir)
  nonhuman_kreport_data  <- process_kreport_folder(nonhuman_kreports_dir)
}

# --- 3.3: Process Bracken Abundance Files ---
if (!args$no_bracken) {
  unaligned_bracken_data <- process_bracken_folder(unaligned_bracken_dir)
  nonhuman_bracken_data  <- process_bracken_folder(nonhuman_bracken_dir)
} else {
  cat("Skipping Bracken processing as per --no-bracken flag.\n")
  unaligned_bracken_data <- tibble()
  nonhuman_bracken_data  <- tibble()
}

# --- 3.4: Process Metadata, Runtime, and Read Statistics ---
sample_metadata <- parse_sample_metadata(args$metadata_dir)
runtime_data    <- parse_runtime_files(args$runtime_dir)
rrstats_data    <- parse_rrstats_files(args$rrstats_dir)


# ==============================================================================
# SECTION 4: DATA FILTERING & PRE-PROCESSING
# ==============================================================================
cat("\n\n=======================================\n")
cat("  SECTION 4: DATA FILTERING & PRE-PROCESSING\n")
cat("=======================================\n\n")

# --- 4.1: Apply Filters to Kraken2 Data ---
unaligned_kreport_filtered <- apply_filters(unaligned_kreport_data, args$min_reads, args$exclude_taxids, args$minimizer_ratio, args$minimizer_threshold)
nonhuman_kreport_filtered  <- apply_filters(nonhuman_kreport_data, args$min_reads, args$exclude_taxids, args$minimizer_ratio, args$minimizer_threshold)

# --- 4.2: Reshape K-report Data to Wide Format ---
if (nrow(unaligned_kreport_filtered) > 0) {
  unaligned_wide <- unaligned_kreport_filtered %>%
    pivot_wider(
      id_cols = c(taxID, name),
      names_from = sample,
      values_from = c(cladeReads, taxonReads),
      names_sep = "_"
    )
} else {
  unaligned_wide <- tibble(taxID = integer(), name = character())
}

if (nrow(nonhuman_kreport_filtered) > 0) {
  nonhuman_wide <- nonhuman_kreport_filtered %>%
    pivot_wider(
      id_cols = c(taxID, name),
      names_from = sample,
      values_from = c(cladeReads, taxonReads),
      names_sep = "_"
    )
} else {
  nonhuman_wide <- tibble(taxID = integer(), name = character())
}

# --- 4.3: Reshape and Merge Bracken Data ---
if (nrow(unaligned_bracken_data) > 0) {
  unaligned_bracken_wide <- unaligned_bracken_data %>%
    pivot_wider(
      id_cols = c(taxID, name),
      names_from = sample,
      values_from = c(bracken_reads, bracken_fraction_total),
      names_sep = "_"
    )
  unaligned_merged <- full_join(unaligned_wide, unaligned_bracken_wide, by = c("taxID", "name"))
  cat("Merged Unaligned k-report data with Bracken data.\n")
} else {
  unaligned_merged <- unaligned_wide
  cat("No Unaligned Bracken data to merge.\n")
}

if (nrow(nonhuman_bracken_data) > 0) {
  nonhuman_bracken_wide <- nonhuman_bracken_data %>%
    pivot_wider(
      id_cols = c(taxID, name),
      names_from = sample,
      values_from = c(bracken_reads, bracken_fraction_total),
      names_sep = "_"
    )
  nonhuman_merged <- full_join(nonhuman_wide, nonhuman_bracken_wide, by = c("taxID", "name"))
  cat("Merged Non-human k-report data with Bracken data.\n")
} else {
  nonhuman_merged <- nonhuman_wide
  cat("No Non-human Bracken data to merge.\n")
}

# --- 4.4: Create Master Metadata Table ---
metadata_list <- list(sample_metadata, runtime_data, rrstats_data)
metadata_list <- metadata_list[sapply(metadata_list, function(df) !is.null(df) && nrow(df) > 0)]

if (length(metadata_list) > 0) {
  master_metadata <- reduce(metadata_list, full_join, by = "sample")
  cat("Created a master metadata table with", nrow(master_metadata), "samples and", ncol(master_metadata), "variables.\n")
} else {
  master_metadata <- tibble(sample = character())
  cat("No metadata files were found or processed. Proceeding without metadata.\n")
}

# --- 4.5: Create Final Long-Format Dataframes for Analysis ---
if (ncol(unaligned_merged) > 2) {
  unaligned_long <- unaligned_merged %>%
    pivot_longer(
      cols = -c(taxID, name),
      names_to = c(".value", "sample"),
      names_pattern = "([a-zA-Z_]+)_(.+)",
      values_drop_na = TRUE
    ) %>%
    left_join(master_metadata, by = "sample")
} else {
  unaligned_long <- tibble()
}

if (ncol(nonhuman_merged) > 2) {
  nonhuman_long <- nonhuman_merged %>%
    pivot_longer(
      cols = -c(taxID, name),
      names_to = c(".value", "sample"),
      names_pattern = "([a-zA-Z_]+)_(.+)",
      values_drop_na = TRUE
    ) %>%
    left_join(master_metadata, by = "sample")
} else {
  nonhuman_long <- tibble()
}

cat("Created final long-format dataframes for analysis.\n")
cat("  - Unaligned dataset:", nrow(unaligned_long), "rows\n")
cat("  - Non-human dataset:", nrow(nonhuman_long), "rows\n")


# ==============================================================================
# SECTION 5: DATA ANALYSIS
# ==============================================================================
cat("\n\n=======================================\n")
cat("  SECTION 5: DATA ANALYSIS\n")
cat("=======================================\n\n")

analysis_results <- list()

#' Run Correlation Analysis Between Two Datasets
#'
#' @param data1 Wide-format dataframe for the first dataset.
#' @param data2 Wide-format dataframe for the second dataset.
#' @param name1 Name of the first dataset (for plot titles).
#' @param name2 Name of the second dataset (for plot titles).
#' @param measure1 The column suffix to correlate from data1.
#' @param measure2 The column suffix to correlate from data2.
#' @return A list containing a summary tibble of correlations and a list of ggplot objects.
run_correlation_analysis <- function(data1, data2, name1, name2, measure1, measure2) {
  cat(paste0("\nRunning correlation analysis: ", name1, " vs ", name2, " on '", measure1, "' vs '", measure2, "'\n"))
  
  cols1 <- grep(paste0("^", measure1, "_"), colnames(data1), value = TRUE)
  cols2 <- grep(paste0("^", measure2, "_"), colnames(data2), value = TRUE)
  
  samples1 <- gsub(paste0("^", measure1, "_"), "", cols1)
  samples2 <- gsub(paste0("^", measure2, "_"), "", cols2)
  
  common_samples <- intersect(samples1, samples2)
  
  if (length(common_samples) == 0) {
    warning("No common samples found for correlation between ", name1, " and ", name2)
    return(list(summary = tibble(), plots = list()))
  }
  cat("Found", length(common_samples), "common samples for correlation.\n")
  
  results <- map(common_samples, ~{
    sample_id <- .x
    col1 <- paste0(measure1, "_", sample_id)
    col2 <- paste0(measure2, "_", sample_id)
    
    comparison_df <- inner_join(
      data1 %>% select(taxID, name, reads1 = all_of(col1)),
      data2 %>% select(taxID, name, reads2 = all_of(col2)),
      by = c("taxID", "name")
    ) %>%
    filter(reads1 > 0 | reads2 > 0)
    
    if (nrow(comparison_df) < 2) return(NULL)
    
    cor_test <- cor.test(~ log10(reads1 + 1) + log10(reads2 + 1), data = comparison_df)
    
    plot_df <- comparison_df %>%
      mutate(
        log_reads1 = log10(reads1 + 1),
        log_reads2 = log10(reads2 + 1),
        label = if_else(
          abs(log_reads1 - log_reads2) > 2 | log_reads1 > 4 | log_reads2 > 4,
          name, 
          NA_character_
        )
      )
      
    p <- ggplot(plot_df, aes(x = log_reads1, y = log_reads2)) +
      geom_point(alpha = 0.6, color = "blue") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      geom_text_repel(aes(label = label), na.rm = TRUE, size = 2.5, max.overlaps = 15) +
      labs(
        title = paste0("Correlation: ", name1, " vs ", name2, " (Sample: ", sample_id, ")"),
        subtitle = sprintf("Pearson's r = %.3f, p-value = %.2e, n = %d",
                           cor_test$estimate, cor_test$p.value, nrow(plot_df)),
        x = paste("log10(", name1, " ", measure1, " + 1)"),
        y = paste("log10(", name2, " ", measure2, " + 1)")
      ) +
      my_ggplot_theme +
      coord_fixed()
      
    list(
      summary = tibble(
        sample = sample_id,
        pearson_r = cor_test$estimate,
        p_value = cor_test$p.value,
        n_species = nrow(plot_df)
      ),
      plot = p
    )
  })
  
  results <- results[!sapply(results, is.null)]
  all_summaries <- map_dfr(results, "summary")
  all_plots <- map(results, "plot")
  names(all_plots) <- all_summaries$sample
  
  cat("Correlation analysis complete.\n")
  return(list(summary = all_summaries, plots = all_plots))
}

#' Calculate Species Frequency in Top N
#'
#' @param data A wide-format dataframe with species counts per sample.
#' @param top_n The number of top species to consider from each sample.
#' @param measure The column prefix to use for ranking (e.g., "cladeReads").
#' @return A tibble with species and their frequency of appearing in the top N.
calculate_species_frequency <- function(data, top_n, measure = "cladeReads") {
  cat("\nCalculating species frequency in Top", top_n, "across all samples...\n")
  
  count_cols <- grep(paste0("^", measure, "_"), colnames(data), value = TRUE)
  if (length(count_cols) == 0) {
    warning("No count columns found with measure: ", measure)
    return(tibble())
  }
  
  top_species_per_sample <- map_dfr(count_cols, ~{
    sample_col <- .x
    data %>%
      select(taxID, name, count = all_of(sample_col)) %>%
      filter(!is.na(count) & count > 0) %>%
      arrange(desc(count)) %>%
      slice_head(n = top_n) %>%
      select(taxID, name)
  })
  
  frequency_table <- top_species_per_sample %>%
    count(taxID, name, name = "topN_freq") %>%
    arrange(desc(topN_freq))
    
  cat("Found", nrow(frequency_table), "unique species in the Top", top_n, "of at least one sample.\n")
  return(frequency_table)
}

# --- 5.1: Correlation Analysis ---
if (!args$no_corr) {
  analysis_results$unaligned_vs_nonhuman_corr <- run_correlation_analysis(
    unaligned_merged, nonhuman_merged, "Unaligned", "Nonhuman", "cladeReads", "cladeReads"
  )
  
  if (nrow(unaligned_bracken_data) > 0) {
    analysis_results$kraken_vs_bracken_unaligned_corr <- run_correlation_analysis(
      unaligned_merged, unaligned_merged, "Kraken", "Bracken", "cladeReads", "bracken_reads"
    )
  }
  
  if (nrow(nonhuman_bracken_data) > 0) {
    analysis_results$kraken_vs_bracken_nonhuman_corr <- run_correlation_analysis(
      nonhuman_merged, nonhuman_merged, "Kraken", "Bracken", "cladeReads", "bracken_reads"
    )
  }
} else {
  cat("Skipping correlation analysis as per --no-corr flag.\n")
}

# --- 5.2: Species Frequency Analysis ---
all_merged_data <- full_join(
  unaligned_merged,
  nonhuman_merged,
  by = c("taxID", "name"),
  suffix = c("_unaligned", "_nonhuman")
)

analysis_results$species_frequency <- calculate_species_frequency(
  all_merged_data, args$top_n, "cladeReads"
)

# --- 5.3: Generate a list of unique species for annotation ---
unique_species_for_annotation <- all_merged_data %>%
  select(taxID, name) %>%
  filter(!is.na(taxID) & !is.na(name)) %>%
  distinct(taxID, .keep_all = TRUE)

cat("\nGenerated a list of", nrow(unique_species_for_annotation), "unique species for annotation.\n")


# ==============================================================================
# SECTION 6: SPECIES ANNOTATION
# ==============================================================================
cat("\n\n=======================================\n")
cat("  SECTION 6: SPECIES ANNOTATION\n")
cat("=======================================\n\n")

#' Annotate a dataframe of species with risk group and HOMD information.
#'
#' @param species_df A dataframe with `taxID` and `name` columns.
#' @param db_dir Path to the directory containing annotation database files.
#' @return An annotated tibble with additional columns for risk group and HOMD info.
annotate_species <- function(species_df, db_dir) {
  cat("Starting species annotation...\n")
  
  epathogen_db_path <- file.path(db_dir, "epathogen-2025-07-15-result.csv")
  homd_db_path <- file.path(db_dir, "HOMD_taxon_table2025-07-15_1752607154.xls")
  
  if (!file.exists(epathogen_db_path) || !file.exists(homd_db_path)) {
    stop("Annotation database files not found in: ", db_dir)
  }
  
  epathogen_db <- read.csv(epathogen_db_path, header = TRUE, stringsAsFactors = FALSE)
  homd_db <- read.table(homd_db_path, header = TRUE, sep = "\t", fill = TRUE, quote = "", stringsAsFactors = FALSE)
  
  cat("Annotation databases loaded successfully.\n")
  
  annotated_df <- species_df %>%
    left_join(
      epathogen_db %>% select(taxID, RiskGroup = Human.classification),
      by = "taxID"
    ) %>%
    mutate(RiskGroup = coalesce(RiskGroup, "NotAnnotated"))
  
  cat("  - Risk Group annotation complete.\n")
  
  homd_lookup <- homd_db %>%
    select(taxID = NCBI.Taxon.ID, HOMD_body_site = Body.Site.s.) %>%
    mutate(taxID = as.integer(taxID)) %>%
    filter(!is.na(taxID)) %>%
    distinct(taxID, .keep_all = TRUE)
    
  annotated_df <- annotated_df %>%
    left_join(homd_lookup, by = "taxID") %>%
    mutate(
      HOMD.Category = case_when(
        is.na(HOMD_body_site) | HOMD_body_site == "" ~ "NotAnnotated",
        grepl("opportunistic pathogen", HOMD_body_site, ignore.case = TRUE) ~ "Opportunist",
        grepl("pathogen", HOMD_body_site, ignore.case = TRUE) ~ "Pathogen",
        TRUE ~ "Microbiome"
      ),
      HOMD_body_site = coalesce(HOMD_body_site, "NotAnnotated")
    )
    
  cat("  - HOMD annotation complete.\n")
  
  annotated_df <- annotated_df %>%
    mutate(
      kingdom = case_when(
        grepl("bacterium|bacterial", name, ignore.case = TRUE) ~ "Bacteria",
        grepl("virus", name, ignore.case = TRUE) ~ "Viruses",
        grepl("fungus|fungi", name, ignore.case = TRUE) ~ "Fungi",
        TRUE ~ "Eukaryota"
      )
    )
  cat("  - Kingdom annotation complete.\n")
  
  return(annotated_df)
}

annotated_species_data <- annotate_species(unique_species_for_annotation, args$db_dir)

# Define expected columns for the final annotated dataframes
# This ensures that even if no data is processed, downstream functions
# that expect these columns will not fail.
expected_cols <- unique(c(
  names(master_metadata),
  "taxID", "name", "cladeReads", "taxonReads", "bracken_reads", "bracken_fraction_total",
  names(annotated_species_data)
))

# Create empty tibbles with the correct structure to avoid errors
unaligned_long_annotated <- tibble::tibble(!!!setNames(rep(list(logical()), length(expected_cols)), expected_cols))
nonhuman_long_annotated <- tibble::tibble(!!!setNames(rep(list(logical()), length(expected_cols)), expected_cols))


if (nrow(unaligned_long) > 0) {
  unaligned_long_annotated <- unaligned_long %>%
    left_join(annotated_species_data, by = c("taxID", "name"))
}

if (nrow(nonhuman_long) > 0) {
  nonhuman_long_annotated <- nonhuman_long %>%
    left_join(annotated_species_data, by = c("taxID", "name"))
}


cat("Annotation data merged back into main dataframes.\n")


# ==============================================================================
# SECTION 7: PLOTTING
# ==============================================================================
cat("\n\n=======================================\n")
cat("  SECTION 7: PLOTTING\n")
cat("=======================================\n\n")

plots <- list()

plots$risk_group_summary <- create_annotation_bar_chart(
  annotated_species_data, "RiskGroup", "Species Count by Risk Group"
)
plots$homd_summary <- create_annotation_bar_chart(
  annotated_species_data, "HOMD.Category", "Species Count by HOMD Category"
)
plots$kingdom_summary <- create_annotation_bar_chart(
  annotated_species_data, "kingdom", "Species Count by Kingdom"
)
cat("Generated annotation summary plots.\n")

all_long_annotated <- bind_rows(
  unaligned_long_annotated %>% mutate(dataset = "Unaligned"),
  nonhuman_long_annotated %>% mutate(dataset = "Non-human")
)

plots$top_species_by_cladeReads <- create_top_species_plot(
  all_long_annotated, "cladeReads", args$top_n, "RiskGroup", "Top Species by Mean Clade Reads"
)
if ("bracken_reads" %in% names(all_long_annotated)) {
  plots$top_species_by_bracken_reads <- create_top_species_plot(
    all_long_annotated, "bracken_reads", args$top_n, "RiskGroup", "Top Species by Mean Bracken Reads"
  )
}
cat("Generated top species plots.\n")

if (nrow(master_metadata) > 0 && "num_input_reads" %in% names(master_metadata)) {
  plots$input_reads_vs_mapped <- ggplot(master_metadata, aes(x = num_input_reads, y = uniquely_mapped_percent)) +
    geom_point(alpha = 0.7, color = "purple") +
    geom_smooth(method = "lm", se = FALSE, color = "orange") +
    scale_x_log10() +
    labs(
      title = "Input Reads vs. Uniquely Mapped Percentage",
      x = "Number of Input Reads (log10 scale)",
      y = "Uniquely Mapped (%)"
    ) +
    my_ggplot_theme
  cat("Generated metadata plots.\n")
}


# ==============================================================================
# SECTION 8: DATA & REPORT SAVING
# ==============================================================================
cat("\n\n=======================================\n")
cat("  SECTION 8: DATA & REPORT SAVING\n")
cat("=======================================\n\n")

save_data(unaligned_merged, "unaligned_merged_wide", args$run_name, output_dir, args$output_prefix)
save_data(nonhuman_merged, "nonhuman_merged_wide", args$run_name, output_dir, args$output_prefix)
save_data(master_metadata, "master_metadata", args$run_name, output_dir, args$output_prefix)
save_data(annotated_species_data, "annotated_species_list", args$run_name, output_dir, args$output_prefix)
save_data(analysis_results$species_frequency, "species_topN_frequency", args$run_name, output_dir, args$output_prefix)

if (!args$no_corr) {
  if (exists("unaligned_vs_nonhuman_corr", where = analysis_results)) {
    save_data(analysis_results$unaligned_vs_nonhuman_corr$summary, "correlation_summary_unaligned_vs_nonhuman", args$run_name, output_dir, args$output_prefix)
    save_plots_to_pdf(analysis_results$unaligned_vs_nonhuman_corr$plots, "correlation_plots_unaligned_vs_nonhuman.pdf", reports_dir, args$output_prefix)
  }
  if (exists("kraken_vs_bracken_unaligned_corr", where = analysis_results)) {
    save_data(analysis_results$kraken_vs_bracken_unaligned_corr$summary, "correlation_summary_kraken_vs_bracken_unaligned", args$run_name, output_dir, args$output_prefix)
    save_plots_to_pdf(analysis_results$kraken_vs_bracken_unaligned_corr$plots, "correlation_plots_kraken_vs_bracken_unaligned.pdf", reports_dir, args$output_prefix)
  }
  if (exists("kraken_vs_bracken_nonhuman_corr", where = analysis_results)) {
    save_data(analysis_results$kraken_vs_bracken_nonhuman_corr$summary, "correlation_summary_kraken_vs_bracken_nonhuman", args$run_name, output_dir, args$output_prefix)
    save_plots_to_pdf(analysis_results$kraken_vs_bracken_nonhuman_corr$plots, "correlation_plots_kraken_vs_bracken_nonhuman.pdf", reports_dir, args$output_prefix)
  }
}

save_plots_to_pdf(plots, "main_analysis_plots.pdf", reports_dir, args$output_prefix)

cat("\nAll data and reports have been saved.\n")


# ==============================================================================
# SECTION 9: BATCH REPORT GENERATION
# ==============================================================================
if (args$generate_report) {
  cat("\n\n=======================================\n")
  cat("  SECTION 9: BATCH REPORT GENERATION\n")
  cat("=======================================\n\n")

  generate_batch_report_pdf <- function(config, master_metadata, annotated_species, plots, correlation_results, output_dir, output_prefix = NULL) {
    cat("\nGenerating final PDF batch report...\n")
    
    final_prefix <- if (!is.null(output_prefix)) paste0(output_prefix, "_") else ""
    report_path <- file.path(output_dir, paste0(final_prefix, "batch_report_", config$run_name, "_", format(Sys.time(), "%Y%m%d"), ".pdf"))
    
    pdf(report_path, width = 11, height = 8.5, onefile = TRUE)
    
    # Page 1: Title Page & Summary
    plot.new()
    grid.text(
      paste("SPARKEN Batch Report:", config$run_name),
      x = 0.5, y = 0.9,
      gp = gpar(fontsize = 20, fontface = "bold")
    )
    
    summary_text <- paste(
      "Part 1: Processing Summary",
      paste("  Run Name:", config$run_name),
      paste("  Processing Date:", format(Sys.time(), "%Y-%m-%d")),
      paste("  Input Directory:", config$input_dir),
      paste("  Total Samples Processed:", if(is.null(master_metadata)) "N/A" else nrow(master_metadata)),
      paste("  Unique Species Identified:", if(is.null(annotated_species)) "N/A" else nrow(annotated_species)),
      "",
      "Part 2: Analysis Configuration",
      paste("  Min. Clade Reads:", config$min_reads),
      paste("  Top N Species Analyzed:", config$top_n),
      paste("  Bracken Processing:", !config$no_bracken),
      paste("  Correlation Analysis:", !config$no_corr),
      paste("  Subspecies Included:", config$subspecies),
      sep = "\n"
    )
    grid.text(summary_text, x = 0.1, y = 0.7, just = "left", gp = gpar(fontsize = 10, fontfamily = "mono"))
    
    # Page 2: Read Statistics Summary
    if (!is.null(master_metadata) && "num_input_reads" %in% names(master_metadata)) {
      plot.new()
      grid.text("Read Statistics Summary", y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
      
      stats_summary_df <- master_metadata %>%
        select(where(is.numeric)) %>%
        select(any_of(c("num_input_reads", "avg_input_read_length", "uniquely_mapped_reads", "percent_mapped_reads", "uniquely_mapped_percent"))) %>%
        summary() %>%
        as.data.frame.matrix() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Statistic")

      if(nrow(stats_summary_df) > 0) {
          grid.table(stats_summary_df, rows = NULL, theme = ttheme_minimal(
              core=list(bg_params = list(fill = brewer.pal(4, "Blues")[1:2], col=NA)),
              colhead=list(fg_params=list(col="navyblue"))
          ))
      }
    }
    
    # Page 3: Runtime Summary
    if (!is.null(master_metadata) && any(grepl("_rt$", names(master_metadata)))) {
        plot.new()
        grid.text("Processing Runtime Summary (seconds)", y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
        
        runtime_summary_df <- master_metadata %>%
          select(ends_with("_rt")) %>%
          summary() %>%
          as.data.frame.matrix() %>%
          t() %>%
          as.data.frame() %>%
          tibble::rownames_to_column("Statistic")
          
        if(nrow(runtime_summary_df) > 0) {
          grid.table(runtime_summary_df, rows = NULL, theme = ttheme_minimal(
              core=list(bg_params = list(fill = brewer.pal(4, "Blues")[1:2], col=NA)),
              colhead=list(fg_params=list(col="navyblue"))
          ))
        }
    }

    # Subsequent Pages: Plots
    for (plot_name in names(plots)) {
      if(!is.null(plots[[plot_name]])) print(plots[[plot_name]])
    }
    
    # Correlation plots
    if (!is.null(correlation_results)) {
      # Filter for actual correlation results which are lists containing a 'plots' list
      corr_analyses_to_plot <- Filter(function(x) is.list(x) && "plots" %in% names(x), correlation_results)

      for (corr_analysis_name in names(corr_analyses_to_plot)) {
        corr_analysis <- corr_analyses_to_plot[[corr_analysis_name]]
        if(length(corr_analysis$plots) > 0) {
          plot.new()
          friendly_name <- gsub("_corr", "", corr_analysis_name) %>% gsub("_", " ", .) %>% stringr::str_to_title()
          grid.text(paste("Correlation Plots:", friendly_name), y = 0.9, gp = gpar(fontsize = 16, fontface = "bold"))
          
          if(nrow(corr_analysis$summary) > 0) {
              grid.table(head(corr_analysis$summary), rows = NULL, theme = ttheme_minimal())
          }

          for (p in corr_analysis$plots) {
            print(p)
          }
        }
      }
    }
    
    dev.off()
    cat("Final batch report saved to:", basename(report_path), "\n")
  }

  generate_batch_report_pdf(
    config = args,
    master_metadata = master_metadata,
    annotated_species = annotated_species_data,
    plots = plots,
    correlation_results = analysis_results,
    output_dir = reports_dir,
    output_prefix = args$output_prefix
  )
}

# Stop the cluster
if (cores > 1) {
  stopCluster(cl)
  cat("\nParallel cluster stopped.\n")
}

main_end_time <- Sys.time()
total_runtime <- difftime(main_end_time, main_start_time, units = "mins")
cat("\n\n=======================================\n")
cat("  Total Pipeline Runtime:", round(total_runtime, 2), "minutes\n")
cat("  Analysis complete for run:", args$run_name, "\n")
cat("  Outputs are located in:", output_dir, "\n")
cat("=======================================\n")
