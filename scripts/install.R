## Non-interactive install script for cluster execution
## Creates an R project file, initializes renv, and installs required packages

## Ensure a CRAN mirror is set for non-interactive installs
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Configure a non-interactive user library path
user_lib <- Sys.getenv("R_LIBS_USER")
if (identical(user_lib, "")) {
  user_lib <- file.path(Sys.getenv("HOME"), "R", paste0(R.version$platform, "-", getRversion()))
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(user_lib)
cat("Using R library path:", user_lib, "\n")

# Ensure renv is installed and load it
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", dependencies = TRUE)
}
library(renv)

# Get the current working directory (should be Ranalysis/)
project_dir <- getwd()
cat("Setting up R project in:", project_dir, "\n")

# Create .Rproj file for the project if it doesn't exist
rproj_file <- file.path(project_dir, "RNASeqAnalysis.Rproj")
if (!file.exists(rproj_file)) {
  rproj_content <- paste0(
    "Version: 1.0\n\n",
    "RestoreWorkspace: Default\n",
    "SaveWorkspace: Default\n",
    "AlwaysSaveHistory: Default\n\n",
    "EnableCodeIndexing: Yes\n",
    "UseSpacesForTab: Yes\n",
    "NumSpacesForTab: 2\n",
    "Encoding: UTF-8\n\n",
    "RnwWeave: Sweave\n",
    "LaTeX: pdfLaTeX\n\n",
    "AutoAppendNewline: Yes\n",
    "StripTrailingWhitespace: Yes\n\n",
    "BuildType: Package\n",
    "PackageUseDevtools: Yes\n",
    "PackageInstallArgs: --no-multiarch --with-keep.source\n"
  )
  writeLines(rproj_content, rproj_file)
  cat("Created R project file:", rproj_file, "\n")
} else {
  cat("R project file already exists:", rproj_file, "\n")
}

# Provide consent to renv operations to avoid prompting
tryCatch({
  if ("consent" %in% ls("package:renv")) {
    renv::consent(provided = TRUE)
  }
}, error = function(e) {
  message("renv::consent not available or failed: ", e$message)
})

# Initialize renv in non-interactive mode if not already initialized
if (!dir.exists(file.path(project_dir, "renv"))) {
  cat("Initializing renv (non-interactive)...\n")
  tryCatch({
    renv::init(bare = TRUE, restart = FALSE)
  }, error = function(e) {
    message("renv::init error (continuing): ", e$message)
  })
} else {
  cat("renv already present, skipping init.\n")
}

# Clean/ensure consistent renv state (best-effort)
cat("Cleaning renv state (best-effort)...\n")
tryCatch({
  renv::clean()
}, error = function(e) {
  message("renv::clean error: ", e$message)
})

cat("Checking renv status...\n")
tryCatch({
  renv::status()
}, error = function(e) {
  message("renv::status error: ", e$message)
})

# Packages to install
packages_to_install <- c(
  "here",
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "tibble",
  "ggrepel",
  "kableExtra",
  "gridExtra",
  "openxlsx",
  "foreach",
  "doParallel",
  "httr2"  # Added for web requests in update_databases
)

# Function to check if package is installed in user library
is_package_installed <- function(pkg, lib_path = user_lib) {
  installed_pkgs <- installed.packages(lib.loc = lib_path)
  return(pkg %in% rownames(installed_pkgs))
}

cat("Installing required packages into:", user_lib, "\n")
for (pkg in packages_to_install) {
  if (!is_package_installed(pkg, user_lib)) {
    cat("Installing package:", pkg, "\n")
    tryCatch({
      install.packages(pkg, dependencies = TRUE, lib = user_lib)
    }, error = function(e) {
      message("Failed to install ", pkg, ": ", e$message)
    })
  } else {
    cat(pkg, "already installed in user library, skipping.\n")
  }
}

# Snapshot and restore renv environment
cat("Creating renv snapshot...\n")
tryCatch({
  renv::snapshot()
}, error = function(e) {
  message("renv::snapshot error: ", e$message)
})

cat("Restoring renv environment...\n")
tryCatch({
  renv::restore()
}, error = function(e) {
  message("renv::restore error: ", e$message)
})

cat("Running renv diagnostics...\n")
tryCatch({
  renv::diagnostics(project = NULL)
}, error = function(e) {
  message("renv::diagnostics error: ", e$message)
})

cat("R project setup complete!\n")



