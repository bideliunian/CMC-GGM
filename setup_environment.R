#!/usr/bin/env Rscript

# CMC-GGM Project Environment Setup Script
# This script installs all required R packages for the graph estimation project
# Combines comprehensive package list with robust error handling

cat("Setting up R environment for CMC-GGM project...\n\n")

# Enhanced function to install packages with better error handling
install_if_missing <- function(packages, category = "packages") {
  failed_packages <- c()
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("Installing", category, "package:", pkg, "\n"))
      tryCatch({
        install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org/")
        library(pkg, character.only = TRUE, quietly = TRUE)
        cat(paste("✓ Successfully installed", pkg, "\n"))
      }, error = function(e) {
        cat(paste("✗ Failed to install", pkg, ":", e$message, "\n"))
        failed_packages <<- c(failed_packages, pkg)
      })
    } else {
      cat(paste("✓ Package", pkg, "already installed\n"))
    }
  }
  
  return(failed_packages)
}

# Function to install Bioconductor packages
install_bioc_if_missing <- function(packages) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  failed_packages <- c()
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("Installing Bioconductor package:", pkg, "\n"))
      tryCatch({
        BiocManager::install(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE, quietly = TRUE)
        cat(paste("✓ Successfully installed", pkg, "\n"))
      }, error = function(e) {
        cat(paste("✗ Failed to install", pkg, ":", e$message, "\n"))
        failed_packages <<- c(failed_packages, pkg)
      })
    } else {
      cat(paste("✓ Bioconductor package", pkg, "already installed\n"))
    }
  }
  
  return(failed_packages)
}

# Core packages that are essential and stable
core_packages <- c(
  # Matrix operations and linear algebra
  "Matrix",
  "matrixcalc",
  "pracma",
  
  # Graph theory and network analysis
  "igraph",
  "huge",
  
  # Statistical modeling and optimization
  "BB",
  "mvtnorm",
  "gglasso",
  
  # Data manipulation and utilities
  "dplyr",
  "tidyr",
  "clue",
  "pdist",
  
  # Parallel computing
  "doParallel",
  "foreach",
  
  # Random number generation
  "randtoolbox",
  
  # Visualization
  "ggplot2"
)

# Optional packages with potential dependency issues
optional_packages <- c(
  # Functional data analysis
  "fda",
  "fdapace"
)

# Packages that may have system dependencies (install separately)
system_dependent_packages <- c(
  # EEG analysis (requires X11 on some systems)
  "eegkit",
  "eegkitdata"
)

# Install core packages (essential)
cat("=== Installing Core Packages ===\n")
failed_core <- install_if_missing(core_packages, "core")

# Install optional packages (nice to have)
cat("\n=== Installing Optional Packages ===\n")
failed_optional <- install_if_missing(optional_packages, "optional")

# Install system-dependent packages (may fail on some systems)
cat("\n=== Installing System-Dependent Packages ===\n")
cat("Note: These packages may fail if system dependencies (like X11) are missing\n")
failed_system <- install_if_missing(system_dependent_packages, "system-dependent")

# Bioconductor packages (if any needed in the future)
# bioc_packages <- c()
# if(length(bioc_packages) > 0) {
#   cat("\n=== Installing Bioconductor Packages ===\n")
#   failed_bioc <- install_bioc_if_missing(bioc_packages)
# }

cat("\n=== Environment Setup Complete ===\n")

# Verify installation by loading key packages
cat("\nVerifying installation of key packages...\n")
key_packages <- c("Matrix", "igraph", "huge", "gglasso", "doParallel", "mvtnorm")

all_loaded <- TRUE
for (pkg in key_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("✓", pkg, "loaded successfully\n"))
  } else {
    cat(paste("✗ Error loading", pkg, "\n"))
    all_loaded <- FALSE
  }
}

# Display R and package version information
cat("\n=== System Information ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")

# Check if parallel processing is available
cat("Available cores for parallel processing:", parallel::detectCores(), "\n")

# Summary of installation results
cat("\n=== Installation Summary ===\n")
if (all_loaded) {
  cat("✓ All essential packages loaded successfully!\n")
} else {
  cat("⚠ Some essential packages failed to load. Check error messages above.\n")
}

# Report failed packages
all_failed <- c(failed_core, failed_optional, failed_system)
if (length(all_failed) > 0) {
  cat("\n⚠ Failed to install the following packages:\n")
  for (pkg in all_failed) {
    cat(paste("  -", pkg, "\n"))
  }
  
  cat("\nTroubleshooting tips:\n")
  if ("eegkit" %in% all_failed || "eegkitdata" %in% all_failed) {
    cat("• EEG packages: Install X11 (XQuartz on macOS) or skip EEG analysis\n")
  }
  if ("fdapace" %in% all_failed) {
    cat("• fdapace: May have complex dependencies, try installing manually\n")
  }
  cat("• For other packages: Check internet connection and R version compatibility\n")
} else {
  cat("✓ All packages installed successfully!\n")
}

cat("\n=== Usage Notes ===\n")
cat("• Core packages are essential for graph estimation\n")
cat("• Optional packages add functionality but are not required\n")
cat("• System-dependent packages may not work on all systems\n")
cat("• You can still run most experiments even if some packages failed\n")

cat("\n=== Next Steps ===\n")
cat("1. Run: source('config_paths.R')\n")
cat("2. Run: source_all_functions()\n")
cat("3. Navigate to Model_* directories to run experiments\n")
cat("4. For complete setup: source('startup.R')\n")

cat("\nEnvironment setup completed!\n") 