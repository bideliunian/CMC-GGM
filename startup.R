# CMC-GGM Project Startup Script
# Source this file to set up the complete environment

cat("Initializing CMC-GGM project environment...\n")

# 1. Load required packages
required_packages <- c("Matrix", "igraph", "huge", "gglasso", "doParallel", "mvtnorm", "dplyr", "ggplot2")

cat("Loading required packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Error: Package", pkg, "not found. Please run setup_environment_minimal.R first.\n"))
  }
}

# 2. Set up paths
cat("Setting up project paths...\n")
source("config_paths.R")

# 3. Load all functions
cat("Loading all utility functions...\n")
source_all_functions()

cat("\n=== CMC-GGM Environment Ready ===\n")
cat("Available functions loaded from Functions/ directory\n")
cat("Path variables configured (e.g., model_a_proj_path)\n")
cat("Parallel processing available with", parallel::detectCores(), "cores\n")

cat("\n=== Quick Start Examples ===\n")
cat("# Generate synthetic data:\n")
cat("data <- gen_data(n=100, p=20, m=3, model='model1', run.ind=1, trans.type='copula')\n\n")

cat("# Run a basic model:\n")
cat("setwd(model_a_proj_path)\n")
cat("# Then run any model script\n\n")

cat("# Set up parallel processing:\n")
cat("cl <- makeCluster(4)\n")
cat("registerDoParallel(cl)\n")
cat("# Your parallel code here\n")
cat("stopCluster(cl)\n\n")

cat("Environment initialization complete!\n") 