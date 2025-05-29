# CMC-GGM Project Path Configuration
# This script helps configure the paths used throughout the project

# Get the current working directory (should be the project root)
project_root <- getwd()

cat("Setting up project paths...\n")
cat("Project root:", project_root, "\n")

# Define all the paths used in the project
paths <- list(
  # Main project directories
  project_root = project_root,
  functions = file.path(project_root, "Functions"),
  
  # Model directories
  model_a = file.path(project_root, "Model_a"),
  model_a_proj = file.path(project_root, "Model_a_proj"),
  model_a_nblasso = file.path(project_root, "Model_a_nblasso"),
  model_a_rank = file.path(project_root, "Model_a_rank"),
  
  model_b = file.path(project_root, "Model_b"),
  model_b_nblasso = file.path(project_root, "Model_b_nblasso"),
  model_b_rank = file.path(project_root, "Model_b_rank"),
  
  model_c = file.path(project_root, "Model_c"),
  model_c_nblasso = file.path(project_root, "Model_c_nblasso"),
  model_c_rank = file.path(project_root, "Model_c_rank"),
  
  # Data and analysis directories
  eeg = file.path(project_root, "EEG"),
  gene_protein = file.path(project_root, "Gene_Protein"),
  fgm = file.path(project_root, "fgm"),
  
  # Output directories
  plot = file.path(project_root, "Plot"),
  textures = file.path(project_root, "Textures"),
  
  # Python projection directories
  proj_cmc = file.path(project_root, "proj_cmc"),
  proj_cmc_k2 = file.path(project_root, "proj_cmc_k2")
)

# Create Results directories if they don't exist
for (model_dir in c("Model_a", "Model_a_proj", "Model_a_nblasso", "Model_a_rank",
                    "Model_b", "Model_b_nblasso", "Model_b_rank",
                    "Model_c", "Model_c_nblasso", "Model_c_rank")) {
  results_dir <- file.path(project_root, model_dir, "Results")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
    cat("Created directory:", results_dir, "\n")
  }
}

# Function to source all function files
source_all_functions <- function() {
  function_files <- list.files(paths$functions, 
                              pattern = "*.R$", 
                              full.names = TRUE, 
                              ignore.case = TRUE)
  
  cat("Sourcing function files:\n")
  for (file in function_files) {
    cat("  -", basename(file), "\n")
    source(file, .GlobalEnv)
  }
  cat("All functions loaded successfully!\n")
}

# Export paths to global environment for easy access
for (name in names(paths)) {
  assign(paste0(name, "_path"), paths[[name]], envir = .GlobalEnv)
}

cat("\n=== Available Path Variables ===\n")
for (name in names(paths)) {
  cat(paste0(name, "_path"), "=", paths[[name]], "\n")
}

cat("\n=== Usage ===\n")
cat("1. Run source_all_functions() to load all utility functions\n")
cat("2. Use the *_path variables in your scripts instead of hardcoded paths\n")
cat("3. Example: setwd(model_a_proj_path)\n")

# Verify key directories exist
cat("\n=== Directory Check ===\n")
key_dirs <- c("Functions", "Model_a_proj", "EEG")
for (dir in key_dirs) {
  full_path <- file.path(project_root, dir)
  if (dir.exists(full_path)) {
    cat("✓", dir, "directory exists\n")
  } else {
    cat("✗", dir, "directory missing\n")
  }
}

cat("\nPath configuration complete!\n") 