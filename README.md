# CMC-GGM: Graph Estimation Computing Project

This R project implements methods for graph estimation computing, focusing on Gaussian graphical models with various transformation techniques.

## Project Structure

```
CMC-GGM/
├── Functions/          # Core utility functions and algorithms
├── Model_a*/          # Model A experiments and variations
├── Model_b*/          # Model B experiments and variations  
├── Model_c*/          # Model C experiments and variations
├── EEG/               # EEG data analysis
├── Gene_Protein/      # Gene-protein network analysis
├── fgm/               # Functional graphical models
├── Plot/              # Plotting and visualization scripts
├── Textures/          # Generated plots and textures
├── proj_cmc*/         # Python projection methods
└── setup files       # Environment setup scripts
```

## Quick Setup

### Prerequisites
- R (version 3.6 or higher)
- RStudio (recommended)
- Git

### Installation

1. **Clone the repository** (if not already done):
   ```bash
   git clone <repository-url>
   cd CMC-GGM
   ```

2. **Set up the R environment**:
   ```r
   # In R or RStudio, run:
   source("setup_environment.R")
   ```

3. **Configure paths**:
   ```r
   # Load path configuration
   source("config_paths.R")
   
   # Load all utility functions
   source_all_functions()
   ```

## Key Features

### Core Algorithms
- **ADMM-based optimization** for group lasso problems
- **Thresholding methods** for precision matrix estimation  
- **Multiple transformation techniques**:
  - Copula transformations
  - Parametric copula transformations
  - Linear transformations

### Model Types
- **Model A**: Basic graphical model estimation
- **Model B**: Extended models with additional constraints
- **Model C**: Complex models with rank constraints

### Applications
- **EEG Analysis**: Functional connectivity analysis
- **Gene-Protein Networks**: Biological network reconstruction
- **Simulation Studies**: Performance evaluation across different scenarios

## Required R Packages

The setup script will automatically install these packages:

### Core Packages
- `Matrix`, `matrixcalc`, `pracma` - Matrix operations
- `igraph`, `huge` - Graph theory and network analysis
- `gglasso`, `BB` - Statistical modeling and optimization
- `mvtnorm` - Multivariate normal distributions

### Data Analysis
- `dplyr`, `tidyr` - Data manipulation
- `ggplot2` - Visualization
- `doParallel`, `foreach` - Parallel computing

### Specialized Packages
- `fda`, `fdapace` - Functional data analysis
- `eegkit`, `eegkitdata` - EEG analysis
- `randtoolbox` - Random number generation

## Usage Examples

### Basic Graph Estimation
```r
# Load configuration and functions
source("config_paths.R")
source_all_functions()

# Set working directory to a model folder
setwd(model_a_proj_path)

# Run a basic experiment
source("model_a_d03_exp.R")
```

### Custom Analysis
```r
# Generate synthetic data
data <- gen_data(n=200, p=50, m=3, model='model1', run.ind=1)

# Apply transformations
data_transformed <- mult.trans(X=data$X, p=50, fun='copula')

# Estimate precision matrix
S <- cov(data_transformed)
result <- vggm.thr(X=S, p=50, rho=1e-3, ridge=TRUE)
```

### Parallel Computing Setup
```r
library(doParallel)

# Set up cluster for parallel processing
cl <- makeCluster(detectCores() - 1)  # Use all cores except one
registerDoParallel(cl)

# Your parallel computation here...

# Don't forget to stop the cluster
stopCluster(cl)
```

## File Organization

### Functions Directory
Contains core algorithms:
- `admm.R` - ADMM optimization routines
- `helper.R` - Utility functions
- `generate_data.R` - Data simulation functions
- `evaluation.R` - Performance evaluation metrics
- `trans.R` - Transformation functions

### Model Directories
Each model directory contains:
- Experiment scripts (e.g., `model_a_d03_exp.R`)
- Power analysis scripts (e.g., `model_a_d03_power.R`)
- Results subdirectory for output files

## Running Experiments

### Single Experiment
```r
# Navigate to model directory
setwd(model_a_proj_path)

# Run specific experiment
source("model_a_d03_exp.R")
```

### Batch Processing
```r
# Run multiple seeds for statistical significance
for(seed in 1:100) {
  system(paste("Rscript model_a_d03_exp.R", seed))
}
```