# DIMCAT Prevalence Analysis

Spatial epidemiological analysis of African Animal Trypanosomosis (AAT) prevalence across Africa using INLA (Integrated Nested Laplace Approximation) for the accompanying paper "Mapping bovine trypanosomosis under diagnostic uncertainty: spatial prevalence estimates in Nigeria and Ethiopia" by Kaye et al.

## Overview

This repository contains the code infrastructure for analyzing bovine trypanosomiasis prevalence using diagnostic test data, environmental covariates, and Bayesian spatial modeling. The analysis produces prevalence maps and infected cattle assessments across Nigeria and Ethiopia.

**Note**: This repository contains the analysis code, some raw data files are excluded for size considerations (see [Data Sources](#data-sources) section below).

## Repository Structure

### Core Analysis (`Code/Prevalence/`)

#### Bovine BCT and PCR Analysis
**Main Directory**: `Code/Prevalence/Bovine BCT and PCR/`

##### Country-Specific Analysis
- **`Analysis_ETH/`** - Ethiopia-specific analysis
- **`Analysis_NGA/`** - Nigeria-specific analysis  

**Key Scripts in Each Country Directory**:
- `Prevalence_plots_optimised.R` - Main prevalence mapping and visualization
- `Prevalence_plots_optimised_fine.R` - High-resolution prevalence analysis
- `Cattle_at_risk.R` - Cattle-at-risk calculations and mapping
- `Cattle_uncertainty_summary.R` - Uncertainty analysis for cattle estimates
- `LGA_choropleth.R` - Administrative unit choropleth mapping
- `WAIC_Analysis.R` - Model comparison using WAIC

##### Model Selection Scripts
- `Model_selection_ETH.r` - Ethiopia model selection and validation
- `Model_selection_NGA.r` - Nigeria model selection and validation


### Supporting Analysis

#### Test Performance (`Code/TestSensSpec/`)
Bayesian latent class models for test sensitivity and specificity:
- `SensSpec.r` - Main sensitivity/specificity analysis
- `Case_adjustment.R` - Prevalence adjustment for test performance

#### Correlation Analysis (`Code/Correlations/`)
Inter-species and spatial correlation analysis:
- `Correlation.r` - Correlation analysis between species/locations

### Data Structure

**Important**: Some raw data files are **not included** in this repository due to size considerations.

#### Required Data Sources (Not Included)

**Environmental Covariates** (Large files - obtain from public sources):
- **Climate**: WorldClim v2.1 data (auto-downloaded via `geodata::worldclim_country()`)
- **Elevation**: SRTM 30s resolution (via `geodata::elevation_30s()`)
- **Livestock**: GLW4 livestock density rasters (https://www.fao.org/livestock-systems/global-distributions/en/)
- **Land Use**: ESA WorldCover 2021 (https://worldcover2021.esa.int/)
- **Tsetse Distribution**: FAO tsetse distribution map (https://openknowledge.fao.org/items/956f7aad-64e2-4bff-af3b-623b2215587c)
- **Population**: GPWv4 population density (https://cran.r-project.org/web/packages/geodata/geodata.pdf)

### Utility Scripts

#### Package Management
- `install_packages.R` - Install all required dependencies
- `test_packages.R` - Test package installations and versions

## Key Dependencies

### R Packages
- **INLA** - Bayesian modeling framework
- **terra/raster** - Spatial data handling  
- **sf/sp** - Spatial data structures
- **geodata** - Access to WorldClim and elevation data
- **afrilearndata** - African administrative boundaries
- **ggplot2** - Visualization
- **dplyr** - Data manipulation

### Data Requirements

#### Setting Up the Analysis Environment

**Step 1: Install Dependencies**
```r
source("install_packages.R")  # Install all required R packages
source("test_packages.R")     # Verify installations
```

**Step 2: Obtain Required Data**

Due to file size considerations, users must obtain some data separately:

**Environmental Covariates** (Automatically Downloaded):
- Climate and elevation data downloaded automatically by analysis scripts
- Livestock and land use data - follow source links above

**Step 3: Expected Directory Structure**
```
Data/
├── ContAtlas_v2/Bovine data/    # AAT data (obtain separately)
├── ContAtlas_v3/                # Updated AAT data (obtain separately)
└── Covariates/                  # Auto-downloaded environmental data
    ├── climate/
    ├── elevation/
    ├── livestock/
    ├── landuse/
    ├── population/
    └── tsenumbspec/
```

## Analysis Workflow

### 1. Data Preparation
```r
# Install dependencies
source("install_packages.R")

# Load and clean diagnostic data
# Extract environmental covariates
# Validate coordinates and sample sizes
```

### 2. Diagnostic test performance inference
```r
# MCMC for inference of sensitivity and specificity
# Case readjustment
```

### 3. Model Fitting
```r
# INLA Bayesian spatial modeling
# PC priors for hyperparameters
# Mesh creation for spatial effects
# Model comparison using WAIC
```

### 4. Prediction and Mapping
```r
# Generate prevalence surfaces
# Quantify uncertainty (confidence intervals)  
# Create publication-quality maps
# Administrative unit aggregation
```

### 5. Cattle-at-Risk Assessment
```r
# Combine prevalence with cattle density
# Apply tsetse distribution masks
# Calculate infected cattle numbers
# Generate risk choropleth maps
```


## Key Outputs

### Prevalence Maps
- Mean prevalence estimates
- Confidence intervals (lower/upper bounds)
- Uncertainty quantification
- Administrative unit aggregation

### Cattle-at-Risk Analysis
- Infected cattle numbers by administrative unit
- Risk category classification
- Burden concentration (Lorenz curves)
- Priority targeting lists

### Model Validation
- WAIC model comparison
- Residual analysis
- Cross-validation assessments


## File Naming Conventions

- **Scripts**: Descriptive names (e.g., `Prevalence_plots_optimised.R`)
- **Models**: Country-specific prefixes (e.g., `Model_selection_ETH.r`)
- **Results**: Systematic naming (e.g., `Projections_model_1.csv`)
- **Plots**: Descriptive with country codes (e.g., `nga_prevalence_mean.png`)

## Coordinate System

All spatial data uses **WGS84 (EPSG:4326)** geographic coordinates.

## Important Notes

### Exclusion of Some Data
Some raw data files are excluded due to file size constraints.

### Data Access Requirements
- **Environmental Data**: Publicly available from international organizations
- **Automated Downloads**: Climate/elevation data downloaded automatically by scripts

### Technical Requirements

### Data Validation
- Always validate coordinate ranges for geographic plausibility
- Check sample size logic (positive cases ≤ total sample)
- Handle missing covariate data appropriately

### Tsetse Distribution
- Binary classification (0/1) from species richness data
- Use nearest neighbor interpolation for binary rasters
- Apply tsetse masks for biologically realistic cattle-at-risk estimates

### Memory Considerations
- INLA models with fine spatial meshes can be memory intensive
- Consider computational resources for continental-scale analysis
