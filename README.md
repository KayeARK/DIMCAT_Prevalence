# DIMCAT Prevalence Analysis

Spatial epidemiological analysis of African Animal Trypanosomiasis (AAT) prevalence across Africa using INLA (Integrated Nested Laplace Approximation) Bayesian modeling.

## Overview

This repository contains R code for analyzing bovine trypanosomiasis prevalence using diagnostic test data, environmental covariates, and Bayesian spatial modeling. The analysis produces prevalence maps, cattle-at-risk assessments, and cost-effectiveness evaluations for diagnostic strategies across African countries.

## Key Features

- **Spatial Bayesian Modeling**: INLA framework with spatial random effects
- **Multi-country Analysis**: Continental Africa and country-specific (Ethiopia, Nigeria) studies
- **Diagnostic Test Integration**: BCT (rapid test) and PCR (laboratory) data
- **Uncertainty Quantification**: Full posterior distributions and confidence intervals
- **Economic Analysis**: Cost-effectiveness and Expected Value of Sample Information (EVSI)
- **Administrative Mapping**: LGA/zone-level prevalence aggregation

## Repository Structure

### Core Analysis (`Code/Prevalence/`)

#### Bovine BCT and PCR Analysis
**Main Directory**: `Code/Prevalence/Bovine BCT and PCR/`

##### Country-Specific Analysis
- **`Analysis_ETH/`** - Ethiopia-specific analysis
- **`Analysis_NGA/`** - Nigeria-specific analysis  
- **`Analysis_ETH_all_PCR/`** - Ethiopia with all PCR data
- **`Analysis_NGA_all_PCR/`** - Nigeria with all PCR data

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
- `Model_selection_ETH_non_linear.r` - Non-linear covariate modeling for Ethiopia
- `Model_selection_NGA_all_PCR.r` - Nigeria PCR-specific model selection

##### Specialized Analysis
- **`Harriet_Layer/`** - Alternative analysis approaches
  - `Analysis_ETH_max/`, `Analysis_ETH_mean/`, `Analysis_NGA_max/` - Different aggregation methods
- **`Cattle at risk distribution inequality/`** - Burden concentration analysis
  - `burden_concentration_analysis.R` - Lorenz curves and Gini coefficients
  - `priority_targeting_list.R` - Priority area identification

#### Individual Test Type Analysis
- **`Bovine BCT/`** - BCT-specific analysis
  - `Continental/` - Africa-wide BCT analysis
  - `INLA_Prevalence_ETH.r`, `INLA_Prevalence_NGA.r` - Country-specific BCT models
- **`Bovine PCR/`** - PCR-specific analysis
  - `INLA_Prevalence.r`, `Model_selection.r` - PCR modeling

### Supporting Analysis

#### Test Performance (`Code/TestSensSpec/`)
Bayesian latent class models for test sensitivity and specificity:
- `SensSpec.r` - Main sensitivity/specificity analysis
- `Case_adjustment.R` - Prevalence adjustment for test performance
- `Results.r` - Test performance results compilation

#### Correlation Analysis (`Code/Correlations/`)
Inter-species and spatial correlation analysis:
- `Correlation.r` - Correlation analysis between species/locations
- `Correlation continental.r` - Continental-scale correlation patterns

### Data Structure (`Data/`)

#### Raw Data
- **`ContAtlas_v2/Bovine data/`** - Diagnostic test data (Excel format)
- **`ContAtlas_v3/`** - Updated data versions

#### Spatial Covariates (`Data/Covariates/`)
- **`livestock/`** - Livestock density rasters by species and year
- **`tsenumbspec/`** - Tsetse fly distribution (continental Africa)
- **`tsetse_Harriet/`** - Alternative tsetse distribution maps
- Climate data (downloaded via `worldclim_country()`)
- Elevation data (SRTM 30s resolution)

### Utility Scripts

#### Package Management
- `install_packages.R` - Install all required dependencies
- `test_packages.R` - Test package installations and versions

#### Continental Analysis
- `Code/INLAContinental.r` - Continental-scale INLA modeling
- `Code/INLACountry.r` - Country-scale INLA setup
- `Code/prelim_analysis_2.R` - Preliminary data exploration

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
- Continental Atlas diagnostic test data (Excel files)
- WorldClim climate data (automatically downloaded)
- Livestock density rasters (GLW4 project)
- Tsetse fly distribution maps
- Administrative boundary data (GADM)

## Analysis Workflow

### 1. Data Preparation
```r
# Install dependencies
source("install_packages.R")

# Load and clean diagnostic data
# Extract environmental covariates
# Validate coordinates and sample sizes
```

### 2. Model Fitting
```r
# INLA Bayesian spatial modeling
# PC priors for hyperparameters
# Mesh creation for spatial effects
# Model comparison using WAIC
```

### 3. Prediction and Mapping
```r
# Generate prevalence surfaces
# Quantify uncertainty (confidence intervals)  
# Create publication-quality maps
# Administrative unit aggregation
```

### 4. Cattle-at-Risk Assessment
```r
# Combine prevalence with cattle density
# Apply tsetse distribution masks
# Calculate infected cattle numbers
# Generate risk choropleth maps
```

## Model Structure

### Response Variable
Binomial distribution: (positive cases, total sample size)

### Fixed Effects
Standardized environmental covariates:
- Temperature (annual mean, seasonality)
- Precipitation (annual, seasonality)  
- Elevation
- Livestock density
- Land cover

### Random Effects
Spatial random field using SPDE approach with Matérn covariance

### Prior Specification
PC (Penalized Complexity) priors for hyperparameters

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

## Usage Examples

### Basic Prevalence Analysis
```r
# Set target country
countries_to_infer <- c("Nigeria")  # or "Ethiopia"

# Run main analysis
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Prevalence_plots_optimised.R")
```

### Cattle-at-Risk Assessment  
```r
# Generate cattle-at-risk maps
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk.R")
```

### Model Selection
```r
# Compare model specifications
source("Code/Prevalence/Bovine BCT and PCR/Model_selection_NGA.r")
```

## File Naming Conventions

- **Scripts**: Descriptive names (e.g., `Prevalence_plots_optimised.R`)
- **Models**: Country-specific prefixes (e.g., `Model_selection_ETH.r`)
- **Results**: Systematic naming (e.g., `Projections_model_1.csv`)
- **Plots**: Descriptive with country codes (e.g., `nga_prevalence_mean.png`)

## Coordinate System

All spatial data uses **WGS84 (EPSG:4326)** geographic coordinates.

## Important Notes

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

## Citation

If you use this code, please cite the relevant publications and acknowledge the DIMCAT project.

## License

Please refer to the project license for usage terms.

## Contact

For questions about the analysis methods or code implementation, please refer to the project documentation or contact the development team.