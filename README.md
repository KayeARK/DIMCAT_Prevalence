# DIMCAT Prevalence Analysis

Spatial epidemiological analysis of African Animal Trypanosomosis (AAT) prevalence across Africa using INLA (Integrated Nested Laplace Approximation) for the accompanying paper "Mapping bovine trypanosomosis under diagnostic uncertainty: spatial prevalence estimates in Nigeria and Ethiopia" by Kaye et al.

## Overview

This repository contains the code infrastructure for analyzing bovine trypanosomiasis prevalence using diagnostic test data, environmental covariates, and Bayesian spatial modeling. The analysis produces prevalence maps and infected cattle assessments across Nigeria and Ethiopia.

**Note**: This repository contains the analysis code, some raw data files are excluded for size considerations (see below).

---

## 1. Software Dependencies and Operating Systems

### Supported Operating Systems
- **macOS**: 10.14 (Mojave) or higher
- **Windows**: 10 or higher (64-bit)
- **Linux**: Ubuntu 18.04+, CentOS 7+, or equivalent distributions

### Core Software Requirements

#### R Environment
- **R**: Version 4.3.0 or higher (tested on 4.3.3)
- **RStudio**: Version 2023.06.0 or higher (recommended IDE)

#### System Dependencies (Auto-handled by R packages)
- **PROJ**: Version 8.0+ (coordinate system transformations)
- **GDAL**: Version 3.4+ (geospatial data abstraction)
- **GEOS**: Version 3.9+ (geometric operations)

### Essential R Packages with Version Numbers

#### Bayesian Modeling
- **INLA**: Version 23.09.09 (Bayesian spatial modeling framework)
- **rstan**: Version 2.32.3 (MCMC for test sensitivity/specificity)

#### Spatial Data Processing
- **terra**: Version 1.7-46 (modern spatial raster handling)
- **raster**: Version 3.6-26 (legacy spatial operations)
- **sf**: Version 1.0-14 (simple features for vector data)
- **sp**: Version 2.0-0 (spatial classes and methods)

#### Data Acquisition
- **geodata**: Version 0.5-8 (WorldClim, elevation data access)
- **afrilearndata**: Development version from GitHub (African boundaries)

#### Visualization
- **ggplot2**: Version 3.4.4 (publication-quality plots)
- **viridis**: Version 0.6.4 (color scales for maps)
- **ggrepel**: Version 0.9.4 (smart label positioning)
- **cowplot**: Version 1.1.1 (plot arrangements)
- **gridExtra**: Version 2.3 (grid layouts)

#### Data Manipulation
- **dplyr**: Version 1.1.3 (data transformation)
- **tidyr**: Version 1.3.0 (data tidying)
- **readxl**: Version 1.4.3 (Excel file reading)

#### Additional Dependencies
- **concaveman**: Version 1.1.0 (boundary creation)
- **devtools**: Version 2.4.5 (GitHub package installation)

---

## 2. Installation Guide

### Step 1: R and RStudio Installation (15-30 minutes)

#### Option A: macOS
```bash
# Install R from CRAN
# Download from: https://cran.r-project.org/bin/macosx/
# Install RStudio from: https://posit.co/download/rstudio-desktop/
```

#### Option B: Windows
```bash
# Download R from: https://cran.r-project.org/bin/windows/base/
# Download RStudio from: https://posit.co/download/rstudio-desktop/
# Run installers with administrator privileges
```

#### Option C: Linux (Ubuntu/Debian)
```bash
# Update system and install R
sudo apt update
sudo apt install -y r-base r-base-dev

# Install system dependencies for spatial packages
sudo apt install -y libproj-dev libgdal-dev libgeos-dev libudunits2-dev

# Download and install RStudio
wget https://download1.rstudio.org/desktop/jammy/amd64/rstudio-2023.06.0-421-amd64.deb
sudo dpkg -i rstudio-2023.06.0-421-amd64.deb
```

### Step 2: R Package Installation (30-60 minutes)

**Estimated Installation Time**: 45 minutes (varies by internet speed and system)

```r
# Navigate to project directory and run installation script
setwd("/path/to/DIMCAT_Prevalence")
source("install_packages.R")

# Verify installations
source("test_packages.R")
```

**Expected Output**: All packages should show ✓ status. Any ✗ indicates installation issues.

### Step 3: Data Setup (10-15 minutes for automated data)

```r
# Climate and elevation data will download automatically during first run
# Manual data downloads required for:
# - Livestock density rasters (GLW4)
# - Land use data (ESA WorldCover)
# - Tsetse distribution maps
```

### Installation Troubleshooting

#### Common Issues:

**INLA Installation Failure**:
```r
# If INLA fails to install
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```

**Spatial Package Compilation Issues (Linux)**:
```bash
# Install additional development libraries
sudo apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev
```

**Memory Issues During Installation**:
```r
# Increase memory limit (Windows)
memory.limit(size=8000)
```

---

## 3. Demo Instructions

### Quick Start Demo (Runtime: ~15-20 minutes)

This demonstration runs the complete analysis workflow starting with test performance inference, then prevalence analysis for Nigeria.

#### Step 1: Verify Installation
```r
# Open RStudio and navigate to project directory
# Set working directory to main project folder
setwd("/path/to/DIMCAT_Prevalence")

# Test all packages are working
source("test_packages.R")
```
**Expected Output**: All packages show ✓ status (green checkmarks)

#### Step 2: Run Test Sensitivity/Specificity Analysis (MCMC)
```r
# Run Bayesian inference for diagnostic test performance
source("Code/TestSensSpec/SensSpec.r")
```

**Expected Runtime**: 5-8 minutes  
**Expected Outputs**:
- Console output showing MCMC progress
- Posterior samples for test sensitivity and specificity
- Model convergence diagnostics

#### Step 3: Run Nigeria Prevalence Analysis
```r
# Run fine-scale prevalence analysis (uses test performance results)
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Prevalence_plots_optimised_fine.R")
```

**Expected Runtime**: 8-12 minutes  
**Expected Outputs**:
- Console output showing INLA model fitting progress
- Prevalence maps saved as PNG files
- Cattle totals printed to console:
  ```
  Total infected cattle (mean): XXXXXX
  Total cattle in tsetse zone: XXXXXXX
  ```

#### Step 4: Generate Cattle-at-Risk Maps
```r
# Create choropleth maps for administrative units
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine.R")
```

**Expected Runtime**: 3-5 minutes  
**Expected Outputs**:
- PDF files: `nga_cattle_at_risk_choropleth_mean.pdf`
- Analysis plots: `nga_cattle_at_risk_mean_analysis.pdf`

#### Step 5: View Results
```r
# Check generated files (they will be in the Analysis_NGA directory)
list.files("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/", pattern = "nga_cattle.*\\.pdf")
```

**Expected Files**:
- `nga_cattle_at_risk_choropleth_mean.pdf` - Choropleth map
- `nga_cattle_at_risk_choropleth_mean_log.pdf` - Log-scale version
- `nga_cattle_at_risk_mean_analysis.pdf` - Distribution analysis

### Demo Validation

**Successful Demo Indicators**:
1. ✅ No error messages during script execution
2. ✅ PDF maps generated showing Nigerian states with color-coded cattle risk
3. ✅ Console output showing realistic cattle numbers (thousands to hundreds of thousands)
4. ✅ Maps display proper geographic boundaries and legend

**Common Demo Issues**:
- **Long download times**: Climate data downloads ~100MB on first run
- **Memory warnings**: Normal for spatial analysis, increase R memory if needed
- **Missing boundaries**: Ensure internet connection for administrative boundary downloads

---

## 4. Instructions for Use

### Project Structure Overview

```
DIMCAT_Prevalence/
├── Code/
│   ├── Prevalence/Bovine BCT and PCR/
│   │   ├── Analysis_NGA/          # Nigeria-specific analysis
│   │   ├── Analysis_ETH/          # Ethiopia-specific analysis
│   │   └── Combined_prevalence_comparison.R
│   ├── TestSensSpec/              # Test performance analysis
│   ├── Correlations/              # Spatial correlation analysis
│   └── Cost effectiveness analysis/
├── Data/
│   ├── ContAtlas_v2/Bovine data/  # Raw diagnostic data
│   └── Covariates/                # Environmental data (auto-downloaded)
└── README.md
```

### Core Analysis Workflows

#### Workflow 1: Complete Analysis Pipeline (Recommended)

**Purpose**: Full analysis pipeline from test performance inference to cattle-at-risk assessment

**Scripts to Run** (in order, all from main directory):
```r
# Set working directory to main project folder
setwd("/path/to/DIMCAT_Prevalence")

# 1. Test Performance Analysis (FIRST STEP)
source("Code/TestSensSpec/SensSpec.r")
source("Code/TestSensSpec/Case_adjustment.R")

# 2. Model Selection and Validation
source("Code/Prevalence/Bovine BCT and PCR/Model_selection_NGA.r")  # or Model_selection_ETH.r

# 3. Main Prevalence Analysis
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Prevalence_plots_optimised_fine.R")  # or Analysis_ETH

# 4. Cattle-at-Risk Assessment
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine.R")  # or Analysis_ETH

# 5. Administrative Unit Analysis
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/LGA_choropleth.R")  # or Analysis_ETH
```

**Expected Runtime**: 60-120 minutes per country  
**Outputs**: Test performance parameters, prevalence maps, cattle estimates, administrative unit summaries

#### Workflow 2: Country-Specific Prevalence Analysis Only

**Purpose**: Generate prevalence maps and cattle-at-risk estimates for individual countries (assumes test performance already estimated)

**Scripts to Run** (in order, from main directory):
```r
# Prevalence analysis for Nigeria
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Prevalence_plots_optimised_fine.R")
source("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine.R")

# OR for Ethiopia
source("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Prevalence_plots_optimised_fine.R")
source("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine.R")
```

**Expected Runtime**: 45-90 minutes per country  
**Outputs**: Prevalence maps, cattle estimates, administrative unit summaries

#### Workflow 3: Cross-Country Comparison

**Purpose**: Compare prevalence patterns and cattle burden between Nigeria and Ethiopia

**Script to Run** (from main directory):
```r
source("Code/Prevalence/Bovine BCT and PCR/Combined_prevalence_comparison.R")
```

**Expected Runtime**: 15-25 minutes  
**Outputs**: Side-by-side comparison maps, correlation analysis

#### Workflow 4: Test Performance Analysis Only

**Purpose**: Estimate diagnostic test sensitivity and specificity using Bayesian methods

**Scripts to Run** (from main directory):
```r
source("Code/TestSensSpec/SensSpec.r")
source("Code/TestSensSpec/Case_adjustment.R")
```

**Expected Runtime**: 20-30 minutes  
**Outputs**: Test performance parameters, adjusted prevalence estimates

### Key Parameters and Customization

#### Geographic Scope Modification
```r
# In any prevalence script, modify countries to analyze
countries_to_infer <- c("Nigeria")  # Single country
countries_to_infer <- c("Nigeria", "Ethiopia")  # Multiple countries
countries_to_infer <- c("Africa")  # Continental analysis
```

#### Spatial Resolution Settings
```r
# Fine-scale analysis (country-level)
mesh <- inla.mesh.2d(loc = coo, offset = c(50, 100), cutoff = 3, max.edge = c(6, 15))

# Coarse-scale analysis (continental)
mesh <- inla.mesh.2d(loc = coo, offset = c(50, 100), cutoff = 1, max.edge = c(30, 60))
```

#### Model Configuration
```r
# PC priors (recommended for new analysis)
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE)

# Matérn covariance (legacy compatibility)
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
```

### Data Requirements and Setup

#### Automatic Data Downloads
**Climate and Elevation**: Downloaded automatically via `geodata` package
- WorldClim v2.1 climate data (~100MB per country)
- SRTM elevation data (~50MB per country)

#### Manual Data Requirements
**Livestock Density**: Download GLW4 cattle density rasters
- Source: https://www.fao.org/livestock-systems/global-distributions/en/
- Place in: `Data/Covariates/livestock/`

**Tsetse Distribution**: Download FAO tsetse maps
- Source: https://openknowledge.fao.org/items/956f7aad-64e2-4bff-af3b-623b2215587c
- Place in: `Data/Covariates/tsenumbspec/`

#### Diagnostic Data Format
Expected Excel format for new countries:
```
Columns Required:
- latitude, longitude (decimal degrees, WGS84)
- sample_size (number of animals tested)
- positive (number of positive results)
- Additional covariates as available
```

### Output Interpretation

#### Prevalence Maps
- **Color Scale**: Viridis (purple = low, yellow = high prevalence)
- **Units**: Proportion (0-1) or percentage (0-100%)
- **Uncertainty**: Displayed via confidence interval maps

#### Cattle-at-Risk Maps
- **Color Scale**: Plasma (dark = low risk, bright = high risk)
- **Units**: Number of cattle at risk per administrative unit
- **Categories**: 
  - `cattle_at_risk`: Areas with both cattle and tsetse
  - `no_cattle_but_risk`: Tsetse areas without cattle
  - `no_risk_or_cattle`: Safe areas

#### Administrative Unit Outputs
- **LGA/Zone Level**: Second-level administrative divisions
- **Metrics**: Mean prevalence, cattle at risk, confidence intervals
- **Ranking**: Priority order for intervention targeting

### Performance Optimization

#### Memory Management
```r
# Increase R memory limit (Windows)
memory.limit(size = 8000)

# Monitor memory usage during analysis
print(paste("Memory usage:", round(memory.size()/1024, 1), "GB"))
```

#### Computational Efficiency
- **Mesh Resolution**: Coarser meshes for exploratory analysis
- **Covariate Selection**: Remove redundant environmental variables
- **Parallel Processing**: INLA automatically uses multiple cores

### Troubleshooting Common Issues

#### Coordinate System Problems
```r
# Validate coordinate ranges before analysis
summary(data$latitude)   # Should be within country bounds
summary(data$longitude)  # Check for coordinate swapping
```

#### Model Convergence Issues
- **Check mesh quality**: Ensure adequate spatial coverage
- **Verify priors**: Use PC priors for better convergence
- **Inspect warnings**: INLA provides diagnostic messages

#### Memory/Performance Issues
- **Reduce mesh resolution** for initial exploration
- **Subset data** for testing new modifications
- **Monitor system resources** during long computations

### Adding New Countries

To extend analysis to additional African countries:

1. **Data Preparation**: Add country diagnostic data to `Data/ContAtlas_v2/Bovine data/`
2. **Geographic Setup**: Determine country bounds and ISO codes
3. **Script Modification**: Update `countries_to_infer` parameter
4. **Validation**: Check coordinate ranges and covariate extraction
5. **Mesh Configuration**: Adjust spatial resolution for country size

### Citation and Usage

When using this analysis framework, please cite:
- **Primary Paper**: Kaye et al. "Mapping bovine trypanosomosis under diagnostic uncertainty: spatial prevalence estimates in Nigeria and Ethiopia"
- **INLA Package**: Rue, H., et al. "Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations"

---

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
- `Cattle_at_risk_fine.R` - Fine-scale cattle-at-risk analysis
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
