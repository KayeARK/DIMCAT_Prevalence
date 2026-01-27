# DIMCAT Prevalence Analysis - Copilot Instructions

## Project Overview
This is a spatial epidemiological research project analyzing African Animal Trypanosomiasis (AAT) prevalence across Africa using INLA (Integrated Nested Laplace Approximation) Bayesian modeling. The codebase processes diagnostic test data, environmental covariates, and produces prevalence maps, cattle-at-risk assessments, and cost-effectiveness analyses for diagnostic strategies.

## Core Architecture

### Data Flow Pattern
1. **Raw Data**: `Data/ContAtlas_v2/Bovine data/` contains diagnostic test data; `Data/ContAtlas_v3/` for newer data versions
2. **Covariate Extraction**: Scripts extract climate, elevation, livestock density, and land use data for each location
3. **Model Building**: INLA models with spatial random effects and environmental covariates
4. **Prediction**: Generate prevalence surfaces and uncertainty estimates across Africa
5. **Post-Processing**: Cattle-at-risk calculations and administrative unit aggregation
6. **Economic Analysis**: Cost-effectiveness and Expected Value of Sample Information (EVSI) calculations

### Directory Structure Logic
- `Code/Prevalence/`: Main analysis divided by diagnostic test type (BCT, PCR, combined)
  - `Continental/`: Africa-wide analysis (excludes Madagascar automatically)
  - `Analysis_[COUNTRY]/`: Country-specific fine-scale analysis with LGA/zone-level outputs
  - `Harriet_Layer/`: Alternative analysis approaches using different model configurations
- `Code/AbsencePresence/`: Presence/absence modeling (simpler than prevalence) 
- `Code/TestSensSpec/`: Bayesian latent class models for test sensitivity/specificity
- `Code/Cost effectiveness analysis/`: Economic decision analysis and EVSI computations
- `Code/Correlations/`: Inter-species correlation analysis
- `Data/Covariates/`: Spatial raster data organized by type (climate, livestock, etc.)

## Key Patterns & Conventions

### Country-Specific Analysis Pattern
Most scripts use this pattern for country-specific analysis:
```r
countries_to_infer=c("Nigeria")  # Set target country
# Data loading uses country code for file paths
r_elv <- elevation_30s(country = "NGA", path = "Data/Covariates")
```
Always use ISO 3-letter country codes (NGA, ZWE, KEN, ETH, etc.) in file paths and geodata functions.

### Continental vs Country Analysis Workflow
**Continental Analysis** (`countries_to_infer=c("Africa")`):
- Uses `Code/Prevalence/[TestType]/Continental/` directory structure
- Automatically excludes Madagascar: `border[!(border[,1] > 40 & border[,1] < 60 & border[,2] > -30 & border[,2] < -10),]`
- Larger spatial meshes with coarser resolution for computational efficiency

**Country-Specific Analysis**:
- Uses `concaveman` package to create analysis boundaries around data points
- Finer spatial resolution and denser meshes
- Country boundaries determined by `afrilearndata::africountries` with filtering

### Covariate Extraction Workflow
Standard pattern across all prevalence scripts:
1. Load base data from Excel files in `Data/ContAtlas_v2/Bovine data/`
2. Clean data: remove NA coordinates, fix sample sizes vs positive cases
3. **Coordinate Validation**: Check coordinate ranges and swap if necessary (common issue)
```r
# Ethiopia example: check if coordinates need swapping
if (any(data$latitude < 3 | data$latitude > 15)) {
  # Implement coordinate swap logic based on geographic bounds
}
```
4. Extract environmental covariates using `terra::extract()` with coordinate pairs
5. Calculate annual averages for monthly climate data (precipitation, temperature)
6. Handle missing values: set livestock densities to 0, interpolate climate data

### INLA Modeling Workflow
**Standard INLA Pipeline** (consistent across all prevalence scripts):
```r
# 1. Mesh creation (country-specific mesh parameters)
mesh <- inla.mesh.2d(loc = coo, offset = c(50, 100), cutoff = 3, max.edge = c(6, 15))  # Fine mesh for countries
mesh <- inla.mesh.2d(loc = coo, offset = c(50, 100), cutoff = 1, max.edge = c(30, 60))  # Coarse mesh for continental

# 2. SPDE setup with PC priors (preferred) or Matérn
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE)  # PC priors (newer scripts)
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)    # Matérn (legacy scripts)
indexs <- inla.spde.make.index("s", spde$n.spde)

# 3. Stacks for estimation and prediction
stk.e <- inla.stack(tag = "est", data = list(y = positive, numtrials = sample_size), ...)
stk.p <- inla.stack(tag = "pred", data = list(y = NA, numtrials = NA), ...)

# 4. Model fitting with standardized controls
res <- inla(formula, data = inla.stack.data(stk.full), family = "binomial", 
           control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE))
```

### Data Validation Patterns
- **Geographic Bounds**: Each country has expected coordinate ranges - validate before analysis
- **Sample Size Logic**: Ensure positive cases ≤ total sample size
- **Coordinate Systems**: Maintain WGS84 (EPSG:4326) throughout analysis
- **Missing Data**: Different strategies for different covariate types (zero-fill vs interpolation)

### INLA Model Structure
- **Response**: Binomial (positive cases, sample size)
- **Fixed Effects**: Standardized environmental covariates
- **Random Effects**: Spatial field using SPDE approach with Matérn covariance
- **Priors**: PC priors for hyperparameters (range, sigma)
- **Model Selection**: WAIC-based comparison stored in `Model_selection.r` files

### File Naming Conventions
- Scripts: Use descriptive names like `INLA_Prevalence.r`, `Model_selection.r`
- Results: Country-specific folders under `Results/` (e.g., `Results/Nigeria/`)
- Outputs: Systematic naming like `Projections_model_1.csv` for ensemble outputs
- Plots: Descriptive names with country codes and plot type

### Visualization Patterns
- Use `ggplot2` with `viridis` color scales for prevalence maps
- **Cattle-at-Risk Maps**: Use `plasma` color palette with specific risk categories:
  ```r
  # Risk categories for tsetse zone classification
  risk_category = case_when(
    cattle_at_risk > 0 ~ "cattle_at_risk",
    cattle_at_risk == 0 & tsetse_zone == TRUE ~ "no_cattle_but_risk", 
    TRUE ~ "no_risk_or_cattle"
  )
  ```
- **Choropleth Maps**: LGA/zone-level administrative boundaries with `gadm()` level 2 data
- Ridge plots with `ggridges` for uncertainty visualization (see `Cattle_at_risk.R`)
- **Combined Plots**: Histogram + boxplot combinations using `grid.arrange()` for distribution analysis
- **Inset Maps**: Continental Africa context using `cowplot::draw_plot()` for country-specific maps
- Correlation heatmaps saved as PNG files with country codes
- PDFs for publication-quality maps and plots

## Critical Dependencies
- **INLA**: Core Bayesian modeling framework
- **terra/raster**: Spatial data handling and covariate extraction
- **geodata**: Access to WorldClim and elevation data
- **sf/sp**: Spatial data structures and operations
- **afrilearndata**: African country boundaries

## Data Dependencies
- Climate data downloaded via `worldclim_country()` to `Data/Covariates/`
- Livestock density rasters in `Data/Covariates/livestock/` by species and year
- **Tsetse distribution**: Binary raster at `Data/Covariates/tsenumbspec/` (values > 1 converted to 1)
- Country-specific diagnostic data in `Data/ContAtlas_v2/Bovine data/`

## Tsetse Zone Classification Patterns
**Critical for cattle-at-risk analysis**:
```r
# Load and process tsetse raster
tsetse_raster <- raster("Data/Covariates/tsenumbspec")
tsetse_raster[tsetse_raster > 1] <- 1  # Convert to binary

# Area-based extraction for LGA/zone classification
lga_tsetse_stats <- raster::extract(tsetse_raster, admin_polygons, fun = mean, na.rm = TRUE)

# Combine raster data with geographical knowledge
tsetse_states <- c("Abia", "Akwa Ibom", "Anambra", "Bayelsa", "Benue", ...)  # Known tsetse belt states
tsetse_zone <- (tsetse_mean > 0) | (admin_unit %in% tsetse_states)
```

## Common Debugging Points
- **Memory Issues**: INLA models with fine spatial meshes can be memory intensive
- **Covariate Extraction**: Check coordinate system consistency between data and rasters
- **Model Convergence**: Monitor INLA warning messages about mesh quality and prior specification
- **File Paths**: Ensure country codes match between script settings and data file structures

## Cost-Effectiveness Analysis Patterns

### EVSI (Expected Value of Sample Information) Workflow
The `Code/Cost effectiveness analysis/` directory contains sophisticated economic analysis:
```r
# Standard EVSI calculation pattern
sample_all_parameters <- function(n_samples) {
  chain_indices <- sample(1:length(se_hct_samples), n_samples, replace = TRUE)
  # Sample from joint posterior for test performance and costs
}
# Compare expected costs under current vs perfect information
```

### Monte Carlo Uncertainty Propagation
- **Test Performance**: Use fitted Bayesian latent class models from `Code/TestSensSpec/latent_class_fit.rds`
- **Cost Parameters**: Gamma distributions for treatment costs and test costs
- **Prevalence**: INLA posterior samples with uncertainty from model ensemble
- **Sample Synchronization**: Use same `chain_indices` across parameter types to preserve correlation

### Economic Analysis Components
- **Decision Strategies**: NO_TEST, HCT (rapid test), PCR (lab test)
- **Cost Structure**: Testing costs, treatment costs, false positive costs, untreated infection costs
- **EVSI Calculation**: Current expected cost minus perfect information expected cost
- **Results Organization**: `results/`, `evsi_results/`, `sampling_design/` subdirectories

## Package Management & Setup

### Dependency Installation Pattern
```r
# Use install_packages.R for project setup
packages_to_test <- c("INLA", "terra", "geodata", "afrilearndata", ...)
# GitHub packages: devtools::install_github("afrimapr/afrilearndata")
```

### Package Testing Workflow
- `test_packages.R` provides systematic package validation
- Critical dependencies: INLA (Bayesian modeling), terra/geodata (spatial data), rstan (MCMC)
- Always test afrilearndata installation from GitHub as it's not on CRAN

## Results Management Patterns

### Output Directory Structure
- Country-specific results: `Code/[AnalysisType]/Results/[CountryCode]/`
- Cross-analysis results: `Code/Cost effectiveness analysis/results/`
- Systematic naming: `[Country][DataType].csv`, `[Country]Predictions.png`

### File Output Conventions
- **Model Results**: WAIC tables, covariate files, projection CSVs
- **Visualizations**: PNG for exploration, PDF for publication
- **Economic Analysis**: CSV files with uncertainty quantification
- **Intermediate Files**: RDS files for fitted models and posterior samples

## Testing & Validation
- Use WAIC for model comparison (lower is better)
- Cross-validation patterns in model selection scripts
- Sensitivity analysis through different mesh configurations
- Visual inspection of predicted surfaces for biological plausibility
- Monte Carlo convergence checks for EVSI calculations (typically 500-1000 samples)

## Key Development Workflows

### Cattle-at-Risk Analysis Pattern
1. Load INLA model predictions and uncertainty bounds
2. Apply tsetse zone classification (raster + geographical knowledge)
3. Create risk categories: `cattle_at_risk`, `no_cattle_but_risk`, `no_risk_or_cattle`
4. Generate choropleth maps with `plasma` palette for cattle counts
5. Create combined plots (histogram + boxplot) for distribution analysis
6. Save both regular and log-scale versions for different data ranges

### Coordinate Validation Workflow
**Critical Step**: Always validate coordinates before analysis
```r
# Check coordinate ranges for geographic plausibility
if (any(data$latitude < min_lat | data$latitude > max_lat)) {
  # Implement coordinate swap logic if needed
  # Ethiopia: lat 3-15°N, lon 33-48°E
  # Nigeria: lat 4-14°N, lon 3-15°E
}
```