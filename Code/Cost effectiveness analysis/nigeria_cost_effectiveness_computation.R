#!/usr/bin/env R
# ==============================================================================
# NIGERIA COST-EFFECTIVENESS ANALYSIS - COMPUTATION
# ==============================================================================
# Bayesian decision-theoretic framework for optimal test selection
# Analyzes NO_TEST, HCT, and PCR strategies across all Nigeria locations
# Uses Expected Value of Sample Information (EVSI) for decision making
# ==============================================================================

library(tidyverse)
library(parallel)

# Set working directory
setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")

cat("=== NIGERIA COST-EFFECTIVENESS ANALYSIS ===\n")
cat("Starting Bayesian framework computation...\n\n")

# ==============================================================================
# PARAMETERS AND TEST CHARACTERISTICS
# ==============================================================================

# ==============================================================================
# LOAD FITTED SENSITIVITY/SPECIFICITY FROM BAYESIAN LATENT CLASS MODEL
# ==============================================================================

cat("Loading fitted test performance parameters from TestSensSpec analysis...\n")

# Load the fitted latent class model
fitted_model <- readRDS("Code/TestSensSpec/latent_class_fit.rds")

# Extract posterior samples for test performance
library(rstan)
se_hct_samples <- extract(fitted_model, 'Se_hct')[[1]]
sp_hct_samples <- extract(fitted_model, 'Sp_hct')[[1]]
se_pcr_samples <- extract(fitted_model, 'Se_pcr')[[1]]
sp_pcr_samples <- extract(fitted_model, 'Sp_pcr')[[1]]

cat("Extracted", length(se_hct_samples), "posterior samples for test performance\n")

# Test performance parameters (from fitted Bayesian model)
test_params <- list(
  hct = list(
    sensitivity_samples = se_hct_samples,     # HCT sensitivity posterior samples
    specificity_samples = sp_hct_samples,     # HCT specificity posterior samples
    sensitivity_mean = mean(se_hct_samples),  # Mean: 0.298
    specificity_mean = mean(sp_hct_samples),  # Mean: 0.996
    cost_shape = 4.0,      # Gamma shape parameter for cost uncertainty
    cost_rate = 8.0        # Gamma rate parameter (shape/rate = mean = 0.50)
  ),
  pcr = list(
    sensitivity_samples = se_pcr_samples,     # PCR sensitivity posterior samples  
    specificity_samples = sp_pcr_samples,     # PCR specificity posterior samples
    sensitivity_mean = mean(se_pcr_samples),  # Mean: 0.689
    specificity_mean = mean(sp_pcr_samples),  # Mean: 0.997
    cost_shape = 9.0,      # Gamma shape parameter for cost uncertainty
    cost_rate = 0.6        # Gamma rate parameter (shape/rate = mean = 15.0)
  )
)

# Treatment and outcome costs with uncertainty (Gamma distributions)
cost_params <- list(
  treatment_cost_shape = 6.25,      # Gamma shape for treatment cost
  treatment_cost_rate = 2.5,        # Gamma rate (mean = 2.50)
  
  cost_untreated_shape = 25.0,      # Gamma shape for untreated cost  
  cost_untreated_rate = 0.5         # Gamma rate (mean = 50.0)
  
  # Note: False positive cost = treatment cost (same sampled value each draw)
)

# Number of Monte Carlo samples for uncertainty propagation
n_uncertainty_samples <- 1000

cat("Test Parameters (with correlated uncertainty):\n")
cat("HCT: Se =", round(test_params$hct$sensitivity_mean, 3), 
    " ± ", round(sd(test_params$hct$sensitivity_samples), 3), 
    ", Sp =", round(test_params$hct$specificity_mean, 3),
    " ± ", round(sd(test_params$hct$specificity_samples), 3), 
    ", Cost ~ Gamma(", test_params$hct$cost_shape, ",", test_params$hct$cost_rate, ")\n")
cat("PCR: Se =", round(test_params$pcr$sensitivity_mean, 3),
    " ± ", round(sd(test_params$pcr$sensitivity_samples), 3),
    ", Sp =", round(test_params$pcr$specificity_mean, 3),
    " ± ", round(sd(test_params$pcr$specificity_samples), 3),
    ", Cost ~ Gamma(", test_params$pcr$cost_shape, ",", test_params$pcr$cost_rate, ")\n")
cat("Note: Se/Sp pairs sampled from same MCMC chain position (preserves joint posterior correlation)\n")
cat("Treatment cost ~ Gamma(", cost_params$treatment_cost_shape, ",", cost_params$treatment_cost_rate, 
    ") [mean = $", round(cost_params$treatment_cost_shape/cost_params$treatment_cost_rate, 2), "]\n")
cat("Cost of untreated infection ~ Gamma(", cost_params$cost_untreated_shape, ",", cost_params$cost_untreated_rate,
    ") [mean = $", round(cost_params$cost_untreated_shape/cost_params$cost_untreated_rate, 2), "]\n")
cat("False positive cost = treatment cost (same value per Monte Carlo draw)\n")
cat("Monte Carlo samples for uncertainty:", n_uncertainty_samples, "\n\n")

# ==============================================================================
# LOAD NIGERIA INLA PROJECTIONS
# ==============================================================================

cat("Loading Nigeria INLA projection data...\n")

# Load all projection files from Nigeria analysis
projection_files <- list.files(
  "Code/Prevalence/Bovine BCT and PCR/Projections_NGA/", 
  pattern = "Projections_model_.*\\.csv$", 
  full.names = TRUE
)

if (length(projection_files) == 0) {
  stop("No projection files found! Please ensure INLA analysis has been completed.")
}

cat("Found", length(projection_files), "projection files\n")

# Load and combine all projections with full uncertainty
all_projections <- map_dfr(projection_files, function(file) {
  model_num <- str_extract(basename(file), "\\d+")
  read_csv(file, show_col_types = FALSE) %>%
    select(Latitude, Longitude, value, variable) %>%
    rename(latitude = Latitude, longitude = Longitude) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    rename(mean = Mean, 
           q025 = `2.5th percentile`, 
           q975 = `97.5th percentile`) %>%
    mutate(model = as.numeric(model_num))
})

cat("Loaded", nrow(all_projections), "total projections\n")
cat("Unique locations:", length(unique(paste(all_projections$longitude, all_projections$latitude))), "\n")
cat("Models per location:", length(unique(all_projections$model)), "\n")

# Check prevalence uncertainty (both between-model and within-model)
prevalence_range_mean <- range(all_projections$mean, na.rm = TRUE)
prevalence_range_q025 <- range(all_projections$q025, na.rm = TRUE)
prevalence_range_q975 <- range(all_projections$q975, na.rm = TRUE)

cat("Mean prevalence range:", round(prevalence_range_mean[1] * 100, 2), "% to", 
    round(prevalence_range_mean[2] * 100, 2), "%\n")
cat("2.5th percentile range:", round(prevalence_range_q025[1] * 100, 2), "% to", 
    round(prevalence_range_q025[2] * 100, 2), "%\n")
cat("97.5th percentile range:", round(prevalence_range_q975[1] * 100, 2), "% to", 
    round(prevalence_range_q975[2] * 100, 2), "%\n")

# Calculate average credible interval width
mean_ci_width <- mean((all_projections$q975 - all_projections$q025) * 100, na.rm = TRUE)
cat("Mean 95% credible interval width:", round(mean_ci_width, 2), "percentage points\n")
cat("Both between-model and within-model uncertainty will be propagated\n\n")

# ==============================================================================
# BAYESIAN DECISION FUNCTIONS WITH UNCERTAINTY PROPAGATION
# ==============================================================================
# 
# ENHANCED UNCERTAINTY PROPAGATION: Both between-model and within-model uncertainty
# 
# BETWEEN-MODEL UNCERTAINTY: 
# - Each location has ~200 different INLA models with different prevalence estimates
# - Monte Carlo sampling randomly selects which model to use for each draw
# 
# WITHIN-MODEL UNCERTAINTY:
# - Each INLA model provides mean + 95% credible interval (2.5%, 97.5%)
# - For each selected model, we sample from approximate normal distribution
# - SD estimated from quantiles: SD ≈ (q975 - q025) / 3.92
# - Samples truncated to [0,1] to ensure valid prevalence values
#
# This captures the FULL uncertainty structure from INLA spatial modeling

# Sample costs from gamma distributions
sample_costs <- function(n_samples) {
  # Sample treatment cost once
  treatment_cost <- rgamma(n_samples, shape = cost_params$treatment_cost_shape, 
                          rate = cost_params$treatment_cost_rate)
  
  list(
    treatment_cost = treatment_cost,
    cost_untreated = rgamma(n_samples, shape = cost_params$cost_untreated_shape, 
                           rate = cost_params$cost_untreated_rate),
    cost_false_positive = treatment_cost  # Same as treatment cost!
  )
}

# Sample test performance parameters
sample_test_params <- function(n_samples) {
  # Sample from posterior distributions - use SAME chain position for joint estimation
  # HCT and PCR parameters were jointly estimated, so preserve correlation structure
  chain_indices <- sample(1:length(test_params$hct$sensitivity_samples), n_samples, replace = TRUE)
  
  list(
    hct_sensitivity = test_params$hct$sensitivity_samples[chain_indices],
    hct_specificity = test_params$hct$specificity_samples[chain_indices],
    hct_cost = rgamma(n_samples, shape = test_params$hct$cost_shape, rate = test_params$hct$cost_rate),
    
    pcr_sensitivity = test_params$pcr$sensitivity_samples[chain_indices],  # Same indices!
    pcr_specificity = test_params$pcr$specificity_samples[chain_indices],  # Same indices!
    pcr_cost = rgamma(n_samples, shape = test_params$pcr$cost_shape, rate = test_params$pcr$cost_rate)
  )
}

# Expected cost for NO_TEST strategy with uncertainty
expected_cost_no_test_uncertain <- function(prevalence, cost_samples) {
  prevalence * cost_samples$cost_untreated
}

# Expected cost for testing strategy with uncertainty
expected_cost_test_uncertain <- function(prevalence, sensitivity, specificity, test_cost, cost_samples) {
  # Probability of positive test result
  prob_positive <- sensitivity * prevalence + (1 - specificity) * (1 - prevalence)
  
  # Expected costs
  cost_testing <- test_cost
  cost_treatment <- prob_positive * cost_samples$treatment_cost
  cost_false_negative <- (1 - sensitivity) * prevalence * cost_samples$cost_untreated
  cost_false_positive <- (1 - specificity) * (1 - prevalence) * cost_samples$cost_false_positive
  
  return(cost_testing + cost_treatment + cost_false_negative + cost_false_positive)
}

# Sample prevalence incorporating both between-model and within-model uncertainty
# Calculate expected costs across uncertainty for one location  
calculate_location_costs <- function(prevalence_means, prevalence_q025, prevalence_q975, n_samples = n_uncertainty_samples) {
  # Create prevalence samples incorporating both between-model and within-model uncertainty
  prevalence_samples <- numeric(n_samples)
  n_models <- length(prevalence_means)
  
  for (i in 1:n_samples) {
    # First, sample which INLA model to use (between-model uncertainty)
    model_idx <- sample(n_models, 1)
    
    # Then sample from within-model uncertainty using approximate normal
    # Estimate SD from quantiles: SD ≈ (q975 - q025) / 3.92
    estimated_sd <- (prevalence_q975[model_idx] - prevalence_q025[model_idx]) / 3.92
    
    # Sample from truncated normal (truncated at 0 and 1)
    prevalence_samples[i] <- pmax(0, pmin(1, 
      rnorm(1, mean = prevalence_means[model_idx], sd = estimated_sd)
    ))
  }
  
  # Sample uncertain parameters
  cost_samples <- sample_costs(n_samples)
  test_samples <- sample_test_params(n_samples)
  
  # Calculate costs for each strategy across uncertainty samples
  # Now prevalence varies across samples too
  costs_no_test <- expected_cost_no_test_uncertain(prevalence_samples, cost_samples)
  
  costs_hct <- expected_cost_test_uncertain(
    prevalence_samples, test_samples$hct_sensitivity, test_samples$hct_specificity, 
    test_samples$hct_cost, cost_samples
  )
  
  costs_pcr <- expected_cost_test_uncertain(
    prevalence_samples, test_samples$pcr_sensitivity, test_samples$pcr_specificity,
    test_samples$pcr_cost, cost_samples
  )
  
  # Return summary statistics
  list(
    cost_no_test_mean = mean(costs_no_test),
    cost_no_test_sd = sd(costs_no_test),
    
    cost_hct_mean = mean(costs_hct),
    cost_hct_sd = sd(costs_hct),
    evsi_hct_mean = mean(pmax(0, costs_no_test - costs_hct)),
    evsi_hct_sd = sd(pmax(0, costs_no_test - costs_hct)),
    
    cost_pcr_mean = mean(costs_pcr),
    cost_pcr_sd = sd(costs_pcr),
    evsi_pcr_mean = mean(pmax(0, costs_no_test - costs_pcr)),
    evsi_pcr_sd = sd(pmax(0, costs_no_test - costs_pcr)),
    
    # Probability each strategy is optimal
    prob_no_test_optimal = mean((costs_no_test <= costs_hct) & (costs_no_test <= costs_pcr)),
    prob_hct_optimal = mean((costs_hct <= costs_no_test) & (costs_hct <= costs_pcr)),
    prob_pcr_optimal = mean((costs_pcr <= costs_no_test) & (costs_pcr <= costs_hct))
  )
}

cat("Defined Bayesian decision functions with full prevalence uncertainty propagation\n")

# ==============================================================================
# BATCH PROCESSING FUNCTION
# ==============================================================================

process_batch <- function(batch_data) {
  # Group by location and collect all prevalence data (mean + quantiles)
  location_results <- batch_data %>%
    group_by(longitude, latitude) %>%
    summarise(
      n_models = n(),
      prevalence_means = list(mean),
      prevalence_q025_vals = list(q025),
      prevalence_q975_vals = list(q975),
      mean_prevalence = mean(mean),
      sd_prevalence = sd(mean),
      q025_prevalence = quantile(mean, 0.025),
      q975_prevalence = quantile(mean, 0.975),
      .groups = "drop"
    )
  
  # Calculate costs with uncertainty for each location
  cat("  Computing costs with full INLA uncertainty propagation...")
  
  cost_results <- pmap_dfr(location_results, function(longitude, latitude, n_models, prevalence_means, 
                                                     prevalence_q025_vals, prevalence_q975_vals,
                                                     mean_prevalence, sd_prevalence, q025_prevalence, q975_prevalence) {
    location_costs <- calculate_location_costs(prevalence_means[[1]], prevalence_q025_vals[[1]], prevalence_q975_vals[[1]])
    
    # Determine optimal strategy based on mean costs
    min_cost_strategy <- which.min(c(location_costs$cost_no_test_mean, 
                                    location_costs$cost_hct_mean, 
                                    location_costs$cost_pcr_mean))
    
    optimal_strategy <- c("NO_TEST", "HCT", "PCR")[min_cost_strategy]
    min_cost <- c(location_costs$cost_no_test_mean, 
                  location_costs$cost_hct_mean, 
                  location_costs$cost_pcr_mean)[min_cost_strategy]
    
    # Calculate maximum EVSI
    max_evsi <- max(location_costs$evsi_hct_mean, location_costs$evsi_pcr_mean, 0)
    
    tibble(
      longitude = longitude,
      latitude = latitude,
      n_models = n_models,
      mean_prevalence = mean_prevalence,
      sd_prevalence = sd_prevalence,
      q025_prevalence = q025_prevalence,
      q975_prevalence = q975_prevalence,
      
      # Cost means and uncertainties
      cost_no_test_mean = location_costs$cost_no_test_mean,
      cost_no_test_sd = location_costs$cost_no_test_sd,
      
      cost_hct_mean = location_costs$cost_hct_mean,
      cost_hct_sd = location_costs$cost_hct_sd,
      evsi_hct_mean = location_costs$evsi_hct_mean,
      evsi_hct_sd = location_costs$evsi_hct_sd,
      
      cost_pcr_mean = location_costs$cost_pcr_mean,
      cost_pcr_sd = location_costs$cost_pcr_sd,
      evsi_pcr_mean = location_costs$evsi_pcr_mean,
      evsi_pcr_sd = location_costs$evsi_pcr_sd,
      
      # Strategy probabilities
      prob_no_test_optimal = location_costs$prob_no_test_optimal,
      prob_hct_optimal = location_costs$prob_hct_optimal,
      prob_pcr_optimal = location_costs$prob_pcr_optimal,
      
      # Optimal strategy (based on mean costs)
      optimal_strategy = optimal_strategy,
      min_cost_mean = min_cost,
      max_evsi_mean = max_evsi
    )
  })
  
  cat(" completed\n")
  return(cost_results)
}

# ==============================================================================
# MAIN COMPUTATION
# ==============================================================================

cat("Starting main computation across all Nigeria locations...\n")

# Create batches for memory efficiency
batch_size <- 1000
unique_locations <- all_projections %>%
  select(longitude, latitude) %>%
  distinct() %>%
  mutate(batch_id = ((row_number() - 1) %/% batch_size) + 1)

n_batches <- max(unique_locations$batch_id)
cat("Processing", nrow(unique_locations), "locations in", n_batches, "batches\n")

# Process each batch
all_results <- tibble()

for (batch_num in 1:n_batches) {
  cat("Processing batch", batch_num, "of", n_batches, "...")
  
  # Get locations for this batch
  batch_locations <- unique_locations %>%
    filter(batch_id == batch_num) %>%
    select(-batch_id)
  
  # Get data for these locations
  batch_data <- all_projections %>%
    inner_join(batch_locations, by = c("longitude", "latitude"))
  
  # Process batch
  batch_results <- process_batch(batch_data)
  
  # Combine results
  all_results <- bind_rows(all_results, batch_results)
  
  cat(" completed (", nrow(batch_results), "locations)\n")
}

cat("\nComputation completed!\n")
cat("Total locations analyzed:", nrow(all_results), "\n")

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

cat("\n=== SUMMARY RESULTS ===\n")

# Strategy distribution
strategy_summary <- all_results %>%
  count(optimal_strategy) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

cat("\nOptimal Strategy Distribution (based on mean costs):\n")
for (i in seq_len(nrow(strategy_summary))) {
  cat("•", strategy_summary$optimal_strategy[i], ":", 
      strategy_summary$percentage[i], "% (", strategy_summary$n[i], "locations)\n")
}

# Uncertainty in optimal strategies
cat("\nStrategy Probability Statistics (accounting for uncertainty):\n")
cat("• Mean probability NO_TEST optimal:", round(mean(all_results$prob_no_test_optimal) * 100, 1), "%\n")
cat("• Mean probability HCT optimal:", round(mean(all_results$prob_hct_optimal) * 100, 1), "%\n") 
cat("• Mean probability PCR optimal:", round(mean(all_results$prob_pcr_optimal) * 100, 1), "%\n")

# Test performance statistics
cat("\nFitted Test Performance (from Bayesian latent class model):\n")
cat("• HCT Se: mean =", round(test_params$hct$sensitivity_mean, 3), 
    ", SD =", round(sd(test_params$hct$sensitivity_samples), 3), "\n")
cat("• HCT Sp: mean =", round(test_params$hct$specificity_mean, 3),
    ", SD =", round(sd(test_params$hct$specificity_samples), 3), "\n")
cat("• PCR Se: mean =", round(test_params$pcr$sensitivity_mean, 3),
    ", SD =", round(sd(test_params$pcr$sensitivity_samples), 3), "\n") 
cat("• PCR Sp: mean =", round(test_params$pcr$specificity_mean, 3),
    ", SD =", round(sd(test_params$pcr$specificity_samples), 3), "\n")

# Prevalence and cost statistics
cat("\nPrevalence Statistics:\n")
cat("• Mean prevalence:", round(mean(all_results$mean_prevalence) * 100, 2), "%\n")
cat("• Prevalence range:", round(min(all_results$mean_prevalence) * 100, 2), "% to", 
    round(max(all_results$mean_prevalence) * 100, 2), "%\n")

cat("\nCost Statistics (with uncertainty):\n")
cat("• Mean cost (optimal strategy): $", round(mean(all_results$min_cost_mean), 2), 
    " ± $", round(mean(sqrt((all_results$cost_no_test_sd^2 + all_results$cost_hct_sd^2 + all_results$cost_pcr_sd^2)/3)), 2), "\n")
cat("• Total EVSI across Nigeria: $", round(sum(all_results$max_evsi_mean)), "\n")
cat("• Mean EVSI uncertainty: ± $", round(mean(sqrt((all_results$evsi_hct_sd^2 + all_results$evsi_pcr_sd^2)/2)), 2), "\n")

# Geographic distribution
cat("\nGeographic Coverage:\n")
cat("• Longitude range:", round(min(all_results$longitude), 2), "° to", 
    round(max(all_results$longitude), 2), "°E\n")
cat("• Latitude range:", round(min(all_results$latitude), 2), "° to", 
    round(max(all_results$latitude), 2), "°N\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\nSaving results...\n")

# Create output directory
dir.create("Code/Cost effectiveness analysis/results", showWarnings = FALSE)

# Save main results
write_csv(all_results, "Code/Cost effectiveness analysis/results/nigeria_cost_effectiveness_results.csv")

# Save summary statistics
summary_stats <- list(
  strategy_distribution = strategy_summary,
  prevalence_stats = list(
    mean = mean(all_results$mean_prevalence),
    min = min(all_results$mean_prevalence),
    max = max(all_results$mean_prevalence),
    sd = sd(all_results$mean_prevalence)
  ),
  cost_stats = list(
    mean_optimal_cost = mean(all_results$min_cost_mean),
    total_evsi = sum(all_results$max_evsi_mean),
    mean_evsi = mean(all_results$max_evsi_mean),
    mean_cost_uncertainty = mean(sqrt((all_results$cost_no_test_sd^2 + all_results$cost_hct_sd^2 + all_results$cost_pcr_sd^2)/3)),
    mean_evsi_uncertainty = mean(sqrt((all_results$evsi_hct_sd^2 + all_results$evsi_pcr_sd^2)/2))
  ),
  strategy_probabilities = list(
    mean_prob_no_test = mean(all_results$prob_no_test_optimal),
    mean_prob_hct = mean(all_results$prob_hct_optimal), 
    mean_prob_pcr = mean(all_results$prob_pcr_optimal)
  ),
  geographic_extent = list(
    lon_min = min(all_results$longitude),
    lon_max = max(all_results$longitude), 
    lat_min = min(all_results$latitude),
    lat_max = max(all_results$latitude)
  )
)

# Save as RDS for easy loading in visualization script
saveRDS(summary_stats, "Code/Cost effectiveness analysis/results/nigeria_summary_statistics.rds")

cat("Results saved to 'Code/Cost effectiveness analysis/results/' directory\n")
cat("• nigeria_cost_effectiveness_results.csv (", nrow(all_results), "locations)\n")
cat("• nigeria_summary_statistics.rds (summary statistics)\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Ready for visualization!\n")