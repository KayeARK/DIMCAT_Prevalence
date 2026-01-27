#!/usr/bin/env R
# ==============================================================================
# SIMPLIFIED TRUE EVSI CALCULATION
# ==============================================================================
# Focuses on the Expected Value of Perfect Information about Test Parameters
# This is computationally tractable and addresses the main uncertainty source
# ==============================================================================

library(tidyverse)
library(rstan)

setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")

cat("=== SIMPLIFIED TRUE EVSI CALCULATION ===\n")
cat("Expected Value of Perfect Information about Test Performance Parameters\n\n")

# ==============================================================================
# SETUP
# ==============================================================================

# Load existing results
existing_results <- read_csv("Code/Cost effectiveness analysis/results/nigeria_cost_effectiveness_results.csv", 
                           show_col_types = FALSE)

# Load test parameters
fitted_model <- readRDS("Code/TestSensSpec/latent_class_fit.rds")
se_hct_samples <- extract(fitted_model, 'Se_hct')[[1]]
sp_hct_samples <- extract(fitted_model, 'Sp_hct')[[1]]
se_pcr_samples <- extract(fitted_model, 'Se_pcr')[[1]]  
sp_pcr_samples <- extract(fitted_model, 'Sp_pcr')[[1]]

# Cost parameters (fixed for simplicity)
treatment_cost <- 2.5
cost_untreated <- 50.0
hct_cost <- 0.5
pcr_cost <- 15.0

n_samples <- 1000

cat("Using", length(se_hct_samples), "test parameter samples\n")
cat("Analyzing", nrow(existing_results), "locations\n\n")

# ==============================================================================
# SIMPLIFIED EVSI CALCULATION
# ==============================================================================

calculate_expected_cost <- function(prevalence, sensitivity, specificity, test_cost, strategy) {
  
  if (strategy == "NO_TEST") {
    return(prevalence * cost_untreated)
  }
  
  # Testing strategy
  prob_positive <- sensitivity * prevalence + (1 - specificity) * (1 - prevalence)
  
  cost_testing <- test_cost
  cost_treatment <- prob_positive * treatment_cost
  cost_false_negative <- (1 - sensitivity) * prevalence * cost_untreated
  cost_false_positive <- (1 - specificity) * (1 - prevalence) * treatment_cost
  
  return(cost_testing + cost_treatment + cost_false_negative + cost_false_positive)
}

# EVSI calculation for one location
calculate_location_evsi <- function(location_data, n_samples = 1000) {
  
  prevalence <- location_data$mean_prevalence
  
  # Sample test parameters from joint posterior
  sample_indices <- sample(seq_along(se_hct_samples), n_samples, replace = TRUE)
  
  # STEP 1: Current Information Value (uncertain test parameters)
  # For each sample of test parameters, find optimal strategy
  current_costs <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    idx <- sample_indices[i]
    
    # Get parameter values for this sample
    se_hct <- se_hct_samples[idx]
    sp_hct <- sp_hct_samples[idx]
    se_pcr <- se_pcr_samples[idx]
    sp_pcr <- sp_pcr_samples[idx]
    
    # Calculate costs for each strategy
    cost_no_test <- calculate_expected_cost(prevalence, NA, NA, 0, "NO_TEST")
    cost_hct <- calculate_expected_cost(prevalence, se_hct, sp_hct, hct_cost, "HCT")
    cost_pcr <- calculate_expected_cost(prevalence, se_pcr, sp_pcr, pcr_cost, "PCR")
    
    # Choose optimal strategy for this parameter sample
    current_costs[i] <- min(cost_no_test, cost_hct, cost_pcr)
  }
  
  current_expected_cost <- mean(current_costs)
  
  # STEP 2: Perfect Information Value (known test parameters)
  # Use expected values of test parameters
  se_hct_mean <- mean(se_hct_samples)
  sp_hct_mean <- mean(sp_hct_samples)
  se_pcr_mean <- mean(se_pcr_samples)
  sp_pcr_mean <- mean(sp_pcr_samples)
  
  # Calculate costs with known (mean) parameters
  cost_no_test_perfect <- calculate_expected_cost(prevalence, NA, NA, 0, "NO_TEST")
  cost_hct_perfect <- calculate_expected_cost(prevalence, se_hct_mean, sp_hct_mean, hct_cost, "HCT")
  cost_pcr_perfect <- calculate_expected_cost(prevalence, se_pcr_mean, sp_pcr_mean, pcr_cost, "PCR")
  
  perfect_expected_cost <- min(cost_no_test_perfect, cost_hct_perfect, cost_pcr_perfect)
  
  # STEP 3: Calculate EVSI
  evsi <- current_expected_cost - perfect_expected_cost
  
  # Additional metrics
  optimal_strategy_perfect <- which.min(c(cost_no_test_perfect, cost_hct_perfect, cost_pcr_perfect))
  strategy_names <- c("NO_TEST", "HCT", "PCR")
  
  return(list(
    longitude = location_data$longitude,
    latitude = location_data$latitude,
    prevalence = prevalence,
    evsi = evsi,
    current_cost = current_expected_cost,
    perfect_cost = perfect_expected_cost,
    optimal_strategy_perfect = strategy_names[optimal_strategy_perfect],
    cost_no_test = cost_no_test_perfect,
    cost_hct = cost_hct_perfect,
    cost_pcr = cost_pcr_perfect
  ))
}

# ==============================================================================
# BATCH PROCESSING
# ==============================================================================

cat("Processing locations in batches...\n")

# Process subset for testing
test_indices <- seq(1, nrow(existing_results), by = 100)  # Every 100th location
cat("Testing on", length(test_indices), "representative locations\n")

# Process test locations
evsi_results <- map_dfr(test_indices, function(idx) {
  if (idx %% 1000 == 1) cat("Processing location", idx, "\n")
  
  location_data <- existing_results[idx, ]
  result <- calculate_location_evsi(location_data, n_samples = n_samples)
  
  data.frame(
    location_idx = idx,
    longitude = result$longitude,
    latitude = result$latitude,
    prevalence = result$prevalence,
    evsi = result$evsi,
    current_cost = result$current_cost,
    perfect_cost = result$perfect_cost,
    optimal_strategy = result$optimal_strategy_perfect,
    cost_no_test = result$cost_no_test,
    cost_hct = result$cost_hct,
    cost_pcr = result$cost_pcr
  )
})

# ==============================================================================
# RESULTS ANALYSIS
# ==============================================================================

cat("\n=== TRUE EVSI RESULTS ===\n")
cat("Locations analyzed:", nrow(evsi_results), "\n")
cat("Mean EVSI: $", round(mean(evsi_results$evsi), 4), "\n")
cat("Median EVSI: $", round(median(evsi_results$evsi), 4), "\n")
cat("EVSI range: $", round(min(evsi_results$evsi), 4), " to $", round(max(evsi_results$evsi), 4), "\n")
cat("Positive EVSI locations:", sum(evsi_results$evsi > 0), "out of", nrow(evsi_results), "\n")

if (sum(evsi_results$evsi > 0) > 0) {
  cat("Mean positive EVSI: $", round(mean(evsi_results$evsi[evsi_results$evsi > 0]), 4), "\n")
}

# Strategy distribution under perfect information
cat("\nOptimal strategies under perfect information:\n")
strategy_dist <- table(evsi_results$optimal_strategy)
for (i in seq_along(strategy_dist)) {
  cat("•", names(strategy_dist)[i], ":", strategy_dist[i], "locations\n")
}

# Compare with current analysis results
comparison_subset <- existing_results[test_indices, ] %>%
  select(longitude, latitude, optimal_strategy, max_evsi_mean) %>%
  left_join(evsi_results, by = c("longitude", "latitude"), suffix = c("_old", "_new"))

cat("\nComparison with existing analysis:\n")
cat("Strategy agreement:", 
    sum(comparison_subset$optimal_strategy_old == comparison_subset$optimal_strategy_new), 
    "out of", nrow(comparison_subset), "locations\n")

# Correlation between old "EVSI" and true EVSI
if (var(evsi_results$evsi) > 0) {
  correlation <- cor(comparison_subset$max_evsi_mean, evsi_results$evsi)
  cat("Correlation between old 'EVSI' and true EVSI:", round(correlation, 3), "\n")
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

dir.create("Code/Cost effectiveness analysis/evsi_results", showWarnings = FALSE)

# Save detailed results
write_csv(evsi_results, "Code/Cost effectiveness analysis/evsi_results/simplified_true_evsi.csv")

# Save comparison
comparison_data <- evsi_results %>%
  left_join(existing_results[test_indices, ], by = c("longitude", "latitude")) %>%
  select(longitude, latitude, prevalence, evsi, current_cost, perfect_cost, 
         optimal_strategy.x, max_evsi_mean, evsi_hct_mean, evsi_pcr_mean,
         prob_hct_optimal, prob_pcr_optimal, prob_no_test_optimal) %>%
  rename(true_evsi = evsi,
         optimal_strategy_perfect_info = optimal_strategy.x,
         old_max_evsi = max_evsi_mean,
         old_hct_evsi = evsi_hct_mean,
         old_pcr_evsi = evsi_pcr_mean)

write_csv(comparison_data, "Code/Cost effectiveness analysis/evsi_results/evsi_detailed_comparison.csv")

cat("\nResults saved:\n")
cat("• simplified_true_evsi.csv - True EVSI calculations\n")
cat("• evsi_detailed_comparison.csv - Comparison with existing analysis\n")

# Quick visualization
if (nrow(evsi_results) > 10) {
  library(ggplot2)
  
  p1 <- ggplot(evsi_results, aes(x = prevalence * 100, y = evsi)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE) +
    labs(title = "True EVSI vs Prevalence",
         subtitle = "Expected Value of Perfect Information about Test Parameters",
         x = "Prevalence (%)", y = "True EVSI ($)") +
    theme_minimal()
  
  ggsave("Code/Cost effectiveness analysis/evsi_results/true_evsi_vs_prevalence.png", 
         p1, width = 10, height = 8, dpi = 300)
  
  cat("• true_evsi_vs_prevalence.png - Visualization\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("True EVSI calculation completed successfully!\n")

# Summary insights
cat("\nKey Insights:\n")
if (mean(evsi_results$evsi) > 0) {
  cat("• Positive EVSI indicates value in resolving test parameter uncertainty\n")
  cat("• Mean EVSI of $", round(mean(evsi_results$evsi), 2), " per location suggests research value\n")
} else {
  cat("• Low/negative EVSI suggests test parameter uncertainty has minimal impact\n")
  cat("• Current parameter estimates appear sufficient for decision making\n")
}

if (var(evsi_results$evsi) > 0) {
  high_evsi_locations <- sum(evsi_results$evsi > quantile(evsi_results$evsi, 0.9))
  cat("• ", high_evsi_locations, " locations have particularly high EVSI (top 10%)\n")
  cat("• These represent priority locations for test validation studies\n")
}