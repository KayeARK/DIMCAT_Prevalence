#!/usr/bin/env R
# ==============================================================================
# DECISION-FOCUSED OPTIMAL SAMPLING DESIGN
# ==============================================================================
# Identifies locations where prevalence uncertainty affects testing strategy decisions
# Focuses on Expected Value of Perfect Information for DECISIONS rather than costs
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(viridis)

setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")

cat("=== DECISION-FOCUSED SAMPLING DESIGN ===\n")
cat("Identifying locations with decision uncertainty due to prevalence uncertainty...\n\n")

# Load results
results <- read_csv("Code/Cost effectiveness analysis/results/nigeria_cost_effectiveness_results.csv", 
                   show_col_types = FALSE)

cat("Loaded", nrow(results), "locations from existing analysis\n")

# ==============================================================================
# DECISION UNCERTAINTY ANALYSIS
# ==============================================================================

# Locations with decision uncertainty are those where:
# 1. No strategy has >95% probability of being optimal
# 2. At least two strategies have >5% probability of being optimal
# 3. The credible interval of prevalence spans decision boundaries

decision_uncertainty <- results %>%
  mutate(
    # Decision uncertainty metrics
    max_strategy_prob = pmax(prob_no_test_optimal, prob_hct_optimal, prob_pcr_optimal),
    decision_entropy = -((prob_no_test_optimal * log(prob_no_test_optimal + 1e-10)) +
                        (prob_hct_optimal * log(prob_hct_optimal + 1e-10)) +
                        (prob_pcr_optimal * log(prob_pcr_optimal + 1e-10))),
    
    # Prevalence uncertainty (coefficient of variation)
    prevalence_cv = sd_prevalence / mean_prevalence,
    
    # High uncertainty flags
    high_decision_uncertainty = max_strategy_prob < 0.95,
    moderate_decision_uncertainty = max_strategy_prob < 0.90,
    
    # Expected value of perfect information for decisions
    # Based on improvement in decision confidence
    decision_evpi = case_when(
      max_strategy_prob < 0.80 ~ max_evsi_mean * (1 - max_strategy_prob) * 2,
      max_strategy_prob < 0.90 ~ max_evsi_mean * (1 - max_strategy_prob) * 1.5,
      max_strategy_prob < 0.95 ~ max_evsi_mean * (1 - max_strategy_prob) * 1,
      TRUE ~ 0
    ),
    
    # Simplified distance to major cities
    distance_to_city = pmin(
      sqrt((latitude - 6.5244)^2 + (longitude - 3.3792)^2),  # Lagos
      sqrt((latitude - 9.0765)^2 + (longitude - 7.3986)^2),  # Abuja
      sqrt((latitude - 12.0022)^2 + (longitude - 8.5920)^2)  # Kano
    ) * 111  # Convert to km
  )

cat("Decision Uncertainty Statistics:\n")
cat("• Locations with high decision uncertainty (max prob <95%):", 
    sum(decision_uncertainty$high_decision_uncertainty), "\n")
cat("• Locations with moderate decision uncertainty (max prob <90%):", 
    sum(decision_uncertainty$moderate_decision_uncertainty), "\n")
cat("• Mean decision entropy:", round(mean(decision_uncertainty$decision_entropy), 3), "\n")
cat("• Max decision EVPI: $", round(max(decision_uncertainty$decision_evpi), 2), "\n\n")

# ==============================================================================
# SAMPLING VALUE CALCULATION
# ==============================================================================

# Sampling costs
sampling_cost_per_location <- function(sample_size, distance_km) {
  cost_per_sample <- 5.0
  travel_cost <- 50.0 + 0.5 * distance_km
  return(cost_per_sample * sample_size + travel_cost)
}

# Calculate optimal sampling for uncertain locations
uncertain_locations <- decision_uncertainty %>%
  filter(decision_evpi > 1.0) %>%  # Only locations with meaningful decision EVPI
  arrange(desc(decision_evpi)) %>%
  slice_head(n = 1000)  # Focus on top 1000 for computational efficiency

cat("Analyzing", nrow(uncertain_locations), "locations with meaningful decision uncertainty...\n")

# Calculate sampling value for different sample sizes
if(nrow(uncertain_locations) > 0) {
  sample_sizes <- c(10, 20, 30, 50, 100)
  
  sampling_analysis <- uncertain_locations %>%
    slice_head(n = 200) %>%  # Limit for speed
    crossing(sample_size = sample_sizes) %>%
    mutate(
      # Estimate how much sampling would reduce decision uncertainty
      # Assume current uncertainty is like having effective sample size based on CV
      current_effective_n = pmin(100, 1 / (prevalence_cv^2)),
      new_effective_n = current_effective_n + sample_size,
      
      # New uncertainty after sampling
      new_cv = sqrt(1 / new_effective_n),
      uncertainty_reduction = 1 - (new_cv / prevalence_cv),
      
      # Sampling improves decision confidence proportional to uncertainty reduction
      decision_improvement = decision_evpi * uncertainty_reduction * 0.5,
      
      # Sampling costs
      sampling_cost = map2_dbl(sample_size, distance_to_city, sampling_cost_per_location),
      
      # Net value
      net_sampling_value = decision_improvement - sampling_cost
    )
  
  # Find optimal sample size per location
  optimal_sampling <- sampling_analysis %>%
    group_by(longitude, latitude) %>%
    filter(net_sampling_value == max(net_sampling_value)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    filter(net_sampling_value > 0)  # Only profitable sampling
  
  cat("Found", nrow(optimal_sampling), "locations where sampling is cost-effective\n")
} else {
  optimal_sampling <- tibble()
  cat("No locations found with sufficient decision uncertainty\n")
}

# ==============================================================================
# ALTERNATIVE SAMPLING PRIORITIES
# ==============================================================================

# Even if formal EVPI is low, identify locations that are scientifically interesting
# for sampling based on other criteria

scientific_priority <- decision_uncertainty %>%
  mutate(
    # Locations near decision boundaries (close to switching strategies)
    near_hct_pcr_boundary = (prob_hct_optimal > 0.2 & prob_pcr_optimal > 0.2),
    near_notest_boundary = (prob_no_test_optimal > 0.1),
    
    # High prevalence uncertainty
    high_prevalence_uncertainty = prevalence_cv > quantile(prevalence_cv, 0.9, na.rm = TRUE),
    
    # Extreme prevalence values (high or low)
    extreme_prevalence = (mean_prevalence > quantile(mean_prevalence, 0.95, na.rm = TRUE)) |
                        (mean_prevalence < quantile(mean_prevalence, 0.05, na.rm = TRUE)),
    
    # Scientific priority score
    scientific_score = (near_hct_pcr_boundary * 3) + 
                      (near_notest_boundary * 2) +
                      (high_prevalence_uncertainty * 2) +
                      (extreme_prevalence * 1) +
                      (decision_entropy > median(decision_entropy, na.rm = TRUE)) * 1,
    
    # Accessibility (closer to cities = easier to sample)
    accessibility_score = case_when(
      distance_to_city < 50 ~ 3,
      distance_to_city < 100 ~ 2,
      distance_to_city < 200 ~ 1,
      TRUE ~ 0
    ),
    
    # Combined priority score
    sampling_priority = scientific_score + accessibility_score
  ) %>%
  filter(sampling_priority >= 4) %>%  # High priority locations
  arrange(desc(sampling_priority))

cat("\nScientific Sampling Priorities:\n")
cat("• Locations near HCT-PCR decision boundary:", sum(scientific_priority$near_hct_pcr_boundary), "\n")
cat("• Locations near NO_TEST boundary:", sum(scientific_priority$near_notest_boundary), "\n") 
cat("• High priority sampling locations:", nrow(scientific_priority), "\n")

# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

cat("\nCreating visualizations...\n")

# Map of decision uncertainty
p_uncertainty <- ggplot(decision_uncertainty, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = decision_entropy, size = decision_evpi), alpha = 0.7) +
  scale_color_viridis_c(name = "Decision\nEntropy", option = "plasma") +
  scale_size_continuous(name = "Decision\nEVPI ($)", range = c(0.5, 3)) +
  labs(title = "Decision Uncertainty Across Nigeria",
       subtitle = "Larger points show higher value of resolving prevalence uncertainty",
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_minimal() +
  coord_fixed(ratio = 1/1.238)

# Scientific priority map
p_priority <- ggplot(scientific_priority, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = scientific_score, size = accessibility_score), alpha = 0.8) +
  scale_color_viridis_c(name = "Scientific\nScore", option = "viridis") +
  scale_size_continuous(name = "Accessibility\nScore", range = c(1, 4)) +
  labs(title = "Scientific Sampling Priorities",
       subtitle = paste0("Top ", nrow(scientific_priority), " locations for additional prevalence sampling"),
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_minimal() +
  coord_fixed(ratio = 1/1.238)

# Strategy probability uncertainty
p_strategy_uncertainty <- decision_uncertainty %>%
  filter(max_strategy_prob < 0.98) %>%
  ggplot(aes(x = latitude, y = longitude)) +
  geom_point(aes(color = optimal_strategy, alpha = 1 - max_strategy_prob), size = 1) +
  scale_color_manual(values = c("NO_TEST" = "#E31A1C", "HCT" = "#1F78B4", "PCR" = "#33A02C"),
                    name = "Most Likely\nStrategy") +
  scale_alpha_continuous(name = "Decision\nUncertainty", range = c(0.3, 1)) +
  labs(title = "Locations with Strategy Decision Uncertainty",
       subtitle = "More transparent points indicate higher decision uncertainty",
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_minimal() +
  coord_fixed(ratio = 1/1.238)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("Saving results...\n")

dir.create("Code/Cost effectiveness analysis/sampling_design", showWarnings = FALSE)

# Save data
write_csv(decision_uncertainty, 
          "Code/Cost effectiveness analysis/sampling_design/decision_uncertainty_analysis.csv")

if(nrow(optimal_sampling) > 0) {
  write_csv(optimal_sampling,
            "Code/Cost effectiveness analysis/sampling_design/optimal_sampling_economic.csv")
}

write_csv(scientific_priority,
          "Code/Cost effectiveness analysis/sampling_design/scientific_sampling_priorities.csv")

# Save plots
ggsave("Code/Cost effectiveness analysis/sampling_design/decision_uncertainty_map.png",
       p_uncertainty, width = 12, height = 10, dpi = 300)

ggsave("Code/Cost effectiveness analysis/sampling_design/scientific_priority_map.png",
       p_priority, width = 12, height = 10, dpi = 300)

ggsave("Code/Cost effectiveness analysis/sampling_design/strategy_uncertainty_map.png",
       p_strategy_uncertainty, width = 12, height = 10, dpi = 300)

# ==============================================================================
# FINAL REPORT
# ==============================================================================

cat("\n=== DECISION-FOCUSED SAMPLING DESIGN COMPLETE ===\n")

cat("\nKey Findings:\n")

if(nrow(optimal_sampling) > 0) {
  cat("ECONOMIC SAMPLING RECOMMENDATIONS:\n")
  cat("•", nrow(optimal_sampling), "locations where sampling has positive economic value\n")
  cat("• Mean net value: $", round(mean(optimal_sampling$net_sampling_value), 2), "\n")
  cat("• Total potential value: $", round(sum(optimal_sampling$net_sampling_value)), "\n")
} else {
  cat("ECONOMIC RESULT: No locations where additional sampling is economically justified\n")
  cat("• This suggests current prevalence estimates are sufficient for decision-making\n")
  cat("• The diagnostic testing strategies are robust to prevalence uncertainty\n")
}

cat("\nSCIENTIFIC SAMPLING RECOMMENDATIONS:\n")
cat("•", nrow(scientific_priority), "high-priority locations for scientific sampling\n")
cat("• Focus on locations near decision boundaries or with high uncertainty\n")
cat("• Mean scientific score:", round(mean(scientific_priority$scientific_score), 1), "\n")

cat("\nDECISION ROBUSTNESS:\n")
cat("• Locations with decision uncertainty (<95% confidence):", 
    sum(decision_uncertainty$high_decision_uncertainty), "of", nrow(decision_uncertainty), "\n")
cat("• This represents", round(100 * mean(decision_uncertainty$high_decision_uncertainty), 1), 
    "% of all locations\n")

cat("\nRECOMMENDATION:\n")
if(nrow(scientific_priority) > 0) {
  cat("• Focus additional sampling on the", min(50, nrow(scientific_priority)), 
      "highest scientific priority locations\n")
  cat("• These provide the most learning about decision boundaries\n")
  cat("• Even if not economically justified, they improve model validation\n")
} else {
  cat("• Current prevalence estimates appear sufficient for policy decisions\n")
  cat("• The cost-effectiveness framework is robust to existing uncertainty\n")
}

cat("\nFiles saved to 'Code/Cost effectiveness analysis/sampling_design/' directory!\n")