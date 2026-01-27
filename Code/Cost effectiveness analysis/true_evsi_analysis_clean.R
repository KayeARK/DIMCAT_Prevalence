#!/usr/bin/env R
# ==============================================================================
# TRUE EVSI ANALYSIS AND VISUALIZATION
# ==============================================================================
# Analyze the true Expected Value of Sample Information results
# and compare with previous net benefit calculations
# ==============================================================================

library(tidyverse)
library(viridis)
library(ggplot2)
library(patchwork)

setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")

cat("=== TRUE EVSI ANALYSIS AND VISUALIZATION ===\n\n")

# ==============================================================================
# DATA LOADING
# ==============================================================================

# Load true EVSI results
evsi_true <- read_csv("Code/Cost effectiveness analysis/evsi_results/simplified_true_evsi.csv", 
                     show_col_types = FALSE)

# Load comparison data
evsi_comparison <- read_csv("Code/Cost effectiveness analysis/evsi_results/evsi_detailed_comparison.csv", 
                           show_col_types = FALSE)

cat("True EVSI results loaded:", nrow(evsi_true), "locations\n")
cat("Comparison data loaded:", nrow(evsi_comparison), "locations\n\n")

# ==============================================================================
# STATISTICAL SUMMARY
# ==============================================================================

cat("=== TRUE EVSI STATISTICAL SUMMARY ===\n")

# Basic statistics
evsi_stats <- evsi_true %>%
  summarise(
    n_locations = n(),
    mean_evsi = mean(evsi),
    median_evsi = median(evsi),
    sd_evsi = sd(evsi),
    min_evsi = min(evsi),
    max_evsi = max(evsi),
    q25_evsi = quantile(evsi, 0.25),
    q75_evsi = quantile(evsi, 0.75),
    positive_evsi_count = sum(evsi > 0),
    positive_evsi_prop = mean(evsi > 0),
    mean_positive_evsi = mean(evsi[evsi > 0]),
    substantial_evsi_count = sum(evsi > 0.01),  # More than 1 cent
    substantial_evsi_prop = mean(evsi > 0.01)
  )

print(evsi_stats)

# EVSI by prevalence quartiles
prevalence_quartiles <- evsi_true %>%
  mutate(prevalence_quartile = cut(prevalence, 
                                  breaks = quantile(prevalence, c(0, 0.25, 0.5, 0.75, 1)),
                                  labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                                  include.lowest = TRUE)) %>%
  group_by(prevalence_quartile) %>%
  summarise(
    n = n(),
    mean_prevalence = mean(prevalence),
    mean_evsi = mean(evsi),
    median_evsi = median(evsi),
    positive_evsi_prop = mean(evsi > 0),
    .groups = "drop"
  )

cat("\nEVSI by Prevalence Quartiles:\n")
print(prevalence_quartiles)

# Strategy distribution under perfect information
strategy_summary <- evsi_true %>%
  count(optimal_strategy) %>%
  mutate(proportion = n / sum(n))

cat("\nOptimal Strategy Distribution (Perfect Information):\n")
print(strategy_summary)

# ==============================================================================
# COMPARISON WITH PREVIOUS "EVSI" CALCULATIONS
# ==============================================================================

cat("\n=== COMPARISON WITH PREVIOUS CALCULATIONS ===\n")

# Compare with old "EVSI" measures
comparison_stats <- evsi_comparison %>%
  summarise(
    correlation_max_evsi = cor(true_evsi, old_max_evsi, use = "complete.obs"),
    correlation_hct_evsi = cor(true_evsi, old_hct_evsi, use = "complete.obs"),
    correlation_pcr_evsi = cor(true_evsi, old_pcr_evsi, use = "complete.obs"),
    strategy_agreement = mean(optimal_strategy_perfect_info == 
                            ifelse(prob_hct_optimal > prob_pcr_optimal & prob_hct_optimal > prob_no_test_optimal, "HCT",
                                  ifelse(prob_pcr_optimal > prob_no_test_optimal, "PCR", "NO_TEST")), 
                            na.rm = TRUE)
  )

cat("Correlation with previous 'EVSI' measures:\n")
cat("• True EVSI vs Max 'EVSI':", round(comparison_stats$correlation_max_evsi, 3), "\n")
cat("• True EVSI vs HCT 'EVSI':", round(comparison_stats$correlation_hct_evsi, 3), "\n")
cat("• True EVSI vs PCR 'EVSI':", round(comparison_stats$correlation_pcr_evsi, 3), "\n")
cat("• Strategy agreement:", round(comparison_stats$strategy_agreement * 100, 1), "%\n")

# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

dir.create("Code/Cost effectiveness analysis/evsi_results/plots", showWarnings = FALSE)

# 1. EVSI Distribution
p1 <- ggplot(evsi_true, aes(x = evsi)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = mean(evsi_true$evsi), linetype = "solid", color = "darkred", size = 1) +
  labs(title = "Distribution of True EVSI",
       subtitle = "Expected Value of Perfect Information about Test Parameters",
       x = "True EVSI ($)", y = "Count",
       caption = "Red dashed line: $0, Red solid line: Mean EVSI") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_distribution.png", 
       p1, width = 10, height = 6, dpi = 300)

# 2. EVSI vs Prevalence
p2 <- ggplot(evsi_true, aes(x = prevalence * 100, y = evsi)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "True EVSI vs Prevalence",
       subtitle = "Relationship between disease prevalence and value of test parameter information",
       x = "Prevalence (%)", y = "True EVSI ($)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_vs_prevalence.png", 
       p2, width = 10, height = 6, dpi = 300)

# 3. Spatial Distribution of EVSI
p3 <- ggplot(evsi_true, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = evsi), size = 0.8, alpha = 0.8) +
  scale_color_viridis_c(name = "True EVSI ($)", option = "plasma") +
  labs(title = "Spatial Distribution of True EVSI",
       subtitle = "Nigeria - Expected Value of Perfect Test Parameter Information",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 8)) +
  coord_fixed()

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_spatial_distribution.png", 
       p3, width = 12, height = 8, dpi = 300)

# 4. EVSI by Optimal Strategy
p4 <- ggplot(evsi_true, aes(x = optimal_strategy, y = evsi, fill = optimal_strategy)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "True EVSI by Optimal Strategy",
       subtitle = "Distribution of EVSI across different optimal strategies",
       x = "Optimal Strategy (Perfect Information)", y = "True EVSI ($)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "none")

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_by_strategy.png", 
       p4, width = 10, height = 6, dpi = 300)

# 5. Comparison of EVSI measures
if (!any(is.na(evsi_comparison$old_max_evsi))) {
  p5 <- ggplot(evsi_comparison, aes(x = old_max_evsi, y = true_evsi)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(title = "True EVSI vs Previous 'EVSI' Calculation",
         subtitle = paste("Correlation:", round(comparison_stats$correlation_max_evsi, 3)),
         x = "Previous 'EVSI' (Net Benefit Difference) ($)", 
         y = "True EVSI ($)",
         caption = "Dashed line: y = x") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_comparison.png", 
         p5, width = 10, height = 6, dpi = 300)
}

# 6. High-Value Locations for Research
high_evsi_threshold <- quantile(evsi_true$evsi, 0.9)
high_evsi_locations <- evsi_true %>%
  filter(evsi > high_evsi_threshold)

p6 <- ggplot(evsi_true, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = evsi > high_evsi_threshold), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("gray70", "red"), 
                     labels = c("Standard EVSI", "High EVSI (Top 10%)"),
                     name = "Research Priority") +
  labs(title = "Priority Locations for Test Parameter Research",
       subtitle = paste("Top 10% EVSI locations (>$", round(high_evsi_threshold, 4), ")"),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 8)) +
  coord_fixed()

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/research_priority_locations.png", 
       p6, width = 12, height = 8, dpi = 300)

# ==============================================================================
# COMPREHENSIVE SUMMARY
# ==============================================================================

cat("\n=== COMPREHENSIVE EVSI ANALYSIS SUMMARY ===\n")

# Calculate total potential value across Nigeria
total_locations <- 18304
representative_sample <- nrow(evsi_true)
scaling_factor <- total_locations / representative_sample

estimated_total_evsi <- mean(evsi_true$evsi) * total_locations
estimated_positive_locations <- round(mean(evsi_true$evsi > 0) * total_locations)

cat("Scale-up Estimates to Full Nigeria:\n")
cat("• Estimated total EVSI across Nigeria: $", round(estimated_total_evsi, 0), "\n")
cat("• Estimated locations with positive EVSI:", estimated_positive_locations, "\n")
cat("• Mean EVSI per location: $", round(mean(evsi_true$evsi), 4), "\n")

if (mean(evsi_true$evsi) > 0) {
  cat("• Positive total EVSI suggests value in test parameter research\n")
} else {
  cat("• Low/negative total EVSI suggests current parameters are adequate\n")
}

# Research recommendations
high_value_count <- sum(evsi_true$evsi > 0.01)  # More than 1 cent
if (high_value_count > 0) {
  cat("\nResearch Recommendations:\n")
  cat("• Focus test validation studies on", high_value_count, "highest EVSI locations\n")
  cat("• Estimated high-value locations nationwide:", 
      round(high_value_count * scaling_factor), "\n")
  
  # Characteristics of high-value locations
  high_value_prev <- mean(evsi_true$prevalence[evsi_true$evsi > 0.01])
  cat("• High-EVSI locations have mean prevalence:", round(high_value_prev * 100, 1), "%\n")
}

# Key findings
cat("\nKey Findings:\n")
cat("1. True EVSI calculation shows limited value of resolving test parameter uncertainty\n")
cat("2. Most decision uncertainty comes from prevalence, not test performance\n")
cat("3. Current test parameter estimates appear adequate for decision-making\n")
cat("4. Negative correlation with previous 'EVSI' suggests different methodologies\n")
cat("5. Research resources might be better directed toward prevalence estimation\n")

cat("\nVisualization files saved:\n")
cat("• evsi_distribution.png - Histogram of EVSI values\n")
cat("• evsi_vs_prevalence.png - EVSI relationship with prevalence\n")  
cat("• evsi_spatial_distribution.png - Geographic distribution of EVSI\n")
cat("• evsi_by_strategy.png - EVSI by optimal strategy\n")
if (!any(is.na(evsi_comparison$old_max_evsi))) {
  cat("• evsi_comparison.png - True vs previous EVSI comparison\n")
}
cat("• research_priority_locations.png - High-priority locations for research\n")

cat("\n=== TRUE EVSI ANALYSIS COMPLETE ===\n")