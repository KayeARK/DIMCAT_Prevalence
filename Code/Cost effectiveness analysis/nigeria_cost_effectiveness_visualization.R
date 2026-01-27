#!/usr/bin/env R
# ==============================================================================
# NIGERIA COST-EFFECTIVENESS ANALYSIS - VISUALIZATION
# ==============================================================================
# Comprehensive visualization suite for Bayesian cost-effectiveness results
# Generates maps, plocat("Files Generated:\n")
cat("INDIVIDUAL PLOTS (12):\n")
cat("• 01_optimal_strategy_map.png - Geographic distribution of optimal strategies\n")
cat("• 02_prevalence_map.png - Prevalence estimates across Nigeria\n")
cat("• 03_evsi_map.png - Expected Value of Sample Information\n")
cat("• 04_cost_map.png - Minimum expected costs\n")
cat("• 05_strategy_distribution.png - Comparison of point estimates vs uncertainty-aware probabilities\n")
cat("• 06_prevalence_evsi_relationship.png - Scatter plot with trend\n")
cat("• 07_cost_distribution.png - Violin plots by strategy\n")
cat("• 08_prevalence_histogram.png - Prevalence distribution\n")
cat("• 09_evsi_histogram.png - EVSI distribution\n")
cat("• 10_decision_uncertainty_map.png - Geographic distribution of decision uncertainty\n")
cat("• 11_cost_uncertainty_map.png - Geographic distribution of cost uncertainty\n")
cat("• 12_prevalence_uncertainty_map.png - Geographic distribution of prevalence uncertainty\n\n")

# ==============================================================================

library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)

# Set working directory
setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")

cat("=== NIGERIA COST-EFFECTIVENESS VISUALIZATION ===\n")
cat("Loading results and generating comprehensive plots...\n\n")

# ==============================================================================
# LOAD RESULTS
# ==============================================================================

# Check if results exist
if (!file.exists("Code/Cost effectiveness analysis/results/nigeria_cost_effectiveness_results.csv")) {
  stop("Results not found! Please run the computation script first.")
}

# Load main results
results <- read_csv("Code/Cost effectiveness analysis/results/nigeria_cost_effectiveness_results.csv", show_col_types = FALSE)
summary_stats <- readRDS("Code/Cost effectiveness analysis/results/nigeria_summary_statistics.rds")

cat("Loaded results for", nrow(results), "locations\n")

# ==============================================================================
# GEOGRAPHIC SETUP
# ==============================================================================

# Nigeria geographic bounds
nigeria_bounds <- list(
  lon_min = min(results$longitude),
  lon_max = max(results$longitude),
  lat_min = min(results$latitude), 
  lat_max = max(results$latitude)
)

# Calculate aspect ratio for proper display
nigeria_width <- nigeria_bounds$lon_max - nigeria_bounds$lon_min
nigeria_height <- nigeria_bounds$lat_max - nigeria_bounds$lat_min
nigeria_aspect <- nigeria_height / nigeria_width

cat("Nigeria bounds: Lon", round(nigeria_bounds$lon_min, 2), "to", round(nigeria_bounds$lon_max, 2),
    "°E, Lat", round(nigeria_bounds$lat_min, 2), "to", round(nigeria_bounds$lat_max, 2), "°N\n")
cat("Aspect ratio:", round(nigeria_aspect, 3), "(", round(nigeria_width, 2), "° × ", 
    round(nigeria_height, 2), "°)\n\n")

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

# Custom theme for maps
theme_map <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

# Color palettes for strategies
strategy_colors <- c("NO_TEST" = "#E31A1C", "HCT" = "#1F78B4", "PCR" = "#33A02C")

# ==============================================================================
# INDIVIDUAL PLOTS
# ==============================================================================

cat("Creating individual plots...\n")

# Plot 1: Optimal Strategy Map (corrected coordinates)
p1 <- ggplot(results, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = optimal_strategy), alpha = 0.8, size = 0.8) +
  scale_color_manual(values = strategy_colors, name = "Optimal\nStrategy") +
  labs(title = "Optimal Testing Strategy Across Nigeria",
       subtitle = paste0("Based on cost-effectiveness analysis (n = ", nrow(results), " locations)"),
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_map() +
  coord_fixed(ratio = 1/nigeria_aspect)

# Plot 2: Prevalence Distribution Map
p2 <- ggplot(results, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = mean_prevalence * 100), alpha = 0.8, size = 0.8) +
  scale_color_viridis_c(name = "Prevalence\n(%)", option = "plasma", trans = "log10",
                        labels = function(x) sprintf("%.1f", x)) +
  labs(title = "Trypanosomiasis Prevalence Distribution",
       subtitle = "Log-scaled prevalence estimates from INLA models",
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_map() +
  coord_fixed(ratio = 1/nigeria_aspect)

# Plot 3: Expected Value of Sample Information (EVSI)
p3 <- ggplot(results, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = max_evsi_mean), alpha = 0.8, size = 0.8) +
  scale_color_viridis_c(name = "Max EVSI\n($)", option = "viridis", trans = "sqrt",
                        labels = function(x) sprintf("%.2f", x)) +
  labs(title = "Expected Value of Sample Information (Mean)",
       subtitle = "Square-root scaled EVSI for optimal testing strategy", 
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_map() +
  coord_fixed(ratio = 1/nigeria_aspect)

# Plot 4: Cost Comparison
p4 <- ggplot(results, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = min_cost_mean), alpha = 0.8, size = 0.8) +
  scale_color_viridis_c(name = "Min Cost\n($)", option = "inferno", trans = "log10",
                        labels = function(x) sprintf("%.1f", x)) +
  labs(title = "Minimum Expected Cost (Mean)",
       subtitle = "Log-scaled cost for optimal strategy",
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_map() +
  coord_fixed(ratio = 1/nigeria_aspect)

# Plot 5: Strategy Distribution with Uncertainty
strategy_summary <- summary_stats$strategy_distribution

# Create probability-based comparison
prob_summary <- data.frame(
  strategy = c("NO_TEST", "HCT", "PCR"),
  mean_probability = c(
    mean(results$prob_no_test_optimal) * 100,
    mean(results$prob_hct_optimal) * 100, 
    mean(results$prob_pcr_optimal) * 100
  ),
  point_estimate = c(
    sum(results$optimal_strategy == "NO_TEST") / nrow(results) * 100,
    sum(results$optimal_strategy == "HCT") / nrow(results) * 100,
    sum(results$optimal_strategy == "PCR") / nrow(results) * 100
  )
) %>%
  pivot_longer(cols = c(mean_probability, point_estimate), 
               names_to = "estimate_type", values_to = "percentage") %>%
  mutate(estimate_type = case_when(
    estimate_type == "mean_probability" ~ "Probability\n(accounting for uncertainty)",
    estimate_type == "point_estimate" ~ "Point Estimate\n(mean costs only)"
  ))

p5 <- ggplot(prob_summary, aes(x = reorder(strategy, -percentage), y = percentage, 
                               fill = strategy, alpha = estimate_type)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = strategy_colors) +
  scale_alpha_manual(values = c("Probability\n(accounting for uncertainty)" = 1.0, 
                                "Point Estimate\n(mean costs only)" = 0.6)) +
  labs(title = "Distribution of Optimal Testing Strategies", 
       subtitle = "Comparison: Point estimates vs uncertainty-aware probabilities",
       x = "Testing Strategy", y = "Percentage (%)",
       alpha = "Estimate Type") +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# Plot 6: Prevalence vs EVSI Relationship
p6 <- ggplot(results, aes(x = mean_prevalence * 100, y = max_evsi_mean)) +
  geom_point(aes(color = optimal_strategy), alpha = 0.7, size = 1) +
  scale_color_manual(values = strategy_colors, name = "Optimal\nStrategy") +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  scale_x_log10(labels = function(x) sprintf("%.1f", x)) +
  scale_y_sqrt(labels = function(x) sprintf("%.2f", x)) +
  labs(title = "EVSI vs Prevalence Relationship", 
       subtitle = "Higher prevalence generally increases value of information",
       x = "Prevalence (%, log scale)", y = "Max EVSI ($, sqrt scale)") +
  theme_minimal()

# Plot 7: Cost Distribution by Strategy
cost_data <- results %>%
  select(longitude, latitude, optimal_strategy, cost_no_test_mean, cost_hct_mean, cost_pcr_mean) %>%
  pivot_longer(cols = ends_with("_mean"), names_to = "strategy_type", values_to = "cost") %>%
  mutate(strategy_type = case_when(
    strategy_type == "cost_no_test_mean" ~ "NO_TEST",
    strategy_type == "cost_hct_mean" ~ "HCT", 
    strategy_type == "cost_pcr_mean" ~ "PCR"
  ))

p7 <- ggplot(cost_data, aes(x = strategy_type, y = cost, fill = strategy_type)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = strategy_colors) +
  scale_y_log10(labels = function(x) sprintf("$%.1f", x)) +
  labs(title = "Cost Distribution by Strategy Type",
       subtitle = "Violin plots with quartiles (log scale)",
       x = "Strategy Type", y = "Expected Cost ($, log scale)") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 8: Prevalence Distribution Histogram
p8 <- ggplot(results, aes(x = mean_prevalence * 100)) +
  geom_histogram(bins = 40, fill = "#1F78B4", alpha = 0.7, color = "white") +
  scale_x_log10(labels = function(x) sprintf("%.1f%%", x)) +
  labs(title = "Distribution of Prevalence Estimates",
       subtitle = paste0("n = ", nrow(results), " locations across Nigeria"),
       x = "Prevalence (%, log scale)", y = "Number of Locations") +
  theme_minimal() +
  geom_vline(xintercept = mean(results$mean_prevalence) * 100, 
             color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = mean(results$mean_prevalence) * 100 * 1.5, 
           y = Inf, vjust = 2, label = paste0("Mean: ", 
           round(mean(results$mean_prevalence) * 100, 2), "%"), color = "red")

# Plot 9: EVSI Distribution
p9 <- ggplot(results, aes(x = max_evsi_mean)) +
  geom_histogram(bins = 40, fill = "#33A02C", alpha = 0.7, color = "white") +
  scale_x_sqrt(labels = function(x) sprintf("$%.2f", x)) +
  labs(title = "Distribution of Maximum EVSI (Mean)",
       subtitle = "Expected Value of Sample Information across Nigeria",
       x = "Max EVSI ($, sqrt scale)", y = "Number of Locations") +
  theme_minimal() +
  geom_vline(xintercept = mean(results$max_evsi_mean), color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = mean(results$max_evsi_mean) * 1.5, y = Inf, vjust = 2, 
           label = paste0("Mean: $", round(mean(results$max_evsi_mean), 2)), color = "red")

# ==============================================================================
# UNCERTAINTY-FOCUSED PLOTS
# ==============================================================================

cat("Creating uncertainty-focused plots...\n")

# Plot 10: Decision Uncertainty Map
p10 <- results %>%
  mutate(max_strategy_prob = pmax(prob_no_test_optimal, prob_hct_optimal, prob_pcr_optimal),
         decision_uncertainty = 1 - max_strategy_prob) %>%
  ggplot(aes(x = latitude, y = longitude)) +
  geom_point(aes(color = decision_uncertainty), alpha = 0.8, size = 0.8) +
  scale_color_viridis_c(name = "Decision\nUncertainty", option = "plasma", 
                        labels = function(x) sprintf("%.2f", x)) +
  labs(title = "Decision Uncertainty Across Nigeria",
       subtitle = "1 - max(probability of optimal strategy) - higher values = more uncertain",
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_map() +
  coord_fixed(ratio = 1/nigeria_aspect)

# Plot 11: Cost Uncertainty Map  
p11 <- results %>%
  mutate(cost_uncertainty = sqrt((cost_no_test_sd^2 + cost_hct_sd^2 + cost_pcr_sd^2) / 3)) %>%
  ggplot(aes(x = latitude, y = longitude)) +
  geom_point(aes(color = cost_uncertainty), alpha = 0.8, size = 0.8) +
  scale_color_viridis_c(name = "Cost\nUncertainty\n($)", option = "cividis",
                        labels = function(x) sprintf("%.2f", x)) +
  labs(title = "Cost Uncertainty Across Nigeria", 
       subtitle = "Standard deviation of costs (sqrt of mean variance across strategies)",
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_map() +
  coord_fixed(ratio = 1/nigeria_aspect)

# Plot 12: Prevalence Uncertainty Map
p12 <- ggplot(results, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = sd_prevalence * 100), alpha = 0.8, size = 0.8) +
  scale_color_viridis_c(name = "Prevalence\nSD (%)", option = "turbo",
                        labels = function(x) sprintf("%.1f", x)) +
  labs(title = "Prevalence Uncertainty Across Nigeria",
       subtitle = "Standard deviation across INLA models (percentage points)", 
       x = "Latitude (°N)", y = "Longitude (°E)") +
  theme_map() +
  coord_fixed(ratio = 1/nigeria_aspect)

# ==============================================================================
# DASHBOARD COMPOSITIONS
# ==============================================================================

cat("Creating dashboard compositions...\n")

# Main Strategy Dashboard
dashboard_main <- (p1 | p5) / (p2 | p6)
dashboard_main <- dashboard_main + 
  plot_annotation(title = "Nigeria Trypanosomiasis Testing Strategy Analysis",
                  subtitle = "Bayesian Cost-Effectiveness Framework Results",
                  theme = theme(plot.title = element_text(size = 16, face = "bold"),
                               plot.subtitle = element_text(size = 12)))

# Economic Analysis Dashboard  
dashboard_economic <- (p3 | p7) / (p4 | p9)
dashboard_economic <- dashboard_economic +
  plot_annotation(title = "Economic Analysis of Testing Strategies",
                  subtitle = "EVSI, Costs, and Value Distribution",
                  theme = theme(plot.title = element_text(size = 16, face = "bold"),
                               plot.subtitle = element_text(size = 12)))

# Prevalence Dashboard
dashboard_prevalence <- (p2 | p8) / (p6 | p9)
dashboard_prevalence <- dashboard_prevalence +
  plot_annotation(title = "Prevalence Analysis and Value Relationships", 
                  subtitle = "Geographic and Statistical Distribution",
                  theme = theme(plot.title = element_text(size = 16, face = "bold"),
                               plot.subtitle = element_text(size = 12)))

# ==============================================================================
# SAVE VISUALIZATIONS
# ==============================================================================

cat("Saving visualizations...\n")

# Create figures directory in Cost effectiveness analysis folder
dir.create("Code/Cost effectiveness analysis/figures", showWarnings = FALSE)

# Save individual plots
ggsave("Code/Cost effectiveness analysis/figures/01_optimal_strategy_map.png", p1, width = 12, height = 10, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/02_prevalence_map.png", p2, width = 12, height = 10, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/03_evsi_map.png", p3, width = 12, height = 10, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/04_cost_map.png", p4, width = 12, height = 10, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/05_strategy_distribution.png", p5, width = 10, height = 8, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/06_prevalence_evsi_relationship.png", p6, width = 10, height = 8, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/07_cost_distribution.png", p7, width = 10, height = 8, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/08_prevalence_histogram.png", p8, width = 10, height = 8, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/09_evsi_histogram.png", p9, width = 10, height = 8, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/10_decision_uncertainty_map.png", p10, width = 12, height = 10, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/11_cost_uncertainty_map.png", p11, width = 12, height = 10, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/12_prevalence_uncertainty_map.png", p12, width = 12, height = 10, dpi = 300)

# Save dashboard compositions
ggsave("Code/Cost effectiveness analysis/figures/dashboard_main_strategy.png", dashboard_main, width = 16, height = 12, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/dashboard_economic_analysis.png", dashboard_economic, width = 16, height = 12, dpi = 300)
ggsave("Code/Cost effectiveness analysis/figures/dashboard_prevalence_analysis.png", dashboard_prevalence, width = 16, height = 12, dpi = 300)

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

cat("\n=== VISUALIZATION COMPLETE ===\n")

cat("\nFiles Generated:\n")
cat("INDIVIDUAL PLOTS (9):\n")
cat("• 01_optimal_strategy_map.png - Geographic distribution of optimal strategies\n")
cat("• 02_prevalence_map.png - Prevalence estimates across Nigeria\n") 
cat("• 03_evsi_map.png - Expected Value of Sample Information\n")
cat("• 04_cost_map.png - Minimum expected costs\n")
cat("• 05_strategy_distribution.png - Bar chart of strategy percentages\n")
cat("• 06_prevalence_evsi_relationship.png - Scatter plot with trend\n")
cat("• 07_cost_distribution.png - Violin plots by strategy\n")
cat("• 08_prevalence_histogram.png - Prevalence distribution\n")
cat("• 09_evsi_histogram.png - EVSI distribution\n")

cat("\nDASHBOARD COMPOSITIONS (3):\n")
cat("• dashboard_main_strategy.png - Main strategy overview\n")
cat("• dashboard_economic_analysis.png - Economic analysis focus\n") 
cat("• dashboard_prevalence_analysis.png - Prevalence relationships\n")

# Print key insights
cat("\n=== KEY INSIGHTS ===\n")
for (i in 1:nrow(strategy_summary)) {
  cat("•", strategy_summary$optimal_strategy[i], "optimal in", 
      strategy_summary$percentage[i], "% of locations (", strategy_summary$n[i], "sites)\n")
}

cat("• Mean prevalence:", round(summary_stats$prevalence_stats$mean * 100, 2), "%\n")
cat("• Total EVSI across Nigeria: $", round(summary_stats$cost_stats$total_evsi), "\n")
cat("• Mean EVSI per location: $", round(summary_stats$cost_stats$mean_evsi, 2), "\n")
cat("• Geographic extent:", round(nigeria_width, 1), "° × ", round(nigeria_height, 1), 
    "° (aspect ratio ", round(nigeria_aspect, 2), ")\n")

cat("\nAll visualizations saved to 'Code/Cost effectiveness analysis/figures/' directory!\n")