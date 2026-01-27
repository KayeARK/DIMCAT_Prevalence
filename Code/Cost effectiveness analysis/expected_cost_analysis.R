#!/usr/bin/env R
# ==============================================================================
# EXPECTED COST ANALYSIS - NIGERIA TRYPANOSOMIASIS DIAGNOSTICS
# ==============================================================================
# Analysis of expected costs under current information (first part of EVSI)
# Shows the actual economic burden at each location given current uncertainty
# ==============================================================================

library(tidyverse)
library(sf)
library(geodata)
library(viridis)
library(gridExtra)
library(ggspatial)
library(cowplot)

setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")

cat("=== EXPECTED COST ANALYSIS ===\n")
cat("Analyzing expected costs under current information (no perfect info correction)\n\n")

# ==============================================================================
# LOAD RESULTS
# ==============================================================================

# Load cost-effectiveness results
results <- read_csv("Code/Cost effectiveness analysis/results/nigeria_cost_effectiveness_results.csv", 
                   show_col_types = FALSE)

cat("Loaded results for", nrow(results), "locations\n")

# ==============================================================================
# EXPECTED COST CALCULATIONS
# ==============================================================================

# Calculate expected cost under optimal strategy (current information)
results <- results %>%
  mutate(
    # Expected cost under current information = cost of optimal strategy
    expected_cost_current = case_when(
      optimal_strategy == "NO_TEST" ~ cost_no_test_mean,
      optimal_strategy == "HCT" ~ cost_hct_mean,
      optimal_strategy == "PCR" ~ cost_pcr_mean,
      TRUE ~ min_cost_mean  # Should be the same as min_cost_mean
    ),
    
    # Expected cost uncertainty (SD of optimal strategy)
    expected_cost_sd = case_when(
      optimal_strategy == "NO_TEST" ~ cost_no_test_sd,
      optimal_strategy == "HCT" ~ cost_hct_sd,
      optimal_strategy == "PCR" ~ cost_pcr_sd,
      TRUE ~ NA_real_
    ),
    
    # Alternative: risk-weighted expected cost (probability-weighted across strategies)
    expected_cost_weighted = prob_no_test_optimal * cost_no_test_mean + 
                             prob_hct_optimal * cost_hct_mean + 
                             prob_pcr_optimal * cost_pcr_mean,
    
    # Cost per prevalence unit (expected cost divided by prevalence)
    cost_per_prevalence = expected_cost_current / pmax(mean_prevalence, 0.001),  # Avoid division by 0
    
    # Cost difference from cheapest strategy (NO_TEST)
    cost_above_no_test = expected_cost_current - cost_no_test_mean,
    
    # Categorize cost levels
    cost_category = case_when(
      expected_cost_current < 5 ~ "Very Low (<$5)",
      expected_cost_current < 15 ~ "Low ($5-15)",
      expected_cost_current < 25 ~ "Moderate ($15-25)", 
      expected_cost_current < 35 ~ "High ($25-35)",
      TRUE ~ "Very High (>$35)"
    ),
    cost_category = factor(cost_category, 
                          levels = c("Very Low (<$5)", "Low ($5-15)", "Moderate ($15-25)", 
                                   "High ($25-35)", "Very High (>$35)"))
  )

# Verify expected_cost_current matches min_cost_mean
cost_match_check <- all(abs(results$expected_cost_current - results$min_cost_mean) < 1e-10, na.rm = TRUE)
cat("Expected cost calculation verified:", cost_match_check, "\n")

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

cat("\n=== EXPECTED COST SUMMARY ===\n")
cat("Expected cost under current information:\n")
cat("  Mean: $", round(mean(results$expected_cost_current), 2), "\n")
cat("  Median: $", round(median(results$expected_cost_current), 2), "\n")
cat("  Range: $", round(min(results$expected_cost_current), 2), " to $", 
    round(max(results$expected_cost_current), 2), "\n")
cat("  Standard deviation: $", round(sd(results$expected_cost_current), 2), "\n")

cat("\n=== OPTIMAL TESTING STRATEGIES (COST-MINIMIZING) ===\n")
optimal_strategy_summary <- results %>%
  group_by(optimal_strategy) %>%
  summarise(
    n_locations = n(),
    percentage = round(n() / nrow(results) * 100, 1),
    mean_cost = round(mean(expected_cost_current), 2),
    median_cost = round(median(expected_cost_current), 2),
    min_cost = round(min(expected_cost_current), 2),
    max_cost = round(max(expected_cost_current), 2),
    mean_prevalence = round(mean(mean_prevalence) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_locations))

print(optimal_strategy_summary)

cat("\nStrategy recommendations:\n")
for(i in seq_len(nrow(optimal_strategy_summary))) {
  strategy <- optimal_strategy_summary$optimal_strategy[i]
  n_loc <- optimal_strategy_summary$n_locations[i]
  pct <- optimal_strategy_summary$percentage[i]
  mean_cost <- optimal_strategy_summary$mean_cost[i]
  mean_prev <- optimal_strategy_summary$mean_prevalence[i]
  
  cat("• ", strategy, ": Optimal for ", n_loc, " locations (", pct, "%) - ", 
      "Mean cost: $", mean_cost, " at ", mean_prev, "% prevalence\n", sep = "")
}

cat("\nCost category distribution:\n")
cost_dist <- table(results$cost_category)
for(i in seq_along(cost_dist)) {
  cat("  ", names(cost_dist)[i], ": ", cost_dist[i], " locations (", 
      round(cost_dist[i]/sum(cost_dist)*100, 1), "%)\n", sep="")
}

cat("\nRelationship with prevalence:\n")
prev_cost_cor <- cor(results$mean_prevalence, results$expected_cost_current)
cat("  Correlation with prevalence: ", round(prev_cost_cor, 3), "\n")

# Cost by prevalence quintiles
results$prev_quintile <- cut(results$mean_prevalence, 
                            breaks = quantile(results$mean_prevalence, probs = 0:5/5),
                            labels = paste0("Q", 1:5), include.lowest = TRUE)

quintile_costs <- results %>%
  group_by(prev_quintile) %>%
  summarise(
    mean_prevalence = mean(mean_prevalence),
    mean_cost = mean(expected_cost_current),
    median_cost = median(expected_cost_current),
    .groups = "drop"
  )
cat("  Cost by prevalence quintile:\n")
print(quintile_costs)

# ==============================================================================
# LOAD GEOGRAPHIC DATA
# ==============================================================================

cat("\nLoading geographic boundaries...\n")

# Get Nigerian boundaries
nigeria_states <- gadm(country = "NGA", level = 1, path = tempdir()) 
nigeria_states_sf <- st_as_sf(nigeria_states)

# Get neighboring countries for context
neighboring_countries <- c("BEN", "NER", "TCD", "CMR")
neighbor_sf_list <- list()

for (country in neighboring_countries) {
  tryCatch({
    country_data <- gadm(country = country, level = 0, path = tempdir())
    neighbor_sf_list[[country]] <- st_as_sf(country_data)
  }, error = function(e) {
    # Silently skip if can't load
  })
}

if(length(neighbor_sf_list) > 0) {
  neighboring_countries_sf <- do.call(rbind, neighbor_sf_list)
} else {
  neighboring_countries_sf <- NULL
}

# ==============================================================================
# VISUALIZATION - HEATMAP STYLE
# ==============================================================================

# Define strategy colors
strategy_colors <- c(
  "NO_TEST" = "#e31a1c",    # Red for no testing
  "HCT" = "#1f78b4",        # Blue for HCT  
  "PCR" = "#33a02c"         # Green for PCR
)

# Expected cost heatmap using stat_summary_2d for smooth surfaces
p_cost_map <- ggplot(data = results, aes(x = latitude, y = longitude)) +
  stat_summary_2d(aes(z = log10(min_cost_mean)), 
                  bins = 30, alpha = 0.9, fun = mean) +
  scale_fill_viridis_c(
    name = "Expected\nCost ($)",
    option = "plasma",
    labels = function(x) scales::dollar_format()(10^x),
    na.value = "transparent"
  ) +
  labs(
    title = "Expected Cost Under Current Information",
    subtitle = "Cost per animal given current test parameter and prevalence uncertainty",
    caption = "Log scale heatmap showing expected cost of optimal diagnostic strategy"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10, color = "grey60"),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    aspect.ratio = 0.7,
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# Add country borders if available
if(!is.null(neighboring_countries_sf)) {
  p_cost_map <- p_cost_map + 
    geom_sf(data = neighboring_countries_sf, fill = NA, color = "grey60", lwd = 0.3, inherit.aes = FALSE)
}

# Add Nigeria state borders
p_cost_map <- p_cost_map + 
  geom_sf(data = nigeria_states_sf, fill = NA, color = "white", lwd = 1.0, inherit.aes = FALSE) +
  coord_sf(xlim = c(3, 15), ylim = c(2, 15.5), expand = TRUE, datum = sf::st_crs(4326))

# Strategy heatmap - use dense points with better styling
p_strategy_map <- ggplot(data = results, aes(x = latitude, y = longitude)) +
  geom_point(aes(color = optimal_strategy), alpha = 0.8, size = 1.0, stroke = 0) +
  scale_color_manual(
    name = "Optimal\nStrategy",
    values = strategy_colors,
    labels = c(
      "HCT" = "HCT (Rapid test)",
      "NO_TEST" = "No testing", 
      "PCR" = "PCR (Molecular)"
    )
  ) +
  labs(
    title = "Cost-Minimizing Diagnostic Strategy by Location",
    subtitle = "Nigeria - Optimal testing approach based on prevalence and costs",
    caption = "Each point represents a location with optimal strategy recommendation"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10, color = "grey60"),
    legend.position = "right",
    legend.key.size = unit(1.2, "cm"),
    aspect.ratio = 0.7,
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# Add country borders if available  
if(!is.null(neighboring_countries_sf)) {
  p_strategy_map <- p_strategy_map + 
    geom_sf(data = neighboring_countries_sf, fill = "grey98", color = "grey70", lwd = 0.3, inherit.aes = FALSE)
}

# Add Nigeria state borders on top
p_strategy_map <- p_strategy_map + 
  geom_sf(data = nigeria_states_sf, fill = NA, color = "white", lwd = 1.0, inherit.aes = FALSE) +
  coord_sf(xlim = c(3, 15), ylim = c(2, 15.5), expand = TRUE, datum = sf::st_crs(4326))

# Expected cost vs prevalence scatter plot
p_cost_vs_prev <- ggplot(results, aes(x = mean_prevalence * 100, y = min_cost_mean)) +
  geom_point(aes(color = optimal_strategy), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    name = "Optimal Strategy",
    values = strategy_colors,
    labels = c("HCT" = "HCT (Rapid test)", "NO_TEST" = "No testing", "PCR" = "PCR (Molecular)")
  ) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = scales::dollar_format()) +
  labs(
    title = "Expected Cost vs Prevalence by Optimal Strategy", 
    subtitle = "PCR becomes cost-effective at high prevalence due to reduced treatment costs",
    x = "Mean Prevalence (%)",
    y = "Expected Cost ($)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# Cost uncertainty explanation
cat("\n=== COST UNCERTAINTY EXPLANATION ===\n")
cat("Cost uncertainty represents the standard deviation of costs across:\n")
cat("• Different prevalence estimates (Monte Carlo samples)\n") 
cat("• Different test performance parameters (sensitivity/specificity)\n")
cat("• Spatial variation in disease prevalence\n")
cat("• Treatment cost variations\n")
cat("Higher uncertainty indicates locations where additional sampling\n")
cat("could provide greater value for reducing decision uncertainty.\n")

# Expected cost distribution by strategy
p_cost_dist <- ggplot(results, aes(x = min_cost_mean, fill = optimal_strategy)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "white") +
  scale_fill_manual(
    name = "Optimal Strategy",
    values = strategy_colors,
    labels = c("HCT" = "HCT (Rapid test)", "NO_TEST" = "No testing", "PCR" = "PCR (Molecular)")
  ) +
  scale_x_continuous(labels = scales::dollar_format()) +
  labs(
    title = "Distribution of Expected Costs by Strategy",
    subtitle = "PCR optimal at moderate-high prevalence, HCT dominates most areas",
    x = "Expected Cost ($)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# Cost uncertainty vs expected cost
p_cost_uncertainty <- ggplot(results, aes(x = min_cost_mean, y = cost_hct_sd)) +
  geom_point(aes(color = optimal_strategy), alpha = 0.6, size = 0.8) +
  scale_color_manual(
    name = "Optimal Strategy", 
    values = strategy_colors,
    labels = c("HCT" = "HCT (Rapid test)", "NO_TEST" = "No testing", "PCR" = "PCR (Molecular)")
  ) +
  scale_x_continuous(labels = scales::dollar_format()) +
  scale_y_continuous(labels = scales::dollar_format()) +
  labs(
    title = "Cost Uncertainty vs Expected Cost",
    subtitle = "Uncertainty represents variability across prevalence and test parameters",
    x = "Expected Cost ($)",
    y = "Cost Uncertainty (SD, $)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

# Create output directory
dir.create("Code/Cost effectiveness analysis/expected_cost_analysis", showWarnings = FALSE)

# Save plots with better dimensions
ggsave("Code/Cost effectiveness analysis/expected_cost_analysis/expected_cost_map.pdf", 
       plot = p_cost_map, width = 16, height = 10, dpi = 300)

ggsave("Code/Cost effectiveness analysis/expected_cost_analysis/strategy_map.pdf", 
       plot = p_strategy_map, width = 16, height = 10, dpi = 300)

ggsave("Code/Cost effectiveness analysis/expected_cost_analysis/cost_vs_prevalence.pdf",
       plot = p_cost_vs_prev, width = 10, height = 8, dpi = 300)

ggsave("Code/Cost effectiveness analysis/expected_cost_analysis/cost_distribution.pdf", 
       plot = p_cost_dist, width = 10, height = 6, dpi = 300)

ggsave("Code/Cost effectiveness analysis/expected_cost_analysis/cost_uncertainty.pdf",
       plot = p_cost_uncertainty, width = 10, height = 8, dpi = 300)

# Save summary data
expected_cost_summary <- results %>%
  select(longitude, latitude, mean_prevalence, expected_cost_current, expected_cost_sd,
         expected_cost_weighted, cost_per_prevalence, cost_above_no_test, 
         optimal_strategy, cost_category, prev_quintile) %>%
  arrange(desc(expected_cost_current))

write_csv(expected_cost_summary, 
          "Code/Cost effectiveness analysis/expected_cost_analysis/expected_cost_summary.csv")

# Save detailed statistics
sink("Code/Cost effectiveness analysis/expected_cost_analysis/expected_cost_statistics.txt")
cat("=== EXPECTED COST ANALYSIS STATISTICS ===\n")
cat("Generated on:", as.character(Sys.time()), "\n\n")

cat("OVERALL STATISTICS:\n")
cat("Total locations:", nrow(results), "\n")
cat("Mean expected cost: $", round(mean(results$expected_cost_current), 2), "\n")
cat("Median expected cost: $", round(median(results$expected_cost_current), 2), "\n")
cat("Cost range: $", round(min(results$expected_cost_current), 2), " to $", 
    round(max(results$expected_cost_current), 2), "\n")
cat("Total economic burden (all locations): $", 
    round(sum(results$expected_cost_current)), "\n\n")

cat("COST BY STRATEGY:\n")
print(optimal_strategy_summary)
cat("\n")

cat("COST CATEGORY DISTRIBUTION:\n")
print(as.data.frame(table(results$cost_category)))
cat("\n")

cat("PREVALENCE-COST RELATIONSHIP:\n")
cat("Correlation: ", round(prev_cost_cor, 3), "\n")
print(quintile_costs)

sink()

# Display plots
print(p_cost_map)
print(p_strategy_map)
print(p_cost_vs_prev)
print(p_cost_dist)
print(p_cost_uncertainty)

# ==============================================================================
# STRATEGY-FOCUSED VISUALIZATIONS  
# ==============================================================================

# Strategy distribution bar plot
p_strategy_count <- results %>%
  count(optimal_strategy) %>%
  mutate(
    percentage = n / sum(n) * 100,
    strategy_label = case_when(
      optimal_strategy == "NO_TEST" ~ "No testing",
      optimal_strategy == "HCT" ~ "HCT (Rapid test)", 
      optimal_strategy == "PCR" ~ "PCR (Molecular)",
      TRUE ~ optimal_strategy
    )
  ) %>%
  ggplot(aes(x = reorder(strategy_label, -n), y = n, fill = optimal_strategy)) +
  geom_col(alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(n, " locations\n(", round(percentage, 1), "%)")), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = strategy_colors,
    guide = "none"
  ) +
  labs(
    title = "Distribution of Cost-Minimizing Testing Strategies",
    subtitle = "Number of locations where each strategy minimizes expected cost",
    x = "Testing Strategy",
    y = "Number of Locations"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

# Strategy vs prevalence ranges
p_strategy_prev_ranges <- results %>%
  mutate(
    strategy_label = case_when(
      optimal_strategy == "NO_TEST" ~ "No testing",
      optimal_strategy == "HCT" ~ "HCT (Rapid test)", 
      optimal_strategy == "PCR" ~ "PCR (Molecular)",
      TRUE ~ optimal_strategy
    )
  ) %>%
  ggplot(aes(x = strategy_label, y = mean_prevalence * 100, fill = optimal_strategy)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_manual(
    values = strategy_colors,
    guide = "none"
  ) +
  labs(
    title = "Prevalence Ranges by Optimal Testing Strategy",
    subtitle = "Distribution of prevalence levels where each strategy minimizes cost",
    x = "Optimal Testing Strategy",
    y = "Mean Prevalence (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_strategy_count)
print(p_strategy_prev_ranges)

# Save the new strategy plots
ggsave("Code/Cost effectiveness analysis/expected_cost_analysis/strategy_distribution.pdf", 
       plot = p_strategy_count, width = 10, height = 8, dpi = 300)

ggsave("Code/Cost effectiveness analysis/expected_cost_analysis/strategy_prevalence_ranges.pdf",
       plot = p_strategy_prev_ranges, width = 10, height = 8, dpi = 300)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Expected cost analysis completed!\n")
cat("Files saved to: Code/Cost effectiveness analysis/expected_cost_analysis/\n")
cat("- expected_cost_map.pdf - Geographic distribution of expected costs\n")
cat("- strategy_map.pdf - Geographic map of optimal testing strategies\n") 
cat("- cost_vs_prevalence.pdf - Cost-prevalence relationship\n") 
cat("- cost_distribution.pdf - Distribution of expected costs\n")
cat("- cost_uncertainty.pdf - Cost uncertainty analysis\n")
cat("- strategy_distribution.pdf - Bar chart of strategy distribution\n")
cat("- strategy_prevalence_ranges.pdf - Strategy vs prevalence ranges\n")
cat("- expected_cost_summary.csv - Summary data by location\n")
cat("- expected_cost_statistics.txt - Detailed statistics\n\n")

cat("KEY INSIGHTS:\n")
cat("• Mean expected cost per animal: $", round(mean(results$min_cost_mean), 2), "\n")
cat("• Cost range: $", round(min(results$min_cost_mean), 2), 
    " to $", round(max(results$min_cost_mean), 2), "\n")
cat("• Total economic burden: $", round(sum(results$min_cost_mean)), 
    " across all analyzed locations\n")

# Strategy insights
strategy_summary <- results %>%
  count(optimal_strategy, sort = TRUE) %>%
  mutate(percentage = n / sum(n) * 100)

most_common_strategy <- strategy_summary$optimal_strategy[1]
most_common_pct <- round(strategy_summary$percentage[1], 1)

cat("\nOPTIMAL TESTING STRATEGIES:\n")
cat("• Most common optimal strategy: ", most_common_strategy, " (", most_common_pct, "% of locations)\n", sep = "")

for(i in seq_len(nrow(strategy_summary))) {
  strategy <- strategy_summary$optimal_strategy[i]
  n_loc <- strategy_summary$n[i] 
  pct <- round(strategy_summary$percentage[i], 1)
  cat("• ", strategy, ": ", n_loc, " locations (", pct, "%)\n", sep = "")
}

# Prevalence insights by strategy
strategy_prev_insights <- results %>%
  group_by(optimal_strategy) %>%
  summarise(
    mean_prev = mean(mean_prevalence) * 100,
    min_prev = min(mean_prevalence) * 100,
    max_prev = max(mean_prevalence) * 100,
    .groups = "drop"
  )

cat("\nPREVALENCE PATTERNS:\n")
for(i in seq_len(nrow(strategy_prev_insights))) {
  strategy <- strategy_prev_insights$optimal_strategy[i]
  mean_prev <- round(strategy_prev_insights$mean_prev[i], 1)
  min_prev <- round(strategy_prev_insights$min_prev[i], 1)
  max_prev <- round(strategy_prev_insights$max_prev[i], 1)
  cat("• ", strategy, ": ", mean_prev, "% avg prevalence (range: ", min_prev, "%-", max_prev, "%)\n", sep = "")
}

cat("\nSTRATEGY RATIONALE:\n")
cat("• HCT dominates at moderate prevalence (3-94%) - good balance of accuracy vs cost\n")
cat("• PCR becomes optimal at high prevalence (27-99%) - accurate testing reduces treatment costs\n")
cat("• No testing optimal only at very low prevalence (3-14%) - treatment cheaper than testing\n")
cat("• PCR justified when avoiding false negatives saves more than extra test costs\n")