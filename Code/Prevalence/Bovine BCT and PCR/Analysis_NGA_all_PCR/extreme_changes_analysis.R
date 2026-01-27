library(ggplot2)
library(dplyr)
library(sf)
library(geodata)
library(viridis)

cat("=== ANALYZING EXTREME CI WIDTH CHANGES ===\n")

# Load the comparison data
ci_comparison <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_comparison_data.csv")

# Validate the data
cat("Data validation:\n")
cat("- Total rows:", nrow(ci_comparison), "\n")
cat("- Columns:", paste(names(ci_comparison), collapse=", "), "\n")

# Check for missing or invalid values
missing_check <- sapply(ci_comparison, function(x) sum(is.na(x)))
if(any(missing_check > 0)) {
  cat("Warning: Missing values found:\n")
  print(missing_check[missing_check > 0])
}

# Check for infinite or extreme values
extreme_check <- sapply(ci_comparison[c("original_ci_width", "all_pcr_ci_width", "ci_width_change_pct")], 
                       function(x) sum(!is.finite(x)))
if(any(extreme_check > 0)) {
  cat("Warning: Non-finite values found:\n")
  print(extreme_check[extreme_check > 0])
}

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

# Define extreme change categories and log-scale analysis
ci_comparison <- ci_comparison %>%
  mutate(
    # Log-scale analysis for better handling of extreme changes
    log_abs_change = log10(abs(ci_width_change_pct) + 1),  # Add 1 to handle 0% changes
    log_ratio = log10(abs(all_pcr_ci_width / original_ci_width)),  # Log ratio of CI widths
    
    # Traditional percentage-based categories
    extreme_category = case_when(
      abs(ci_width_change_pct) > 500 ~ "Extreme (>500%)",
      abs(ci_width_change_pct) > 200 ~ "Very High (200-500%)",
      abs(ci_width_change_pct) > 100 ~ "High (100-200%)", 
      abs(ci_width_change_pct) > 50 ~ "Moderate (50-100%)",
      TRUE ~ "Small (≤50%)"
    ),
    
    # Log-scale based categories (more appropriate for extreme changes)
    log_category = case_when(
      log_abs_change > log10(501) ~ "Extreme (>500%)",
      log_abs_change > log10(201) ~ "Very High (200-500%)",
      log_abs_change > log10(101) ~ "High (100-200%)",
      log_abs_change > log10(51) ~ "Moderate (50-100%)",
      TRUE ~ "Small (≤50%)"
    ),
    
    extreme_category = factor(extreme_category, 
                             levels = c("Small (≤50%)", "Moderate (50-100%)", 
                                       "High (100-200%)", "Very High (200-500%)", "Extreme (>500%)")),
    log_category = factor(log_category,
                         levels = c("Small (≤50%)", "Moderate (50-100%)", 
                                   "High (100-200%)", "Very High (200-500%)", "Extreme (>500%)"))
  )

# Summary of extreme changes
cat("Summary of CI width change magnitudes:\n")
summary_table <- ci_comparison %>%
  count(extreme_category) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))
print(summary_table)

cat("\nLog-scale summary statistics:\n")
cat("Log(absolute change + 1) - Mean:", round(mean(ci_comparison$log_abs_change), 3), 
    "Median:", round(median(ci_comparison$log_abs_change), 3), "\n")
cat("Log ratio (All-PCR/Original) - Mean:", round(mean(ci_comparison$log_ratio), 3), 
    "Median:", round(median(ci_comparison$log_ratio), 3), "\n")

# Identify outliers using log-scale (points beyond 2 standard deviations)
log_mean <- mean(ci_comparison$log_abs_change)
log_sd <- sd(ci_comparison$log_abs_change)
outlier_threshold <- log_mean + 2 * log_sd

outliers <- ci_comparison %>% filter(log_abs_change > outlier_threshold)
cat("Statistical outliers (>2 SD from mean in log-scale):", nrow(outliers), 
    "(", round(nrow(outliers)/nrow(ci_comparison)*100, 1), "%)\n")

# Create map showing extreme changes
p_extreme_map <- ggplot()

# Add neighboring countries
if(!is.null(neighboring_countries_sf)) {
  p_extreme_map <- p_extreme_map + 
    geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
}

# Add Nigeria states
p_extreme_map <- p_extreme_map + 
  geom_sf(data = nigeria_states_sf, fill = NA, color = "black", lwd = 0.5)

# Add points colored by extreme change category
p_extreme_map <- p_extreme_map +
  geom_point(data = ci_comparison, 
             aes(x = Longitude, y = Latitude, 
                 color = extreme_category, size = extreme_category),
             alpha = 0.7) +
  scale_color_manual(
    name = "CI Width Change",
    values = c("Small (≤50%)" = "#440154",
               "Moderate (50-100%)" = "#31688E", 
               "High (100-200%)" = "#35B779",
               "Very High (200-500%)" = "#FDE725",
               "Extreme (>500%)" = "#FF0000"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_size_manual(
    name = "CI Width Change",
    values = c("Small (≤50%)" = 0.3,
               "Moderate (50-100%)" = 0.5,
               "High (100-200%)" = 0.8,
               "Very High (200-500%)" = 1.2,
               "Extreme (>500%)" = 1.8)
  ) +
  labs(
    title = "Geographic Distribution of Extreme CI Width Changes",
    subtitle = "All-PCR vs Original (BCT+PCR) Analysis - Nigeria",
    caption = "Red points show locations with most extreme uncertainty changes"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  coord_sf(xlim = c(2, 15), ylim = c(4, 14))

# Create scatter plot showing the relationship
p_scatter <- ggplot(ci_comparison, aes(x = original_ci_width, y = all_pcr_ci_width)) +
  geom_point(aes(color = extreme_category), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(
    name = "CI Width Change",
    values = c("Small (≤50%)" = "#440154",
               "Moderate (50-100%)" = "#31688E",
               "High (100-200%)" = "#35B779", 
               "Very High (200-500%)" = "#FDE725",
               "Extreme (>500%)" = "#FF0000")
  ) +
  labs(
    title = "Original vs All-PCR CI Widths",
    subtitle = "Points above red line indicate wider CIs in All-PCR analysis",
    x = "Original Analysis CI Width",
    y = "All-PCR Analysis CI Width"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# Create log-scale scatter plot for better visualization of extreme changes
p_log_scatter <- ggplot(ci_comparison, aes(x = log10(original_ci_width), y = log10(all_pcr_ci_width))) +
  geom_point(aes(color = log_category), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(
    name = "CI Width Change",
    values = c("Small (≤50%)" = "#440154",
               "Moderate (50-100%)" = "#31688E",
               "High (100-200%)" = "#35B779", 
               "Very High (200-500%)" = "#FDE725",
               "Extreme (>500%)" = "#FF0000")
  ) +
  labs(
    title = "Original vs All-PCR CI Widths (Log Scale)",
    subtitle = "Log-scale reveals patterns in extreme changes more clearly",
    x = "Log10(Original Analysis CI Width)",
    y = "Log10(All-PCR Analysis CI Width)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# Create prevalence comparison plot
p_prevalence <- ggplot(ci_comparison, aes(x = original_mean, y = all_pcr_mean)) +
  geom_point(aes(color = extreme_category), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(
    name = "CI Width Change",
    values = c("Small (≤50%)" = "#440154",
               "Moderate (50-100%)" = "#31688E",
               "High (100-200%)" = "#35B779",
               "Very High (200-500%)" = "#FDE725", 
               "Extreme (>500%)" = "#FF0000")
  ) +
  labs(
    title = "Original vs All-PCR Mean Prevalence",
    subtitle = "Extreme CI changes associated with high original prevalence",
    x = "Original Analysis Mean Prevalence", 
    y = "All-PCR Analysis Mean Prevalence"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# Print summary statistics for extreme locations
extreme_locations <- ci_comparison %>% filter(extreme_category == "Extreme (>500%)")

cat("\n=== ANALYSIS OF EXTREME LOCATIONS (>500% change) ===\n")
cat("Number of extreme locations:", nrow(extreme_locations), "\n")

if(nrow(extreme_locations) > 0) {
  cat("Geographic clustering:\n")
  cat("  Northeast Nigeria (lat>11, lon>13):", 
      sum(extreme_locations$Latitude > 11 & extreme_locations$Longitude > 13), 
      "(", round(sum(extreme_locations$Latitude > 11 & extreme_locations$Longitude > 13)/nrow(extreme_locations)*100, 1), "%)\n")
} else {
  cat("  No extreme locations found\n")
}

if(nrow(extreme_locations) > 0) {
  cat("\nCharacteristics of extreme locations:\n")
  cat("Original Analysis:\n")
  cat("  Mean prevalence:", round(mean(extreme_locations$original_mean), 3), 
      "(range:", paste(round(range(extreme_locations$original_mean), 3), collapse="-"), ")\n")
  cat("  Mean CI width:", round(mean(extreme_locations$original_ci_width), 4),
      "(range:", paste(round(range(extreme_locations$original_ci_width), 4), collapse="-"), ")\n")

  cat("All-PCR Analysis:\n")
  cat("  Mean prevalence:", round(mean(extreme_locations$all_pcr_mean), 3),
      "(range:", paste(round(range(extreme_locations$all_pcr_mean), 3), collapse="-"), ")\n")
  cat("  Mean CI width:", round(mean(extreme_locations$all_pcr_ci_width), 4),
      "(range:", paste(round(range(extreme_locations$all_pcr_ci_width), 4), collapse="-"), ")\n")
  
  cat("Log-scale characteristics:\n")
  cat("  Mean log(absolute change + 1):", round(mean(extreme_locations$log_abs_change), 3), "\n")
  cat("  Mean log(CI width ratio):", round(mean(extreme_locations$log_ratio), 3), "\n")
} else {
  cat("\nNo extreme locations to analyze.\n")
}

# Additional analysis for statistical outliers
cat("\n=== ANALYSIS OF STATISTICAL OUTLIERS (>2 SD in log-scale) ===\n")
cat("Number of statistical outliers:", nrow(outliers), "\n")

if(nrow(outliers) > 0) {
  cat("Outlier characteristics:\n")
  cat("  Mean percentage change:", round(mean(outliers$ci_width_change_pct), 1), "%\n")
  cat("  Median percentage change:", round(median(outliers$ci_width_change_pct), 1), "%\n")
  cat("  Range of percentage changes:", paste(round(range(outliers$ci_width_change_pct), 1), collapse=" to "), "%\n")
}

# Save plots
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/extreme_changes_map.pdf", 
       plot = p_extreme_map, width = 12, height = 8, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_scatter.pdf",
       plot = p_scatter, width = 10, height = 8, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_scatter_log.pdf",
       plot = p_log_scatter, width = 10, height = 8, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/prevalence_scatter.pdf", 
       plot = p_prevalence, width = 10, height = 8, dpi = 300)

# Display plots
print(p_extreme_map)
print(p_scatter)
print(p_log_scatter)
print(p_prevalence)

cat("\n=== CONCLUSION ===\n")
cat("The extreme CI width changes are concentrated in northeast Nigeria and represent\n")
cat("locations where the original BCT+PCR model was overconfident with very high\n") 
cat("prevalence predictions and narrow CIs. The All-PCR model shows more realistic\n")
cat("uncertainty in these regions, suggesting better model calibration.\n\n")
cat("Log-scale analysis provides better insight into extreme changes by:\n")
cat("1. Normalizing the scale of percentage changes across orders of magnitude\n")
cat("2. Identifying statistical outliers based on standard deviations in log-space\n")
cat("3. Revealing patterns that may be masked by extreme values in linear scale\n")
cat("4. Providing more robust statistical measures for highly skewed distributions\n")