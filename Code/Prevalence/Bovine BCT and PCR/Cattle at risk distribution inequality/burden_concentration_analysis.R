library(ggplot2)
library(dplyr)
library(readr)
library(sf)
library(geodata)
library(viridis)
library(gridExtra)
library(grid)
library(ggridges)
library(scales)

cat("=== BURDEN CONCENTRATION ANALYSIS ===\n")
cat("Analyzing how concentrated the cattle-at-risk burden is across administrative units\n\n")

# ==============================================================================
# NIGERIA BURDEN CONCENTRATION ANALYSIS
# ==============================================================================

cat("1. NIGERIA - LGA BURDEN CONCENTRATION:\n")
cat("=====================================\n")

# Load Nigeria cattle data
nigeria_cattle_data <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Nigeria_cattle_uncertainty_fine.csv")

# Load tsetse data for proper classification
library(terra)
library(raster)
tsetse_raster <- raster("Data/Covariates/tsenumbspec")
tsetse_raster[tsetse_raster > 1] <- 1

# Get Nigerian LGAs
nigeria_lgas <- gadm(country = "NGA", level = 2, path = tempdir())
nigeria_lgas_sf <- st_as_sf(nigeria_lgas)

# Extract tsetse for Nigeria
nga_tsetse_stats <- raster::extract(tsetse_raster, nigeria_lgas_sf, fun = mean, na.rm = TRUE, df = TRUE)
nga_tsetse_data <- data.frame(
  NAME_2 = nigeria_lgas_sf$NAME_2,
  NAME_1 = nigeria_lgas_sf$NAME_1,
  tsetse_mean = nga_tsetse_stats[,2]
)
nga_tsetse_data$tsetse_mean[is.na(nga_tsetse_data$tsetse_mean)] <- 0

# Apply tsetse zone classification
nga_tsetse_states <- c("Abia", "Akwa Ibom", "Anambra", "Bayelsa", "Benue", "Cross River", 
                       "Delta", "Ebonyi", "Edo", "Ekiti", "Enugu", "Imo", "Lagos", "Ogun", 
                       "Ondo", "Osun", "Oyo", "Rivers", "Taraba", "Adamawa", "Plateau", 
                       "Nassarawa", "Kwara", "Kogi", "Federal Capital Territory")

nga_tsetse_data$tsetse_zone <- (nga_tsetse_data$tsetse_mean > 0) | 
                               (nga_tsetse_data$NAME_1 %in% nga_tsetse_states)

# Merge cattle data with tsetse data
nga_burden_analysis <- nigeria_cattle_data %>%
  dplyr::filter(state != "Total") %>%  # Exclude the total row
  dplyr::select(NAME_2 = state, 
                mean_cattle_at_risk = mean_cattle_mean,
                lower_cattle_at_risk = data_uncertainty_lower,
                upper_cattle_at_risk = data_uncertainty_upper) %>%
  right_join(nga_tsetse_data, by = "NAME_2") %>%
  dplyr::mutate(
    # Handle missing data
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    lower_cattle_at_risk = ifelse(is.na(lower_cattle_at_risk), 0, lower_cattle_at_risk),
    upper_cattle_at_risk = ifelse(is.na(upper_cattle_at_risk), 0, upper_cattle_at_risk)
  ) %>%
  # Only include LGAs with cattle at risk > 0 for concentration analysis
  dplyr::filter(mean_cattle_at_risk > 0) %>%
  dplyr::arrange(desc(mean_cattle_at_risk))

cat("Nigeria analysis includes", nrow(nga_burden_analysis), "LGAs with cattle at risk > 0\n")

# Calculate concentration metrics for Nigeria
nga_total_cattle <- sum(nga_burden_analysis$mean_cattle_at_risk, na.rm = TRUE)
nga_burden_analysis$cumulative_cattle <- cumsum(nga_burden_analysis$mean_cattle_at_risk)
nga_burden_analysis$cumulative_percent <- nga_burden_analysis$cumulative_cattle / nga_total_cattle * 100
nga_burden_analysis$lga_rank <- seq_len(nrow(nga_burden_analysis))
nga_burden_analysis$lga_percentile <- nga_burden_analysis$lga_rank / nrow(nga_burden_analysis) * 100

# Key concentration statistics for Nigeria
nga_top_10_pct_count <- ceiling(nrow(nga_burden_analysis) * 0.10)
nga_top_20_pct_count <- ceiling(nrow(nga_burden_analysis) * 0.20)
nga_top_50_pct_count <- ceiling(nrow(nga_burden_analysis) * 0.50)

nga_top_10_pct_burden <- sum(nga_burden_analysis$mean_cattle_at_risk[1:nga_top_10_pct_count], na.rm = TRUE)
nga_top_20_pct_burden <- sum(nga_burden_analysis$mean_cattle_at_risk[1:nga_top_20_pct_count], na.rm = TRUE)
nga_top_50_pct_burden <- sum(nga_burden_analysis$mean_cattle_at_risk[1:nga_top_50_pct_count], na.rm = TRUE)

nga_top_10_pct_share <- nga_top_10_pct_burden / nga_total_cattle * 100
nga_top_20_pct_share <- nga_top_20_pct_burden / nga_total_cattle * 100
nga_top_50_pct_share <- nga_top_50_pct_burden / nga_total_cattle * 100

cat("\nNIGERIA CONCENTRATION RESULTS:\n")
cat("Total cattle at risk:", format(nga_total_cattle, big.mark = ","), "\n")
cat("Top 10% of LGAs (", nga_top_10_pct_count, " LGAs) account for ", 
    round(nga_top_10_pct_share, 1), "% of burden\n", sep = "")
cat("Top 20% of LGAs (", nga_top_20_pct_count, " LGAs) account for ", 
    round(nga_top_20_pct_share, 1), "% of burden\n", sep = "")
cat("Top 50% of LGAs (", nga_top_50_pct_count, " LGAs) account for ", 
    round(nga_top_50_pct_share, 1), "% of burden\n", sep = "")

# Show top 10 LGAs
cat("\nTop 10 highest burden LGAs in Nigeria:\n")
nga_top_10 <- nga_burden_analysis %>% 
  slice_head(n = 10) %>%
  dplyr::select(NAME_2, NAME_1, mean_cattle_at_risk, cumulative_percent)
print(as.data.frame(nga_top_10))

# ==============================================================================
# ETHIOPIA BURDEN CONCENTRATION ANALYSIS
# ==============================================================================

cat("\n2. ETHIOPIA - ZONE BURDEN CONCENTRATION:\n")
cat("=======================================\n")

# Load Ethiopia cattle data
ethiopia_cattle_data <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Ethiopia_cattle_uncertainty_fine.csv")

# Get Ethiopian zones
ethiopia_zones <- gadm(country = "ETH", level = 3, path = tempdir())
ethiopia_zones_sf <- st_as_sf(ethiopia_zones)

# Extract tsetse for Ethiopia
eth_tsetse_stats <- raster::extract(tsetse_raster, ethiopia_zones_sf, fun = mean, na.rm = TRUE, df = TRUE)
eth_tsetse_data <- data.frame(
  NAME_3 = ethiopia_zones_sf$NAME_3,
  NAME_2 = ethiopia_zones_sf$NAME_2,
  NAME_1 = ethiopia_zones_sf$NAME_1,
  tsetse_mean = eth_tsetse_stats[,2]
)
eth_tsetse_data$tsetse_mean[is.na(eth_tsetse_data$tsetse_mean)] <- 0

# Apply REFINED Ethiopian tsetse zone classification (no Oromia over-classification)
eth_tsetse_regions <- c("Gambela Peoples", "Benshangul-Gumaz")
eth_tsetse_data$tsetse_zone <- (eth_tsetse_data$tsetse_mean > 0) | 
                               (eth_tsetse_data$NAME_1 %in% eth_tsetse_regions)

# Merge cattle data with tsetse data
eth_burden_analysis <- ethiopia_cattle_data %>%
  dplyr::filter(region != "Total") %>%  # Exclude the total row
  dplyr::select(NAME_3 = region, 
                mean_cattle_at_risk = mean_cattle_mean,
                lower_cattle_at_risk = data_uncertainty_lower,
                upper_cattle_at_risk = data_uncertainty_upper) %>%
  right_join(eth_tsetse_data, by = "NAME_3") %>%
  dplyr::mutate(
    # Handle missing data
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    lower_cattle_at_risk = ifelse(is.na(lower_cattle_at_risk), 0, lower_cattle_at_risk),
    upper_cattle_at_risk = ifelse(is.na(upper_cattle_at_risk), 0, upper_cattle_at_risk)
  ) %>%
  # Only include zones with cattle at risk > 0 for concentration analysis
  dplyr::filter(mean_cattle_at_risk > 0) %>%
  dplyr::arrange(desc(mean_cattle_at_risk))

cat("Ethiopia analysis includes", nrow(eth_burden_analysis), "zones with cattle at risk > 0\n")

# Calculate concentration metrics for Ethiopia
eth_total_cattle <- sum(eth_burden_analysis$mean_cattle_at_risk, na.rm = TRUE)
eth_burden_analysis$cumulative_cattle <- cumsum(eth_burden_analysis$mean_cattle_at_risk)
eth_burden_analysis$cumulative_percent <- eth_burden_analysis$cumulative_cattle / eth_total_cattle * 100
eth_burden_analysis$zone_rank <- seq_len(nrow(eth_burden_analysis))
eth_burden_analysis$zone_percentile <- eth_burden_analysis$zone_rank / nrow(eth_burden_analysis) * 100

# Key concentration statistics for Ethiopia
eth_top_10_pct_count <- ceiling(nrow(eth_burden_analysis) * 0.10)
eth_top_20_pct_count <- ceiling(nrow(eth_burden_analysis) * 0.20)
eth_top_50_pct_count <- ceiling(nrow(eth_burden_analysis) * 0.50)

eth_top_10_pct_burden <- sum(eth_burden_analysis$mean_cattle_at_risk[1:eth_top_10_pct_count], na.rm = TRUE)
eth_top_20_pct_burden <- sum(eth_burden_analysis$mean_cattle_at_risk[1:eth_top_20_pct_count], na.rm = TRUE)
eth_top_50_pct_burden <- sum(eth_burden_analysis$mean_cattle_at_risk[1:eth_top_50_pct_count], na.rm = TRUE)

eth_top_10_pct_share <- eth_top_10_pct_burden / eth_total_cattle * 100
eth_top_20_pct_share <- eth_top_20_pct_burden / eth_total_cattle * 100
eth_top_50_pct_share <- eth_top_50_pct_burden / eth_total_cattle * 100

cat("\nETHIOPIA CONCENTRATION RESULTS:\n")
cat("Total cattle at risk:", format(eth_total_cattle, big.mark = ","), "\n")
cat("Top 10% of zones (", eth_top_10_pct_count, " zones) account for ", 
    round(eth_top_10_pct_share, 1), "% of burden\n", sep = "")
cat("Top 20% of zones (", eth_top_20_pct_count, " zones) account for ", 
    round(eth_top_20_pct_share, 1), "% of burden\n", sep = "")
cat("Top 50% of zones (", eth_top_50_pct_count, " zones) account for ", 
    round(eth_top_50_pct_share, 1), "% of burden\n", sep = "")

# Show top 10 zones
cat("\nTop 10 highest burden zones in Ethiopia:\n")
eth_top_10 <- eth_burden_analysis %>% 
  slice_head(n = 10) %>%
  dplyr::select(NAME_3, NAME_1, mean_cattle_at_risk, cumulative_percent)
print(as.data.frame(eth_top_10))

# ==============================================================================
# COMPARATIVE ANALYSIS & VISUALIZATION
# ==============================================================================

cat("\n3. COMPARATIVE ANALYSIS:\n")
cat("=======================\n")

# Gini coefficient calculation function
calculate_gini <- function(x) {
  x <- sort(x[x > 0])  # Remove zeros and sort
  n <- length(x)
  if (n == 0) return(0)
  
  # Calculate Gini coefficient
  gini <- (2 * sum((1:n) * x)) / (n * sum(x)) - (n + 1) / n
  return(gini)
}

nga_gini <- calculate_gini(nga_burden_analysis$mean_cattle_at_risk)
eth_gini <- calculate_gini(eth_burden_analysis$mean_cattle_at_risk)


# Create comprehensive visualization
cat("\n4. CREATING VISUALIZATIONS...\n")

# Lorenz curve data
nga_lorenz <- nga_burden_analysis %>%
  dplyr::mutate(
    country = "Nigeria",
    unit_percentile = lga_percentile,
    burden_cumulative = cumulative_percent / 100
  ) %>%
  dplyr::select(country, unit_percentile, burden_cumulative)

eth_lorenz <- eth_burden_analysis %>%
  dplyr::mutate(
    country = "Ethiopia", 
    unit_percentile = zone_percentile,
    burden_cumulative = cumulative_percent / 100
  ) %>%
  dplyr::select(country, unit_percentile, burden_cumulative)

# Combine for plotting
lorenz_data <- rbind(nga_lorenz, eth_lorenz) %>%
  rbind(data.frame(country = c("Nigeria", "Ethiopia"), 
                   unit_percentile = c(0, 0), 
                   burden_cumulative = c(0, 0)))

# Create Lorenz curve plot
lorenz_plot <- ggplot(lorenz_data, aes(x = unit_percentile/100, y = burden_cumulative, color = country)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", size = 0.8) +
  scale_color_manual(values = c("Nigeria" = "#E31A1C", "Ethiopia" = "#1F78B4")) +
  scale_x_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.2)) +
  labs(x = "Cumulative % of administrative units (ranked highest to lowest burden)",
    y = "Cumulative % of total infected cattle",
    color = "Country"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

# Top burden areas plot
burden_comparison <- data.frame(
  Country = rep(c("Nigeria", "Ethiopia"), each = 3),
  Percentile = rep(c("Top 10%", "Top 20%", "Top 50%"), 2),
  Share = c(nga_top_10_pct_share, nga_top_20_pct_share, nga_top_50_pct_share,
            eth_top_10_pct_share, eth_top_20_pct_share, eth_top_50_pct_share),
  Count = c(nga_top_10_pct_count, nga_top_20_pct_count, nga_top_50_pct_count,
            eth_top_10_pct_count, eth_top_20_pct_count, eth_top_50_pct_count)
)

burden_plot <- ggplot(burden_comparison, aes(x = Percentile, y = Share, fill = Country)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Nigeria" = "#E31A1C", "Ethiopia" = "#1F78B4")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 100, 20)) +
  labs(x = "Administrative unit percentile",
    y = "Share of total burden (%)",
    fill = "Country"
  ) +
  geom_text(aes(label = paste0(round(Share, 1), "%\n(n=", Count, ")")), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  ylim(0, max(burden_comparison$Share) * 1.15)

# Distribution plots
nga_dist_data <- nga_burden_analysis %>%
  dplyr::mutate(log_burden = log10(mean_cattle_at_risk + 1)) %>%
  dplyr::select(mean_cattle_at_risk, log_burden) %>%
  dplyr::mutate(country = "Nigeria")

eth_dist_data <- eth_burden_analysis %>%
  dplyr::mutate(log_burden = log10(mean_cattle_at_risk + 1)) %>%
  dplyr::select(mean_cattle_at_risk, log_burden) %>%
  dplyr::mutate(country = "Ethiopia")

dist_data <- rbind(nga_dist_data, eth_dist_data)

# Log-scale distribution
log_dist_plot <- ggplot(dist_data, aes(x = log_burden, y = country, fill = country)) +
  geom_density_ridges(alpha = 0.7, scale = 2) +
  scale_fill_manual(values = c("Nigeria" = "#E31A1C", "Ethiopia" = "#1F78B4")) +
  labs(x = expression("Infected cattle (" * log[10] * ")"),
    y = "Country"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Save plots
ggsave("Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/burden_concentration_lorenz.pdf", lorenz_plot, width = 10, height = 8, dpi = 300)
ggsave("Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/burden_concentration_percentiles.pdf", burden_plot, width = 10, height = 6, dpi = 300)
ggsave("Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/burden_distribution_ridges.pdf", log_dist_plot, width = 10, height = 6, dpi = 300)

# Create combined summary plot
combined_plot <- grid.arrange(
  lorenz_plot,
  grid.arrange(burden_plot, log_dist_plot, ncol = 1),
  ncol = 2,
  widths = c(1, 1)
)

ggsave("Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/burden_concentration_summary.pdf", combined_plot, width = 16, height = 10, dpi = 300)
