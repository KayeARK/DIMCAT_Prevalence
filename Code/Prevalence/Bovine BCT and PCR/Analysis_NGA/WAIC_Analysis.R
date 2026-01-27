rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(readr)
library(stringr)

# Set country for analysis
country_code <- "NGA"
country_name <- "Nigeria"

cat("Analyzing WAIC increase data for", country_name, "...\n")

#===============================================================================
# LOAD AND PROCESS WAIC INCREASE DATA
#===============================================================================

# Get list of all WAIC increase files
waic_files <- list.files(paste0("Code/Prevalence/Bovine BCT and PCR/WAIC_increase_", country_code), 
                         pattern = "WAIC_increase_model_.*\\.csv$", 
                         full.names = TRUE)

cat("Found", length(waic_files), "WAIC increase files\n")

# Function to load and process a single WAIC file
load_waic_file <- function(file_path) {
  model_num <- as.numeric(str_extract(basename(file_path), "\\d+"))
  data <- read_csv(file_path, show_col_types = FALSE)
  data$model <- model_num
  return(data)
}

# Load all WAIC data
cat("Loading WAIC data from all models...\n")
all_waic_data <- do.call(rbind, lapply(waic_files, load_waic_file))

# Function to rename covariates for better display
rename_covariates <- function(covariate_names) {
  covariate_names <- case_when(
    covariate_names == "pop_den" ~ "Population density",
    covariate_names == "tmax" ~ "Maximum temperature", 
    covariate_names == "tavg" ~ "Average temperature",
    covariate_names == "tmin" ~ "Minimum temperature",
    covariate_names == "human_fp" ~ "Human footprint",
    covariate_names == "elevation" ~ "Elevation",
    covariate_names == "precipitation" ~ "Precipitation",
    covariate_names == "tree" ~ "Tree cover",
    covariate_names == "grassland" ~ "Grassland",
    covariate_names == "shrub" ~ "Shrubland",
    covariate_names == "cropland" ~ "Cropland",
    covariate_names == "built" ~ "Built area",
    covariate_names == "bare" ~ "Bare ground",
    covariate_names == "water" ~ "Water",
    covariate_names == "wetland" ~ "Wetland",
    covariate_names == "mangrove" ~ "Mangrove",
    covariate_names == "moss" ~ "Moss",
    covariate_names == "cattle" ~ "Cattle density",
    covariate_names == "goat" ~ "Goat density",
    covariate_names == "horse" ~ "Horse density",
    covariate_names == "pig" ~ "Pig density",
    covariate_names == "sheep" ~ "Sheep density",
    covariate_names == "buffalo" ~ "Buffalo density",
    covariate_names == "tsetse" ~ "Tsetse presence",
    covariate_names == "tsetse_habitat" ~ "Tsetse habitat",
    covariate_names == "total_animal" ~ "Total livestock",
    TRUE ~ covariate_names  # Keep original name if no match found
  )
  return(covariate_names)
}

# Check raw covariate names before renaming
cat("Raw covariates before renaming:\n")
print(sort(unique(all_waic_data$covariates)))

# Apply covariate renaming
all_waic_data$covariates <- rename_covariates(all_waic_data$covariates)

# Check for any NA values after renaming
na_count <- sum(is.na(all_waic_data$covariates))
if(na_count > 0) {
  cat("WARNING:", na_count, "covariates became NA after renaming\n")
}

# Basic data summary
cat("Data dimensions:", dim(all_waic_data), "\n")
cat("Unique covariates after renaming:", paste(sort(unique(all_waic_data$covariates)), collapse = ", "), "\n")
cat("Number of models:", length(unique(all_waic_data$model)), "\n")

# Check if mangrove is present
mangrove_present <- "Mangrove" %in% all_waic_data$covariates
cat("Mangrove present in data:", mangrove_present, "\n")
if(mangrove_present) {
  mangrove_data <- all_waic_data[all_waic_data$covariates == "Mangrove", ]
  cat("Mangrove data points:", nrow(mangrove_data), "\n")
  cat("Mangrove WAIC range:", range(mangrove_data$waic_increase, na.rm = TRUE), "\n")
}

#===============================================================================
# DATA PROCESSING AND NORMALIZATION
#===============================================================================

# Calculate normalized WAIC increase within each model
all_waic_data <- all_waic_data %>%
  group_by(model) %>%
  mutate(
    waic_normalized = waic_increase / sum(waic_increase, na.rm = TRUE),
    waic_standardized = scale(waic_increase)[,1],
    waic_rank = rank(-waic_increase)
  ) %>%
  ungroup()

# Calculate aggregate statistics across all models
waic_summary <- all_waic_data %>%
  group_by(covariates) %>%
  summarise(
    mean_waic = mean(waic_increase, na.rm = TRUE),
    median_waic = median(waic_increase, na.rm = TRUE),
    sd_waic = sd(waic_increase, na.rm = TRUE),
    # Calculate 95% credible intervals (2.5% and 97.5% quantiles)
    ci_lower = quantile(waic_increase, 0.025, na.rm = TRUE),
    ci_upper = quantile(waic_increase, 0.975, na.rm = TRUE),
    sum_normalized = sum(waic_normalized, na.rm = TRUE),
    mean_normalized = mean(waic_normalized, na.rm = TRUE),
    mean_rank = mean(waic_rank, na.rm = TRUE),
    count = n(),
    .groups = 'drop'
  ) %>%
  # Re-normalize sum_normalized to sum to 1 across all covariates
  mutate(sum_normalized = sum_normalized / sum(sum_normalized, na.rm = TRUE)) %>%
  arrange(desc(sum_normalized))

print("WAIC Summary Statistics:")
print(waic_summary)

#===============================================================================
# VISUALIZATION 1: AGGREGATE NORMALIZED WAIC IMPORTANCE
#===============================================================================

p1 <- ggplot(waic_summary, aes(x = reorder(covariates, sum_normalized), y = sum_normalized)) +
  geom_col(fill = "#2E86AB", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = round(sum_normalized, 3)), hjust = -0.1, size = 3.5, 
            fontface = "bold", color = "#2E86AB") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "#2E86AB"),
    plot.caption = element_text(size = 9, color = "gray40", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Aggregate normalized WAIC importance -", country_name),
    x = "Covariate",
    y = "Normalized importance score (sum = 1)",
    caption = "Higher values indicate greater importance across all models"
  )

print(p1)
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/WAIC_Aggregate_Importance_", country_code, ".png"), 
       plot = p1, width = 12, height = 8, dpi = 300)

#===============================================================================
# VISUALIZATION 2: MEAN WAIC INCREASE WITH ERROR BARS
#===============================================================================

p2 <- ggplot(waic_summary, aes(x = reorder(covariates, mean_waic), y = mean_waic)) +
  geom_col(fill = "#A23B72", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.3, alpha = 0.9, color = "#732147", size = 0.8) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "#A23B72"),
    plot.caption = element_text(size = 9, color = "gray40", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Mean WAIC increase by covariate -", country_name),
    x = "Covariate", 
    y = "Mean WAIC increase with 95% credible intervals",
    caption = "Error bars show 95% credible intervals across different model configurations"
  )

print(p2)
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/WAIC_Mean_Importance_", country_code, ".png"), 
       plot = p2, width = 12, height = 8, dpi = 300)

#===============================================================================
# VISUALIZATION 3: HEATMAP OF WAIC INCREASES ACROSS MODELS
#===============================================================================

# Create heatmap data (sample subset for visualization if too many models)
if(length(unique(all_waic_data$model)) > 50) {
  sample_models <- sample(unique(all_waic_data$model), 50)
  heatmap_data <- all_waic_data %>% filter(model %in% sample_models)
  subtitle <- "Showing random sample of 50 models"
} else {
  heatmap_data <- all_waic_data
  subtitle <- paste("Showing all", length(unique(all_waic_data$model)), "models")
}

p3 <- ggplot(heatmap_data, aes(x = factor(model), y = reorder(covariates, waic_increase), 
                               fill = waic_normalized)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_viridis_c(name = "Normalized\nWAIC\nincrease", option = "plasma", 
                       trans = "sqrt", labels = scales::number_format(accuracy = 0.01)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "#E16462"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("WAIC increase heatmap -", country_name),
    subtitle = subtitle,
    x = "Model number",
    y = "Covariate"
  )

print(p3)
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/WAIC_Heatmap_", country_code, ".png"), 
       plot = p3, width = 15, height = 8, dpi = 300)

#===============================================================================
# VISUALIZATION 4: DISTRIBUTION OF WAIC INCREASES BY COVARIATE
#===============================================================================

p4 <- ggplot(all_waic_data, aes(x = reorder(covariates, waic_increase, median), 
                                y = waic_increase)) +
  geom_boxplot(fill = "#F18F01", alpha = 0.8, outlier.alpha = 0.6, 
               outlier.color = "#C73E1D", outlier.size = 1.2, 
               color = "#B8860B", size = 0.6) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "#F18F01"),
    plot.caption = element_text(size = 9, color = "gray40", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Distribution of WAIC increases -", country_name),
    x = "Covariate",
    y = "WAIC increase",
    caption = "Boxplots show variability across all model configurations"
  )

print(p4)
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/WAIC_Distribution_", country_code, ".png"), 
       plot = p4, width = 12, height = 8, dpi = 300)

#===============================================================================
# VISUALIZATION 5: RANKING CONSISTENCY ACROSS MODELS
#===============================================================================

rank_consistency <- all_waic_data %>%
  group_by(covariates) %>%
  summarise(
    mean_rank = mean(waic_rank, na.rm = TRUE),
    sd_rank = sd(waic_rank, na.rm = TRUE),
    rank_consistency = 1 / (1 + sd_rank), # Higher = more consistent ranking
    .groups = 'drop'
  ) %>%
  arrange(mean_rank)

p5 <- ggplot(rank_consistency, aes(x = mean_rank, y = rank_consistency)) +
  geom_point(size = 5, alpha = 0.8, color = "#8E44AD", stroke = 1.2, 
             fill = "#BB8FCE", shape = 21) +
  geom_text(aes(label = covariates), hjust = 0.1, vjust = -0.7, size = 3.5, 
            fontface = "bold", color = "#5B2C87") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black", face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "#8E44AD"),
    plot.caption = element_text(size = 9, color = "gray40", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  scale_x_reverse(limits = c(15, 0)) +
  labs(
    title = paste("Ranking consistency vs average rank -", country_name),
    x = "Average rank (lower = more important)",
    y = "Ranking consistency (higher = more stable)",
    caption = "Covariates in top-right quadrant are consistently important across models"
  )

print(p5)
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/WAIC_Ranking_Consistency_", country_code, ".png"), 
       plot = p5, width = 12, height = 8, dpi = 300)

#===============================================================================
# VISUALIZATION 6: TOP COVARIATES RADAR PLOT
#===============================================================================

# Select top 8 covariates by AGGREGATE IMPORTANCE (sum_normalized)
# This uses the cumulative normalized WAIC importance across all models

# Debug: show the ranking to verify selection
cat("\nTop 10 covariates by aggregate importance (sum_normalized):\n")
top_ranking <- waic_summary %>% 
  arrange(desc(sum_normalized)) %>% 
  head(10)
print(top_ranking %>% dplyr::select(covariates, sum_normalized, mean_rank))

top_covariates <- waic_summary %>% 
  top_n(8, sum_normalized) %>% 
  pull(covariates)

cat("\nSelected top 8 covariates for radar plot:\n")
cat(paste(top_covariates, collapse = ", "), "\n")

radar_data <- waic_summary %>%
  filter(covariates %in% top_covariates) %>%
  dplyr::select(covariates, mean_waic, sum_normalized, mean_rank) %>%
  mutate(
    mean_waic_scaled = scales::rescale(mean_waic, to = c(0, 1)),
    sum_normalized_scaled = scales::rescale(sum_normalized, to = c(0, 1)),
    rank_scaled = scales::rescale(-mean_rank, to = c(0, 1)) # Flip rank so high = good
  )

# Debug: verify the radar data before pivoting
cat("\nRadar data before pivoting (ordered by sum_normalized):\n")
print(radar_data %>% arrange(desc(sum_normalized)) %>% dplyr::select(covariates, sum_normalized, mean_rank))

radar_data <- radar_data %>%
  pivot_longer(cols = ends_with("_scaled"), names_to = "Metric", values_to = "value") %>%
  mutate(Metric = recode(Metric, 
                        "mean_waic_scaled" = "Mean WAIC increase",
                        "sum_normalized_scaled" = "Aggregate importance", 
                        "rank_scaled" = "Rank importance"))

p6 <- ggplot(radar_data, aes(x = Metric, y = value, group = covariates, color = covariates)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 2) +
  coord_polar() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("#E31A1C", "#FF7F00", "#1F78B4", "#33A02C", 
                                "#6A3D9A", "#B15928", "#A6CEE3", "#FDBF6F"),
                     name = "Covariate") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(
    title = paste("Top covariates performance radar -", country_name),
    caption = "Scaled metrics (0-1): outer edge = better performance"
  )

print(p6)
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/WAIC_Radar_Plot_", country_code, ".png"), 
       plot = p6, width = 10, height = 10, dpi = 300)

cat("\n=== WAIC Analysis Summary ===\n")
cat("Total models analyzed:", length(unique(all_waic_data$model)), "\n")
cat("Total covariates:", length(unique(all_waic_data$covariates)), "\n")
cat("\nTop 5 most important covariates (by aggregate normalized WAIC):\n")
print(head(waic_summary[c("covariates", "sum_normalized", "mean_rank")], 5))

cat("\nAll visualizations and summary data saved to Analysis_", country_code, "/ folder\n")