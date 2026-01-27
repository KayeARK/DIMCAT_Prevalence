rm(list=ls())

library(dplyr)

# Set the directory containing the Cattle_at_risk files
data_dir <- "Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk/"

# Get list of all projection files
projection_files <- list.files(data_dir, pattern = "Projections_model_.*\\.csv", full.names = TRUE)
cat("Found", length(projection_files), "projection files\n")

# Initialize list to store all data
all_data <- list()

# Read all files
for (i in seq_along(projection_files)) {
  file_path <- projection_files[i]
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Add model identifier
  data$model_id <- i
  
  # Store in list
  all_data[[i]] <- data
  
  if (i %% 50 == 0) {
    cat("Processed", i, "files...\n")
  }
}

# Combine all data into one dataframe
combined_data <- do.call(rbind, all_data)

# Remove the Total rows for regional analysis (we'll compute totals separately)
regional_data <- combined_data[combined_data$region != "Total", ]

# Compute uncertainty metrics by region
uncertainty_summary <- regional_data %>%
  group_by(region) %>%
  summarise(
    # Basic statistics for mean cattle at risk
    mean_cattle_mean = mean(cattle_mean, na.rm = TRUE),
    
    # Data uncertainty: 95% CI for the mean cattle at risk across all models
    data_uncertainty_lower = quantile(cattle_mean, 0.025, na.rm = TRUE),
    data_uncertainty_upper = quantile(cattle_mean, 0.975, na.rm = TRUE),
    data_uncertainty_width = data_uncertainty_upper - data_uncertainty_lower,
    
    # Model uncertainty: Mean of the 95% CIs from each model
    model_uncertainty_lower = mean(cattle_lower, na.rm = TRUE),
    model_uncertainty_upper = mean(cattle_upper, na.rm = TRUE),  
    model_uncertainty_width = model_uncertainty_upper - model_uncertainty_lower,
    
    # Additional statistics
    median_cattle_mean = median(cattle_mean, na.rm = TRUE),
    sd_cattle_mean = sd(cattle_mean, na.rm = TRUE),
    
    # Prevalence statistics
    mean_value_mean = mean(mean_value, na.rm = TRUE),
    prevalence_data_uncertainty_lower = quantile(mean_value, 0.025, na.rm = TRUE),
    prevalence_data_uncertainty_upper = quantile(mean_value, 0.975, na.rm = TRUE),
    prevalence_model_uncertainty_lower = mean(lower_value, na.rm = TRUE),
    prevalence_model_uncertainty_upper = mean(upper_value, na.rm = TRUE),
    
    # Sample size
    n_models = n(),
    .groups = 'drop'
  )

# Compute totals across all regions
total_data <- combined_data[combined_data$region == "Total", ]

total_summary <- total_data %>%
  summarise(
    region = "Total",
    
    # Basic statistics for mean cattle at risk
    mean_cattle_mean = mean(cattle_mean, na.rm = TRUE),
    
    # Data uncertainty: 95% CI for the mean cattle at risk across all models
    data_uncertainty_lower = quantile(cattle_mean, 0.025, na.rm = TRUE),
    data_uncertainty_upper = quantile(cattle_mean, 0.975, na.rm = TRUE),
    data_uncertainty_width = data_uncertainty_upper - data_uncertainty_lower,
    
    # Model uncertainty: Mean of the 95% CIs from each model
    model_uncertainty_lower = mean(cattle_lower, na.rm = TRUE),
    model_uncertainty_upper = mean(cattle_upper, na.rm = TRUE),
    model_uncertainty_width = model_uncertainty_upper - model_uncertainty_lower,
    
    # Additional statistics
    median_cattle_mean = median(cattle_mean, na.rm = TRUE),
    sd_cattle_mean = sd(cattle_mean, na.rm = TRUE),
    
    # Prevalence statistics
    mean_value_mean = mean(mean_value, na.rm = TRUE),
    prevalence_data_uncertainty_lower = quantile(mean_value, 0.025, na.rm = TRUE),
    prevalence_data_uncertainty_upper = quantile(mean_value, 0.975, na.rm = TRUE),
    prevalence_model_uncertainty_lower = mean(lower_value, na.rm = TRUE),
    prevalence_model_uncertainty_upper = mean(upper_value, na.rm = TRUE),
    
    # Sample size
    n_models = n(),
    .groups = 'drop'
  )

# Combine regional and total summaries
final_summary <- rbind(uncertainty_summary, total_summary)

# Round columns appropriately for readability
# Cattle numbers (including uncertainty bounds): round to 0 decimal places
cattle_cols <- grep("cattle|data_uncertainty_lower|data_uncertainty_upper|model_uncertainty_lower|model_uncertainty_upper|data_uncertainty_width|model_uncertainty_width", names(final_summary), value = TRUE)
final_summary[cattle_cols] <- lapply(final_summary[cattle_cols], function(x) round(x, 0))

# Prevalence values: round to 3 decimal places
prevalence_cols <- grep("value|prevalence", names(final_summary), value = TRUE)
prevalence_cols <- setdiff(prevalence_cols, cattle_cols)  # Remove any overlap with cattle cols
final_summary[prevalence_cols] <- lapply(final_summary[prevalence_cols], function(x) round(x, 3))

# Other numeric columns (like sd): round to 2 decimal places
other_numeric_cols <- names(final_summary)[sapply(final_summary, is.numeric)]
other_numeric_cols <- setdiff(other_numeric_cols, c(cattle_cols, prevalence_cols, "n_models"))
final_summary[other_numeric_cols] <- lapply(final_summary[other_numeric_cols], function(x) round(x, 2))

# Add uncertainty ratio (data uncertainty / model uncertainty)
final_summary <- final_summary %>%
  mutate(
    uncertainty_ratio = round(data_uncertainty_width / model_uncertainty_width, 2)
  )



# Create a focused summary with key columns
key_summary <- final_summary %>%
  dplyr::select(region, mean_value_mean, mean_cattle_mean, 
         data_uncertainty_lower, data_uncertainty_upper, data_uncertainty_width,
         model_uncertainty_lower, model_uncertainty_upper, model_uncertainty_width,
         uncertainty_ratio)

# Save the focused summary  
focused_output_file <- "Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Ethiopia_cattle_uncertainty.csv"
write.csv(key_summary, focused_output_file, row.names = FALSE)

cat("Focused summary saved to:", focused_output_file, "\n")

# Print interpretation
cat("\n=== INTERPRETATION ===\n")
cat("Data uncertainty: Variability across ensemble model realizations\n")
cat("Model uncertainty: Average prediction interval from individual models\n")
cat("Uncertainty ratio > 1: Data uncertainty dominates\n")
cat("Uncertainty ratio < 1: Model uncertainty dominates\n")