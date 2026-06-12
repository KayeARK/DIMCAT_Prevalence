rm(list=ls())

library(dplyr)
library(sf)
library(geodata)

# Load state-level correction factors
correction_factors <- read.csv(
  "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_correction_factors.csv",
  stringsAsFactors = FALSE
) %>%
  dplyr::select(state, correction_factor)

# Remove national total row if present
correction_factors <- correction_factors %>%
  filter(state != "Total")


#### LOOK UP TABLE FOR LGA TO STATE MAPPING ####
# Nigeria admin level 2 (LGAs)
nga_lga <- gadm(
  country = "NGA",
  level = 2,
  path = tempdir()
)

nga_lga <- st_as_sf(nga_lga)

lga_lookup <- nga_lga %>%
  st_drop_geometry() %>%
  dplyr::select(
    lga = NAME_2,
    state = NAME_1
  ) %>%
  distinct()

# Set the directory containing the Cattle_at_risk files
data_dir <- "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine/"

# Get list of all projection files
projection_files <- list.files(data_dir, pattern = "Projections_model_.*\\.csv", full.names = TRUE)
cat("Found", length(projection_files), "projection files\n")

# Initialize list to store all data
all_data <- list()

# Read all files
for (i in seq_along(projection_files)) {
  file_path <- projection_files[i]
  data <- read.csv(file_path, stringsAsFactors = FALSE)

# Separate Total row
total_row_original <- data[data$state == "Total", ]
data <- data[data$state != "Total", ]

# Rename LGA column
data <- data %>%
  rename(lga = state)

# Attach parent state
data <- data %>%
  left_join(lga_lookup, by = "lga")

if(any(is.na(data$state))){
  stop(
    "Missing state assignment for: ",
    paste(unique(data$lga[is.na(data$state)]),
          collapse = ", ")
  )
}

# Attach correction factor
data <- data %>%
  left_join(correction_factors, by = "state")

if(any(is.na(data$correction_factor))){
  stop(
    "Missing correction factor for: ",
    paste(unique(data$state[is.na(data$correction_factor)]),
          collapse = ", ")
  )
}

# Apply correction
data <- data %>%
  mutate(
    cattle_mean  = cattle_mean  * correction_factor,
    cattle_lower = cattle_lower * correction_factor,
    cattle_upper = cattle_upper * correction_factor
  )

# Remove temporary columns
data <- data %>%
  dplyr::select(-state, -correction_factor) %>%
  rename(state = lga)
  
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

# Remove the Total rows for state analysis (we'll compute totals separately)
state_data <- combined_data[combined_data$state != "Total", ]

# Compute uncertainty metrics by state
uncertainty_summary <- state_data %>%
  group_by(state) %>%
  summarise(
    # Basic statistics for mean cattle at risk
    mean_value_mean = mean(mean_value, na.rm = TRUE),
    mean_cattle_mean = mean(cattle_mean, na.rm = TRUE),
    
    # Data uncertainty: 95% CI for the mean cattle at risk across all models
    data_uncertainty_lower = quantile(cattle_mean, 0.025, na.rm = TRUE),
    data_uncertainty_upper = quantile(cattle_mean, 0.975, na.rm = TRUE),
    data_uncertainty_width = data_uncertainty_upper - data_uncertainty_lower,
    
    # Model uncertainty: Mean of the 95% CIs from each model
    model_uncertainty_lower = mean(cattle_lower, na.rm = TRUE),
    model_uncertainty_upper = mean(cattle_upper, na.rm = TRUE),  
    model_uncertainty_width = model_uncertainty_upper - model_uncertainty_lower,
    
    
    # Sample size
    n_models = n(),
    .groups = 'drop'
  )

# Compute corrected national totals from LGA estimates

total_summary <- combined_data %>%
  group_by(model_id) %>%
  summarise(
    mean_value = mean(mean_value, na.rm = TRUE),
    cattle_mean = sum(cattle_mean, na.rm = TRUE),
    cattle_lower = sum(cattle_lower, na.rm = TRUE),
    cattle_upper = sum(cattle_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  summarise(
    state = "Total",

    mean_value_mean = mean(mean_value, na.rm = TRUE),
    mean_cattle_mean = mean(cattle_mean, na.rm = TRUE),

    data_uncertainty_lower = quantile(cattle_mean, 0.025, na.rm = TRUE),
    data_uncertainty_upper = quantile(cattle_mean, 0.975, na.rm = TRUE),
    data_uncertainty_width = data_uncertainty_upper - data_uncertainty_lower,

    model_uncertainty_lower = mean(cattle_lower, na.rm = TRUE),
    model_uncertainty_upper = mean(cattle_upper, na.rm = TRUE),
    model_uncertainty_width = model_uncertainty_upper - model_uncertainty_lower,

    n_models = n()
  )

# Combine state and total summaries
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
  dplyr::select(state, mean_value_mean, mean_cattle_mean, 
         data_uncertainty_lower, data_uncertainty_upper, data_uncertainty_width,
         model_uncertainty_lower, model_uncertainty_upper, model_uncertainty_width,
         uncertainty_ratio)

# Save the focused summary  
focused_output_file <- "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Nigeria_cattle_uncertainty_fine.csv"
write.csv(key_summary, focused_output_file, row.names = FALSE)

cat("Focused summary saved to:", focused_output_file, "\n")

# Print interpretation
cat("\n=== INTERPRETATION ===\n")
cat("Data uncertainty: Variability across ensemble model realizations\n")
cat("Model uncertainty: Average prediction interval from individual models\n")
cat("Uncertainty ratio > 1: Data uncertainty dominates\n")
cat("Uncertainty ratio < 1: Model uncertainty dominates\n")