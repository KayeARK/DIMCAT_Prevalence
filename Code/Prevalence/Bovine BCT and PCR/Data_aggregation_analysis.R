rm(list=ls())

library(readxl)
library(dplyr)

cat("=== DATA AGGREGATION ANALYSIS ===\n")
cat("Aggregating data by location and reporting datapoints by test type\n\n")

#===============================================================================
# NIGERIA DATA ANALYSIS
#===============================================================================

cat("NIGERIA ANALYSIS:\n")
cat("-----------------\n")

# Load Nigeria BCT data
cat("Loading Nigeria BCT data...\n")
nga_bct_data <- read_excel("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")

# Check column names and structure
cat("BCT data columns:", paste(colnames(nga_bct_data), collapse = ", "), "\n")

# Filter for Nigeria data only
nga_bct_filtered <- nga_bct_data %>% 
  filter(Country == "Nigeria" | Country == "NG" | grepl("Nigeria", Country, ignore.case = TRUE))

cat("Nigeria BCT raw datapoints:", nrow(nga_bct_filtered), "\n")

# Aggregate Nigeria BCT data by location
if(nrow(nga_bct_filtered) > 0) {
  # Check what location columns are available
  location_cols <- grep("location|lat|lon|coord|site|place", colnames(nga_bct_filtered), ignore.case = TRUE, value = TRUE)
  cat("Available location columns:", paste(location_cols, collapse = ", "), "\n")
  
  # Check for sample size and positive columns
  sample_cols <- grep("sample|tested|examined|n_", colnames(nga_bct_filtered), ignore.case = TRUE, value = TRUE)
  positive_cols <- grep("positive|infected|cases", colnames(nga_bct_filtered), ignore.case = TRUE, value = TRUE)
  
  cat("Sample size columns:", paste(sample_cols, collapse = ", "), "\n")
  cat("Positive columns:", paste(positive_cols, collapse = ", "), "\n")
  
  # Show structure
  cat("\nNigeria BCT data structure:\n")
  print(str(nga_bct_filtered))
}

# Load Nigeria PCR data  
cat("\nLoading Nigeria PCR data...\n")
nga_pcr_data <- read_excel("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")

# Filter for Nigeria data only
nga_pcr_filtered <- nga_pcr_data %>% 
  filter(Country == "Nigeria" | Country == "NG" | grepl("Nigeria", Country, ignore.case = TRUE))

cat("Nigeria PCR raw datapoints:", nrow(nga_pcr_filtered), "\n")

# Aggregate Nigeria PCR data by location
if(nrow(nga_pcr_filtered) > 0) {
  # Check what location columns are available
  location_cols_pcr <- grep("location|lat|lon|coord|site|place", colnames(nga_pcr_filtered), ignore.case = TRUE, value = TRUE)
  cat("Available location columns:", paste(location_cols_pcr, collapse = ", "), "\n")
  
  # Check for sample size and positive columns
  sample_cols_pcr <- grep("sample|tested|examined|n_", colnames(nga_pcr_filtered), ignore.case = TRUE, value = TRUE)
  positive_cols_pcr <- grep("positive|infected|cases", colnames(nga_pcr_filtered), ignore.case = TRUE, value = TRUE)
  
  cat("Sample size columns:", paste(sample_cols_pcr, collapse = ", "), "\n")
  cat("Positive columns:", paste(positive_cols_pcr, collapse = ", "), "\n")
  
  # Show structure
  cat("\nNigeria PCR data structure:\n")
  print(str(nga_pcr_filtered))
}

#===============================================================================
# ETHIOPIA DATA ANALYSIS
#===============================================================================

cat("\n\nETHIOPIA ANALYSIS:\n")
cat("------------------\n")

# Filter for Ethiopia data from both datasets
eth_bct_filtered <- nga_bct_data %>% 
  filter(Country == "Ethiopia" | Country == "ET" | grepl("Ethiopia", Country, ignore.case = TRUE))

cat("Ethiopia BCT raw datapoints:", nrow(eth_bct_filtered), "\n")

eth_pcr_filtered <- nga_pcr_data %>% 
  filter(Country == "Ethiopia" | Country == "ET" | grepl("Ethiopia", Country, ignore.case = TRUE))

cat("Ethiopia PCR raw datapoints:", nrow(eth_pcr_filtered), "\n")

# Show structure if data exists
if(nrow(eth_bct_filtered) > 0) {
  cat("\nEthiopia BCT data structure:\n")
  print(str(eth_bct_filtered))
}

if(nrow(eth_pcr_filtered) > 0) {
  cat("\nEthiopia PCR data structure:\n")
  print(str(eth_pcr_filtered))
}

cat("\n=== RAW DATA SUMMARY ===\n")
cat("Nigeria BCT raw datapoints:", nrow(nga_bct_filtered), "\n")
cat("Nigeria PCR raw datapoints:", nrow(nga_pcr_filtered), "\n")
cat("Ethiopia BCT raw datapoints:", nrow(eth_bct_filtered), "\n") 
cat("Ethiopia PCR raw datapoints:", nrow(eth_pcr_filtered), "\n")

#===============================================================================
# PERFORM LOCATION-BASED AGGREGATION
#===============================================================================

cat("\n=== LOCATION-BASED AGGREGATION ===\n")

# Function to perform aggregation by location
aggregate_by_location <- function(data, test_type, country) {
  cat("\n", country, test_type, "aggregation:\n")
  
  if(nrow(data) == 0) {
    cat("No data available\n")
    return(NULL)
  }
  
  # Check for duplicate locations
  unique_locations <- unique(data$Location)
  cat("Unique locations before aggregation:", length(unique_locations), "\n")
  cat("Total datapoints before aggregation:", nrow(data), "\n")
  
  # Find locations with multiple datapoints
  location_counts <- table(data$Location)
  multiple_datapoints <- location_counts[location_counts > 1]
  
  if(length(multiple_datapoints) > 0) {
    cat("Locations with multiple datapoints:\n")
    print(multiple_datapoints)
  } else {
    cat("No locations with multiple datapoints found\n")
  }
  
  # Perform aggregation
  aggregated_data <- data %>%
    group_by(Location, Latitude, Longitude) %>%
    summarise(
      total_tested = sum(Number_of_animal_tested, na.rm = TRUE),
      total_positive = sum(Number_of_infections, na.rm = TRUE),
      n_studies = n(),
      prevalence = total_positive / total_tested * 100,
      .groups = 'drop'
    ) %>%
    filter(total_tested > 0)  # Remove locations with no valid tests
  
  cat("Unique locations after aggregation:", nrow(aggregated_data), "\n")
  cat("Total tests after aggregation:", sum(aggregated_data$total_tested), "\n")
  cat("Total positives after aggregation:", sum(aggregated_data$total_positive), "\n")
  
  # Check if any aggregation actually occurred
  studies_combined <- sum(aggregated_data$n_studies > 1)
  cat("Locations where studies were combined:", studies_combined, "\n")
  
  return(aggregated_data)
}

# Aggregate Nigeria data
nga_bct_agg <- aggregate_by_location(nga_bct_filtered, "BCT", "Nigeria")
nga_pcr_agg <- aggregate_by_location(nga_pcr_filtered, "PCR", "Nigeria")

# Aggregate Ethiopia data  
eth_bct_agg <- aggregate_by_location(eth_bct_filtered, "BCT", "Ethiopia")
eth_pcr_agg <- aggregate_by_location(eth_pcr_filtered, "PCR", "Ethiopia")

#===============================================================================
# FINAL SUMMARY
#===============================================================================

cat("\n=== FINAL SUMMARY ===\n")
cat("NIGERIA:\n")
cat("BCT: ", nrow(nga_bct_filtered), " raw datapoints ->", ifelse(is.null(nga_bct_agg), 0, nrow(nga_bct_agg)), "aggregated locations\n")
cat("PCR: ", nrow(nga_pcr_filtered), " raw datapoints ->", ifelse(is.null(nga_pcr_agg), 0, nrow(nga_pcr_agg)), "aggregated locations\n")

cat("\nETHIOPIA:\n")
cat("BCT: ", nrow(eth_bct_filtered), " raw datapoints ->", ifelse(is.null(eth_bct_agg), 0, nrow(eth_bct_agg)), "aggregated locations\n")
cat("PCR: ", nrow(eth_pcr_filtered), " raw datapoints ->", ifelse(is.null(eth_pcr_agg), 0, nrow(eth_pcr_agg)), "aggregated locations\n")

# Check total samples consistency
if(!is.null(nga_bct_agg)) {
  cat("\nNigeria BCT total samples: Raw =", sum(nga_bct_filtered$Number_of_animal_tested, na.rm = TRUE), 
      ", Aggregated =", sum(nga_bct_agg$total_tested), "\n")
}
if(!is.null(nga_pcr_agg)) {
  cat("Nigeria PCR total samples: Raw =", sum(nga_pcr_filtered$Number_of_animal_tested, na.rm = TRUE), 
      ", Aggregated =", sum(nga_pcr_agg$total_tested), "\n")
}
if(!is.null(eth_bct_agg)) {
  cat("Ethiopia BCT total samples: Raw =", sum(eth_bct_filtered$Number_of_animal_tested, na.rm = TRUE), 
      ", Aggregated =", sum(eth_bct_agg$total_tested), "\n")
}
if(!is.null(eth_pcr_agg)) {
  cat("Ethiopia PCR total samples: Raw =", sum(eth_pcr_filtered$Number_of_animal_tested, na.rm = TRUE), 
      ", Aggregated =", sum(eth_pcr_agg$total_tested), "\n")
}