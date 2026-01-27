#===============================================================================
# DATE RANGE ANALYSIS FOR SURVEY DATA
#===============================================================================
# Simple script to analyze survey date ranges for Nigeria and Ethiopia data

library(readxl)
library(dplyr)

cat("=== DATE RANGE ANALYSIS ===\n")
cat("Loading survey data to examine temporal coverage\n\n")

# Load Nigeria BCT data
cat("Loading Nigeria BCT data...\n")
nga_bct <- read_excel("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
cat("BCT data dimensions:", dim(nga_bct), "\n")

# Load Nigeria PCR data
cat("Loading Nigeria PCR data...\n")
nga_pcr <- read_excel("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")
cat("PCR data dimensions:", dim(nga_pcr), "\n")

# Check column names for date-related columns
cat("\nAvailable columns in BCT data:\n")
print(colnames(nga_bct))

cat("\nAvailable columns in PCR data:\n")
print(colnames(nga_pcr))

# Filter for Nigeria and Ethiopia data
cat("\n=== COUNTRY FILTERING ===\n")

# Check unique countries
cat("Unique countries in BCT data:\n")
print(unique(nga_bct$Country))

cat("\nUnique countries in PCR data:\n")
print(unique(nga_pcr$Country))

# Filter Nigeria data
nga_bct_filtered <- nga_bct %>% 
  filter(Country == "Nigeria" | Country == "NG" | grepl("Nigeria", Country, ignore.case = TRUE))

nga_pcr_filtered <- nga_pcr %>% 
  filter(Country == "Nigeria" | Country == "NG" | grepl("Nigeria", Country, ignore.case = TRUE))

# Filter Ethiopia data
eth_bct_filtered <- nga_bct %>% 
  filter(Country == "Ethiopia" | Country == "ET" | grepl("Ethiopia", Country, ignore.case = TRUE))

eth_pcr_filtered <- nga_pcr %>% 
  filter(Country == "Ethiopia" | Country == "ET" | grepl("Ethiopia", Country, ignore.case = TRUE))

cat("\nCountry-specific data counts:\n")
cat("Nigeria BCT:", nrow(nga_bct_filtered), "records\n")
cat("Nigeria PCR:", nrow(nga_pcr_filtered), "records\n")
cat("Ethiopia BCT:", nrow(eth_bct_filtered), "records\n")
cat("Ethiopia PCR:", nrow(eth_pcr_filtered), "records\n")

#===============================================================================
# DATE ANALYSIS FUNCTION
#===============================================================================

analyze_date_range <- function(data, country_name, test_type) {
  cat("\n", country_name, test_type, "Date Range Analysis:\n")
  cat("===============================================\n")
  
  if(nrow(data) == 0) {
    cat("No data available\n")
    return(NULL)
  }
  
  # Check for date columns
  date_cols <- grep("year|month|date", colnames(data), ignore.case = TRUE, value = TRUE)
  cat("Available date columns:", paste(date_cols, collapse = ", "), "\n")
  
  # Analyze start years
  if("Survey_start_year" %in% colnames(data)) {
    start_years <- data$Survey_start_year[!is.na(data$Survey_start_year)]
    if(length(start_years) > 0) {
      cat("Survey start years range:", min(start_years), "to", max(start_years), "\n")
      cat("Number of records with start year:", length(start_years), "out of", nrow(data), "\n")
    } else {
      cat("No valid start years found\n")
    }
  } else {
    cat("Survey_start_year column not found\n")
  }
  
  # Analyze end years  
  if("Survey_end_year" %in% colnames(data)) {
    end_years <- data$Survey_end_year[!is.na(data$Survey_end_year)]
    if(length(end_years) > 0) {
      cat("Survey end years range:", min(end_years), "to", max(end_years), "\n")
      cat("Number of records with end year:", length(end_years), "out of", nrow(data), "\n")
    } else {
      cat("No valid end years found\n")
    }
  } else {
    cat("Survey_end_year column not found\n")
  }
  
  # Create a combined date range analysis
  if("Survey_start_year" %in% colnames(data) && "Survey_end_year" %in% colnames(data)) {
    # Get all years mentioned (both start and end)
    all_years <- c(data$Survey_start_year, data$Survey_end_year)
    all_years <- all_years[!is.na(all_years)]
    
    if(length(all_years) > 0) {
      cat("Overall temporal coverage:", min(all_years), "to", max(all_years), "\n")
      cat("Temporal span:", max(all_years) - min(all_years), "years\n")
      
      # Show distribution by decade
      decades <- floor(all_years / 10) * 10
      decade_counts <- table(decades)
      cat("Studies by decade:\n")
      for(decade in names(decade_counts)) {
        cat("  ", decade, "s:", decade_counts[decade], "year mentions\n")
      }
    }
  } else {
    cat("Cannot perform combined temporal analysis - missing year columns\n")
  }
  
  # Analyze months if available
  if("Survey_start_month" %in% colnames(data)) {
    start_months <- data$Survey_start_month[!is.na(data$Survey_start_month)]
    if(length(start_months) > 0) {
      month_counts <- table(start_months)
      cat("Survey start months distribution:\n")
      print(month_counts)
    }
  }
  
  return(list(
    start_years = if("Survey_start_year" %in% colnames(data)) data$Survey_start_year[!is.na(data$Survey_start_year)] else NULL,
    end_years = if("Survey_end_year" %in% colnames(data)) data$Survey_end_year[!is.na(data$Survey_end_year)] else NULL,
    all_years = if("Survey_start_year" %in% colnames(data) && "Survey_end_year" %in% colnames(data)) {
      all_years <- c(data$Survey_start_year, data$Survey_end_year)
      all_years[!is.na(all_years)]
    } else NULL
  ))
}

#===============================================================================
# AREA SIZE ANALYSIS FUNCTION
#===============================================================================

analyze_area_size_range <- function(data, country_name, test_type) {
  cat("\n", country_name, test_type, "Area Size Range Analysis:\n")
  cat("===============================================\n")
  
  if(nrow(data) == 0) {
    cat("No data available\n")
    return(NULL)
  }
  
  # Check for area size columns
  area_cols <- grep("area|size|Area_size", colnames(data), ignore.case = TRUE, value = TRUE)
  cat("Available area-related columns:", paste(area_cols, collapse = ", "), "\n")
  
  # Analyze Area_size specifically
  if("Area_size" %in% colnames(data)) {
    # Work with raw data as-is
    area_raw <- data$Area_size
    area_values <- area_raw[!is.na(area_raw) & area_raw != ""]
    
    cat("Total Area_size entries:", length(area_raw), "\n")
    cat("Non-empty Area_size entries:", length(area_values), "\n")
    
    if(length(area_values) > 0) {
      # Show unique values and their frequencies
      area_table <- table(area_values)
      area_table_sorted <- sort(area_table, decreasing = TRUE)
      
      cat("Unique Area_size values found (", length(area_table_sorted), " different values):\n", sep = "")
      
      # Show ALL unique values, regardless of how many
      for(i in 1:length(area_table_sorted)) {
        cat("  '", names(area_table_sorted)[i], "': ", area_table_sorted[i], " studies\n", sep = "")
      }
      
      # Summary statistics
      cat("\nArea_size summary:\n")
      cat("  Total studies with area data:", length(area_values), "\n")
      cat("  Number of unique area specifications:", length(unique(area_values)), "\n")
      cat("  Most common area specification: '", names(area_table_sorted)[1], "' (", 
          area_table_sorted[1], " studies)\n", sep = "")
      
    } else {
      cat("No Area_size values found\n")
    }
  } else {
    cat("Area_size column not found\n")
  }
  
  return(list(
    area_values = if("Area_size" %in% colnames(data)) {
      area_raw <- data$Area_size
      area_values <- area_raw[!is.na(area_raw) & area_raw != ""]
      area_values
    } else NULL,
    summary_stats = if("Area_size" %in% colnames(data)) {
      area_raw <- data$Area_size
      area_values <- area_raw[!is.na(area_raw) & area_raw != ""]
      if(length(area_values) > 0) {
        area_table <- table(area_values)
        list(
          total_entries = length(area_values),
          unique_values = length(unique(area_values)),
          most_common = names(sort(area_table, decreasing = TRUE))[1],
          most_common_count = max(area_table),
          value_frequencies = area_table
        )
      } else NULL
    } else NULL
  ))
}

#===============================================================================
# ANALYZE DATE RANGES
#===============================================================================

# Analyze date ranges for Nigeria
nga_bct_dates <- analyze_date_range(nga_bct_filtered, "Nigeria", "BCT")
nga_pcr_dates <- analyze_date_range(nga_pcr_filtered, "Nigeria", "PCR")

# Analyze date ranges for Ethiopia  
eth_bct_dates <- analyze_date_range(eth_bct_filtered, "Ethiopia", "BCT")
eth_pcr_dates <- analyze_date_range(eth_pcr_filtered, "Ethiopia", "PCR")

#===============================================================================
# ANALYZE AREA SIZE RANGES
#===============================================================================

# Analyze area size ranges for Nigeria
nga_bct_areas <- analyze_area_size_range(nga_bct_filtered, "Nigeria", "BCT")
nga_pcr_areas <- analyze_area_size_range(nga_pcr_filtered, "Nigeria", "PCR")

# Analyze area size ranges for Ethiopia
eth_bct_areas <- analyze_area_size_range(eth_bct_filtered, "Ethiopia", "BCT")
eth_pcr_areas <- analyze_area_size_range(eth_pcr_filtered, "Ethiopia", "PCR")

#===============================================================================
# COMBINED DATE RANGE SUMMARY
#===============================================================================

cat("\n\n=== COMBINED DATE RANGE SUMMARY ===\n")
cat("=====================================\n")

# Combine all years from all datasets
all_nigeria_years <- c()
all_ethiopia_years <- c()

if(!is.null(nga_bct_dates$all_years)) all_nigeria_years <- c(all_nigeria_years, nga_bct_dates$all_years)
if(!is.null(nga_pcr_dates$all_years)) all_nigeria_years <- c(all_nigeria_years, nga_pcr_dates$all_years)
if(!is.null(eth_bct_dates$all_years)) all_ethiopia_years <- c(all_ethiopia_years, eth_bct_dates$all_years)
if(!is.null(eth_pcr_dates$all_years)) all_ethiopia_years <- c(all_ethiopia_years, eth_pcr_dates$all_years)

if(length(all_nigeria_years) > 0) {
  cat("NIGERIA OVERALL TEMPORAL COVERAGE:", min(all_nigeria_years), "to", max(all_nigeria_years), "\n")
  cat("  Span:", max(all_nigeria_years) - min(all_nigeria_years), "years\n")
} else {
  cat("NIGERIA: No temporal data found\n")
}

if(length(all_ethiopia_years) > 0) {
  cat("ETHIOPIA OVERALL TEMPORAL COVERAGE:", min(all_ethiopia_years), "to", max(all_ethiopia_years), "\n") 
  cat("  Span:", max(all_ethiopia_years) - min(all_ethiopia_years), "years\n")
} else {
  cat("ETHIOPIA: No temporal data found\n")
}

# Overall dataset coverage
all_years_combined <- c(all_nigeria_years, all_ethiopia_years)
if(length(all_years_combined) > 0) {
  cat("COMBINED DATASET TEMPORAL COVERAGE:", min(all_years_combined), "to", max(all_years_combined), "\n")
  cat("  Total span:", max(all_years_combined) - min(all_years_combined), "years\n")
} else {
  cat("No temporal data found in any dataset\n")
}

cat("\nDate range analysis complete.\n")

#===============================================================================
# COMBINED AREA SIZE SUMMARY
#===============================================================================

cat("\n\n=== COMBINED AREA SIZE SUMMARY ===\n")
cat("===================================\n")

# Combine all area values from all datasets (as raw text)
all_nigeria_areas <- c()
all_ethiopia_areas <- c()

if(!is.null(nga_bct_areas$area_values)) all_nigeria_areas <- c(all_nigeria_areas, nga_bct_areas$area_values)
if(!is.null(nga_pcr_areas$area_values)) all_nigeria_areas <- c(all_nigeria_areas, nga_pcr_areas$area_values)
if(!is.null(eth_bct_areas$area_values)) all_ethiopia_areas <- c(all_ethiopia_areas, eth_bct_areas$area_values)
if(!is.null(eth_pcr_areas$area_values)) all_ethiopia_areas <- c(all_ethiopia_areas, eth_pcr_areas$area_values)

if(length(all_nigeria_areas) > 0) {
  nigeria_area_table <- table(all_nigeria_areas)
  nigeria_area_sorted <- sort(nigeria_area_table, decreasing = TRUE)
  cat("NIGERIA AREA SIZE DATA:\n")
  cat("  Total studies with area data:", length(all_nigeria_areas), "\n")
  cat("  Unique area specifications:", length(unique(all_nigeria_areas)), "\n")
  cat("  Most common area specification: '", names(nigeria_area_sorted)[1], "' (", 
      nigeria_area_sorted[1], " studies)\n", sep = "")
  cat("  Top 3 area specifications:\n")
  for(i in 1:min(3, length(nigeria_area_sorted))) {
    cat("    '", names(nigeria_area_sorted)[i], "': ", nigeria_area_sorted[i], " studies\n", sep = "")
  }
} else {
  cat("NIGERIA: No area size data found\n")
}

if(length(all_ethiopia_areas) > 0) {
  ethiopia_area_table <- table(all_ethiopia_areas)
  ethiopia_area_sorted <- sort(ethiopia_area_table, decreasing = TRUE)
  cat("ETHIOPIA AREA SIZE DATA:\n")
  cat("  Total studies with area data:", length(all_ethiopia_areas), "\n")
  cat("  Unique area specifications:", length(unique(all_ethiopia_areas)), "\n")
  cat("  Most common area specification: '", names(ethiopia_area_sorted)[1], "' (", 
      ethiopia_area_sorted[1], " studies)\n", sep = "")
  cat("  Top 3 area specifications:\n")
  for(i in 1:min(3, length(ethiopia_area_sorted))) {
    cat("    '", names(ethiopia_area_sorted)[i], "': ", ethiopia_area_sorted[i], " studies\n", sep = "")
  }
} else {
  cat("ETHIOPIA: No area size data found\n")
}

# Overall dataset area coverage (raw text analysis)
all_areas_combined <- c(all_nigeria_areas, all_ethiopia_areas)
if(length(all_areas_combined) > 0) {
  combined_area_table <- table(all_areas_combined)
  combined_area_sorted <- sort(combined_area_table, decreasing = TRUE)
  cat("COMBINED DATASET AREA SIZE DATA:\n")
  cat("  Total studies with area data:", length(all_areas_combined), "\n")
  cat("  Unique area specifications across all datasets:", length(unique(all_areas_combined)), "\n")
  cat("  Most common area specification overall: '", names(combined_area_sorted)[1], "' (", 
      combined_area_sorted[1], " studies)\n", sep = "")
  
  cat("  Top 5 area specifications overall:\n")
  for(i in 1:min(5, length(combined_area_sorted))) {
    cat("    '", names(combined_area_sorted)[i], "': ", combined_area_sorted[i], " studies\n", sep = "")
  }
} else {
  cat("No area size data found in any dataset\n")
}

cat("\nArea size analysis complete.\n")
cat("\nAll analyses completed successfully.\n")