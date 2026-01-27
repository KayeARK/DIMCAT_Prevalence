library(ggplot2)
library(dplyr)
library(readr)
library(sf)
library(geodata)
library(viridis)
library(gridExtra)
library(grid)
library(gtable)
library(readxl)
library(ggspatial)
library(cowplot)

# ==============================================================================
# CI WIDTH COMPARISON: Analysis_NGA vs Analysis_NGA_all_PCR
# ==============================================================================
# This script compares the 95% CI widths between the original Nigeria analysis
# (BCT + PCR) and the all-PCR analysis to show percentage changes in uncertainty
# ==============================================================================

cat("=== CI WIDTH COMPARISON ANALYSIS ===\n")
cat("Comparing uncertainty between Analysis_NGA and Analysis_NGA_all_PCR\n\n")

# Get Nigerian administrative boundaries at LGA level (level 2)
nigeria_lgas <- gadm(country = "NGA", level = 2, path = tempdir())
nigeria_lgas_sf <- st_as_sf(nigeria_lgas)

# Also get state boundaries for context
nigeria_states <- gadm(country = "NGA", level = 1, path = tempdir()) 
nigeria_states_sf <- st_as_sf(nigeria_states)

# Get neighboring countries for geographic context
cat("Loading neighboring countries for geographic context...\n")
neighboring_countries <- c("BEN", "NER", "TCD", "CMR", "CAF")
neighbor_sf_list <- list()
neighboring_countries_sf <- NULL  # Initialize variable

for (country in neighboring_countries) {
  tryCatch({
    country_data <- gadm(country = country, level = 0, path = tempdir())
    neighbor_sf_list[[country]] <- st_as_sf(country_data)
    cat("Loaded", country, "\n")
  }, error = function(e) {
    cat("Failed to load", country, ":", e$message, "\n")
  })
}

# Combine all neighboring countries into one sf object
if(length(neighbor_sf_list) > 0) {
  neighboring_countries_full <- do.call(rbind, neighbor_sf_list)
  
  # Validate and repair geometries to avoid intersection errors
  cat("Validating and repairing geometries...\n")
  neighboring_countries_full <- st_make_valid(neighboring_countries_full)
  nigeria_states_valid <- st_make_valid(nigeria_states_sf)
  
  # Create a buffer around Nigeria to crop neighboring countries to border regions only
  nigeria_buffer <- st_buffer(st_union(nigeria_states_valid), dist = 1.5)  # ~150km buffer
  nigeria_buffer <- st_make_valid(nigeria_buffer)
  
  # Use st_crop instead of st_intersection for more robust clipping
  tryCatch({
    neighboring_countries_sf <- st_crop(neighboring_countries_full, nigeria_buffer)
    cat("Successfully loaded and cropped", nrow(neighboring_countries_sf), "neighboring country border regions\n")
  }, error = function(e) {
    cat("st_crop failed, trying alternative approach:", e$message, "\n")
    # Fallback: use a simpler bbox-based crop
    nigeria_bbox <- st_bbox(nigeria_buffer)
    neighboring_countries_sf <<- st_crop(neighboring_countries_full, nigeria_bbox)
    cat("Successfully loaded and cropped using bbox approach\n")
  })
  
} else {
  neighboring_countries_sf <- NULL
  cat("No neighboring countries loaded\n")
}

# Load continental Africa for inset map
cat("Loading continental Africa for inset map...\n")
africa_countries_sf <- NULL
nigeria_country_sf <- NULL
tryCatch({
  # Get all African countries
  africa_iso3 <- c("DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CMR", "CPV", "CAF", 
                   "TCD", "COM", "COG", "COD", "CIV", "DJI", "EGY", "GNQ", "ERI", 
                   "ETH", "GAB", "GMB", "GHA", "GIN", "GNB", "KEN", "LSO", "LBR", 
                   "LBY", "MDG", "MWI", "MLI", "MRT", "MUS", "MAR", "MOZ", "NAM", 
                   "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", 
                   "SSD", "SDN", "SWZ", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE")
  
  africa_list <- list()
  for(country_code in africa_iso3) {
    tryCatch({
      country_data <- gadm(country = country_code, level = 0, path = tempdir())
      africa_list[[country_code]] <- st_as_sf(country_data)
    }, error = function(e) {
      # Silently skip countries that fail to load
    })
  }
  
  if(length(africa_list) > 0) {
    africa_countries_sf <- do.call(rbind, africa_list)
    # Separate Nigeria for highlighting
    nigeria_country_sf <- africa_countries_sf %>% filter(GID_0 == "NGA")
    cat("Successfully loaded continental Africa with", nrow(africa_countries_sf), "countries\n")
  }
}, error = function(e) {
  cat("Failed to load continental Africa:", e$message, "\n")
})

# ==============================================================================
# LOAD AND PROCESS PROJECTION DATA FROM BOTH ANALYSES
# ==============================================================================

# Function to load and process projection files from a directory
load_ensemble_projections <- function(projections_dir, analysis_name) {
  cat("\n=== Loading", analysis_name, "projections ===\n")
  projection_files <- list.files(projections_dir, pattern = "Projections_model_.*\\.csv", full.names = TRUE)
  cat("Found", length(projection_files), "projection files\n")
  
  # Function to load and process a single projection file
  load_projection_file <- function(file_path) {
    tryCatch({
      data <- read.csv(file_path)
      data$model_id <- basename(file_path)
      return(data)
    }, error = function(e) {
      cat("Error loading", file_path, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Load all projection files
  all_projections <- do.call(rbind, lapply(projection_files, load_projection_file))
  all_projections <- all_projections[!is.null(all_projections), ]
  
  cat("Total rows loaded:", nrow(all_projections), "\n")
  cat("Unique models:", length(unique(all_projections$model_id)), "\n")
  
  # Average across all models for each coordinate and variable combination
  ensemble_projections <- all_projections %>%
    dplyr::group_by(.data$Latitude, .data$Longitude, .data$variable) %>%
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      n_models = n(),
      .groups = 'drop'
    )
  
  cat("Ensemble projections created with", nrow(ensemble_projections), "unique coordinate-variable combinations\n")
  
  return(ensemble_projections)
}

# Load projections from both analyses
ensemble_original <- load_ensemble_projections("Code/Prevalence/Bovine BCT and PCR/Projections_NGA/", "Original (BCT+PCR)")
ensemble_all_pcr <- load_ensemble_projections("Code/Prevalence/Bovine BCT and PCR/Projections_NGA_all_PCR/", "All-PCR")

# ==============================================================================
# CALCULATE CI WIDTHS FOR BOTH ANALYSES
# ==============================================================================

# Function to calculate CI width from ensemble projections
calculate_ci_widths <- function(ensemble_data, analysis_name) {
  cat("\n=== Calculating CI widths for", analysis_name, "===\n")
  
  # Extract 2.5th and 97.5th percentiles
  lower <- ensemble_data %>%
    filter(.data$variable == "2.5th percentile") %>%
    dplyr::select(.data$Longitude, .data$Latitude, lower_value = .data$value)
  
  upper <- ensemble_data %>%
    filter(.data$variable == "97.5th percentile") %>%
    dplyr::select(.data$Longitude, .data$Latitude, upper_value = .data$value)
  
  mean_val <- ensemble_data %>%
    filter(.data$variable == "Mean") %>%
    dplyr::select(.data$Longitude, .data$Latitude, mean_value = .data$value)
  
  # Join and calculate CI width
  ci_data <- lower %>%
    left_join(upper, by = c("Longitude", "Latitude")) %>%
    left_join(mean_val, by = c("Longitude", "Latitude")) %>%
    mutate(
      ci_width = .data$upper_value - .data$lower_value,
      relative_ci_width = .data$ci_width / .data$mean_value  # Relative to mean for comparison
    ) %>%
    filter(!is.na(.data$ci_width), !is.na(.data$mean_value), .data$mean_value > 0)
  
  cat("CI widths calculated for", nrow(ci_data), "locations\n")
  cat("Mean CI width:", round(mean(ci_data$ci_width), 4), "\n")
  cat("Mean relative CI width:", round(mean(ci_data$relative_ci_width), 4), "\n")
  
  return(ci_data)
}

# Calculate CI widths for both analyses
ci_original <- calculate_ci_widths(ensemble_original, "Original (BCT+PCR)")
ci_all_pcr <- calculate_ci_widths(ensemble_all_pcr, "All-PCR")

# ==============================================================================
# COMPARE CI WIDTHS AND CALCULATE PERCENTAGE CHANGES
# ==============================================================================

cat("\n=== Comparing CI widths between analyses ===\n")

# Function to process projections data (swap coordinates if needed and validate)
process_projections_coords <- function(proj_data, analysis_name = "") {
  cat("Processing coordinates for", analysis_name, "...\n")
  
  # Check coordinate ranges before any processing
  lat_range <- range(proj_data$Latitude, na.rm = TRUE)
  lon_range <- range(proj_data$Longitude, na.rm = TRUE)
  
  cat("  Original - Latitude range:", round(lat_range, 2), "\n")
  cat("  Original - Longitude range:", round(lon_range, 2), "\n")
  
  # Nigeria should have roughly: Latitude 4-14°N, Longitude 3-15°E
  # If coordinates are outside reasonable bounds, they might be swapped
  lat_valid <- lat_range[1] >= 3 && lat_range[2] <= 15
  lon_valid <- lon_range[1] >= 2 && lon_range[2] <= 16
  
  if(!lat_valid || !lon_valid) {
    cat("  Coordinates appear to be outside Nigeria bounds, attempting to correct...\n")
    
    # Try swapping
    temp <- proj_data$Longitude
    proj_data$Longitude <- proj_data$Latitude
    proj_data$Latitude <- temp
    
    # Check new ranges
    new_lat_range <- range(proj_data$Latitude, na.rm = TRUE)
    new_lon_range <- range(proj_data$Longitude, na.rm = TRUE)
    
    cat("  After swap - Latitude range:", round(new_lat_range, 2), "\n")
    cat("  After swap - Longitude range:", round(new_lon_range, 2), "\n")
    
    # Verify the swap improved things
    new_lat_valid <- new_lat_range[1] >= 3 && new_lat_range[2] <= 15
    new_lon_valid <- new_lon_range[1] >= 2 && new_lon_range[2] <= 16
    
    if(!new_lat_valid || !new_lon_valid) {
      cat("  Warning: Coordinates still outside expected Nigeria bounds after swap!\n")
    } else {
      cat("  Coordinate swap successful - now within Nigeria bounds\n")
    }
  } else {
    cat("  Coordinates appear correct - no swap needed\n")
  }
  
  return(proj_data)
}

# Process coordinates for both datasets
ci_original <- process_projections_coords(ci_original, "Original (BCT+PCR)")
ci_all_pcr <- process_projections_coords(ci_all_pcr, "All-PCR")

# Join the two datasets to compare CI widths at matching locations
ci_comparison <- ci_original %>%
  dplyr::select(.data$Longitude, .data$Latitude, 
                original_ci_width = .data$ci_width,
                original_relative_ci_width = .data$relative_ci_width,
                original_mean = .data$mean_value) %>%
  inner_join(
    ci_all_pcr %>% 
      dplyr::select(.data$Longitude, .data$Latitude, 
                    all_pcr_ci_width = .data$ci_width,
                    all_pcr_relative_ci_width = .data$relative_ci_width,
                    all_pcr_mean = .data$mean_value),
    by = c("Longitude", "Latitude")
  ) %>%
  mutate(
    # Calculate percentage change in CI width
    ci_width_change_pct = ((.data$all_pcr_ci_width - .data$original_ci_width) / .data$original_ci_width) * 100,
    ci_width_change_abs = .data$all_pcr_ci_width - .data$original_ci_width,
    
    # Calculate percentage change in relative CI width
    relative_ci_change_pct = ((.data$all_pcr_relative_ci_width - .data$original_relative_ci_width) / .data$original_relative_ci_width) * 100,
    
    # Calculate mean value change
    mean_change_pct = ((.data$all_pcr_mean - .data$original_mean) / .data$original_mean) * 100,
    
    # Log-scale analysis for better handling of extreme changes
    log_abs_change = log10(abs(.data$ci_width_change_pct) + 1),  # Add 1 to handle 0% changes
    log_ratio = log10(abs(.data$all_pcr_ci_width / .data$original_ci_width)),  # Log ratio of CI widths
    
    # Categorize changes
    change_category = case_when(
      .data$ci_width_change_pct < -20 ~ "Large Decrease (>20%)",
      .data$ci_width_change_pct < -10 ~ "Moderate Decrease (10-20%)",
      .data$ci_width_change_pct < -5 ~ "Small Decrease (5-10%)",
      .data$ci_width_change_pct >= -5 & .data$ci_width_change_pct <= 5 ~ "No Change (±5%)",
      .data$ci_width_change_pct > 5 & .data$ci_width_change_pct <= 10 ~ "Small Increase (5-10%)",
      .data$ci_width_change_pct > 10 & .data$ci_width_change_pct <= 20 ~ "Moderate Increase (10-20%)",
      .data$ci_width_change_pct > 20 ~ "Large Increase (>20%)",
      TRUE ~ "Other"
    ),
    
    # Log-scale based categories
    log_change_category = case_when(
      abs(.data$ci_width_change_pct) > 100 ~ "Extreme (>100%)",
      abs(.data$ci_width_change_pct) > 50 ~ "Very High (50-100%)",
      abs(.data$ci_width_change_pct) > 20 ~ "High (20-50%)",
      abs(.data$ci_width_change_pct) > 10 ~ "Moderate (10-20%)",
      TRUE ~ "Small (≤10%)"
    ),
    
    log_change_category = factor(.data$log_change_category,
                                levels = c("Small (≤10%)", "Moderate (10-20%)", 
                                          "High (20-50%)", "Very High (50-100%)", "Extreme (>100%)"))
  )

cat("Comparison dataset created with", nrow(ci_comparison), "matching locations\n")

# Validate that we have reasonable CI width values
original_ci_summary <- summary(ci_comparison$original_ci_width)
all_pcr_ci_summary <- summary(ci_comparison$all_pcr_ci_width)

cat("\nCI Width Validation:\n")
cat("Original Analysis CI widths:\n")
print(original_ci_summary)
cat("\nAll-PCR Analysis CI widths:\n")
print(all_pcr_ci_summary)

# Check for any extreme values that might indicate errors
extreme_changes <- ci_comparison %>%
  filter(abs(.data$ci_width_change_pct) > 500)  # More than 500% change

if(nrow(extreme_changes) > 0) {
  cat("\nWarning:", nrow(extreme_changes), "locations have extreme CI width changes (>500%)\n")
  cat("This might indicate data quality issues\n")
}

cat("\nCI Width Change Summary:\n")
cat("Mean CI width change:", round(mean(ci_comparison$ci_width_change_pct), 2), "%\n")
cat("Median CI width change:", round(median(ci_comparison$ci_width_change_pct), 2), "%\n")
cat("Standard deviation of changes:", round(sd(ci_comparison$ci_width_change_pct), 2), "%\n")

cat("\nLog-scale analysis:\n")
cat("Mean log(absolute change + 1):", round(mean(ci_comparison$log_abs_change), 3), "\n")
cat("Median log(absolute change + 1):", round(median(ci_comparison$log_abs_change), 3), "\n")
cat("Mean log(CI width ratio):", round(mean(ci_comparison$log_ratio), 3), "\n")

# Identify statistical outliers using log-scale
log_mean <- mean(ci_comparison$log_abs_change)
log_sd <- sd(ci_comparison$log_abs_change)
outlier_threshold <- log_mean + 2 * log_sd

statistical_outliers <- ci_comparison %>% filter(.data$log_abs_change > outlier_threshold)
cat("Statistical outliers (>2 SD in log-scale):", nrow(statistical_outliers), 
    "(", round(nrow(statistical_outliers)/nrow(ci_comparison)*100, 1), "%)\n")

# Print summary of changes
change_summary <- ci_comparison %>%
  count(.data$change_category) %>%
  mutate(percentage = round(.data$n / sum(.data$n) * 100, 1)) %>%
  arrange(desc(.data$n))

cat("\nCI Width Change Summary:\n")
print(change_summary)

# ==============================================================================
# AGGREGATE BY LGA FOR CHOROPLETH MAPPING
# ==============================================================================

cat("\n=== Aggregating by LGA ===\n")

# Function to aggregate data by LGA using spatial join
aggregate_by_lga <- function(comparison_data, value_column, analysis_name = "") {
  cat("Aggregating", analysis_name, "data by LGA...\n")
  
  # Convert to spatial points
  comparison_sf <- st_as_sf(comparison_data, 
                           coords = c("Longitude", "Latitude"), 
                           crs = 4326)
  
  # Spatial join with LGA polygons
  joined_data <- st_join(comparison_sf, nigeria_lgas_sf)
  
  # Aggregate by LGA
  lga_summary <- joined_data %>%
    st_drop_geometry() %>%
    dplyr::filter(!is.na(.data$GID_2)) %>%
    dplyr::group_by(.data$GID_2, .data$NAME_2, .data$NAME_1) %>%
    dplyr::summarise(
      lga_value = mean(.data[[value_column]], na.rm = TRUE),
      median_value = median(.data[[value_column]], na.rm = TRUE),
      min_value = min(.data[[value_column]], na.rm = TRUE),
      max_value = max(.data[[value_column]], na.rm = TRUE),
      n_points = n(),
      .groups = 'drop'
    )
  
  # Rename the lga_value column to match the metric
  col_name <- paste0("lga_", gsub("_", "", value_column))
  names(lga_summary)[names(lga_summary) == "lga_value"] <- col_name
  
  # Handle LGAs without data using spatial interpolation
  all_lgas <- nigeria_lgas_sf %>%
    st_drop_geometry() %>%
    dplyr::select(.data$GID_2, .data$NAME_2, .data$NAME_1)
  
  missing_lgas <- all_lgas %>%
    dplyr::anti_join(lga_summary, by = "GID_2")
  
  cat("Found", nrow(missing_lgas), "LGAs without data. Filling using spatial interpolation...\n")
  
  if(nrow(missing_lgas) > 0) {
    for(i in seq_len(nrow(missing_lgas))) {
      lga_id <- missing_lgas$GID_2[i]
      
      # Get centroid of the missing LGA
      lga_centroid <- nigeria_lgas_sf %>%
        dplyr::filter(.data$GID_2 == lga_id) %>%
        st_centroid(of_largest_polygon = TRUE)
      
      # Find distances to all comparison points
      distances <- st_distance(lga_centroid, comparison_sf)
      distances_numeric <- as.numeric(distances)
      
      # Get indices of 5 closest points
      n_neighbors <- min(5, nrow(comparison_sf))
      closest_indices <- order(distances_numeric)[1:n_neighbors]
      
      # Calculate inverse distance weighted average
      closest_data <- comparison_data[closest_indices, ]
      weights <- 1 / (distances_numeric[closest_indices] + 1e-10)
      weights <- weights / sum(weights)
      
      interpolated_value <- sum(closest_data[[value_column]] * weights)
      
      # Add to summary
      col_name <- paste0("lga_", gsub("_", "", value_column))
      new_row <- data.frame(
        GID_2 = lga_id,
        NAME_2 = missing_lgas$NAME_2[i],
        NAME_1 = missing_lgas$NAME_1[i],
        stringsAsFactors = FALSE
      )
      new_row[[col_name]] <- interpolated_value
      new_row[["median_value"]] <- interpolated_value
      new_row[["min_value"]] <- interpolated_value
      new_row[["max_value"]] <- interpolated_value
      new_row[["n_points"]] <- 0
      
      lga_summary <- rbind(lga_summary, new_row)
    }
  }
  
  cat("LGA aggregation complete:", nrow(lga_summary), "LGAs\n")
  return(lga_summary)
}

# Aggregate different metrics by LGA
lga_ci_change_pct <- aggregate_by_lga(ci_comparison, "ci_width_change_pct", "CI width % change")
lga_ci_change_abs <- aggregate_by_lga(ci_comparison, "ci_width_change_abs", "CI width absolute change")
lga_relative_change <- aggregate_by_lga(ci_comparison, "relative_ci_change_pct", "Relative CI width % change")

# ==============================================================================
# CREATE CHOROPLETH MAPS
# ==============================================================================

# Function to create choropleth map
create_comparison_choropleth <- function(lga_data, value_column, title_text, 
                                       subtitle_text = "", 
                                       legend_title = "Change (%)",
                                       color_scale = "RdBu",
                                       add_inset = TRUE,
                                       center_scale = TRUE,
                                       use_log_scale = FALSE) {
  
  # Join LGA-level data back to LGA polygons
  nigeria_lgas_final <- nigeria_lgas_sf %>%
    left_join(lga_data, by = "GID_2")
  
  # Create log-scale transformation if requested
  if(use_log_scale) {
    nigeria_lgas_final <- nigeria_lgas_final %>%
      mutate(
        log_value = sign(.data[[value_column]]) * log10(abs(.data[[value_column]]) + 1)
      )
    plot_column <- "log_value"
    
    # Adjust legend title for log scale
    if(grepl("%", legend_title)) {
      legend_title <- paste("Log10(|", gsub("\\n", " ", legend_title), "| + 1)")
    }
  } else {
    plot_column <- value_column
  }
  
  # Determine color scale limits
  data_range <- range(nigeria_lgas_final[[plot_column]], na.rm = TRUE)
  
  if(center_scale) {
    # Center the scale around 0 for percentage changes
    max_abs <- max(abs(data_range), na.rm = TRUE)
    scale_limits <- c(-max_abs, max_abs)
    scale_breaks <- pretty(scale_limits, n = 5)
  } else {
    scale_limits <- data_range
    scale_breaks <- pretty(data_range, n = 5)
  }
  
  # Create choropleth map
  p <- ggplot(nigeria_lgas_final)
  
  # Add neighboring countries as background context
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Add Nigeria LGAs with data
  p <- p + 
    geom_sf(aes(fill = !!sym(plot_column)), lwd = 0.1, color = "white") +
    geom_sf(data = nigeria_states_sf, fill = NA, color = "black", lwd = 0.5)
  
  # Choose color scale
  if(color_scale == "RdBu") {
    p <- p + scale_fill_distiller(
      name = legend_title,
      type = "div",
      palette = "RdBu",
      direction = 1,
      limits = scale_limits,
      breaks = scale_breaks,
      na.value = "grey90",
      labels = if(use_log_scale) {
        function(x) round(x, 2)
      } else if(grepl("%", legend_title)) {
        function(x) paste0(round(x, 1), "%")
      } else {
        function(x) round(x, 3)
      }
    )
  } else if(color_scale == "viridis") {
    p <- p + scale_fill_viridis_c(
      name = legend_title,
      option = "viridis",
      direction = 1,
      limits = scale_limits,
      breaks = scale_breaks,
      na.value = "grey90",
      labels = if(use_log_scale) {
        function(x) round(x, 2)
      } else if(grepl("%", legend_title)) {
        function(x) paste0(round(x, 1), "%")
      } else {
        function(x) round(x, 3)
      }
    )
  }
  
  p <- p +
    # Add scale bar and north arrow
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                    bar_cols = c("black", "white"), line_width = 1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                          pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                          style = north_arrow_fancy_orienteering(text_size = 8)) +
    theme_void() +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      caption = "Data: INLA prevalence projections; Geometry: GADM"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0.5, size = 9, color = "grey60"),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    )
  
  # Add inset map if requested
  if(add_inset && !is.null(africa_countries_sf) && !is.null(nigeria_country_sf)) {
    inset_plot <- ggplot() +
      geom_sf(data = africa_countries_sf, fill = "grey90", color = "white", size = 0.1) +
      geom_sf(data = nigeria_country_sf, fill = "#440154", color = "white", size = 0.3) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 0.5),
        plot.margin = margin(2.1, 2.1, 2.1, 2.1)
      ) +
      coord_sf(expand = FALSE)
    
    # Combine with inset
    p <- ggdraw(p) +
      draw_plot(inset_plot, x = 0.70, y = 0.02, width = 0.25, height = 0.25)
  }
  
  return(p)
}

# Create comparison choropleth maps
p_ci_change_pct <- create_comparison_choropleth(
  lga_ci_change_pct, "lga_ciwidthchangepct",
  "Change in 95% CI Width: All-PCR vs Original (BCT+PCR)",
  "Percentage change in uncertainty at LGA level",
  "CI Width\nChange (%)",
  color_scale = "RdBu",
  center_scale = TRUE,
  use_log_scale = FALSE
)

p_ci_change_pct_log <- create_comparison_choropleth(
  lga_ci_change_pct, "lga_ciwidthchangepct",
  "Change in 95% CI Width: All-PCR vs Original (Log Scale)",
  "Log-scale percentage change in uncertainty at LGA level",
  "CI Width\nChange (%)",
  color_scale = "RdBu",
  center_scale = TRUE,
  use_log_scale = TRUE
)

p_ci_change_abs <- create_comparison_choropleth(
  lga_ci_change_abs, "lga_ciwidthchangeabs", 
  "Absolute Change in 95% CI Width: All-PCR vs Original",
  "Absolute change in prevalence uncertainty (percentage points)",
  "CI Width\nChange (pp)",
  color_scale = "RdBu",
  center_scale = TRUE,
  use_log_scale = FALSE
)

p_relative_change <- create_comparison_choropleth(
  lga_relative_change, "lga_relativecichangepct",
  "Change in Relative 95% CI Width: All-PCR vs Original", 
  "Percentage change in CI width relative to mean prevalence",
  "Relative CI\nChange (%)",
  color_scale = "RdBu",
  center_scale = TRUE,
  use_log_scale = FALSE
)

p_relative_change_log <- create_comparison_choropleth(
  lga_relative_change, "lga_relativecichangepct",
  "Change in Relative 95% CI Width: All-PCR vs Original (Log Scale)", 
  "Log-scale percentage change in CI width relative to mean prevalence",
  "Relative CI\nChange (%)",
  color_scale = "RdBu",
  center_scale = TRUE,
  use_log_scale = TRUE
)

# ==============================================================================
# CREATE SUMMARY STATISTICS AND DISTRIBUTIONS
# ==============================================================================

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Original Analysis (BCT+PCR):\n")
cat("  Mean CI width:", round(mean(ci_original$ci_width), 4), "\n")
cat("  Mean relative CI width:", round(mean(ci_original$relative_ci_width), 4), "\n")

cat("\nAll-PCR Analysis:\n")
cat("  Mean CI width:", round(mean(ci_all_pcr$ci_width), 4), "\n")
cat("  Mean relative CI width:", round(mean(ci_all_pcr$relative_ci_width), 4), "\n")

cat("\nComparison:\n")
cat("  Mean CI width change:", round(mean(ci_comparison$ci_width_change_pct), 2), "%\n")
cat("  Median CI width change:", round(median(ci_comparison$ci_width_change_pct), 2), "%\n")
cat("  CI width change range:", round(range(ci_comparison$ci_width_change_pct), 2), "%\n")

# Create distribution comparison plot
p_distributions <- ggplot(ci_comparison) +
  geom_histogram(aes(x = ci_width_change_pct, fill = "CI Width Change"), 
                 bins = 50, alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = mean(ci_comparison$ci_width_change_pct), 
             linetype = "solid", color = "darkred", size = 1) +
  scale_fill_viridis_d(option = "viridis", name = "") +
  labs(
    title = "Distribution of CI Width Changes",
    subtitle = "All-PCR vs Original (BCT+PCR) Analysis",
    x = "CI Width Change (%)",
    y = "Count",
    caption = "Red dashed line: 0% change, Red solid line: Mean change"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none"
  )

# Create log-scale distribution plot
p_log_distributions <- ggplot(ci_comparison) +
  geom_histogram(aes(x = log_abs_change, fill = "Log Change"), 
                 bins = 50, alpha = 0.7, color = "white") +
  geom_vline(xintercept = log_mean, linetype = "solid", color = "darkred", size = 1) +
  geom_vline(xintercept = outlier_threshold, linetype = "dashed", color = "red", size = 1) +
  scale_fill_viridis_d(option = "plasma", name = "") +
  labs(
    title = "Distribution of CI Width Changes (Log Scale)",
    subtitle = "Log10(|change%| + 1) - Better for extreme values",
    x = "Log10(|CI Width Change %| + 1)",
    y = "Count",
    caption = "Red solid line: Mean, Red dashed line: 2 SD threshold"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none"
  )

# Create log-scale scatter plot for CI width comparison
p_log_scatter <- ggplot(ci_comparison, aes(x = log10(original_ci_width), y = log10(all_pcr_ci_width))) +
  geom_point(aes(color = log_change_category), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(
    name = "Change Category",
    values = c("Small (≤10%)" = "#440154",
               "Moderate (10-20%)" = "#31688E",
               "High (20-50%)" = "#35B779",
               "Very High (50-100%)" = "#FDE725", 
               "Extreme (>100%)" = "#FF0000")
  ) +
  labs(
    title = "CI Width Comparison (Log Scale)",
    subtitle = "Log-scale reveals patterns masked by extreme values",
    x = "Log10(Original CI Width)",
    y = "Log10(All-PCR CI Width)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# Create simple CI width change plot (absolute difference)
p_ci_width_change <- ggplot(ci_comparison, aes(x = original_ci_width, y = ci_width_change_abs)) +
  geom_point(aes(color = ifelse(ci_width_change_abs > 0, "Increased", "Decreased")), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  scale_color_manual(
    name = "CI Width Change",
    values = c("Increased" = "#E31A1C", "Decreased" = "#1F78B4")
  ) +
  labs(
    title = "Change in 95% CI Width: All-PCR vs Original",
    subtitle = "Simple difference in credible interval width (no percentages or logs)",
    x = "Original CI Width",
    y = "CI Width Change (All-PCR - Original)",
    caption = "Points above red line: CI width increased; below: CI width decreased"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    plot.caption = element_text(hjust = 0.5, size = 9, color = "grey60"),
    legend.position = "bottom"
  )

# Create histogram for simple CI width change
p_ci_width_change_hist <- ggplot(ci_comparison, aes(x = ci_width_change_abs)) +
  geom_histogram(aes(fill = ifelse(ci_width_change_abs > 0, "Increased", "Decreased")), 
                 bins = 50, alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = mean(ci_comparison$ci_width_change_abs), 
             linetype = "solid", color = "darkred", size = 1) +
  scale_fill_manual(
    name = "CI Width Change",
    values = c("Increased" = "#1F78B4", "Decreased" = "#E31A1C")
  ) +
  labs(
    title = "Distribution of CI Width Changes",
    subtitle = "Simple difference in credible interval width (All-PCR - Original)",
    x = "CI Width Change (All-PCR - Original)",
    y = "Count",
    caption = "Red dashed line: no change (0), Red solid line: mean change"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.caption = element_text(hjust = 0.5, size = 9, color = "grey60"),
    legend.position = "bottom"
  )

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

# Create output directory
dir.create("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison", showWarnings = FALSE)

# Save choropleth maps
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_change_percentage.pdf", 
       plot = p_ci_change_pct, width = 12, height = 10, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_change_percentage_log.pdf", 
       plot = p_ci_change_pct_log, width = 12, height = 10, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_change_absolute.pdf", 
       plot = p_ci_change_abs, width = 12, height = 10, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/relative_ci_change.pdf", 
       plot = p_relative_change, width = 12, height = 10, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/relative_ci_change_log.pdf", 
       plot = p_relative_change_log, width = 12, height = 10, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_change_distribution.pdf", 
       plot = p_distributions, width = 10, height = 6, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_change_distribution_log.pdf", 
       plot = p_log_distributions, width = 10, height = 6, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_scatter_log.pdf", 
       plot = p_log_scatter, width = 10, height = 8, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_change_simple.pdf", 
       plot = p_ci_width_change, width = 10, height = 8, dpi = 300)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_change_histogram.pdf", 
       plot = p_ci_width_change_hist, width = 10, height = 6, dpi = 300)

# Save comparison data
write.csv(ci_comparison, "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/ci_width_comparison_data.csv", row.names = FALSE)
write.csv(lga_ci_change_pct, "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/lga_ci_change_percentage.csv", row.names = FALSE)

# Display plots
print(p_ci_change_pct)
print(p_ci_change_pct_log)
print(p_ci_change_abs) 
print(p_relative_change)
print(p_relative_change_log)
print(p_distributions)
print(p_log_distributions)
print(p_log_scatter)
print(p_ci_width_change)
print(p_ci_width_change_hist)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("CI width comparison analysis completed successfully!\n")
cat("Files saved to: Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/CI_comparison/\n")
cat("- ci_width_change_percentage.pdf - Percentage change choropleth (linear scale)\n")
cat("- ci_width_change_percentage_log.pdf - Percentage change choropleth (log scale)\n")
cat("- ci_width_change_absolute.pdf - Absolute change choropleth\n") 
cat("- relative_ci_change.pdf - Relative CI change choropleth (linear scale)\n")
cat("- relative_ci_change_log.pdf - Relative CI change choropleth (log scale)\n")
cat("- ci_change_distribution.pdf - Distribution of changes (linear scale)\n")
cat("- ci_change_distribution_log.pdf - Distribution of changes (log scale)\n")
cat("- ci_width_scatter_log.pdf - Log-scale scatter plot\n")
cat("- ci_width_change_simple.pdf - Simple CI width change plot (no percentages/logs)\n")
cat("- ci_width_change_histogram.pdf - Histogram of CI width changes (no percentages/logs)\n")
cat("- ci_width_comparison_data.csv - Point-level comparison data\n")
cat("- lga_ci_change_percentage.csv - LGA-level summary data\n\n")
cat("Log-scale analysis provides:\n")
cat("1. Better handling of extreme percentage changes\n")
cat("2. Statistical outlier detection based on standard deviations\n")
cat("3. Clearer visualization of patterns masked by extreme values\n")
cat("4. More robust statistical measures for skewed distributions\n")