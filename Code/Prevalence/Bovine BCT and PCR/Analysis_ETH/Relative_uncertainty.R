rm(list=ls())

library(ggplot2)
library(see)
library(zoom)
library(ggpattern)
library(ggnewscale)
library(concaveman)
library(dplyr)
library(terra)
library(sf)
library(readxl)
library(geodata)
library(raster)
library(viridis)
library(gridExtra)
library(grid)
library(ggspatial)



# Determine number of units from Ethiopian projection files
n_datasets <- 1000
sample_file <- "Code/Prevalence/Bovine BCT and PCR/Projections_ETH/Projections_model_1.csv"
if (file.exists(sample_file)) {
  sample_data <- read.csv(sample_file)
  n_units <- nrow(sample_data[sample_data$variable == "Mean", ])
  cat("Ethiopian projection files have", n_units, "units per model\n")
} else {
  stop("Cannot find Ethiopian projection files to determine n_units")
}

check_symmetry <- function(dpm_mean, dpm_lower, dpm_upper, tol = 0.1) {
  # dpm_* are numeric vectors of same length
  # tol = allowable asymmetry ratio (10% by default)

  # midpoint of CI
  midpoint <- (dpm_lower + dpm_upper) / 2
  
  # deviation of mean from midpoint
  dev <- abs(dpm_mean - midpoint)
  
  # half-width of CI
  halfwidth <- (dpm_upper - dpm_lower) / 2
  
  # relative asymmetry (0 = perfectly symmetric, 0.5 = mean shifted halfway to one side)
  rel_asym <- dev / halfwidth
  
  # flag units where asymmetry exceeds tolerance
  flagged <- rel_asym > tol
  
  return(list(rel_asym = rel_asym,
              flagged = flagged,
              frac_flagged = mean(flagged, na.rm=TRUE)))
}

logit <- function(p, eps = 1e-6) {
  p <- pmax(pmin(p, 1 - eps), eps)  # avoid Inf
  log(p / (1 - p))
}

# Containers to accumulate
means_matrix <- matrix(NA, nrow=n_units, ncol=n_datasets)
within_var_accum <- matrix(NA, nrow=n_units, ncol=n_datasets)

for (i in 1:n_datasets){
dpm <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Projections_ETH/Projections_model_",i,".csv"))


#extract all data where 4th column is "Mean"
dpm_mean <- logit(dpm[dpm[,4]=="Mean",]$value)
dpm_lower <- logit(dpm[dpm[,4]=="2.5th percentile",]$value)
dpm_upper <- logit(dpm[dpm[,4]=="97.5th percentile",]$value)

means_matrix[, i] <- dpm_mean
within_sd <- (dpm_upper - dpm_lower) / (2*1.96)
within_var_accum[, i] <- within_sd^2

}

# After loop: compute averages and variances across datasets
between_var <- apply(means_matrix, 1, var, na.rm=TRUE)
within_var  <- rowMeans(within_var_accum, na.rm=TRUE)

# Relative importance
total_var <- between_var + within_var
rel_within  <- ifelse(total_var > 0, within_var / total_var, NA)
rel_between <- ifelse(total_var > 0, between_var / total_var, NA)

# Debug: Check point-level variation
cat("Point-level uncertainty ranges:\n")
cat("  rel_between range:", range(rel_between, na.rm=TRUE), "\n")
cat("  rel_within range:", range(rel_within, na.rm=TRUE), "\n")
cat("  between_var range:", range(between_var, na.rm=TRUE), "\n")
cat("  Non-NA rel_between points:", sum(!is.na(rel_between)), "\n")

#make a new dataframe with Latitude, Longitude, between_var, within_var, rel_within, rel_between
dpm_var <- data.frame(Latitude = dpm$Latitude[dpm$variable=="Mean"],
                      Longitude = dpm$Longitude[dpm$variable=="Mean"],
                      between_var = between_var,
                      within_var = within_var,
                      rel_within = rel_within,
                      rel_between = rel_between)

# Get Ethiopian administrative boundaries for choropleth mapping
ethiopia_regions <- gadm(country = "ETH", level = 1, path = tempdir())
ethiopia_regions_sf <- st_as_sf(ethiopia_regions)
ethiopia_zones <- gadm(country = "ETH", level = 3, path = tempdir()) 
ethiopia_zones_sf <- st_as_sf(ethiopia_zones)

# Load neighboring countries for geographic context
cat("Loading neighboring countries for geographic context...\n")
neighbor_countries <- c("Sudan" = "SDN", "South Sudan" = "SSD", "Kenya" = "KEN", 
                       "Somalia" = "SOM", "Djibouti" = "DJI", "Eritrea" = "ERI", "Uganda" = "UGA")

neighboring_countries_list <- list()
for(i in 1:length(neighbor_countries)) {
  country_name <- names(neighbor_countries)[i]
  country_code <- neighbor_countries[i]
  
  tryCatch({
    country_data <- gadm(country = country_code, level = 0, path = tempdir())
    country_sf <- st_as_sf(country_data)
    neighboring_countries_list[[country_name]] <- country_sf
    cat("Loaded", country_name, "(", country_code, ")\n")
  }, error = function(e) {
    cat("Failed to load", country_name, "(", country_code, "):", e$message, "\n")
  })
}

# Combine all neighboring countries into single sf object
if(length(neighboring_countries_list) > 0) {
  neighboring_countries_full <- do.call(rbind, neighboring_countries_list)
  
  # Validate and repair geometries to avoid intersection errors
  cat("Validating and repairing geometries...\n")
  neighboring_countries_full <- st_make_valid(neighboring_countries_full)
  ethiopia_regions_valid <- st_make_valid(ethiopia_regions_sf)
  
  # Create a buffer around Ethiopia to crop neighboring countries to border regions only
  ethiopia_buffer <- st_buffer(st_union(ethiopia_regions_valid), dist = 1.5)  # ~150km buffer
  ethiopia_buffer <- st_make_valid(ethiopia_buffer)
  
  # Use st_crop instead of st_intersection for more robust clipping
  tryCatch({
    neighboring_countries_sf <- st_crop(neighboring_countries_full, ethiopia_buffer)
    cat("Successfully loaded and cropped", nrow(neighboring_countries_sf), "neighboring country border regions\n")
  }, error = function(e) {
    cat("st_crop failed, trying alternative approach:", e$message, "\n")
    # Fallback: use a simpler bbox-based crop
    ethiopia_bbox <- st_bbox(ethiopia_buffer)
    neighboring_countries_sf <<- st_crop(neighboring_countries_full, ethiopia_bbox)
    cat("Successfully loaded and cropped using bbox approach\n")
  })
  
} else {
  neighboring_countries_sf <- NULL
  cat("No neighboring countries loaded\n")
}

# Load bovine data from Excel files
bovine_bct_raw <- read_excel("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
bovine_pcr_raw <- read_excel("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")

# Add test type identifier to each dataset
bovine_bct_raw$Test_Type <- "BCT/HCT"
bovine_pcr_raw$Test_Type <- "PCR"

# Standardize column names - ensure both datasets have the same column structure
# First, get the column names of both datasets
bct_cols <- names(bovine_bct_raw)
pcr_cols <- names(bovine_pcr_raw)

cat("BCT columns:", paste(bct_cols, collapse = ", "), "\n")
cat("PCR columns:", paste(pcr_cols, collapse = ", "), "\n")

# Rename prevalence columns to standardized names
if("TPR" %in% names(bovine_bct_raw)) {
  bovine_bct_raw$Prevalence_Rate <- bovine_bct_raw$TPR
  bovine_bct_raw$TPR <- NULL
}
if("T_ATPR" %in% names(bovine_pcr_raw)) {
  bovine_pcr_raw$Prevalence_Rate <- bovine_pcr_raw$T_ATPR
  bovine_pcr_raw$T_ATPR <- NULL
}

# Find common columns between the two datasets
common_cols <- intersect(names(bovine_bct_raw), names(bovine_pcr_raw))
cat("Common columns:", paste(common_cols, collapse = ", "), "\n")

# Select only common columns from both datasets before combining
bovine_bct_clean <- bovine_bct_raw[, common_cols, drop = FALSE]
bovine_pcr_clean <- bovine_pcr_raw[, common_cols, drop = FALSE]

# Now combine the datasets with matching column structures
bovine_data_raw <- rbind(bovine_bct_clean, bovine_pcr_clean)

# Process bovine data for mapping
process_bovine_data <- function(data) {
  # Clean and filter the data
  bovine_clean <- data %>%
    # Remove rows with missing coordinates
    dplyr::filter(!is.na(Longitude), !is.na(Latitude)) %>%
    # Remove rows with missing or invalid sample sizes
    dplyr::filter(!is.na(Number_of_animal_tested), Number_of_animal_tested > 0) %>%
    # Calculate prevalence if not already calculated
    dplyr::mutate(
      Number_of_infections = ifelse(is.na(Number_of_infections), 0, Number_of_infections),
      Prevalence = ifelse(Number_of_animal_tested > 0, Number_of_infections / Number_of_animal_tested, 0)
    )
  
  # Convert to spatial points (no coordinate swap needed for this dataset)
  bovine_sf <- st_as_sf(bovine_clean, 
                        coords = c("Longitude", "Latitude"), 
                        crs = 4326)
  
  # Filter for points within Ethiopia using spatial intersection
  bovine_ethiopia <- st_filter(bovine_sf, ethiopia_regions_sf)
  
  # Add back coordinate columns for plotting
  coords <- st_coordinates(bovine_ethiopia)
  bovine_ethiopia$lon <- coords[, 1]
  bovine_ethiopia$lat <- coords[, 2]
  
  return(bovine_ethiopia)
}

# Process the bovine data
bovine_ethiopia_sf <- process_bovine_data(bovine_data_raw)
cat("Found", nrow(bovine_ethiopia_sf), "bovine data points within Ethiopia\n")
cat("  - BCT tests:", sum(bovine_ethiopia_sf$Test_Type == "BCT/HCT"), "\n")
cat("  - PCR tests:", sum(bovine_ethiopia_sf$Test_Type == "PCR"), "\n")

# Function to process uncertainty data (create spatial points without coordinate swap)
process_uncertainty_data <- function(uncertainty_data, value_col) {
  # Create a copy and select relevant columns
  proj_data <- uncertainty_data[, c("Longitude", "Latitude", value_col)]
  names(proj_data)[3] <- "value"  # Rename the value column for consistency

  # Swap longitude and latitude
  temp <- proj_data$Longitude
  proj_data$Longitude <- proj_data$Latitude
  proj_data$Latitude <- temp
  
  # Remove any rows with missing coordinates or values
  proj_data <- proj_data[!is.na(proj_data$Longitude) & !is.na(proj_data$Latitude) & !is.na(proj_data$value), ]
  
  # Convert to spatial points
  proj_sf <- st_as_sf(proj_data, 
                      coords = c("Longitude", "Latitude"), 
                      crs = 4326)
  
  return(proj_sf)
}

# Function to aggregate uncertainty data by region using spatial join with gap filling
aggregate_uncertainty_by_zone <- function(uncertainty_sf) {
  joined_data <- st_join(uncertainty_sf, ethiopia_zones_sf)
  
  # First, get zones with direct data points
  zone_uncertainty <- joined_data %>%
    st_drop_geometry() %>%
    dplyr::filter(!is.na(GID_3)) %>%  # Only keep points that fell within zones
    dplyr::group_by(GID_3, NAME_2) %>%
    dplyr::summarise(
      zone_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      n_points = n(),
      .groups = 'drop'
    )
  
  # Find zones without data (grey areas)
  all_zones <- ethiopia_zones_sf %>%
    st_drop_geometry() %>%
    dplyr::select(GID_3, NAME_2)
  
  missing_zones <- all_zones %>%
    dplyr::anti_join(zone_uncertainty, by = "GID_3")
  
  cat("Found", nrow(missing_zones), "zones without uncertainty data. Filling using spatial interpolation...\n")
  
  if(nrow(missing_zones) > 0) {
    # For each missing zone, find closest uncertainty points and interpolate
    for(i in seq_len(nrow(missing_zones))) {
      zone_id <- missing_zones$GID_3[i]
      
      # Get centroid of the missing zone (suppress attribute warning)
      zone_centroid <- ethiopia_zones_sf %>%
        dplyr::filter(GID_3 == zone_id) %>%
        st_centroid(of_largest_polygon = TRUE)
      
      # Find distances to all uncertainty points
      distances <- st_distance(zone_centroid, uncertainty_sf)
      
      # Convert distances to numeric (remove units)
      distances_numeric <- as.numeric(distances)
      
      # Get indices of 5 closest points (or fewer if not enough points available)
      n_neighbors <- min(5, nrow(uncertainty_sf))
      closest_indices <- order(distances_numeric)[1:n_neighbors]
      
      # Calculate inverse distance weighted average
      closest_points <- uncertainty_sf[closest_indices, ]
      weights <- 1 / (distances_numeric[closest_indices] + 1e-10)  # Add small value to avoid division by zero
      weights <- weights / sum(weights)  # Normalize weights
      
      # Calculate weighted average uncertainty
      interpolated_value <- sum(closest_points$value * weights)
      
      # Add to zone_uncertainty data
      new_row <- data.frame(
        GID_3 = zone_id,
        NAME_2 = missing_zones$NAME_2[i],
        zone_value = interpolated_value,
        median_value = interpolated_value,
        min_value = interpolated_value,
        max_value = interpolated_value,
        n_points = 0  # Indicate this is interpolated
      )
      
      zone_uncertainty <- rbind(zone_uncertainty, new_row)
    }
  }
  
  return(zone_uncertainty)
}

# Process uncertainty data for choropleth mapping
rel_between_sf <- process_uncertainty_data(dpm_var, "rel_between")
rel_within_sf <- process_uncertainty_data(dpm_var, "rel_within")
between_var_sf <- process_uncertainty_data(dpm_var, "between_var")

# Aggregate by zone with spatial interpolation
zone_rel_between <- aggregate_uncertainty_by_zone(rel_between_sf)
zone_rel_within <- aggregate_uncertainty_by_zone(rel_within_sf)
zone_between_var <- aggregate_uncertainty_by_zone(between_var_sf)

# Debug: Check zonal aggregation results
cat("\nZonal aggregation results:\n")
cat("  rel_between zonal range:", range(zone_rel_between$zone_value, na.rm=TRUE), "\n")
cat("  rel_within zonal range:", range(zone_rel_within$zone_value, na.rm=TRUE), "\n")
cat("  between_var zonal range:", range(zone_between_var$zone_value, na.rm=TRUE), "\n")
cat("  Number of zones with data:\n")
cat("    rel_between:", nrow(zone_rel_between), "\n")
cat("    rel_within:", nrow(zone_rel_within), "\n")
cat("    between_var:", nrow(zone_between_var), "\n")

# Function to create choropleth map for uncertainty measures
create_uncertainty_choropleth <- function(zone_data, title_suffix, legend_title, color_limits = c(0, 1)) {
  # Join zone-level uncertainty back to zone polygons
  ethiopia_zones_final <- ethiopia_zones_sf %>%
    left_join(zone_data, by = "GID_3")
  
  # Debug: Check what we're about to plot
  cat("Choropleth data for", title_suffix, ":\n")
  cat("  Data range:", range(ethiopia_zones_final$zone_value, na.rm=TRUE), "\n")
  cat("  Non-NA zones:", sum(!is.na(ethiopia_zones_final$zone_value)), "\n")
  cat("  Color limits:", color_limits, "\n")
  
  # Create zone-level choropleth map
  p <- ggplot(ethiopia_zones_final)
  
  # Add neighboring countries as background context (if available)
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Add Ethiopia zones with uncertainty data
  p <- p + 
    geom_sf(aes(fill = zone_value), lwd = 0.05, color = "black") +
    geom_sf(data = ethiopia_regions_sf, fill = NA, color = "black", lwd = 0.5) +  # Add region boundaries
    # Add bovine data points (different colors for BCT and PCR)
    geom_point(data = st_drop_geometry(bovine_ethiopia_sf), 
               aes(x = lon, y = lat, size = Number_of_animal_tested, shape = Test_Type), colour="black", fill="white", stroke=0.8,
               alpha = 0.8) +  # Increased from 0.7 to 0.8 for better visibility
    scale_fill_viridis_c(
      name = legend_title,
      na.value = "grey90",
      direction = 1,
      option = "magma",
      limits = color_limits,
      labels = if(max(color_limits) <= 1) scales::percent_format(accuracy = 0.1) else scales::number_format(accuracy = 0.001)
    ) +
    # Add size legend for sample size
    scale_size_continuous(name = "Sample size", 
                         range = c(1, 10),  # Increased from c(1, 4) to make sizes more distinguishable
                         breaks = c(10, 50, 100, 200),
                         labels = c("10", "50", "100", "200+")) +
    # Add color scale for test types
    scale_color_manual(name = "Test type",
                      values = c("BCT/HCT" = "red", "PCR" = "blue"),
                      guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
                      scale_shape_manual(name = "Test type",labels = c(
    "BCT/HCT" = "BCT",
    "PCR" = "PCR"
  ),values = c("BCT/HCT" = 21,"PCR" = 22)
)+
    # Add scale bar and north arrow
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                    bar_cols = c("black", "white"), line_width = 1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                          pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                          style = north_arrow_fancy_orienteering(text_size = 8)) +
    theme_void() +
    labs(
      title = paste("Uncertainty analysis by district in Ethiopia", title_suffix)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.spacing.y = unit(1, "cm"),  # Add vertical spacing between legends
      legend.box.spacing = unit(1, "cm")  # Add spacing between legend boxes
    )
  
  return(p)
}

# Function to create combined histogram and boxplot for uncertainty measures
create_uncertainty_combined_plot <- function(zone_data, estimate_name, x_label) {
  # Histogram with viridis color scale matching choropleth
  p_hist <- ggplot(zone_data, aes(x = zone_value, fill = after_stat(x))) +
    geom_histogram(binwidth = 0.04, color = "white", linewidth = 0.2,
                   boundary = 0, closed = "left") +
    scale_fill_viridis_c(option = "magma", guide = "none") +  # Same as choropleth, no legend
    scale_x_continuous(labels = if(max(zone_data$zone_value, na.rm = TRUE) <= 1) scales::percent_format() else scales::number_format(accuracy = 0.001),
                       breaks = seq(0, 1, by = 0.2),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title.y = element_blank(),  # Remove y-axis title
      axis.text.y = element_blank(),   # Remove y-axis tick labels
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),  # Remove all grid lines
      plot.margin = margin(t = 5, r = 15, b = 0, l = 15),  # More left/right margins
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(
      title = paste("Distribution of", estimate_name, "by district in Ethiopia"),
      x = ""
    )
  
  # Boxplot with viridis color matching choropleth
  viridis_color <- viridis::viridis(1, begin = 0.5, end = 0.5, option = "magma")  # Mid-range viridis color
  
  p_box <- ggplot(zone_data, aes(x = zone_value, y = 1)) +
    geom_boxplot(fill = viridis_color, alpha = 0.8, width = 0.8) +  # Much wider boxplot, no jitter
    scale_x_continuous(labels = if(max(zone_data$zone_value, na.rm = TRUE) <= 1) scales::percent_format() else scales::number_format(accuracy = 0.001),
                       breaks = seq(0, 1, by = 0.2),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0.2, 1.8), expand = c(0, 0)) +  # More vertical space
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 5)),
      axis.text.x = element_text(size = 9),
      axis.ticks.x = element_line(),
      axis.title.y = element_blank(),  # Remove y-axis title
      axis.text.y = element_blank(),   # Remove y-axis text
      axis.ticks.y = element_blank(),  # Remove y-axis ticks
      panel.grid.major = element_blank(),  # Remove all grid lines
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 15, b = 5, l = 15)  # Same margins as histogram for alignment
    ) +
    labs(
      x = x_label
    )
  
  # Create title
  title <- grid::textGrob(paste("Uncertainty Analysis -", estimate_name), 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  # Combine plots with custom heights
  combined_plot <- grid.arrange(
    title,
    p_hist,
    p_box,
    nrow = 3,
    heights = unit(c(0.8, 3.5, 2.5), "cm")  # Small title, histogram, much taller boxplot
  )
  
  return(combined_plot)
}

# Create choropleth maps
p_rel_between <- create_uncertainty_choropleth(zone_rel_between, "\n (relative variability from diagnostic uncertainty)", "Relative\nimportance", c(0, 1))
p_rel_within <- create_uncertainty_choropleth(zone_rel_within, "(relative variability from coverage uncertainty)", "Relative\nimportance", c(0, 1))
p_between_var <- create_uncertainty_choropleth(zone_between_var, "(Variance from diagnostic uncertainty)", "Variance", c(min(zone_between_var$zone_value, na.rm = TRUE), max(zone_between_var$zone_value, na.rm = TRUE)))

# Create combined histogram and boxplot analyses
p_rel_between_combined <- create_uncertainty_combined_plot(zone_rel_between, "relative variability from diagnostic uncertainty", "Relative importance of variability from diagnostic uncertainty")
p_rel_within_combined <- create_uncertainty_combined_plot(zone_rel_within, "relative variability from coverage uncertainty", "Relative importance of variability from coverage uncertainty")
p_between_var_combined <- create_uncertainty_combined_plot(zone_between_var, "variance from diagnostic uncertainty", "Variance from diagnostic uncertainty")

# Save all choropleth maps
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/region_rel_between_choropleth.pdf", plot = p_rel_between, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/region_rel_within_choropleth.pdf", plot = p_rel_within, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/region_between_var_choropleth.pdf", plot = p_between_var, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each uncertainty measure)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/region_rel_between_analysis.pdf", plot = p_rel_between_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/region_rel_within_analysis.pdf", plot = p_rel_within_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/region_between_var_analysis.pdf", plot = p_between_var_combined, width = 10, height = 6)

cat("=== UNCERTAINTY ESTIMATES SUMMARY ===\n\n")

cat("Relative Between-Dataset Variability:\n")
cat("Total zones with data:", nrow(zone_rel_between), "\n")
if(nrow(zone_rel_between) > 0) {
  print(summary(zone_rel_between$zone_value))
}

cat("\nRelative Within-Dataset Variability:\n")
cat("Total zones with data:", nrow(zone_rel_within), "\n")
if(nrow(zone_rel_within) > 0) {
  print(summary(zone_rel_within$zone_value))
}

cat("\nBetween-Dataset Variance:\n")
cat("Total zones with data:", nrow(zone_between_var), "\n")
if(nrow(zone_between_var) > 0) {
  print(summary(zone_between_var$zone_value))
}



#plot rel_between
ggplot(dpm_var) + geom_tile(aes(Longitude, Latitude, fill = rel_between)) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Relative importance of\nbetween-dataset variability",
    low = "blue", high = "orange",
    limits = c(0, 1),
    guide = guide_colorbar(order = 2)
  ) +
  new_scale_fill() +  # Add a new fill scale
  #geom_tile(data=water_cover, aes(Longitude, Latitude,fill=water_cover,alpha=water_cover),inherit.aes = FALSE)+
  scale_fill_manual(
  name = "",  # Colorbar title
  values = c("High uncertainty" = "grey50", "Water" = "steelblue1"),  # Specific colors for categories
  limits = c("High uncertainty", "Water"),
  guide = guide_legend(order = 1)
)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Relative_uncertainty_between.png",width = 12, height = 8)


#plot histogram of dpm_var$rel_between
ggplot(dpm_var, aes(x = rel_between)) +
  geom_histogram(binwidth = 0.05, fill = "black", color = "black", alpha = 0.7) +
  labs(
       x = "Relative importance of between-dataset variability",
       y = "Frequency")
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Relative_uncertainty_between_histogram.png",width = 8, height = 6)








# ============================================================
# FIGURE 6B SOURCE DATA
#
# Ethiopian uncertainty decomposition maps and associated
# histogram/boxplot panels.
#
# Sheets:
#   README
#   District_values
#   Distribution_values
#   Grid_values
#   Survey_points
#   Plot_settings
#   Summary_statistics
#
# ============================================================

# Explicit namespaces are used because raster may mask
# functions such as select().

# ------------------------------------------------------------
# 1. Output location
# ------------------------------------------------------------

figure6b_output_dir <- paste0(
  "Code/Prevalence/Bovine BCT and PCR/",
  "Analysis_ETH"
)

dir.create(
  figure6b_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Unique Ethiopian level-3 polygon lookup
# ------------------------------------------------------------

figure6b_district_lookup <- ethiopia_zones_sf |>
  sf::st_drop_geometry() |>
  dplyr::transmute(
    GID_3 = as.character(GID_3),
    Region = as.character(NAME_1),
    Zone = as.character(NAME_2),
    District = as.character(NAME_3)
  ) |>
  dplyr::distinct(
    GID_3,
    .keep_all = TRUE
  )

stopifnot(
  nrow(figure6b_district_lookup) ==
    dplyr::n_distinct(
      figure6b_district_lookup$GID_3
    )
)

# ------------------------------------------------------------
# 3. Standardise uncertainty-result objects
# ------------------------------------------------------------

figure6b_rel_between <- zone_rel_between |>
  dplyr::as_tibble() |>
  dplyr::transmute(
    GID_3 = as.character(GID_3),

    Zone_from_result =
      as.character(NAME_2),

    Relative_diagnostic_uncertainty =
      as.numeric(zone_value),

    Median_grid_value =
      as.numeric(median_value),

    Minimum_grid_value =
      as.numeric(min_value),

    Maximum_grid_value =
      as.numeric(max_value),

    Number_of_grid_points =
      as.integer(n_points),

    Value_source = dplyr::if_else(
      Number_of_grid_points > 0L,
      "Direct grid-point aggregation",
      paste(
        "Five-nearest-point",
        "inverse-distance interpolation"
      )
    )
  )

figure6b_rel_within <- zone_rel_within |>
  dplyr::as_tibble() |>
  dplyr::transmute(
    GID_3 = as.character(GID_3),

    Relative_coverage_uncertainty =
      as.numeric(zone_value),

    Coverage_median_grid_value =
      as.numeric(median_value),

    Coverage_minimum_grid_value =
      as.numeric(min_value),

    Coverage_maximum_grid_value =
      as.numeric(max_value),

    Coverage_number_of_grid_points =
      as.integer(n_points),

    Coverage_value_source = dplyr::if_else(
      Coverage_number_of_grid_points > 0L,
      "Direct grid-point aggregation",
      paste(
        "Five-nearest-point",
        "inverse-distance interpolation"
      )
    )
  )

figure6b_between_var <- zone_between_var |>
  dplyr::as_tibble() |>
  dplyr::transmute(
    GID_3 = as.character(GID_3),

    Diagnostic_variance =
      as.numeric(zone_value),

    Diagnostic_variance_median =
      as.numeric(median_value),

    Diagnostic_variance_minimum =
      as.numeric(min_value),

    Diagnostic_variance_maximum =
      as.numeric(max_value),

    Variance_number_of_grid_points =
      as.integer(n_points),

    Variance_value_source = dplyr::if_else(
      Variance_number_of_grid_points > 0L,
      "Direct grid-point aggregation",
      paste(
        "Five-nearest-point",
        "inverse-distance interpolation"
      )
    )
  )

stopifnot(
  nrow(figure6b_rel_between) ==
    dplyr::n_distinct(
      figure6b_rel_between$GID_3
    ),

  nrow(figure6b_rel_within) ==
    dplyr::n_distinct(
      figure6b_rel_within$GID_3
    ),

  nrow(figure6b_between_var) ==
    dplyr::n_distinct(
      figure6b_between_var$GID_3
    )
)

# ------------------------------------------------------------
# 4. Combined district-level table
# ------------------------------------------------------------

figure6b_district_values <- figure6b_district_lookup |>
  dplyr::left_join(
    figure6b_rel_between,
    by = "GID_3",
    relationship = "one-to-one"
  ) |>
  dplyr::left_join(
    figure6b_rel_within,
    by = "GID_3",
    relationship = "one-to-one"
  ) |>
  dplyr::left_join(
    figure6b_between_var,
    by = "GID_3",
    relationship = "one-to-one"
  ) |>
  dplyr::mutate(
    Relative_components_sum =
      Relative_diagnostic_uncertainty +
      Relative_coverage_uncertainty,

    Relative_components_difference_from_one =
      Relative_components_sum - 1,

    Diagnostic_value_interpolated =
      Number_of_grid_points == 0L,

    Coverage_value_interpolated =
      Coverage_number_of_grid_points == 0L,

    Variance_value_interpolated =
      Variance_number_of_grid_points == 0L
  ) |>
  dplyr::arrange(
    Region,
    Zone,
    District,
    GID_3
  )

stopifnot(
  nrow(figure6b_district_values) ==
    dplyr::n_distinct(
      figure6b_district_values$GID_3
    )
)

# ------------------------------------------------------------
# 5. Long-format plotted district values
# ------------------------------------------------------------

figure6b_distribution_values <- dplyr::bind_rows(

  figure6b_district_values |>
    dplyr::transmute(
      GID_3,
      Region,
      Zone,
      District,

      Uncertainty_measure =
        "Relative diagnostic uncertainty",

      Plot_object =
        "zone_rel_between",

      Source_column =
        "zone_value",

      Plotted_value =
        Relative_diagnostic_uncertainty,

      Number_of_grid_points =
        Number_of_grid_points,

      Interpolated =
        Diagnostic_value_interpolated
    ),

  figure6b_district_values |>
    dplyr::transmute(
      GID_3,
      Region,
      Zone,
      District,

      Uncertainty_measure =
        "Relative coverage uncertainty",

      Plot_object =
        "zone_rel_within",

      Source_column =
        "zone_value",

      Plotted_value =
        Relative_coverage_uncertainty,

      Number_of_grid_points =
        Coverage_number_of_grid_points,

      Interpolated =
        Coverage_value_interpolated
    ),

  figure6b_district_values |>
    dplyr::transmute(
      GID_3,
      Region,
      Zone,
      District,

      Uncertainty_measure =
        "Diagnostic variance",

      Plot_object =
        "zone_between_var",

      Source_column =
        "zone_value",

      Plotted_value =
        Diagnostic_variance,

      Number_of_grid_points =
        Variance_number_of_grid_points,

      Interpolated =
        Variance_value_interpolated
    )

) |>
  dplyr::filter(
    !is.na(Plotted_value)
  ) |>
  dplyr::arrange(
    Uncertainty_measure,
    Region,
    Zone,
    District
  )

# ------------------------------------------------------------
# 6. Prediction-grid uncertainty values
# ------------------------------------------------------------

# process_uncertainty_data() swaps the original Longitude and
# Latitude columns before creating the plotting coordinates.
# Both original and corrected coordinate interpretations are
# retained here.

figure6b_grid_values <- dpm_var |>
  dplyr::as_tibble() |>
  dplyr::transmute(
    Grid_cell_ID =
      dplyr::row_number(),

    Original_Latitude_column =
      as.numeric(Latitude),

    Original_Longitude_column =
      as.numeric(Longitude),

    Plot_longitude =
      as.numeric(Latitude),

    Plot_latitude =
      as.numeric(Longitude),

    Diagnostic_variance =
      as.numeric(between_var),

    Coverage_variance =
      as.numeric(within_var),

    Total_variance =
      Diagnostic_variance +
      Coverage_variance,

    Relative_coverage_uncertainty =
      as.numeric(rel_within),

    Relative_diagnostic_uncertainty =
      as.numeric(rel_between),

    Relative_components_sum =
      Relative_coverage_uncertainty +
      Relative_diagnostic_uncertainty
  )

# ------------------------------------------------------------
# 7. Survey-point overlays
# ------------------------------------------------------------

if (
  exists("bovine_ethiopia_sf") &&
  !is.null(bovine_ethiopia_sf) &&
  nrow(bovine_ethiopia_sf) > 0
) {

  figure6b_survey_points <- bovine_ethiopia_sf |>
    sf::st_drop_geometry() |>
    dplyr::transmute(
      Survey_point_ID =
        dplyr::row_number(),

      Longitude =
        as.numeric(lon),

      Latitude =
        as.numeric(lat),

      Sample_size =
        as.integer(Number_of_animal_tested),

      Number_infected =
        as.integer(Number_of_infections),

      Observed_prevalence =
        as.numeric(Prevalence),

      Test_type = dplyr::case_when(
        Test_Type == "BCT/HCT" ~ "HCT/BCT",
        Test_Type == "BCT" ~ "HCT/BCT",
        Test_Type == "HCT/BCT" ~ "HCT/BCT",
        Test_Type == "PCR" ~ "PCR",
        TRUE ~ as.character(Test_Type)
      )
    ) |>
    dplyr::arrange(
      Test_type,
      Longitude,
      Latitude
    )

} else {

  figure6b_survey_points <- data.frame(
    Survey_point_ID = integer(0),
    Longitude = numeric(0),
    Latitude = numeric(0),
    Sample_size = integer(0),
    Number_infected = integer(0),
    Observed_prevalence = numeric(0),
    Test_type = character(0)
  )
}

# ------------------------------------------------------------
# 8. Summary-statistics helper
# ------------------------------------------------------------

figure6b_make_summary <- function(
    values,
    measure_name
) {

  values <- values[
    is.finite(values)
  ]

  data.frame(
    Uncertainty_measure =
      measure_name,

    Number_of_districts =
      length(values),

    Mean =
      mean(
        values,
        na.rm = TRUE
      ),

    Standard_deviation =
      stats::sd(
        values,
        na.rm = TRUE
      ),

    Minimum =
      min(
        values,
        na.rm = TRUE
      ),

    Lower_2.5_percentile =
      as.numeric(
        stats::quantile(
          values,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

    Median =
      stats::median(
        values,
        na.rm = TRUE
      ),

    Upper_97.5_percentile =
      as.numeric(
        stats::quantile(
          values,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

    Maximum =
      max(
        values,
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE
  )
}

figure6b_summary_statistics <- dplyr::bind_rows(

  figure6b_make_summary(
    figure6b_district_values$
      Relative_diagnostic_uncertainty,
    "Relative diagnostic uncertainty"
  ),

  figure6b_make_summary(
    figure6b_district_values$
      Relative_coverage_uncertainty,
    "Relative coverage uncertainty"
  ),

  figure6b_make_summary(
    figure6b_district_values$
      Diagnostic_variance,
    "Diagnostic variance"
  )
)

# ------------------------------------------------------------
# 9. Add interpolation counts
# ------------------------------------------------------------

figure6b_interpolation_summary <- data.frame(
  Uncertainty_measure = c(
    "Relative diagnostic uncertainty",
    "Relative coverage uncertainty",
    "Diagnostic variance"
  ),

  Directly_aggregated_districts = c(
    sum(
      !figure6b_district_values$
        Diagnostic_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      !figure6b_district_values$
        Coverage_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      !figure6b_district_values$
        Variance_value_interpolated,
      na.rm = TRUE
    )
  ),

  Interpolated_districts = c(
    sum(
      figure6b_district_values$
        Diagnostic_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      figure6b_district_values$
        Coverage_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      figure6b_district_values$
        Variance_value_interpolated,
      na.rm = TRUE
    )
  ),

  stringsAsFactors = FALSE
)

figure6b_summary_statistics <-
  figure6b_summary_statistics |>
  dplyr::left_join(
    figure6b_interpolation_summary,
    by = "Uncertainty_measure",
    relationship = "one-to-one"
  )

# ------------------------------------------------------------
# 10. Plot settings
# ------------------------------------------------------------

figure6b_plot_settings <- data.frame(
  Setting = c(
    "Administrative level",
    "Administrative identifier",
    "Prediction datasets",
    "Prediction grid cells",
    "Analysis scale",
    "Diagnostic variance",
    "Coverage variance",
    "Total variance",
    "Relative diagnostic uncertainty",
    "Relative coverage uncertainty",
    "District aggregation",
    "Missing-district interpolation",
    "Interpolation neighbours",
    "Interpolation weighting",
    "Relative-map colour limits",
    "Diagnostic-variance colour limits",
    "Colour palette",
    "Histogram bin width",
    "Histogram boundary",
    "Histogram closure",
    "Histogram x-axis limits",
    "Histogram x-axis breaks",
    "Boxplot width",
    "Survey-point size range",
    "Survey-point size breaks",
    "Survey-point shape",
    "Coordinate reference system"
  ),

  Value = c(
    "GADM level 3 Ethiopian administrative districts",
    "GID_3",
    as.character(n_datasets),
    as.character(n_units),
    "Logit prevalence scale",

    paste(
      "Variance across the 1,000",
      "corrected-prevalence model datasets"
    ),

    paste(
      "Mean squared within-dataset standard deviation,",
      "estimated from the 2.5th and 97.5th percentiles"
    ),

    "Diagnostic variance plus coverage variance",

    paste(
      "Diagnostic variance divided by",
      "total variance"
    ),

    paste(
      "Coverage variance divided by",
      "total variance"
    ),

    paste(
      "Arithmetic mean of prediction-grid values",
      "within each GADM level-3 district"
    ),

    paste(
      "Districts containing no prediction-grid points",
      "were assigned an inverse-distance-weighted value"
    ),

    "Five nearest prediction-grid points",

    paste(
      "Inverse distance, with 1e-10 added",
      "to avoid division by zero"
    ),

    "0 to 1",

    paste0(
      format(
        min(
          zone_between_var$zone_value,
          na.rm = TRUE
        ),
        scientific = TRUE,
        digits = 8
      ),
      " to ",
      format(
        max(
          zone_between_var$zone_value,
          na.rm = TRUE
        ),
        scientific = TRUE,
        digits = 8
      )
    ),

    "viridis magma",
    "0.04",
    "0",
    "Left-closed",
    "0 to 1",
    "0, 0.2, 0.4, 0.6, 0.8 and 1",
    "0.8",
    "1 to 10",
    "10, 50, 100 and 200",
    "HCT/BCT = 21; PCR = 22",
    "WGS84, EPSG:4326"
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 11. README
# ------------------------------------------------------------

figure6b_readme <- data.frame(
  Item = c(
    "Workbook description",
    "Associated figure",
    "District_values sheet",
    "Distribution_values sheet",
    "Grid_values sheet",
    "Survey_points sheet",
    "Plot_settings sheet",
    "Summary_statistics sheet",
    "Relative diagnostic uncertainty",
    "Relative coverage uncertainty",
    "Diagnostic variance",
    "Coverage variance",
    "Value_source",
    "Interpolated values",
    "Grid coordinate columns",
    "Administrative geometry",
    "Units"
  ),

  Description = c(
    paste(
      "Source data underlying Figure 6B, including",
      "Ethiopian district-level uncertainty measures,",
      "prediction-grid uncertainty values, distribution",
      "plot inputs and bovine survey-point overlays."
    ),

    "Figure 6B.",

    paste(
      "One row per unique Ethiopian GADM level-3",
      "administrative polygon. Contains all three mapped",
      "uncertainty measures, within-district summaries",
      "and interpolation indicators."
    ),

    paste(
      "Long-format district values used directly by",
      "the choropleths, histograms and boxplots."
    ),

    paste(
      "Prediction-grid uncertainty decomposition before",
      "aggregation to Ethiopian districts."
    ),

    paste(
      "Bovine survey locations plotted over each",
      "uncertainty choropleth."
    ),

    paste(
      "Transformations, interpolation rules, scale",
      "limits and plotting parameters."
    ),

    paste(
      "Descriptive summaries and counts of directly",
      "aggregated and interpolated districts."
    ),

    paste(
      "Proportion of total logit-scale variance",
      "attributed to variation across corrected-",
      "prevalence datasets."
    ),

    paste(
      "Proportion of total logit-scale variance",
      "attributed to uncertainty within each spatial",
      "model projection."
    ),

    paste(
      "Variance across corrected-prevalence model",
      "datasets on the logit prevalence scale."
    ),

    paste(
      "Mean within-dataset variance estimated from",
      "the 2.5th and 97.5th percentiles."
    ),

    paste(
      "Indicates whether the district value was",
      "calculated from prediction-grid points lying",
      "within the polygon or by interpolation."
    ),

    paste(
      "An interpolated district has",
      "Number_of_grid_points equal to zero and was",
      "assigned a value from the five nearest points."
    ),

    paste(
      "The projection files store geographical",
      "coordinates in reversed Longitude and Latitude",
      "columns. Plot_longitude and Plot_latitude are",
      "the corrected coordinates used for aggregation."
    ),

    paste(
      "Administrative geometry is obtained separately",
      "from GADM and is not embedded in the workbook."
    ),

    paste(
      "Relative uncertainty measures are proportions;",
      "variance measures are on the squared logit-",
      "prevalence scale."
    )
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 12. Create workbook
# ------------------------------------------------------------

figure6b_wb <- openxlsx::createWorkbook(
  creator = "AR Kaye",
  title = "Figure 6B source data",
  subject = paste(
    "Ethiopia uncertainty decomposition maps,",
    "histograms and boxplots"
  ),
  category = "Research data"
)

openxlsx::addWorksheet(
  figure6b_wb,
  "README"
)

openxlsx::addWorksheet(
  figure6b_wb,
  "District_values"
)

openxlsx::addWorksheet(
  figure6b_wb,
  "Distribution_values"
)

openxlsx::addWorksheet(
  figure6b_wb,
  "Grid_values"
)

openxlsx::addWorksheet(
  figure6b_wb,
  "Survey_points"
)

openxlsx::addWorksheet(
  figure6b_wb,
  "Plot_settings"
)

openxlsx::addWorksheet(
  figure6b_wb,
  "Summary_statistics"
)

openxlsx::writeData(
  figure6b_wb,
  "README",
  figure6b_readme
)

openxlsx::writeData(
  figure6b_wb,
  "District_values",
  figure6b_district_values
)

openxlsx::writeData(
  figure6b_wb,
  "Distribution_values",
  figure6b_distribution_values
)

openxlsx::writeData(
  figure6b_wb,
  "Grid_values",
  figure6b_grid_values
)

openxlsx::writeData(
  figure6b_wb,
  "Survey_points",
  figure6b_survey_points
)

openxlsx::writeData(
  figure6b_wb,
  "Plot_settings",
  figure6b_plot_settings
)

openxlsx::writeData(
  figure6b_wb,
  "Summary_statistics",
  figure6b_summary_statistics
)

# ------------------------------------------------------------
# 13. Workbook styles
# ------------------------------------------------------------

figure6b_header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "bottom",
  borderColour = "#FFFFFF"
)

figure6b_decimal_style <- openxlsx::createStyle(
  numFmt = "0.000000"
)

figure6b_variance_style <- openxlsx::createStyle(
  numFmt = "0.0000000000"
)

figure6b_integer_style <- openxlsx::createStyle(
  numFmt = "0"
)

figure6b_percentage_style <- openxlsx::createStyle(
  numFmt = "0.000%"
)

figure6b_wrap_style <- openxlsx::createStyle(
  wrapText = TRUE,
  valign = "top"
)

# ------------------------------------------------------------
# 14. Formatting helper
# ------------------------------------------------------------

figure6b_format_data_sheet <- function(
    workbook,
    sheet_name,
    data_object
) {

  if (ncol(data_object) < 1) {
    return(invisible(NULL))
  }

  openxlsx::addStyle(
    workbook,
    sheet = sheet_name,
    style = figure6b_header_style,
    rows = 1,
    cols = seq_len(ncol(data_object)),
    gridExpand = TRUE
  )

  openxlsx::setColWidths(
    workbook,
    sheet = sheet_name,
    cols = seq_len(ncol(data_object)),
    widths = "auto"
  )

  openxlsx::freezePane(
    workbook,
    sheet = sheet_name,
    firstRow = TRUE
  )

  openxlsx::addFilter(
    workbook,
    sheet = sheet_name,
    row = 1,
    cols = seq_len(ncol(data_object))
  )
}

figure6b_format_data_sheet(
  figure6b_wb,
  "District_values",
  figure6b_district_values
)

figure6b_format_data_sheet(
  figure6b_wb,
  "Distribution_values",
  figure6b_distribution_values
)

figure6b_format_data_sheet(
  figure6b_wb,
  "Grid_values",
  figure6b_grid_values
)

figure6b_format_data_sheet(
  figure6b_wb,
  "Survey_points",
  figure6b_survey_points
)

figure6b_format_data_sheet(
  figure6b_wb,
  "Summary_statistics",
  figure6b_summary_statistics
)

# ------------------------------------------------------------
# 15. README formatting
# ------------------------------------------------------------

openxlsx::addStyle(
  figure6b_wb,
  "README",
  figure6b_header_style,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::addStyle(
  figure6b_wb,
  "README",
  figure6b_wrap_style,
  rows = 2:(nrow(figure6b_readme) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::setColWidths(
  figure6b_wb,
  "README",
  cols = 1,
  widths = 32
)

openxlsx::setColWidths(
  figure6b_wb,
  "README",
  cols = 2,
  widths = 90
)

openxlsx::freezePane(
  figure6b_wb,
  "README",
  firstRow = TRUE
)

# ------------------------------------------------------------
# 16. Plot-settings formatting
# ------------------------------------------------------------

openxlsx::addStyle(
  figure6b_wb,
  "Plot_settings",
  figure6b_header_style,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::addStyle(
  figure6b_wb,
  "Plot_settings",
  figure6b_wrap_style,
  rows = 2:(nrow(figure6b_plot_settings) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::setColWidths(
  figure6b_wb,
  "Plot_settings",
  cols = 1,
  widths = 42
)

openxlsx::setColWidths(
  figure6b_wb,
  "Plot_settings",
  cols = 2,
  widths = 90
)

openxlsx::freezePane(
  figure6b_wb,
  "Plot_settings",
  firstRow = TRUE
)

# ------------------------------------------------------------
# 17. Numeric formatting: District_values
# ------------------------------------------------------------

figure6b_relative_columns <- match(
  c(
    "Relative_diagnostic_uncertainty",
    "Relative_coverage_uncertainty",
    "Relative_components_sum",
    "Relative_components_difference_from_one"
  ),
  names(figure6b_district_values)
)

figure6b_relative_columns <-
  figure6b_relative_columns[
    !is.na(figure6b_relative_columns)
  ]

if (
  nrow(figure6b_district_values) > 0 &&
  length(figure6b_relative_columns) > 0
) {

  openxlsx::addStyle(
    figure6b_wb,
    "District_values",
    figure6b_percentage_style,
    rows = 2:(nrow(figure6b_district_values) + 1),
    cols = figure6b_relative_columns,
    gridExpand = TRUE
  )
}

figure6b_variance_columns <- match(
  c(
    "Diagnostic_variance",
    "Diagnostic_variance_median",
    "Diagnostic_variance_minimum",
    "Diagnostic_variance_maximum"
  ),
  names(figure6b_district_values)
)

figure6b_variance_columns <-
  figure6b_variance_columns[
    !is.na(figure6b_variance_columns)
  ]

if (
  nrow(figure6b_district_values) > 0 &&
  length(figure6b_variance_columns) > 0
) {

  openxlsx::addStyle(
    figure6b_wb,
    "District_values",
    figure6b_variance_style,
    rows = 2:(nrow(figure6b_district_values) + 1),
    cols = figure6b_variance_columns,
    gridExpand = TRUE
  )
}

figure6b_integer_columns <- match(
  c(
    "Number_of_grid_points",
    "Coverage_number_of_grid_points",
    "Variance_number_of_grid_points"
  ),
  names(figure6b_district_values)
)

figure6b_integer_columns <-
  figure6b_integer_columns[
    !is.na(figure6b_integer_columns)
  ]

if (
  nrow(figure6b_district_values) > 0 &&
  length(figure6b_integer_columns) > 0
) {

  openxlsx::addStyle(
    figure6b_wb,
    "District_values",
    figure6b_integer_style,
    rows = 2:(nrow(figure6b_district_values) + 1),
    cols = figure6b_integer_columns,
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 18. Numeric formatting: Distribution_values
# ------------------------------------------------------------

figure6b_distribution_value_column <- match(
  "Plotted_value",
  names(figure6b_distribution_values)
)

if (
  nrow(figure6b_distribution_values) > 0 &&
  !is.na(figure6b_distribution_value_column)
) {

  openxlsx::addStyle(
    figure6b_wb,
    "Distribution_values",
    figure6b_decimal_style,
    rows = 2:(nrow(figure6b_distribution_values) + 1),
    cols = figure6b_distribution_value_column,
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 19. Numeric formatting: Grid_values
# ------------------------------------------------------------

figure6b_grid_coordinate_columns <- match(
  c(
    "Original_Latitude_column",
    "Original_Longitude_column",
    "Plot_longitude",
    "Plot_latitude"
  ),
  names(figure6b_grid_values)
)

figure6b_grid_coordinate_columns <-
  figure6b_grid_coordinate_columns[
    !is.na(figure6b_grid_coordinate_columns)
  ]

if (
  nrow(figure6b_grid_values) > 0 &&
  length(figure6b_grid_coordinate_columns) > 0
) {

  openxlsx::addStyle(
    figure6b_wb,
    "Grid_values",
    figure6b_decimal_style,
    rows = 2:(nrow(figure6b_grid_values) + 1),
    cols = figure6b_grid_coordinate_columns,
    gridExpand = TRUE
  )
}

figure6b_grid_variance_columns <- match(
  c(
    "Diagnostic_variance",
    "Coverage_variance",
    "Total_variance"
  ),
  names(figure6b_grid_values)
)

figure6b_grid_variance_columns <-
  figure6b_grid_variance_columns[
    !is.na(figure6b_grid_variance_columns)
  ]

if (
  nrow(figure6b_grid_values) > 0 &&
  length(figure6b_grid_variance_columns) > 0
) {

  openxlsx::addStyle(
    figure6b_wb,
    "Grid_values",
    figure6b_variance_style,
    rows = 2:(nrow(figure6b_grid_values) + 1),
    cols = figure6b_grid_variance_columns,
    gridExpand = TRUE
  )
}

figure6b_grid_relative_columns <- match(
  c(
    "Relative_coverage_uncertainty",
    "Relative_diagnostic_uncertainty",
    "Relative_components_sum"
  ),
  names(figure6b_grid_values)
)

figure6b_grid_relative_columns <-
  figure6b_grid_relative_columns[
    !is.na(figure6b_grid_relative_columns)
  ]

if (
  nrow(figure6b_grid_values) > 0 &&
  length(figure6b_grid_relative_columns) > 0
) {

  openxlsx::addStyle(
    figure6b_wb,
    "Grid_values",
    figure6b_percentage_style,
    rows = 2:(nrow(figure6b_grid_values) + 1),
    cols = figure6b_grid_relative_columns,
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 20. Numeric formatting: Survey_points
# ------------------------------------------------------------

if (nrow(figure6b_survey_points) > 0) {

  openxlsx::addStyle(
    figure6b_wb,
    "Survey_points",
    figure6b_decimal_style,
    rows = 2:(nrow(figure6b_survey_points) + 1),
    cols = c(2, 3),
    gridExpand = TRUE
  )

  openxlsx::addStyle(
    figure6b_wb,
    "Survey_points",
    figure6b_integer_style,
    rows = 2:(nrow(figure6b_survey_points) + 1),
    cols = c(1, 4, 5),
    gridExpand = TRUE
  )

  openxlsx::addStyle(
    figure6b_wb,
    "Survey_points",
    figure6b_percentage_style,
    rows = 2:(nrow(figure6b_survey_points) + 1),
    cols = 6,
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 21. Save workbook
# ------------------------------------------------------------

figure6b_output_file <- file.path(
  figure6b_output_dir,
  "Figure_6B_source_data.xlsx"
)

openxlsx::saveWorkbook(
  figure6b_wb,
  file = figure6b_output_file,
  overwrite = TRUE
)

message(
  "Figure 6B source-data workbook saved to: ",
  figure6b_output_file
)

message(
  "District rows exported: ",
  nrow(figure6b_district_values)
)

message(
  "Distribution values exported: ",
  nrow(figure6b_distribution_values)
)

message(
  "Prediction-grid rows exported: ",
  nrow(figure6b_grid_values)
)

message(
  "Survey points exported: ",
  nrow(figure6b_survey_points)
)