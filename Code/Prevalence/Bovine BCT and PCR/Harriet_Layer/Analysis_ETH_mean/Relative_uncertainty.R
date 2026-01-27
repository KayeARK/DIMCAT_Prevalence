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
n_datasets <- 200
sample_file <- "Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Projections_ETH_mean/Projections_model_1.csv"
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
dpm <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Projections_ETH_mean/Projections_model_",i,".csv"))


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
ethiopia_zones <- gadm(country = "ETH", level = 2, path = tempdir()) 
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

# Function to process uncertainty data (swap coordinates and create spatial points)
process_uncertainty_data <- function(uncertainty_data, value_col) {
  # Create a copy and select relevant columns
  proj_data <- uncertainty_data[, c("Longitude", "Latitude", value_col)]
  names(proj_data)[3] <- "value"  # Rename the value column for consistency
  
  # Swap longitude and latitude (they are reversed in the data)
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
    dplyr::filter(!is.na(GID_2)) %>%  # Only keep points that fell within zones
    dplyr::group_by(GID_2, NAME_2) %>%
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
    dplyr::select(GID_2, NAME_2)
  
  missing_zones <- all_zones %>%
    dplyr::anti_join(zone_uncertainty, by = "GID_2")
  
  cat("Found", nrow(missing_zones), "zones without uncertainty data. Filling using spatial interpolation...\n")
  
  if(nrow(missing_zones) > 0) {
    # For each missing zone, find closest uncertainty points and interpolate
    for(i in seq_len(nrow(missing_zones))) {
      zone_id <- missing_zones$GID_2[i]
      
      # Get centroid of the missing zone (suppress attribute warning)
      zone_centroid <- ethiopia_zones_sf %>%
        dplyr::filter(GID_2 == zone_id) %>%
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
        GID_2 = zone_id,
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
    left_join(zone_data, by = "GID_2")
  
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
    geom_sf(aes(fill = zone_value), lwd = 0.1, color = "white") +
    geom_sf(data = ethiopia_regions_sf, fill = NA, color = "black", lwd = 0.3) +  # Add region boundaries
    # Add bovine data points (different colors for BCT and PCR)
    geom_point(data = st_drop_geometry(bovine_ethiopia_sf), 
               aes(x = lon, y = lat, size = Number_of_animal_tested, color = Test_Type),
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
                         range = c(1, 6),  # Increased from c(1, 4) to make sizes more distinguishable
                         breaks = c(10, 50, 100, 200),
                         labels = c("10", "50", "100", "200+")) +
    # Add color scale for test types
    scale_color_manual(name = "Test type",
                      values = c("BCT/HCT" = "red", "PCR" = "blue"),
                      guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    # Add scale bar and north arrow
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                    bar_cols = c("black", "white"), line_width = 1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                          pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                          style = north_arrow_fancy_orienteering(text_size = 8)) +
    theme_void() +
    labs(
      title = paste("Uncertainty analysis by zone in Ethiopia", title_suffix)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
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
                       breaks = scales::pretty_breaks(n = 5),
                       limits = c(min(zone_data$zone_value, na.rm = TRUE), max(zone_data$zone_value, na.rm = TRUE)),
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
      title = paste("Distribution of", estimate_name, "by zone"),
      x = ""
    )
  
  # Boxplot with viridis color matching choropleth
  viridis_color <- viridis::viridis(1, begin = 0.5, end = 0.5, option = "magma")  # Mid-range viridis color
  
  p_box <- ggplot(zone_data, aes(x = zone_value, y = 1)) +
    geom_boxplot(fill = viridis_color, alpha = 0.8, width = 0.8) +  # Much wider boxplot, no jitter
    scale_x_continuous(labels = if(max(zone_data$zone_value, na.rm = TRUE) <= 1) scales::percent_format() else scales::number_format(accuracy = 0.001),
                       breaks = scales::pretty_breaks(n = 5),
                       limits = c(min(zone_data$zone_value, na.rm = TRUE), max(zone_data$zone_value, na.rm = TRUE)),
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
p_rel_between <- create_uncertainty_choropleth(zone_rel_between, "(relative variability from data uncertainty)", "Relative\nimportance", c(0, 1))
p_rel_within <- create_uncertainty_choropleth(zone_rel_within, "(relative variability from model uncertainty)", "Relative\nimportance", c(0, 1))
p_between_var <- create_uncertainty_choropleth(zone_between_var, "(Variance from data uncertainty)", "Variance", c(min(zone_between_var$zone_value, na.rm = TRUE), max(zone_between_var$zone_value, na.rm = TRUE)))

# Create combined histogram and boxplot analyses
p_rel_between_combined <- create_uncertainty_combined_plot(zone_rel_between, "relative variability from data uncertainty", "Relative importance of variability from data uncertainty")
p_rel_within_combined <- create_uncertainty_combined_plot(zone_rel_within, "relative variability from model uncertainty", "Relative importance of variability from model uncertainty")
p_between_var_combined <- create_uncertainty_combined_plot(zone_between_var, "variance from data uncertainty", "Variance from data uncertainty")

# Save all choropleth maps
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/region_rel_between_choropleth.pdf", plot = p_rel_between, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/region_rel_within_choropleth.pdf", plot = p_rel_within, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/region_between_var_choropleth.pdf", plot = p_between_var, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each uncertainty measure)
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/region_rel_between_analysis.pdf", plot = p_rel_between_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/region_rel_within_analysis.pdf", plot = p_rel_within_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/region_between_var_analysis.pdf", plot = p_between_var_combined, width = 10, height = 6)

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
ggplot(dpm_var) + geom_tile(aes(Latitude, Longitude, fill = rel_between)) +
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
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/Relative_uncertainty_between.png",width = 12, height = 8)


#plot histogram of dpm_var$rel_between
ggplot(dpm_var, aes(x = rel_between)) +
  geom_histogram(binwidth = 0.05, fill = "black", color = "black", alpha = 0.7) +
  labs(
       x = "Relative importance of between-dataset variability",
       y = "Frequency")
ggsave("Code/Prevalence/Bovine BCT and PCR/Harriet_Layer/Analysis_ETH_mean/Relative_uncertainty_between_histogram.png",width = 8, height = 6)
