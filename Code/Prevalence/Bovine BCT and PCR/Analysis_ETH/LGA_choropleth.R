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

# Get Ethiopian administrative boundaries at zone level (level 2)
# Ethiopia GADM data goes to level 2 (Zones), similar to Nigeria's LGAs
ethiopia_zones <- gadm(country = "ETH", level = 3, path = tempdir())
ethiopia_zones_sf <- st_as_sf(ethiopia_zones)


# Also get region boundaries for context
ethiopia_regions <- gadm(country = "ETH", level = 1, path = tempdir()) 
ethiopia_regions_sf <- st_as_sf(ethiopia_regions)

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

# Load continental Africa for inset map
cat("Loading continental Africa for inset map...\n")
africa_countries_sf <- NULL
ethiopia_country_sf <- NULL
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
    # Separate Ethiopia for highlighting
    ethiopia_country_sf <- africa_countries_sf %>% filter(GID_0 == "ETH")
    cat("Successfully loaded continental Africa with", nrow(africa_countries_sf), "countries\n")
  }
}, error = function(e) {
  cat("Failed to load continental Africa:", e$message, "\n")
})

# Load and average all projection files in Projections_ETH directory
projections_dir <- "Code/Prevalence/Bovine BCT and PCR/Projections_ETH/"
projection_files <- list.files(projections_dir, pattern = "Projections_model_.*\\.csv", full.names = TRUE)

cat("Found", length(projection_files), "projection files to average\n")

# Function to load and process a single projection file
load_projection_file <- function(file_path) {
  tryCatch({
    data <- read.csv(file_path)
    # Add model identifier
    data$model_id <- basename(file_path)
    return(data)
  }, error = function(e) {
    cat("Error loading", file_path, ":", e$message, "\n")
    return(NULL)
  })
}

# Load all projection files
cat("Loading all projection files...\n")
all_projections <- do.call(rbind, lapply(projection_files, load_projection_file))

# Remove any NULL entries (failed loads)
all_projections <- all_projections[!is.null(all_projections), ]

cat("Total rows loaded:", nrow(all_projections), "\n")
cat("Unique models:", length(unique(all_projections$model_id)), "\n")

# Average across all models for each coordinate and variable combination
cat("Computing ensemble averages...\n")
ensemble_projections <- all_projections %>%
  dplyr::group_by(Latitude, Longitude, variable) %>%
  dplyr::summarise(
    value = mean(value, na.rm = TRUE),
    n_models = n(),
    .groups = 'drop'
  )

cat("Ensemble projections created with", nrow(ensemble_projections), "unique coordinate-variable combinations\n")
cat("Variables available:", unique(ensemble_projections$variable), "\n")

# Prepare data for mean, 2.5th and 97.5th percentiles from ensemble
projections_mean <- ensemble_projections %>%
  dplyr::filter(variable == "Mean") %>%
  dplyr::select(Longitude, Latitude, prevalence = value)

projections_lower <- ensemble_projections %>%
  dplyr::filter(variable == "2.5th percentile") %>%
  dplyr::select(Longitude, Latitude, prevalence = value)

projections_upper <- ensemble_projections %>%
  dplyr::filter(variable == "97.5th percentile") %>%
  dplyr::select(Longitude, Latitude, prevalence = value)

# Function to process projections data (check if coordinates need swapping for Ethiopia)
process_projections <- function(proj_data) {
  # Ethiopia coordinate ranges: Latitude 3-15°N, Longitude 33-48°E
  # Check if coordinates are in the correct range
  lat_in_range <- all(proj_data$Latitude >= 3 & proj_data$Latitude <= 15, na.rm = TRUE)
  lon_in_range <- all(proj_data$Longitude >= 33 & proj_data$Longitude <= 48, na.rm = TRUE)
  
  # If coordinates appear swapped, swap them
  if (!lat_in_range || !lon_in_range) {
    cat("Coordinates appear to be swapped, correcting...\n")
    temp <- proj_data$Longitude
    proj_data$Longitude <- proj_data$Latitude
    proj_data$Latitude <- temp
  }
  
  # Remove any rows with missing coordinates
  proj_data <- proj_data[!is.na(proj_data$Longitude) & !is.na(proj_data$Latitude), ]
  
  # Convert to spatial points
  proj_sf <- st_as_sf(proj_data, 
                      coords = c("Longitude", "Latitude"), 
                      crs = 4326)
  
  return(proj_sf)
}

# Process all three datasets
projections_mean_sf <- process_projections(projections_mean)
projections_lower_sf <- process_projections(projections_lower)
projections_upper_sf <- process_projections(projections_upper)

# Function to aggregate prevalence data by zone using spatial join with gap filling
aggregate_by_zone <- function(projections_sf) {
  joined_data <- st_join(projections_sf, ethiopia_zones_sf)
  
  # First, get zones with direct data points
  zone_prevalence <- joined_data %>%
    st_drop_geometry() %>%
    dplyr::filter(!is.na(GID_3)) %>%  # Only keep points that fell within zones
    dplyr::group_by(GID_3, NAME_2, NAME_1) %>%
    dplyr::summarise(
      zone_prevalence = mean(prevalence, na.rm = TRUE),
      median_prevalence = median(prevalence, na.rm = TRUE),
      min_prevalence = min(prevalence, na.rm = TRUE),
      max_prevalence = max(prevalence, na.rm = TRUE),
      n_points = n(),
      .groups = 'drop'
    )
  
  # Find zones without data (grey areas)
  all_zones <- ethiopia_zones_sf %>%
    st_drop_geometry() %>%
    dplyr::select(GID_3, NAME_2, NAME_1)
  
  missing_zones <- all_zones %>%
    dplyr::anti_join(zone_prevalence, by = "GID_3")
  
  cat("Found", nrow(missing_zones), "zones without data. Filling using spatial interpolation...\n")
  
  if(nrow(missing_zones) > 0) {
    # For each missing zone, find closest projection points and interpolate
    for(i in seq_len(nrow(missing_zones))) {
      zone_id <- missing_zones$GID_3[i]
      
      # Get centroid of the missing zone (suppress attribute warning)
      zone_centroid <- ethiopia_zones_sf %>%
        dplyr::filter(GID_3 == zone_id) %>%
        st_centroid(of_largest_polygon = TRUE)
      
      # Find distances to all projection points
      distances <- st_distance(zone_centroid, projections_sf)
      
      # Convert distances to numeric (remove units)
      distances_numeric <- as.numeric(distances)
      
      # Get indices of 5 closest points (or fewer if not enough points available)
      n_neighbors <- min(5, nrow(projections_sf))
      closest_indices <- order(distances_numeric)[1:n_neighbors]
      
      # Calculate inverse distance weighted average
      closest_points <- projections_sf[closest_indices, ]
      weights <- 1 / (distances_numeric[closest_indices] + 1e-10)  # Add small value to avoid division by zero
      weights <- weights / sum(weights)  # Normalize weights
      
      # Calculate weighted average prevalence
      interpolated_prevalence <- sum(closest_points$prevalence * weights)
      
      # Add to zone_prevalence data
      new_row <- data.frame(
        GID_3 = zone_id,
        NAME_2 = missing_zones$NAME_2[i],
        NAME_1 = missing_zones$NAME_1[i],
        zone_prevalence = interpolated_prevalence,
        median_prevalence = interpolated_prevalence,
        min_prevalence = interpolated_prevalence,
        max_prevalence = interpolated_prevalence,
        n_points = 0  # Indicate this is interpolated
      )
      
      zone_prevalence <- rbind(zone_prevalence, new_row)
    }
  }
  
  return(zone_prevalence)
}

# Aggregate all three datasets
zone_prevalence_mean_raw <- aggregate_by_zone(projections_mean_sf)
zone_prevalence_lower_raw <- aggregate_by_zone(projections_lower_sf)
zone_prevalence_upper_raw <- aggregate_by_zone(projections_upper_sf)

# Standardize column names for consistency with plotting functions
zone_prevalence_mean <- zone_prevalence_mean_raw %>%
  dplyr::rename(mean_prevalence = zone_prevalence)

zone_prevalence_lower <- zone_prevalence_lower_raw %>%
  dplyr::mutate(
    lower_prevalence = zone_prevalence,
    mean_prevalence = zone_prevalence
  ) %>%
  dplyr::select(-zone_prevalence)

zone_prevalence_upper <- zone_prevalence_upper_raw %>%
  dplyr::mutate(
    upper_prevalence = zone_prevalence,
    mean_prevalence = zone_prevalence
  ) %>%
  dplyr::select(-zone_prevalence)

# Calculate 95% CI width for uncertainty assessment
zone_uncertainty <- zone_prevalence_lower %>%
  dplyr::select(GID_3, lower = lower_prevalence) %>%
  left_join(zone_prevalence_upper %>% dplyr::select(GID_3, upper = upper_prevalence), by = "GID_3") %>%
  dplyr::mutate(
    ci_width = upper - lower,
    high_uncertainty = ci_width > 0.80
  )

cat("Zones with high uncertainty (CI width > 0.80):", sum(zone_uncertainty$high_uncertainty), "out of", nrow(zone_uncertainty), "\n")

# Load bovine data from Excel files for overlaying on percentile maps
bovine_bct_raw <- readxl::read_excel("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
bovine_pcr_raw <- readxl::read_excel("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")

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
  
  # Convert to spatial points and filter for Ethiopian coordinate ranges
  # Ethiopia coordinate ranges: Latitude 3-15°N, Longitude 33-48°E
  bovine_ethiopia_filtered <- bovine_clean %>%
    dplyr::filter(Longitude >= 33, Longitude <= 48, Latitude >= 3, Latitude <= 15)
  
  if(nrow(bovine_ethiopia_filtered) == 0) {
    cat("No bovine data points found within Ethiopian boundaries\n")
    return(NULL)
  }
  
  bovine_sf <- st_as_sf(bovine_ethiopia_filtered, 
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
if(!is.null(bovine_ethiopia_sf)) {
  cat("Found", nrow(bovine_ethiopia_sf), "bovine data points within Ethiopia\n")
  cat("  - BCT tests:", sum(bovine_ethiopia_sf$Test_Type == "BCT/HCT"), "\n")
  cat("  - PCR tests:", sum(bovine_ethiopia_sf$Test_Type == "PCR"), "\n")
} else {
  cat("No bovine data available for Ethiopia - maps will show choropleth only\n")
}

# Function to create choropleth map (matching Nigerian version exactly)
create_choropleth <- function(lga_data, title_suffix, data_type = "mean", add_bovine = FALSE, show_uncertainty = FALSE, add_inset = FALSE) {
  
  # Handle the different data types
  if(data_type == "mean") {
    fill_column <- "mean_prevalence"
  } else if(data_type == "lower") {
    fill_column <- "lower_prevalence"  
  } else if(data_type == "upper") {
    fill_column <- "upper_prevalence"
  } else {
    fill_column <- data_type
  }
  
  # Join zone-level prevalence back to zone polygons
  ethiopia_zones_final <- ethiopia_zones_sf %>%
    left_join(lga_data %>% st_drop_geometry() %>% dplyr::select(GID_3, all_of(fill_column)), 
              by = "GID_3")
  
  # Add uncertainty information if requested
  if(show_uncertainty) {
    ethiopia_zones_final <- ethiopia_zones_final %>%
      left_join(zone_uncertainty, by = "GID_3")
  }
  
  # Create zone-level choropleth map (matching Nigerian styling exactly)
  p <- ggplot(ethiopia_zones_final)
  
  # Add neighboring countries as background context (if available)
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Add Ethiopia zones with prevalence data
  p <- p + 
    geom_sf(aes(fill = !!sym(fill_column)), lwd = 0.1, color = "white") +
    geom_sf(data = ethiopia_regions_sf, fill = NA, color = "black", lwd = 0.5)
  
  # Add high uncertainty overlay if requested
  if(show_uncertainty && "high_uncertainty" %in% names(ethiopia_zones_final)) {
    high_uncertainty_zones <- ethiopia_zones_final %>% 
      dplyr::filter(high_uncertainty == TRUE)
    
    if(nrow(high_uncertainty_zones) > 0) {
      p <- p + 
        geom_sf(data = high_uncertainty_zones, fill = "white", alpha = 0.4, 
                color = "red", lwd = 0.3)
    }
  }
  
  p <- p + scale_fill_viridis_c(
      name = "AAT\nprevalence",
      na.value = "grey90",
      direction = 1,
      option = "viridis",
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 0.1)
    ) +
    # Add scale bar and north arrow
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                    bar_cols = c("black", "white"), line_width = 1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                          pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                          style = north_arrow_fancy_orienteering(text_size = 8)) +
    theme_void()
  
  p <- p +
    labs(
      title = paste("AAT prevalence by district in Ethiopia", title_suffix)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0.5, size = 9, color = "grey60"),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    )
  
  # Add uncertainty caption if high uncertainty areas exist (like Nigerian version)
  if(show_uncertainty && "high_uncertainty" %in% names(ethiopia_zones_final) && 
     any(ethiopia_zones_final$high_uncertainty, na.rm = TRUE)) {
    p <- p + labs(caption = "Areas with white overlay and red border indicate high uncertainty (CI width > 0.80)")
  }
  
  # Add bovine data points if requested and available (for percentile maps)
  if(add_bovine && !is.null(bovine_ethiopia_sf) && nrow(bovine_ethiopia_sf) > 0) {
    p <- p + geom_point(data = st_drop_geometry(bovine_ethiopia_sf), 
                       aes(x = lon, y = lat, size = Number_of_animal_tested, color = Test_Type),
                       alpha = 0.8)  # Increased alpha for better visibility
    
    # Add additional scales for bovine data
    p <- p +
      # Add size legend for sample size
      scale_size_continuous(name = "Sample size", 
                           range = c(1, 6),  # Increased from c(1, 4) to make sizes more distinguishable
                           breaks = c(10, 50, 100, 200),
                           labels = c("10", "50", "100", "200+")) +
      # Add color scale for test types
      scale_color_manual(name = "Test type",
                        values = c("BCT/HCT" = "red", "PCR" = "blue"),
                        guide = guide_legend(override.aes = list(size = 3, alpha = 0.8)))
  }
  
  # Add inset map if requested and data is available
  if(add_inset && !is.null(africa_countries_sf) && !is.null(ethiopia_country_sf)) {
    # Create inset map showing Ethiopia's location in Africa
    inset_plot <- ggplot() +
      geom_sf(data = africa_countries_sf, fill = "grey90", color = "white", size = 0.1) +
      geom_sf(data = ethiopia_country_sf, fill = "#440154", color = "white", size = 0.3) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 0.5),
        plot.margin = margin(2.1, 2.1, 2.1, 2.1)
      ) +
      coord_sf(expand = FALSE)
    
    # Combine main plot with inset using cowplot
    p <- ggdraw(p) +
      draw_plot(inset_plot, x = 0.70, y = 0.65, width = 0.25, height = 0.25)
  }
  
  return(p)
}

# Create all three maps using zone_prevalence data (renamed to lga_data for consistency)
lga_prevalence_mean <- zone_prevalence_mean
lga_prevalence_lower <- zone_prevalence_lower
lga_prevalence_upper <- zone_prevalence_upper

# Create choropleth maps - mean map first, then percentile maps with bovine data overlay
p_mean <- create_choropleth(zone_prevalence_mean, "(mean)", show_uncertainty = TRUE, add_inset = TRUE)
p_lower <- create_choropleth(zone_prevalence_lower, "(2.5th percentile)", data_type = "lower", add_bovine = TRUE, show_uncertainty = TRUE)
p_upper <- create_choropleth(zone_prevalence_upper, "(97.5th percentile)", data_type = "upper", add_bovine = TRUE, show_uncertainty = TRUE)

# Create combined histogram and boxplot for each estimate type
create_combined_plot <- function(lga_data, estimate_name) {
  # Get zone data without geometry for plotting
  plot_data <- lga_data %>% st_drop_geometry() %>% dplyr::filter(!is.na(mean_prevalence))
  
  # Histogram with viridis color scale matching choropleth
  p_hist <- ggplot(plot_data, aes(x = mean_prevalence, fill = after_stat(x))) +
    geom_histogram(binwidth = 0.04, color = "white", linewidth = 0.2,
                   boundary = 0, closed = "left") +
    scale_fill_viridis_c(option = "viridis", guide = "none") +  # Same as choropleth, no legend
    scale_x_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.25),
                       limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
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
      title = paste("Distribution of", estimate_name, "prevalence by district in Ethiopia"),
      x = ""
    )
  
  # Boxplot with viridis color matching choropleth
  viridis_color <- viridis::viridis(1, begin = 0.5, end = 0.5, option = "viridis")  # Mid-range viridis color
  
  p_box <- ggplot(plot_data, aes(x = mean_prevalence, y = 1)) +  # Use y = 1 instead of y = ""
    geom_boxplot(fill = viridis_color, alpha = 0.8, width = 0.8) +  # Much wider boxplot, no jitter
    scale_x_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.25),
                       limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0.2, 1.8), expand = c(0, 0)) +  # More vertical space
    theme_minimal() +  # Use theme_minimal instead of theme_void for better visibility
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
      x = "AAT prevalence"
    )
  
  # Use a simpler approach with grid.arrange for better stability
  library(gridExtra)
  library(grid)
  
  # Create title
  title <- grid::textGrob(paste("AAT prevalence Analysis -", estimate_name), 
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

# Create the three combined plots
p_mean_combined <- create_combined_plot(lga_prevalence_mean, "mean")
p_lower_combined <- create_combined_plot(lga_prevalence_lower, "2.5th percentile") 
p_upper_combined <- create_combined_plot(lga_prevalence_upper, "97.5th percentile")

# Save all plots to Analysis_ETH directory
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/eth_choropleth_mean.pdf", plot = p_mean, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/eth_choropleth_lower.pdf", plot = p_lower, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/eth_choropleth_upper.pdf", plot = p_upper, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each estimate type) - more compact dimensions
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/eth_mean_analysis.pdf", plot = p_mean_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/eth_lower_analysis.pdf", plot = p_lower_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/eth_upper_analysis.pdf", plot = p_upper_combined, width = 10, height = 6)

# Display all plots
print(p_mean)
print(p_lower) 
print(p_upper)
grid.draw(p_mean_combined)
grid.draw(p_lower_combined)
grid.draw(p_upper_combined)