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

# Get Nigerian administrative boundaries at LGA level (level 2)
# Nigeria GADM data only goes to level 2 (Local Government Areas), not level 3 (wards)
nigeria_lgas <- gadm(country = "NGA", level = 2, path = tempdir())
nigeria_lgas_sf <- st_as_sf(nigeria_lgas)

# Also get state boundaries for context
nigeria_states <- gadm(country = "NGA", level = 1, path = tempdir()) 
nigeria_states_sf <- st_as_sf(nigeria_states)

# Get neighboring countries for geographic context
cat("Loading neighboring countries for geographic context...\n")
neighboring_countries <- c("BEN", "NER", "TCD", "CMR", "CAF")  # Benin, Niger, Chad, Cameroon, Central African Republic
neighbor_sf_list <- list()

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

# Load and average all projection files in Projections_NGA directory
projections_dir <- "Code/Prevalence/Bovine BCT and PCR/Projections_NGA/"
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
  filter(variable == "Mean") %>%
  dplyr::select(Longitude, Latitude, prevalence = value)

projections_lower <- ensemble_projections %>%
  filter(variable == "2.5th percentile") %>%
  dplyr::select(Longitude, Latitude, prevalence = value)

projections_upper <- ensemble_projections %>%
  filter(variable == "97.5th percentile") %>%
  dplyr::select(Longitude, Latitude, prevalence = value)

# Function to process projections data (swap coordinates and create spatial points)
process_projections <- function(proj_data) {
  # Swap longitude and latitude (they are reversed in the data)
  temp <- proj_data$Longitude
  proj_data$Longitude <- proj_data$Latitude
  proj_data$Latitude <- temp
  
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

# Function to aggregate prevalence data by LGA using spatial join with gap filling
aggregate_by_lga <- function(projections_sf) {
  joined_data <- st_join(projections_sf, nigeria_lgas_sf)
  
  # First, get LGAs with direct data points
  lga_prevalence <- joined_data %>%
    st_drop_geometry() %>%
    dplyr::filter(!is.na(GID_2)) %>%  # Only keep points that fell within LGAs
    dplyr::group_by(GID_2, NAME_2, NAME_1) %>%
    dplyr::summarise(
      lga_prevalence = mean(prevalence, na.rm = TRUE),
      median_prevalence = median(prevalence, na.rm = TRUE),
      min_prevalence = min(prevalence, na.rm = TRUE),
      max_prevalence = max(prevalence, na.rm = TRUE),
      n_points = n(),
      .groups = 'drop'
    )
  
  # Find LGAs without data (grey areas)
  all_lgas <- nigeria_lgas_sf %>%
    st_drop_geometry() %>%
    dplyr::select(GID_2, NAME_2, NAME_1)
  
  missing_lgas <- all_lgas %>%
    dplyr::anti_join(lga_prevalence, by = "GID_2")
  
  cat("Found", nrow(missing_lgas), "LGAs without data. Filling using spatial interpolation...\n")
  
  if(nrow(missing_lgas) > 0) {
    # For each missing LGA, find closest projection points and interpolate
    for(i in seq_len(nrow(missing_lgas))) {
      lga_id <- missing_lgas$GID_2[i]
      
      # Get centroid of the missing LGA (suppress attribute warning)
      lga_centroid <- nigeria_lgas_sf %>%
        dplyr::filter(GID_2 == lga_id) %>%
        st_centroid(of_largest_polygon = TRUE)
      
      # Find distances to all projection points
      distances <- st_distance(lga_centroid, projections_sf)
      
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
      
      # Add to lga_prevalence data
      new_row <- data.frame(
        GID_2 = lga_id,
        NAME_2 = missing_lgas$NAME_2[i],
        NAME_1 = missing_lgas$NAME_1[i],
        lga_prevalence = interpolated_prevalence,
        median_prevalence = interpolated_prevalence,
        min_prevalence = interpolated_prevalence,
        max_prevalence = interpolated_prevalence,
        n_points = 0  # Indicate this is interpolated
      )
      
      lga_prevalence <- rbind(lga_prevalence, new_row)
    }
  }
  
  return(lga_prevalence)
}

# Aggregate all three datasets
lga_prevalence_mean <- aggregate_by_lga(projections_mean_sf)
lga_prevalence_lower <- aggregate_by_lga(projections_lower_sf)
lga_prevalence_upper <- aggregate_by_lga(projections_upper_sf)

# Calculate 95% CI width for uncertainty assessment
lga_uncertainty <- lga_prevalence_lower %>%
  dplyr::select(GID_2, lower = lga_prevalence) %>%
  left_join(lga_prevalence_upper %>% dplyr::select(GID_2, upper = lga_prevalence), by = "GID_2") %>%
  dplyr::mutate(
    ci_width = upper - lower,
    high_uncertainty = ci_width > 0.80
  )

cat("LGAs with high uncertainty (CI width > 0.80):", sum(lga_uncertainty$high_uncertainty), "out of", nrow(lga_uncertainty), "\n")

# Load bovine data from Excel files for overlaying on percentile maps
bovine_bct_raw <- readxl::read_excel("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
bovine_pcr_raw <- readxl::read_excel("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")

# Add test type identifier to each dataset
bovine_bct_raw$Test_Type <- "BCT"
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
  
  # Filter for points within Nigeria using spatial intersection
  bovine_nigeria <- st_filter(bovine_sf, nigeria_states_sf)
  
  # Add back coordinate columns for plotting
  coords <- st_coordinates(bovine_nigeria)
  bovine_nigeria$lon <- coords[, 1]
  bovine_nigeria$lat <- coords[, 2]
  
  return(bovine_nigeria)
}

# Process the bovine data
bovine_nigeria_sf <- process_bovine_data(bovine_data_raw)
cat("Found", nrow(bovine_nigeria_sf), "bovine data points within Nigeria\n")
cat("  - BCT tests:", sum(bovine_nigeria_sf$Test_Type == "BCT"), "\n")
cat("  - PCR tests:", sum(bovine_nigeria_sf$Test_Type == "PCR"), "\n")

# Function to create choropleth map
create_choropleth <- function(lga_data, title_suffix, add_bovine_data = FALSE, show_uncertainty = FALSE, add_inset = FALSE) {
  # Join LGA-level prevalence back to LGA polygons
  nigeria_lgas_final <- nigeria_lgas_sf %>%
    left_join(lga_data, by = "GID_2")
  
  # Add uncertainty information if requested
  if(show_uncertainty) {
    nigeria_lgas_final <- nigeria_lgas_final %>%
      left_join(lga_uncertainty, by = "GID_2")
  }
  
  # Create LGA-level choropleth map
  p <- ggplot(nigeria_lgas_final)
  
  # Add neighboring countries as background context (if available)
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Add Nigeria LGAs with prevalence data
  p <- p + 
    geom_sf(aes(fill = lga_prevalence), lwd = 0.1, color = "white") +
    geom_sf(data = nigeria_states_sf, fill = NA, color = "black", lwd = 0.5)
  
  # Add high uncertainty overlay if requested
  if(show_uncertainty && "high_uncertainty" %in% names(nigeria_lgas_final)) {
    high_uncertainty_lgas <- nigeria_lgas_final %>% 
      dplyr::filter(high_uncertainty == TRUE)
    
    if(nrow(high_uncertainty_lgas) > 0) {
      p <- p + 
        geom_sf(data = high_uncertainty_lgas, fill = "white", alpha = 0.4, 
                color = "red", lwd = 0.3)
    }
  }
  
  # Add bovine data points if requested (only for percentile plots)
  if(add_bovine_data) {
    p <- p + geom_point(data = st_drop_geometry(bovine_nigeria_sf), 
                       aes(x = lon, y = lat, size = Number_of_animal_tested, color = Test_Type),
                       alpha = 0.8)  # Increased alpha for better visibility
  }
  
  p <- p +
    scale_fill_viridis_c(
      name = "AAT\nprevalence",
      na.value = "grey90",
      direction = 1,
      option = "viridis",
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 0.1)
    )
  
  # Add uncertainty legend if high uncertainty areas exist
  if(show_uncertainty && "high_uncertainty" %in% names(nigeria_lgas_final) && 
     any(nigeria_lgas_final$high_uncertainty, na.rm = TRUE)) {
    p <- p + 
      labs(caption = "Areas with white overlay and red border indicate high uncertainty (CI width > 0.80)")
  }
  
  # Add additional scales for bovine data if included
  if(add_bovine_data) {
    p <- p +
      # Add size legend for sample size
      scale_size_continuous(name = "Sample size", 
                           range = c(1, 6),  # Increased from c(1, 4) to make sizes more distinguishable
                           breaks = c(10, 50, 100, 200),
                           labels = c("10", "50", "100", "200+")) +
      # Add color scale for test types
      scale_color_manual(name = "Test type",
                        values = c("BCT" = "red", "PCR" = "blue"),
                        guide = guide_legend(override.aes = list(size = 3, alpha = 0.8)))
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
      title = paste("AAT prevalence by LGA in Nigeria", title_suffix)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0.5, size = 9, color = "grey60"),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    )
  
  # Add extra spacing between legends if bovine data is included
  if(add_bovine_data) {
    p <- p + theme(
      legend.spacing.y = unit(1, "cm"),  # Add vertical spacing between legends
      legend.box.spacing = unit(1, "cm")  # Add spacing between legend boxes
    )
  }
  
  # Add inset map if requested and data is available
  if(add_inset && !is.null(africa_countries_sf) && !is.null(nigeria_country_sf)) {
    # Create inset map showing Nigeria's location in Africa
    inset_plot <- ggplot() +
      geom_sf(data = africa_countries_sf, fill = "grey90", color = "white", size = 0.1) +
      geom_sf(data = nigeria_country_sf, fill = "#440154", color = "white", size = 0.3) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 0.5),
        plot.margin = margin(2.1, 2.1, 2.1, 2.1)
      ) +
      coord_sf(expand = FALSE)
    
    # Combine main plot with inset using cowplot
    p <- ggdraw(p) +
      draw_plot(inset_plot, x = 0.70, y = 0.02, width = 0.25, height = 0.25)
  }
  
  return(p)
}

# Create all three maps (only add bovine data to percentile plots, show uncertainty on percentile plots)
p_mean <- create_choropleth(lga_prevalence_mean, "(mean)", add_bovine_data = FALSE, show_uncertainty = TRUE, add_inset = TRUE)
p_lower <- create_choropleth(lga_prevalence_lower, "(2.5th percentile)", add_bovine_data = TRUE, show_uncertainty = TRUE)
p_upper <- create_choropleth(lga_prevalence_upper, "(97.5th percentile)", add_bovine_data = TRUE, show_uncertainty = TRUE)

# Create combined histogram and boxplot for each estimate type
create_combined_plot <- function(lga_data, estimate_name) {
  # Histogram with viridis color scale matching choropleth
  p_hist <- ggplot(lga_data, aes(x = lga_prevalence, fill = after_stat(x))) +
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
      title = paste("Distribution of", estimate_name, "prevalence by LGA in Nigeria"),
      x = ""
    )
  
  # Boxplot with viridis color matching choropleth
  viridis_color <- viridis::viridis(1, begin = 0.5, end = 0.5, option = "viridis")  # Mid-range viridis color
  
  p_box <- ggplot(lga_data, aes(x = lga_prevalence, y = 1)) +  # Use y = 1 instead of y = ""
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

# Save all plots to current working directory
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_choropleth_mean.pdf", plot = p_mean, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_choropleth_lower.pdf", plot = p_lower, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_choropleth_upper.pdf", plot = p_upper, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each estimate type) - more compact dimensions
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_mean_analysis.pdf", plot = p_mean_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_lower_analysis.pdf", plot = p_lower_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_upper_analysis.pdf", plot = p_upper_combined, width = 10, height = 6)

# Display all plots
print(p_mean)
print(p_lower) 
print(p_upper)
grid.draw(p_mean_combined)
grid.draw(p_lower_combined)
grid.draw(p_upper_combined)

# Print summary statistics for all three estimates
cat("=== PREVALENCE ESTIMATES SUMMARY ===\n\n")

cat("Mean Prevalence:\n")
cat("Total LGAs with data:", nrow(lga_prevalence_mean), "\n")
if(nrow(lga_prevalence_mean) > 0) {
  print(summary(lga_prevalence_mean$lga_prevalence))
}

cat("\n2.5th Percentile Prevalence:\n")
cat("Total LGAs with data:", nrow(lga_prevalence_lower), "\n")
if(nrow(lga_prevalence_lower) > 0) {
  print(summary(lga_prevalence_lower$lga_prevalence))
}

cat("\n97.5th Percentile Prevalence:\n")
cat("Total LGAs with data:", nrow(lga_prevalence_upper), "\n")
if(nrow(lga_prevalence_upper) > 0) {
  print(summary(lga_prevalence_upper$lga_prevalence))
}