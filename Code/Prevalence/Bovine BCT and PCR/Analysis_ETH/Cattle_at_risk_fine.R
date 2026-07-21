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

# Get Ethiopian administrative boundaries at zone level (level 3)
# Ethiopia GADM data goes to level 3 (Zones), similar to Nigeria's LGAs
ethiopia_zones <- gadm(country = "ETH", level = 3, path = tempdir())
ethiopia_zones_sf <- st_as_sf(ethiopia_zones)

# Also get region boundaries for context
ethiopia_regions <- gadm(country = "ETH", level = 1, path = tempdir()) 
ethiopia_regions_sf <- st_as_sf(ethiopia_regions)

# Get neighboring countries for geographic context
cat("Loading neighboring countries for geographic context...\n")
neighboring_countries <- c("SDN", "SSD", "KEN", "SOM", "DJI", "ERI")  # Sudan, South Sudan, Kenya, Somalia, Djibouti, Eritrea
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

# Load cattle at risk data from the summary file
cat("Loading cattle at risk data from Ethiopia_cattle_uncertainty_fine.csv...\n")
cattle_summary_data <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Ethiopia_cattle_uncertainty_fine.csv")

cat("Loaded cattle at risk data with", nrow(cattle_summary_data), "regions\n")
cat("Columns available:", paste(names(cattle_summary_data), collapse = ", "), "\n")

# MANUAL EXCLUSION: Remove problem regions that have cattle at risk but no tsetse presence
# These regions are biologically impossible and likely due to spatial processing artifacts
problem_regions_to_exclude <- c("Becho", "Jarso", "Kersa", "Wara Jarso", "Limu")
cat("Manually excluding", length(problem_regions_to_exclude), "problem regions without tsetse presence...\n")

# Zero out cattle values for these problem regions
cattle_summary_data$mean_cattle_mean[cattle_summary_data$region %in% problem_regions_to_exclude] <- 0
cattle_summary_data$data_uncertainty_lower[cattle_summary_data$region %in% problem_regions_to_exclude] <- 0
cattle_summary_data$data_uncertainty_upper[cattle_summary_data$region %in% problem_regions_to_exclude] <- 0

excluded_regions <- cattle_summary_data$region[cattle_summary_data$region %in% problem_regions_to_exclude]
# Load tsetse distribution data to properly classify tsetse zones
cat("Loading tsetse distribution data for proper zone classification...\n")
library(terra)
library(raster)
tsetse_raster <- raster("Data/Covariates/tsenumbspec")
# Convert values > 1 to 1 (binary presence/absence) as done in tsetse plot script
tsetse_raster[tsetse_raster > 1] <- 1

# Use area-based extraction instead of centroid-based to capture zones that partially overlap with tsetse areas
cat("Extracting tsetse presence for each zone using area-based method...\n")

# Extract tsetse values for each zone using area statistics
# This will give us the mean tsetse value within each zone
zone_tsetse_stats <- raster::extract(tsetse_raster, ethiopia_zones_sf, fun = mean, na.rm = TRUE, df = TRUE)

# Create a data frame linking zone names to tsetse presence
zone_tsetse_data <- data.frame(
  NAME_3 = ethiopia_zones_sf$NAME_3,
  NAME_2 = ethiopia_zones_sf$NAME_2,  # Include region names for geographical rules
  NAME_1 = ethiopia_zones_sf$NAME_1,  # Include state names for geographical rules
  tsetse_mean = zone_tsetse_stats[,2]  # Second column contains the mean values
)

# Replace NA values with 0 (no tsetse) and create binary tsetse zone indicator
zone_tsetse_data$tsetse_mean[is.na(zone_tsetse_data$tsetse_mean)] <- 0

# SIMPLIFIED: Use only biological tsetse presence, no geographical exclusions
# A zone is considered tsetse zone if it has any meaningful tsetse presence
zone_tsetse_data$tsetse_zone <- zone_tsetse_data$tsetse_mean > 0

cat("Found", sum(zone_tsetse_data$tsetse_zone, na.rm = TRUE), "zones with tsetse presence out of", nrow(zone_tsetse_data), "total zones\n")
cat("Using simple tsetse presence logic: no geographical exclusions\n")
cat("Zones with any tsetse (>0):", sum(zone_tsetse_data$tsetse_mean > 0, na.rm = TRUE), "\n")
cat("Zones without tsetse (=0):", sum(zone_tsetse_data$tsetse_mean == 0, na.rm = TRUE), "\n")

# Create zone-level data for plotting - use the pre-masked cattle data from Prevalence_plots_optimised_fine.R
# The cattle data is already correctly masked (cattle at risk = 0 where no tsetse)
zone_cattle_at_risk_mean <- cattle_summary_data %>%
  dplyr::select(NAME_3 = region, 
                mean_cattle_at_risk = mean_cattle_mean,
                risk_value = mean_value_mean) %>%
  right_join(zone_tsetse_data, by = "NAME_3") %>%  # Use right_join to keep all zones
  dplyr::mutate(
    # If no cattle data, assume 0 cattle in that zone
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    log_cattle_at_risk = log10(mean_cattle_at_risk + 1),  # Add log transformation
    # Simplified risk categories: trust the pre-masked data from Prevalence_plots_optimised_fine.R
    risk_category = case_when(
      mean_cattle_at_risk > 0 ~ "cattle_at_risk",      # Areas with cattle at risk (already tsetse-masked)
      tsetse_zone == TRUE ~ "no_cattle_but_risk",      # Tsetse zones with no cattle at risk
      TRUE ~ "no_risk_or_cattle"                       # Areas without tsetse (no risk possible)
    )
  )

zone_cattle_at_risk_lower <- cattle_summary_data %>%
  dplyr::select(NAME_3 = region, 
                lower_cattle_at_risk = data_uncertainty_lower,
                mean_cattle_at_risk = mean_cattle_mean,
                risk_value = mean_value_mean) %>%
  right_join(zone_tsetse_data, by = "NAME_3") %>%  # Use right_join to keep all zones
  dplyr::mutate(
    # If no cattle data, assume 0 cattle in that zone
    lower_cattle_at_risk = ifelse(is.na(lower_cattle_at_risk), 0, lower_cattle_at_risk),
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    log_cattle_at_risk = log10(lower_cattle_at_risk + 1),
    # Simplified: trust pre-masked data
    risk_category = case_when(
      lower_cattle_at_risk > 0 ~ "cattle_at_risk",
      tsetse_zone == TRUE ~ "no_cattle_but_risk",
      TRUE ~ "no_risk_or_cattle"
    )
  )

zone_cattle_at_risk_upper <- cattle_summary_data %>%
  dplyr::select(NAME_3 = region, 
                upper_cattle_at_risk = data_uncertainty_upper,
                mean_cattle_at_risk = mean_cattle_mean,
                risk_value = mean_value_mean) %>%
  right_join(zone_tsetse_data, by = "NAME_3") %>%  # Use right_join to keep all zones
  dplyr::mutate(
    # If no cattle data, assume 0 cattle in that zone
    upper_cattle_at_risk = ifelse(is.na(upper_cattle_at_risk), 0, upper_cattle_at_risk),
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    log_cattle_at_risk = log10(upper_cattle_at_risk + 1),
    # Simplified: trust pre-masked data
    risk_category = case_when(
      upper_cattle_at_risk > 0 ~ "cattle_at_risk",
      tsetse_zone == TRUE ~ "no_cattle_but_risk",
      TRUE ~ "no_risk_or_cattle"
    )
  )

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

# Function to create choropleth map for cattle at risk with risk categories
create_choropleth <- function(zone_data, title_suffix, data_type = "mean", add_bovine = FALSE, use_log = FALSE, add_inset = FALSE) {
  
  # Handle the different data types
  if(data_type == "mean") {
    fill_column <- if(use_log) "log_cattle_at_risk" else "mean_cattle_at_risk"
  } else if(data_type == "lower") {
    fill_column <- if(use_log) "log_cattle_at_risk" else "lower_cattle_at_risk"
  } else if(data_type == "upper") {
    fill_column <- if(use_log) "log_cattle_at_risk" else "upper_cattle_at_risk"
  } else {
    fill_column <- data_type
  }
  
  # Join zone-level cattle at risk back to zone polygons
  ethiopia_zones_final <- ethiopia_zones_sf %>%
    left_join(zone_data %>% st_drop_geometry() %>% dplyr::select(NAME_3, all_of(fill_column), risk_category), 
              by = "NAME_3")
  
  # Create zone-level choropleth map
  p <- ggplot(ethiopia_zones_final)
  
  # Add neighboring countries as background context (if available)
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Separate areas with cattle at risk from zero areas
  zones_with_cattle <- ethiopia_zones_final %>% dplyr::filter(risk_category == "cattle_at_risk" | is.na(risk_category))
  zones_no_cattle_but_risk <- ethiopia_zones_final %>% dplyr::filter(risk_category == "no_cattle_but_risk")
  zones_no_risk <- ethiopia_zones_final %>% dplyr::filter(risk_category == "no_risk_or_cattle")
  
  # Add areas with no cattle but risk (tsetse zones without cattle) in light blue
  if(nrow(zones_no_cattle_but_risk) > 0) {
    p <- p + geom_sf(data = zones_no_cattle_but_risk, fill = "lightblue", color = "black", lwd = 0.05)
  }
  
  # Add areas with no risk/cattle in light grey  
  if(nrow(zones_no_risk) > 0) {
    p <- p + geom_sf(data = zones_no_risk, fill = "lightgrey", color = "black", lwd = 0.05)
  }
  
  # Add Ethiopia zones with cattle at risk data (main choropleth)
  if(nrow(zones_with_cattle) > 0) {
    p <- p + geom_sf(data = zones_with_cattle, aes(fill = !!sym(fill_column)), lwd = 0.05, color = "black")
  }
  
  # Add region boundaries
  p <- p + geom_sf(data = ethiopia_regions_sf, fill = NA, color = "black", lwd = 0.5)
  
  # Calculate appropriate limits for cattle at risk - STANDARDIZED across countries
  max_cattle_at_risk <- max(zones_with_cattle[[fill_column]], na.rm = TRUE)
  
  # Set standardized limits for comparison across Ethiopia and Nigeria
  # These values should be adjusted based on the data range across both countries
  standardized_limits <- if(use_log) {
    c(0, 5.6)  # For log scale - adjust based on log10 values
  } else {
    c(0, 500000)  # For regular scale - 500k cattle max, adjust as needed
  }
  
  cat("Ethiopia", data_type, "- Using standardized limits:", paste(standardized_limits, collapse = " to "), 
      "(country max:", round(max_cattle_at_risk), ")\n")
  
  # Create appropriate title and legend based on log transformation
  legend_name <- if(use_log) expression(atop("Infected cattle", "("*log[10] * ")")) else "Cattle\nat risk"
  
  p <- p + scale_fill_viridis_c(
      name = legend_name,
      na.value = "grey90",
      direction = 1,
      option = "plasma",  # Different color palette for cattle at risk
      limits = standardized_limits,  # Use standardized limits for comparison
      labels = if(use_log) {
        function(x) format(x, digits = 2, scientific = FALSE)  # Just show log values as-is
      } else {
        function(x) format(x, big.mark = ",", scientific = FALSE)
      }
    ) +
    # Add scale bar and north arrow
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                    bar_cols = c("black", "white"), line_width = 1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                          pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                          style = north_arrow_fancy_orienteering(text_size = 8)) +
    theme_void()
  
  p <- p +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    )
  
  return(p)
}

# Create all choropleth maps - regular and log versions
# Regular scale maps
p_mean <- create_choropleth(zone_cattle_at_risk_mean, "(mean)", add_inset = TRUE)
p_lower <- create_choropleth(zone_cattle_at_risk_lower, "(2.5th percentile)", data_type = "lower", add_bovine = TRUE)
p_upper <- create_choropleth(zone_cattle_at_risk_upper, "(97.5th percentile)", data_type = "upper", add_bovine = TRUE)

# Log scale maps
p_mean_log <- create_choropleth(zone_cattle_at_risk_mean, "(mean - log scale)", use_log = TRUE, add_inset = TRUE)
p_lower_log <- create_choropleth(zone_cattle_at_risk_lower, "(2.5th percentile - log scale)", data_type = "lower", use_log = TRUE, add_bovine = TRUE)
p_upper_log <- create_choropleth(zone_cattle_at_risk_upper, "(97.5th percentile - log scale)", data_type = "upper", use_log = TRUE, add_bovine = TRUE)

# Create combined histogram and boxplot for each estimate type
create_combined_plot <- function(zone_data, estimate_name, use_log = FALSE) {
  # Get zone data without geometry for plotting
  value_column <- if(use_log) "log_cattle_at_risk" else "mean_cattle_at_risk"
  plot_data <- zone_data %>% st_drop_geometry() %>% dplyr::filter(!is.na(!!sym(value_column)))
  
  # Histogram with plasma color scale matching choropleth
  p_hist <- ggplot(plot_data, aes(x = !!sym(value_column), fill = after_stat(x))) +
    geom_histogram(bins = 30, color = "white", linewidth = 0.2,
                   boundary = 0, closed = "left") +
    scale_fill_viridis_c(option = "plasma", guide = "none") +  # Same as choropleth, no legend
    scale_x_continuous(labels = if(use_log) {
      function(x) format(x, digits = 2, scientific = FALSE)  # Just show log values as-is
    } else {
      function(x) format(x, big.mark = ",", scientific = FALSE)
    }, expand = expansion(mult = c(0.02, 0.02))) +
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
      title = paste("Distribution of", estimate_name, "cattle at risk by district in Ethiopia"),
      x = ""
    )
  
  # Boxplot with plasma color matching choropleth
  plasma_color <- viridis::plasma(1, begin = 0.5, end = 0.5)  # Mid-range plasma color
  
  p_box <- ggplot(plot_data, aes(x = !!sym(value_column), y = 1)) +  # Use y = 1 instead of y = ""
    geom_boxplot(fill = plasma_color, alpha = 0.8, width = 0.8) +  # Much wider boxplot, no jitter
    scale_x_continuous(labels = if(use_log) {
      function(x) format(x, digits = 2, scientific = FALSE)  # Just show log values as-is
    } else {
      function(x) format(x, big.mark = ",", scientific = FALSE)
    }, expand = expansion(mult = c(0.02, 0.02))) +
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
  labs(x = expression("Infected cattle (" * log[10] * ")"))
  
  # Use a simpler approach with grid.arrange for better stability
  library(gridExtra)
  library(grid)
  
  # Create title
  title <- grid::textGrob(paste("Cattle at Risk Analysis -", estimate_name), 
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

# Create the six combined plots (3 regular + 3 log)
p_mean_combined <- create_combined_plot(zone_cattle_at_risk_mean, "mean")
p_lower_combined <- create_combined_plot(zone_cattle_at_risk_lower, "2.5th percentile") 
p_upper_combined <- create_combined_plot(zone_cattle_at_risk_upper, "97.5th percentile")

p_mean_log_combined <- create_combined_plot(zone_cattle_at_risk_mean, "mean (log scale)", use_log = TRUE)
p_lower_log_combined <- create_combined_plot(zone_cattle_at_risk_lower, "2.5th percentile (log scale)", use_log = TRUE) 
p_upper_log_combined <- create_combined_plot(zone_cattle_at_risk_upper, "97.5th percentile (log scale)", use_log = TRUE)

# Save all plots to Analysis_ETH directory
# Regular scale plots
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_choropleth_mean.pdf", plot = p_mean, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_choropleth_lower.pdf", plot = p_lower, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_choropleth_upper.pdf", plot = p_upper, width = 12, height = 10)

# Log scale plots
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_choropleth_mean_log.pdf", plot = p_mean_log, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_choropleth_lower_log.pdf", plot = p_lower_log, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_choropleth_upper_log.pdf", plot = p_upper_log, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each estimate type) - more compact dimensions
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_mean_analysis.pdf", plot = p_mean_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_lower_analysis.pdf", plot = p_lower_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_upper_analysis.pdf", plot = p_upper_combined, width = 10, height = 6)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_mean_log_analysis.pdf", plot = p_mean_log_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_lower_log_analysis.pdf", plot = p_lower_log_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_plots/eth_cattle_at_risk_upper_log_analysis.pdf", plot = p_upper_log_combined, width = 10, height = 6)






# ============================================================
# FIGURE 5B SOURCE DATA
#
# Exports the data required to reconstruct:
#   1. Mean cattle-at-risk choropleth
#   2. Lower 2.5th-percentile choropleth
#   3. Upper 97.5th-percentile choropleth
#   4. Log10-transformed choropleths
#   5. Histograms and boxplots
#   6. Survey-point overlays
#   7. Tsetse-risk categories
#   8. Manual exclusions
#
# This block deliberately follows the existing plotting code.
# ============================================================

library(openxlsx)

# Always use explicit namespace because raster::select can mask
# dplyr::select.

figure5b_output_dir <- paste0(
  "Code/Prevalence/Bovine BCT and PCR/",
  "Analysis_ETH/Cattle_at_risk_fine_plots"
)

dir.create(
  figure5b_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Administrative polygon lookup
# ------------------------------------------------------------

figure5b_polygon_lookup <- ethiopia_zones_sf |>
  sf::st_drop_geometry() |>
  dplyr::transmute(
    GID_3 = as.character(GID_3),
    Region = as.character(NAME_1),
    Zone = as.character(NAME_2),
    Administrative_area = as.character(NAME_3)
  )

stopifnot(
  nrow(figure5b_polygon_lookup) ==
    dplyr::n_distinct(figure5b_polygon_lookup$GID_3)
)

figure5b_repeated_names <- figure5b_polygon_lookup |>
  dplyr::count(
    Administrative_area,
    name = "Number_of_GADM_polygons"
  ) |>
  dplyr::filter(
    Number_of_GADM_polygons > 1
  ) |>
  dplyr::arrange(
    dplyr::desc(Number_of_GADM_polygons),
    Administrative_area
  )

# ------------------------------------------------------------
# 2. Exact source objects used by the plotting functions
# ------------------------------------------------------------

figure5b_mean_source <- zone_cattle_at_risk_mean |>
  sf::st_drop_geometry() |>
  dplyr::transmute(
    Source_row_ID = dplyr::row_number(),
    Administrative_area = as.character(NAME_3),
    Region_from_source = as.character(NAME_1),
    Zone_from_source = as.character(NAME_2),
    Mean_cattle_at_risk =
      as.numeric(mean_cattle_at_risk),
    Mean_log10_cattle_at_risk =
      as.numeric(log_cattle_at_risk),
    Risk_value =
      as.numeric(risk_value),
    Tsetse_mean =
      as.numeric(tsetse_mean),
    Tsetse_zone =
      as.logical(tsetse_zone),
    Mean_risk_category =
      as.character(risk_category)
  )

figure5b_lower_source <- zone_cattle_at_risk_lower |>
  sf::st_drop_geometry() |>
  dplyr::transmute(
    Source_row_ID = dplyr::row_number(),
    Administrative_area = as.character(NAME_3),
    Region_from_source = as.character(NAME_1),
    Zone_from_source = as.character(NAME_2),
    Lower_2.5_cattle_at_risk =
      as.numeric(lower_cattle_at_risk),
    Mean_cattle_at_risk_in_lower_object =
      as.numeric(mean_cattle_at_risk),
    Lower_log10_cattle_at_risk =
      as.numeric(log_cattle_at_risk),
    Risk_value =
      as.numeric(risk_value),
    Tsetse_mean =
      as.numeric(tsetse_mean),
    Tsetse_zone =
      as.logical(tsetse_zone),
    Lower_risk_category =
      as.character(risk_category)
  )

figure5b_upper_source <- zone_cattle_at_risk_upper |>
  sf::st_drop_geometry() |>
  dplyr::transmute(
    Source_row_ID = dplyr::row_number(),
    Administrative_area = as.character(NAME_3),
    Region_from_source = as.character(NAME_1),
    Zone_from_source = as.character(NAME_2),
    Upper_97.5_cattle_at_risk =
      as.numeric(upper_cattle_at_risk),
    Mean_cattle_at_risk_in_upper_object =
      as.numeric(mean_cattle_at_risk),
    Upper_log10_cattle_at_risk =
      as.numeric(log_cattle_at_risk),
    Risk_value =
      as.numeric(risk_value),
    Tsetse_mean =
      as.numeric(tsetse_mean),
    Tsetse_zone =
      as.logical(tsetse_zone),
    Upper_risk_category =
      as.character(risk_category)
  )

# ------------------------------------------------------------
# 3. Recreate the exact joins used by the maps
# ------------------------------------------------------------

# These joins intentionally allow many-to-many matching because
# that is how the existing map code behaves when NAME_3 occurs
# more than once.

figure5b_mean_map <- figure5b_polygon_lookup |>
  dplyr::left_join(
    figure5b_mean_source,
    by = "Administrative_area",
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(
    Map_row_ID = dplyr::row_number()
  ) |>
  dplyr::select(
    Map_row_ID,
    GID_3,
    Region,
    Zone,
    Administrative_area,
    Source_row_ID,
    Region_from_source,
    Zone_from_source,
    Mean_cattle_at_risk,
    Mean_log10_cattle_at_risk,
    Risk_value,
    Tsetse_mean,
    Tsetse_zone,
    Mean_risk_category
  )

figure5b_lower_map <- figure5b_polygon_lookup |>
  dplyr::left_join(
    figure5b_lower_source,
    by = "Administrative_area",
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(
    Map_row_ID = dplyr::row_number()
  ) |>
  dplyr::select(
    Map_row_ID,
    GID_3,
    Region,
    Zone,
    Administrative_area,
    Source_row_ID,
    Region_from_source,
    Zone_from_source,
    Lower_2.5_cattle_at_risk,
    Mean_cattle_at_risk_in_lower_object,
    Lower_log10_cattle_at_risk,
    Risk_value,
    Tsetse_mean,
    Tsetse_zone,
    Lower_risk_category
  )

figure5b_upper_map <- figure5b_polygon_lookup |>
  dplyr::left_join(
    figure5b_upper_source,
    by = "Administrative_area",
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(
    Map_row_ID = dplyr::row_number()
  ) |>
  dplyr::select(
    Map_row_ID,
    GID_3,
    Region,
    Zone,
    Administrative_area,
    Source_row_ID,
    Region_from_source,
    Zone_from_source,
    Upper_97.5_cattle_at_risk,
    Mean_cattle_at_risk_in_upper_object,
    Upper_log10_cattle_at_risk,
    Risk_value,
    Tsetse_mean,
    Tsetse_zone,
    Upper_risk_category
  )

# ------------------------------------------------------------
# 4. Map values in long format
# ------------------------------------------------------------

figure5b_map_values <- dplyr::bind_rows(

  figure5b_mean_map |>
    dplyr::transmute(
      Map_row_ID,
      GID_3,
      Region,
      Zone,
      Administrative_area,
      Source_row_ID,
      Estimate = "Mean",
      Scale = "Raw",
      Plotted_value = Mean_cattle_at_risk,
      Risk_category = Mean_risk_category,
      Tsetse_mean,
      Tsetse_zone
    ),

  figure5b_mean_map |>
    dplyr::transmute(
      Map_row_ID,
      GID_3,
      Region,
      Zone,
      Administrative_area,
      Source_row_ID,
      Estimate = "Mean",
      Scale = "Log10",
      Plotted_value = Mean_log10_cattle_at_risk,
      Risk_category = Mean_risk_category,
      Tsetse_mean,
      Tsetse_zone
    ),

  figure5b_lower_map |>
    dplyr::transmute(
      Map_row_ID,
      GID_3,
      Region,
      Zone,
      Administrative_area,
      Source_row_ID,
      Estimate = "Lower 2.5th percentile",
      Scale = "Raw",
      Plotted_value = Lower_2.5_cattle_at_risk,
      Risk_category = Lower_risk_category,
      Tsetse_mean,
      Tsetse_zone
    ),

  figure5b_lower_map |>
    dplyr::transmute(
      Map_row_ID,
      GID_3,
      Region,
      Zone,
      Administrative_area,
      Source_row_ID,
      Estimate = "Lower 2.5th percentile",
      Scale = "Log10",
      Plotted_value = Lower_log10_cattle_at_risk,
      Risk_category = Lower_risk_category,
      Tsetse_mean,
      Tsetse_zone
    ),

  figure5b_upper_map |>
    dplyr::transmute(
      Map_row_ID,
      GID_3,
      Region,
      Zone,
      Administrative_area,
      Source_row_ID,
      Estimate = "Upper 97.5th percentile",
      Scale = "Raw",
      Plotted_value = Upper_97.5_cattle_at_risk,
      Risk_category = Upper_risk_category,
      Tsetse_mean,
      Tsetse_zone
    ),

  figure5b_upper_map |>
    dplyr::transmute(
      Map_row_ID,
      GID_3,
      Region,
      Zone,
      Administrative_area,
      Source_row_ID,
      Estimate = "Upper 97.5th percentile",
      Scale = "Log10",
      Plotted_value = Upper_log10_cattle_at_risk,
      Risk_category = Upper_risk_category,
      Tsetse_mean,
      Tsetse_zone
    )

) |>
  dplyr::mutate(
    Manually_excluded =
      Administrative_area %in%
        problem_regions_to_exclude
  ) |>
  dplyr::arrange(
    Scale,
    Estimate,
    Region,
    Zone,
    Administrative_area,
    GID_3,
    Map_row_ID
  )

# ------------------------------------------------------------
# 5. Exact values used by histograms and boxplots
# ------------------------------------------------------------

# This follows create_combined_plot() exactly:
# raw plots use mean_cattle_at_risk for all three supplied
# objects; log plots use each object's log_cattle_at_risk.

figure5b_distribution_values <- dplyr::bind_rows(

  figure5b_mean_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Plot_label = "Mean",
      Scale = "Raw",
      Source_object = "zone_cattle_at_risk_mean",
      Source_column = "mean_cattle_at_risk",
      Plotted_value = Mean_cattle_at_risk
    ),

  figure5b_lower_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Plot_label = "Lower 2.5th percentile",
      Scale = "Raw",
      Source_object = "zone_cattle_at_risk_lower",
      Source_column = "mean_cattle_at_risk",
      Plotted_value =
        Mean_cattle_at_risk_in_lower_object
    ),

  figure5b_upper_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Plot_label = "Upper 97.5th percentile",
      Scale = "Raw",
      Source_object = "zone_cattle_at_risk_upper",
      Source_column = "mean_cattle_at_risk",
      Plotted_value =
        Mean_cattle_at_risk_in_upper_object
    ),

  figure5b_mean_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Plot_label = "Mean",
      Scale = "Log10",
      Source_object = "zone_cattle_at_risk_mean",
      Source_column = "log_cattle_at_risk",
      Plotted_value =
        Mean_log10_cattle_at_risk
    ),

  figure5b_lower_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Plot_label = "Lower 2.5th percentile",
      Scale = "Log10",
      Source_object = "zone_cattle_at_risk_lower",
      Source_column = "log_cattle_at_risk",
      Plotted_value =
        Lower_log10_cattle_at_risk
    ),

  figure5b_upper_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Plot_label = "Upper 97.5th percentile",
      Scale = "Log10",
      Source_object = "zone_cattle_at_risk_upper",
      Source_column = "log_cattle_at_risk",
      Plotted_value =
        Upper_log10_cattle_at_risk
    )

) |>
  dplyr::filter(
    !is.na(Plotted_value)
  ) |>
  dplyr::arrange(
    Scale,
    Plot_label,
    Source_row_ID
  )

# ------------------------------------------------------------
# 6. Estimate-specific values
# ------------------------------------------------------------

# This sheet contains the actual mean, lower and upper values,
# irrespective of the current histogram-function column choice.

figure5b_estimate_values <- dplyr::bind_rows(

  figure5b_mean_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Estimate = "Mean",
      Raw_value = Mean_cattle_at_risk,
      Log10_value = Mean_log10_cattle_at_risk
    ),

  figure5b_lower_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Estimate = "Lower 2.5th percentile",
      Raw_value = Lower_2.5_cattle_at_risk,
      Log10_value = Lower_log10_cattle_at_risk
    ),

  figure5b_upper_source |>
    dplyr::transmute(
      Source_row_ID,
      Region = Region_from_source,
      Zone = Zone_from_source,
      Administrative_area,
      Estimate = "Upper 97.5th percentile",
      Raw_value = Upper_97.5_cattle_at_risk,
      Log10_value = Upper_log10_cattle_at_risk
    )

) |>
  dplyr::arrange(
    Estimate,
    Region,
    Zone,
    Administrative_area,
    Source_row_ID
  )

# ------------------------------------------------------------
# 7. Survey points
# ------------------------------------------------------------

if (
  exists("bovine_ethiopia_sf") &&
  !is.null(bovine_ethiopia_sf) &&
  nrow(bovine_ethiopia_sf) > 0
) {

  figure5b_survey_points <- bovine_ethiopia_sf |>
    sf::st_drop_geometry() |>
    dplyr::transmute(
      Survey_point_ID = dplyr::row_number(),
      Longitude = as.numeric(lon),
      Latitude = as.numeric(lat),
      Sample_size =
        as.integer(Number_of_animal_tested),
      Number_infected =
        as.integer(Number_of_infections),
      Observed_prevalence =
        as.numeric(Prevalence),
      Test_type = dplyr::case_when(
        Test_Type == "BCT/HCT" ~ "HCT/BCT",
        Test_Type == "BCT" ~ "HCT/BCT",
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

  figure5b_survey_points <- data.frame(
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
# 8. Manual exclusions
# ------------------------------------------------------------

figure5b_manual_exclusions <- data.frame(
  Administrative_area =
    problem_regions_to_exclude,

  Applied_action =
    paste(
      "Mean, lower and upper cattle-at-risk values",
      "set to zero"
    ),

  Matching_method =
    "Matched using cattle_summary_data$region",

  Reason = paste(
    "Estimated cattle at risk occurred without tsetse",
    "presence and was treated as a spatial-processing",
    "artefact."
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 9. Plot settings
# ------------------------------------------------------------

figure5b_plot_settings <- data.frame(
  Setting = c(
    "Administrative polygons",
    "Polygon identifier",
    "Map join field",
    "Map relationship",
    "Raw map colour limits",
    "Log map colour limits",
    "Log transformation",
    "Colour palette",
    "Histogram bins",
    "Histogram boundary",
    "Histogram closure",
    "Boxplot width",
    "No cattle but tsetse risk",
    "No risk or cattle",
    "Tsetse classification",
    "Survey CRS"
  ),

  Value = c(
    "Ethiopia GADM level 3",
    "GID_3",
    "NAME_3",
    "Many-to-many retained to reproduce plotting code",
    "0 to 500000",
    "0 to 5.6",
    "log10(cattle at risk + 1)",
    "viridis plasma",
    "30",
    "0",
    "Left-closed",
    "0.8",
    "Light blue",
    "Light grey",
    "tsetse_mean > 0",
    "WGS84, EPSG:4326"
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10. README
# ------------------------------------------------------------

figure5b_readme <- data.frame(
  Item = c(
    "Workbook description",
    "Associated figure",
    "Map_values sheet",
    "Distribution_values sheet",
    "Estimate_values sheet",
    "Survey_points sheet",
    "Repeated_names sheet",
    "Manual_exclusions sheet",
    "Plot_settings sheet",
    "Administrative geometry",
    "Repeated administrative names",
    "Units"
  ),

  Description = c(
    paste(
      "Source data underlying Figure 5B, including",
      "Ethiopian cattle-at-risk choropleths,",
      "histograms, boxplots and survey overlays."
    ),

    "Figure 5B.",

    paste(
      "Exact rows produced by joining GADM level-3",
      "polygons to the mean, lower and upper cattle-risk",
      "objects using NAME_3."
    ),

    paste(
      "Exact numerical values passed to the current",
      "histogram and boxplot function."
    ),

    paste(
      "Actual mean, lower and upper cattle-at-risk",
      "estimates in raw and log10 form."
    ),

    "Survey points retained within Ethiopia.",

    paste(
      "Administrative names occurring in more than one",
      "GADM level-3 polygon."
    ),

    "Administrative names manually assigned zero risk.",

    paste(
      "Binning, transformations, colour scales and",
      "other reconstruction settings."
    ),

    paste(
      "Polygon geometry is obtained separately from GADM",
      "and is not embedded in the workbook."
    ),

    paste(
      "Repeated names and corresponding many-to-many",
      "map rows are retained because they are present in",
      "the existing plotting workflow."
    ),

    "Cattle-at-risk values are estimated numbers of cattle."
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 11. Create workbook
# ------------------------------------------------------------

wb <- openxlsx::createWorkbook(
  creator = "AR Kaye",
  title = "Figure 5B source data",
  subject = paste(
    "Ethiopia cattle-at-risk maps, histograms and boxplots"
  ),
  category = "Research data"
)

openxlsx::addWorksheet(wb, "README")
openxlsx::addWorksheet(wb, "Map_values")
openxlsx::addWorksheet(wb, "Distribution_values")
openxlsx::addWorksheet(wb, "Estimate_values")
openxlsx::addWorksheet(wb, "Survey_points")
openxlsx::addWorksheet(wb, "Repeated_names")
openxlsx::addWorksheet(wb, "Manual_exclusions")
openxlsx::addWorksheet(wb, "Plot_settings")

openxlsx::writeData(
  wb,
  "README",
  figure5b_readme
)

openxlsx::writeData(
  wb,
  "Map_values",
  figure5b_map_values
)

openxlsx::writeData(
  wb,
  "Distribution_values",
  figure5b_distribution_values
)

openxlsx::writeData(
  wb,
  "Estimate_values",
  figure5b_estimate_values
)

openxlsx::writeData(
  wb,
  "Survey_points",
  figure5b_survey_points
)

openxlsx::writeData(
  wb,
  "Repeated_names",
  figure5b_repeated_names
)

openxlsx::writeData(
  wb,
  "Manual_exclusions",
  figure5b_manual_exclusions
)

openxlsx::writeData(
  wb,
  "Plot_settings",
  figure5b_plot_settings
)

# ------------------------------------------------------------
# 12. Styling
# ------------------------------------------------------------

header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "bottom",
  borderColour = "#FFFFFF"
)

count_style <- openxlsx::createStyle(
  numFmt = "#,##0.000"
)

decimal_style <- openxlsx::createStyle(
  numFmt = "0.000000"
)

integer_style <- openxlsx::createStyle(
  numFmt = "0"
)

wrap_style <- openxlsx::createStyle(
  wrapText = TRUE,
  valign = "top"
)

format_data_sheet <- function(
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
    style = header_style,
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

format_data_sheet(
  wb,
  "Map_values",
  figure5b_map_values
)

format_data_sheet(
  wb,
  "Distribution_values",
  figure5b_distribution_values
)

format_data_sheet(
  wb,
  "Estimate_values",
  figure5b_estimate_values
)

format_data_sheet(
  wb,
  "Survey_points",
  figure5b_survey_points
)

format_data_sheet(
  wb,
  "Repeated_names",
  figure5b_repeated_names
)

format_data_sheet(
  wb,
  "Manual_exclusions",
  figure5b_manual_exclusions
)

# README
openxlsx::addStyle(
  wb,
  "README",
  header_style,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::addStyle(
  wb,
  "README",
  wrap_style,
  rows = 2:(nrow(figure5b_readme) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::setColWidths(
  wb,
  "README",
  cols = 1,
  widths = 30
)

openxlsx::setColWidths(
  wb,
  "README",
  cols = 2,
  widths = 90
)

openxlsx::freezePane(
  wb,
  "README",
  firstRow = TRUE
)

# Plot settings
openxlsx::addStyle(
  wb,
  "Plot_settings",
  header_style,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::addStyle(
  wb,
  "Plot_settings",
  wrap_style,
  rows = 2:(nrow(figure5b_plot_settings) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::setColWidths(
  wb,
  "Plot_settings",
  cols = 1,
  widths = 36
)

openxlsx::setColWidths(
  wb,
  "Plot_settings",
  cols = 2,
  widths = 85
)

openxlsx::freezePane(
  wb,
  "Plot_settings",
  firstRow = TRUE
)

# ------------------------------------------------------------
# 13. Numeric formatting
# ------------------------------------------------------------

map_value_column <- match(
  "Plotted_value",
  names(figure5b_map_values)
)

if (
  nrow(figure5b_map_values) > 0 &&
  !is.na(map_value_column)
) {
  openxlsx::addStyle(
    wb,
    "Map_values",
    count_style,
    rows = 2:(nrow(figure5b_map_values) + 1),
    cols = map_value_column,
    gridExpand = TRUE
  )
}

distribution_value_column <- match(
  "Plotted_value",
  names(figure5b_distribution_values)
)

if (
  nrow(figure5b_distribution_values) > 0 &&
  !is.na(distribution_value_column)
) {
  openxlsx::addStyle(
    wb,
    "Distribution_values",
    count_style,
    rows = 2:(nrow(figure5b_distribution_values) + 1),
    cols = distribution_value_column,
    gridExpand = TRUE
  )
}

estimate_numeric_columns <- match(
  c("Raw_value", "Log10_value"),
  names(figure5b_estimate_values)
)

estimate_numeric_columns <- estimate_numeric_columns[
  !is.na(estimate_numeric_columns)
]

if (
  nrow(figure5b_estimate_values) > 0 &&
  length(estimate_numeric_columns) > 0
) {
  openxlsx::addStyle(
    wb,
    "Estimate_values",
    count_style,
    rows = 2:(nrow(figure5b_estimate_values) + 1),
    cols = estimate_numeric_columns,
    gridExpand = TRUE
  )
}

if (nrow(figure5b_survey_points) > 0) {

  openxlsx::addStyle(
    wb,
    "Survey_points",
    decimal_style,
    rows = 2:(nrow(figure5b_survey_points) + 1),
    cols = c(2, 3, 6),
    gridExpand = TRUE
  )

  openxlsx::addStyle(
    wb,
    "Survey_points",
    integer_style,
    rows = 2:(nrow(figure5b_survey_points) + 1),
    cols = c(1, 4, 5),
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 14. Save workbook
# ------------------------------------------------------------

figure5b_output_file <- file.path(
  figure5b_output_dir,
  "Figure_5B_source_data.xlsx"
)

openxlsx::saveWorkbook(
  wb,
  file = figure5b_output_file,
  overwrite = TRUE
)

message(
  "Figure 5B source-data workbook saved to: ",
  figure5b_output_file
)

message(
  "Map values exported: ",
  nrow(figure5b_map_values)
)

message(
  "Distribution values exported: ",
  nrow(figure5b_distribution_values)
)

message(
  "Estimate-specific values exported: ",
  nrow(figure5b_estimate_values)
)

message(
  "Survey points exported: ",
  nrow(figure5b_survey_points)
)

message(
  "Repeated administrative names documented: ",
  nrow(figure5b_repeated_names)
)