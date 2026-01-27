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
cattle_summary_data <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Ethiopia_cattle_uncertainty_fine_whole_country.csv")

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
    p <- p + geom_sf(data = zones_no_cattle_but_risk, fill = "lightblue", color = "white", lwd = 0.1)
  }
  
  # Add areas with no risk/cattle in light grey  
  if(nrow(zones_no_risk) > 0) {
    p <- p + geom_sf(data = zones_no_risk, fill = "lightgrey", color = "white", lwd = 0.1)
  }
  
  # Add Ethiopia zones with cattle at risk data (main choropleth)
  if(nrow(zones_with_cattle) > 0) {
    p <- p + geom_sf(data = zones_with_cattle, aes(fill = !!sym(fill_column)), lwd = 0.1, color = "white")
  }
  
  # Add region boundaries
  p <- p + geom_sf(data = ethiopia_regions_sf, fill = NA, color = "black", lwd = 0.5)
  
  # Calculate appropriate limits for cattle at risk
  max_cattle_at_risk <- max(zones_with_cattle[[fill_column]], na.rm = TRUE)
  
  # Create appropriate title and legend based on log transformation
  legend_name <- if(use_log) expression(atop("Infected cattle", "("*log[10] * ")")) else "Cattle\nat risk"
  
  p <- p + scale_fill_viridis_c(
      name = legend_name,
      na.value = "grey90",
      direction = 1,
      option = "plasma",  # Different color palette for cattle at risk
      limits = c(0, max_cattle_at_risk),
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
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_choropleth_mean.pdf", plot = p_mean, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_choropleth_lower.pdf", plot = p_lower, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_choropleth_upper.pdf", plot = p_upper, width = 12, height = 10)

# Log scale plots
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_choropleth_mean_log.pdf", plot = p_mean_log, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_choropleth_lower_log.pdf", plot = p_lower_log, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_choropleth_upper_log.pdf", plot = p_upper_log, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each estimate type) - more compact dimensions
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_mean_analysis.pdf", plot = p_mean_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_lower_analysis.pdf", plot = p_lower_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_upper_analysis.pdf", plot = p_upper_combined, width = 10, height = 6)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_mean_log_analysis.pdf", plot = p_mean_log_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_lower_log_analysis.pdf", plot = p_lower_log_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_fine_whole_country_plots/eth_cattle_at_risk_upper_log_analysis.pdf", plot = p_upper_log_combined, width = 10, height = 6)

