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

# Also get state boundaries (level 1) for context
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

# Load cattle at risk data from the summary file
cat("Loading cattle at risk data from Nigeria_cattle_uncertainty_fine.csv...\n")
cattle_summary_data <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Nigeria_cattle_uncertainty_fine.csv")

cat("Loaded cattle at risk data with", nrow(cattle_summary_data), "states\n")
cat("Columns available:", paste(names(cattle_summary_data), collapse = ", "), "\n")

# Load tsetse distribution data to properly classify tsetse zones
cat("Loading tsetse distribution data for proper zone classification...\n")
library(terra)
library(raster)
tsetse_raster <- raster("Data/Covariates/tsenumbspec")
# Convert values > 1 to 1 (binary presence/absence) as done in tsetse plot script
tsetse_raster[tsetse_raster > 1] <- 1

# Use area-based extraction instead of centroid-based to capture LGAs that partially overlap with tsetse areas
cat("Extracting tsetse presence for each LGA using area-based method...\n")

# Extract tsetse values for each LGA using area statistics
# This will give us the mean tsetse value within each LGA
lga_tsetse_stats <- raster::extract(tsetse_raster, nigeria_lgas_sf, fun = mean, na.rm = TRUE, df = TRUE)

# Create a data frame linking LGA names to tsetse presence
lga_tsetse_data <- data.frame(
  NAME_2 = nigeria_lgas_sf$NAME_2,
  NAME_1 = nigeria_lgas_sf$NAME_1,  # Include state names for geographical rules
  tsetse_mean = lga_tsetse_stats[,2]  # Second column contains the mean values
)

# Replace NA values with 0 (no tsetse) and create binary tsetse zone indicator
lga_tsetse_data$tsetse_mean[is.na(lga_tsetse_data$tsetse_mean)] <- 0

# Apply geographical knowledge: certain states are entirely within tsetse zones
# These states should be considered tsetse zones regardless of raster values
tsetse_states <- c("Abia", "Akwa Ibom", "Anambra", "Bayelsa", "Benue", "Cross River", 
                   "Delta", "Ebonyi", "Edo", "Ekiti", "Enugu", "Imo", "Lagos", "Ogun", 
                   "Ondo", "Osun", "Oyo", "Rivers", "Taraba", "Adamawa", "Plateau", 
                   "Nassarawa", "Kwara", "Kogi", "Federal Capital Territory")

# An LGA is considered tsetse zone if:
# 1. ANY part of it has tsetse in the raster (mean > 0), OR
# 2. It's in a state known to be entirely within the tsetse belt
lga_tsetse_data$tsetse_zone <- (lga_tsetse_data$tsetse_mean > 0) | 
                               (lga_tsetse_data$NAME_1 %in% tsetse_states)

cat("Found", sum(lga_tsetse_data$tsetse_zone, na.rm = TRUE), "LGAs with tsetse presence out of", nrow(lga_tsetse_data), "total LGAs\n")
cat("Using area-based extraction + geographical knowledge: LGAs with partial tsetse coverage are now correctly classified\n")
cat("Raster-based tsetse LGAs:", sum(lga_tsetse_data$tsetse_mean > 0, na.rm = TRUE), "\n")
cat("Geographically-assigned tsetse LGAs:", sum(lga_tsetse_data$NAME_1 %in% tsetse_states & lga_tsetse_data$tsetse_mean == 0, na.rm = TRUE), "\n")

# Create LGA-level data for plotting - match LGA names to GADM LGA names
# The CSV uses 'state' column and 'mean_cattle_mean', 'data_uncertainty_lower', 'data_uncertainty_upper'
# Now use actual tsetse data to distinguish between no cattle vs no risk areas
lga_cattle_at_risk_mean <- cattle_summary_data %>%
  dplyr::select(NAME_2 = state, 
                mean_cattle_at_risk = mean_cattle_mean,
                risk_value = mean_value_mean) %>%
  right_join(lga_tsetse_data, by = "NAME_2") %>%  # Use right_join to keep all LGAs
  dplyr::mutate(
    # If no cattle data, assume 0 cattle in that LGA
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    log_cattle_at_risk = log10(mean_cattle_at_risk + 1),  # Add log transformation
    # Create categories using actual tsetse distribution data
    risk_category = case_when(
      mean_cattle_at_risk > 0 ~ "cattle_at_risk",
      mean_cattle_at_risk == 0 & tsetse_zone == TRUE ~ "no_cattle_but_risk", # Tsetse zones with no cattle
      TRUE ~ "no_risk_or_cattle"  # Areas outside tsetse zones
    )
  )

lga_cattle_at_risk_lower <- cattle_summary_data %>%
  dplyr::select(NAME_2 = state, 
                lower_cattle_at_risk = data_uncertainty_lower,
                mean_cattle_at_risk = mean_cattle_mean,
                risk_value = mean_value_mean) %>%
  right_join(lga_tsetse_data, by = "NAME_2") %>%  # Use right_join to keep all LGAs
  dplyr::mutate(
    # If no cattle data, assume 0 cattle in that LGA
    lower_cattle_at_risk = ifelse(is.na(lower_cattle_at_risk), 0, lower_cattle_at_risk),
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    log_cattle_at_risk = log10(lower_cattle_at_risk + 1),
    risk_category = case_when(
      lower_cattle_at_risk > 0 ~ "cattle_at_risk",
      lower_cattle_at_risk == 0 & tsetse_zone == TRUE ~ "no_cattle_but_risk",
      TRUE ~ "no_risk_or_cattle"
    )
  )

lga_cattle_at_risk_upper <- cattle_summary_data %>%
  dplyr::select(NAME_2 = state, 
                upper_cattle_at_risk = data_uncertainty_upper,
                mean_cattle_at_risk = mean_cattle_mean,
                risk_value = mean_value_mean) %>%
  right_join(lga_tsetse_data, by = "NAME_2") %>%  # Use right_join to keep all LGAs
  dplyr::mutate(
    # If no cattle data, assume 0 cattle in that LGA
    upper_cattle_at_risk = ifelse(is.na(upper_cattle_at_risk), 0, upper_cattle_at_risk),
    mean_cattle_at_risk = ifelse(is.na(mean_cattle_at_risk), 0, mean_cattle_at_risk),
    log_cattle_at_risk = log10(upper_cattle_at_risk + 1),
    risk_category = case_when(
      upper_cattle_at_risk > 0 ~ "cattle_at_risk",
      upper_cattle_at_risk == 0 & tsetse_zone == TRUE ~ "no_cattle_but_risk",
      TRUE ~ "no_risk_or_cattle"
    )
  )

# Calculate 95% CI width for uncertainty assessment
lga_uncertainty <- lga_cattle_at_risk_lower %>%
  dplyr::select(NAME_2, lower = lower_cattle_at_risk) %>%
  left_join(lga_cattle_at_risk_upper %>% dplyr::select(NAME_2, upper = upper_cattle_at_risk), 
            by = "NAME_2", relationship = "many-to-many") %>%
  dplyr::mutate(
    ci_width = upper - lower,
    high_uncertainty = ci_width > quantile(ci_width, 0.9, na.rm = TRUE)  # Use 90th percentile for cattle at risk
  )

cat("LGAs with high uncertainty (CI width > 90th percentile):", sum(lga_uncertainty$high_uncertainty, na.rm = TRUE), "out of", nrow(lga_uncertainty), "\n")

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
  
  # Convert to spatial points and filter for Nigerian coordinate ranges
  # Nigeria coordinate ranges: Latitude 4-14°N, Longitude 3-15°E
  bovine_nigeria_filtered <- bovine_clean %>%
    dplyr::filter(Longitude >= 3, Longitude <= 15, Latitude >= 4, Latitude <= 14)
  
  if(nrow(bovine_nigeria_filtered) == 0) {
    cat("No bovine data points found within Nigerian boundaries\n")
    return(NULL)
  }
  
  bovine_sf <- st_as_sf(bovine_nigeria_filtered, 
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
if(!is.null(bovine_nigeria_sf)) {
  cat("Found", nrow(bovine_nigeria_sf), "bovine data points within Nigeria\n")
  cat("  - BCT tests:", sum(bovine_nigeria_sf$Test_Type == "BCT/HCT"), "\n")
  cat("  - PCR tests:", sum(bovine_nigeria_sf$Test_Type == "PCR"), "\n")
} else {
  cat("No bovine data available for Nigeria - maps will show choropleth only\n")
}

# Function to create choropleth map for cattle at risk
create_choropleth <- function(lga_data, title_suffix, data_type = "mean", add_bovine = FALSE, use_log = FALSE, add_inset = FALSE) {
  
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
  
  # Join LGA-level cattle at risk back to LGA polygons
  nigeria_lgas_final <- nigeria_lgas_sf %>%
    left_join(lga_data %>% st_drop_geometry() %>% dplyr::select(NAME_2, all_of(fill_column), risk_category), 
              by = "NAME_2")
  
  # Create LGA-level choropleth map
  p <- ggplot(nigeria_lgas_final)
  
  # Add neighboring countries as background context (if available)
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Separate areas with cattle at risk from zero areas
  lgas_with_cattle <- nigeria_lgas_final %>% dplyr::filter(risk_category == "cattle_at_risk" | is.na(risk_category))
  lgas_no_cattle_but_risk <- nigeria_lgas_final %>% dplyr::filter(risk_category == "no_cattle_but_risk")
  lgas_no_risk <- nigeria_lgas_final %>% dplyr::filter(risk_category == "no_risk_or_cattle")
  
  # Add areas with no cattle but risk (tsetse zones without cattle) in light blue
  if(nrow(lgas_no_cattle_but_risk) > 0) {
    p <- p + geom_sf(data = lgas_no_cattle_but_risk, fill = "lightblue", color = "white", lwd = 0.1)
  }
  
  # Add areas with no risk/cattle in light grey  
  if(nrow(lgas_no_risk) > 0) {
    p <- p + geom_sf(data = lgas_no_risk, fill = "lightgrey", color = "white", lwd = 0.1)
  }
  
  # Add Nigeria LGAs with cattle at risk data (main choropleth)
  if(nrow(lgas_with_cattle) > 0) {
    p <- p + geom_sf(data = lgas_with_cattle, aes(fill = !!sym(fill_column)), lwd = 0.1, color = "white")
  }
  
  # Add state boundaries for reference
  p <- p + geom_sf(data = nigeria_states_sf, fill = NA, color = "black", lwd = 0.5)
  
  # Calculate appropriate limits for cattle at risk
  max_cattle_at_risk <- max(lgas_with_cattle[[fill_column]], na.rm = TRUE)
  
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
      plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0.5, size = 9, color = "grey60"),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    )

  
  return(p)
}

# Create all three maps using lga_cattle_at_risk data
# Create choropleth maps - mean map first, then percentile maps with bovine data overlay
p_mean <- create_choropleth(lga_cattle_at_risk_mean, "(mean)", add_inset = TRUE)
p_lower <- create_choropleth(lga_cattle_at_risk_lower, "(2.5th percentile)", data_type = "lower", add_bovine = TRUE)
p_upper <- create_choropleth(lga_cattle_at_risk_upper, "(97.5th percentile)", data_type = "upper", add_bovine = TRUE)

# Create combined histogram and boxplot for each estimate type
create_combined_plot <- function(lga_data, estimate_name) {
  # Get LGA data without geometry for plotting
  plot_data <- lga_data %>% st_drop_geometry() %>% dplyr::filter(!is.na(mean_cattle_at_risk))
  
  # Histogram with plasma color scale matching choropleth
  p_hist <- ggplot(plot_data, aes(x = mean_cattle_at_risk, fill = after_stat(x))) +
    geom_histogram(bins = 30, color = "white", linewidth = 0.2,
                   boundary = 0, closed = "left") +
    scale_fill_viridis_c(option = "plasma", guide = "none") +  # Same as choropleth, no legend
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE),
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
      title = paste("Distribution of", estimate_name, "cattle at risk by LGA in Nigeria"),
      x = ""
    )
  
  # Boxplot with plasma color matching choropleth
  plasma_color <- viridis::plasma(1, begin = 0.5, end = 0.5)  # Mid-range plasma color
  
  p_box <- ggplot(plot_data, aes(x = mean_cattle_at_risk, y = 1)) +  # Use y = 1 instead of y = ""
    geom_boxplot(fill = plasma_color, alpha = 0.8, width = 0.8) +  # Much wider boxplot, no jitter
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                       expand = expansion(mult = c(0.02, 0.02))) +
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
      x = "Cattle at risk"
    )
  
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

# Create all choropleth maps - regular and log versions
# Regular scale maps
p_mean <- create_choropleth(lga_cattle_at_risk_mean, "(mean)", add_inset = TRUE)
p_lower <- create_choropleth(lga_cattle_at_risk_lower, "(2.5th percentile)", data_type = "lower", add_bovine = TRUE)
p_upper <- create_choropleth(lga_cattle_at_risk_upper, "(97.5th percentile)", data_type = "upper", add_bovine = TRUE)

# Log scale maps
p_mean_log <- create_choropleth(lga_cattle_at_risk_mean, "(mean - log scale)", use_log = TRUE, add_inset = TRUE)
p_lower_log <- create_choropleth(lga_cattle_at_risk_lower, "(2.5th percentile - log scale)", data_type = "lower", use_log = TRUE, add_bovine = TRUE)
p_upper_log <- create_choropleth(lga_cattle_at_risk_upper, "(97.5th percentile - log scale)", data_type = "upper", use_log = TRUE, add_bovine = TRUE)

# Create combined histogram and boxplot for each estimate type
create_combined_plot <- function(lga_data, estimate_name, use_log = FALSE) {
  # Get LGA data without geometry for plotting
  value_column <- if(use_log) "log_cattle_at_risk" else "mean_cattle_at_risk"
  plot_data <- lga_data %>% st_drop_geometry() %>% dplyr::filter(!is.na(!!sym(value_column)))
  
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
      title = paste("Distribution of", estimate_name, "cattle at risk by LGA in Nigeria"),
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
    labs(
      x = if(use_log) "Log10(Cattle at risk + 1)" else "Cattle at risk"
    )
  
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
p_mean_combined <- create_combined_plot(lga_cattle_at_risk_mean, "mean")
p_lower_combined <- create_combined_plot(lga_cattle_at_risk_lower, "2.5th percentile") 
p_upper_combined <- create_combined_plot(lga_cattle_at_risk_upper, "97.5th percentile")

p_mean_log_combined <- create_combined_plot(lga_cattle_at_risk_mean, "mean (log scale)", use_log = TRUE)
p_lower_log_combined <- create_combined_plot(lga_cattle_at_risk_lower, "2.5th percentile (log scale)", use_log = TRUE) 
p_upper_log_combined <- create_combined_plot(lga_cattle_at_risk_upper, "97.5th percentile (log scale)", use_log = TRUE)

# Save all plots to Analysis_NGA directory
# Regular scale plots
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_choropleth_mean.pdf", plot = p_mean, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_choropleth_lower.pdf", plot = p_lower, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_choropleth_upper.pdf", plot = p_upper, width = 12, height = 10)

# Log scale plots
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_choropleth_mean_log.pdf", plot = p_mean_log, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_choropleth_lower_log.pdf", plot = p_lower_log, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_choropleth_upper_log.pdf", plot = p_upper_log, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each estimate type) - more compact dimensions
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_mean_analysis.pdf", plot = p_mean_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_lower_analysis.pdf", plot = p_lower_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_upper_analysis.pdf", plot = p_upper_combined, width = 10, height = 6)

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_mean_log_analysis.pdf", plot = p_mean_log_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_lower_log_analysis.pdf", plot = p_lower_log_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_fine_plots/nga_cattle_at_risk_upper_log_analysis.pdf", plot = p_upper_log_combined, width = 10, height = 6)
