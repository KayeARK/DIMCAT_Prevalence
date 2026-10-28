# NGA display update: matches ETH relative colour ranges; crops histogram/boxplot views.
# Analysis, coordinate handling, input files and exported numerical values unchanged.
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

# The supplied histogram supports a focused linear 0-25% between-dataset scale.
# Complementary within-dataset values use 75-100%; scales are labelled separately.
# Calculation, aggregation, file selection and coordinate handling are unchanged.
relative_between_limits <- c(0, 0.25)
relative_within_limits <- c(0.75, 1)

make_uncertainty_spec <- function(limits, proportion = TRUE, breaks = NULL) {
  if (is.null(breaks)) breaks <- pretty(limits, n = 5)
  breaks <- sort(unique(c(limits, breaks[breaks >= limits[1] & breaks <= limits[2]])))
  labels <- if (proportion) scales::percent(breaks, accuracy = 1) else
    scales::label_number()(breaks)
  if (proportion && limits[1] > 0) labels[1] <- paste0("<=", labels[1])
  if (proportion && limits[2] < 1) labels[length(labels)] <- paste0(">=", labels[length(labels)])
  list(limits = limits, breaks = breaks, labels = labels,
       proportion = proportion, binwidth = if (proportion) 0.01 else diff(limits)/25)
}
uncertainty_fill_scale <- function(spec, name = NULL, guide = "colourbar") {
  ggplot2::scale_fill_viridis_c(name = name, option = "magma",
    limits = spec$limits, breaks = spec$breaks, labels = spec$labels,
    trans = "identity", oob = scales::squish, na.value = "grey90", guide = guide)
}




#make empty dataframe of size 54912*4
dpm_all <- rep(0,54912)
n_datasets <- 1000
n_units <- 18304

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
dpm <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Projections_NGA/Projections_model_",i,".csv"))


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

#make a new dataframe with Latitude, Longitude, between_var, within_var, rel_within, rel_between
dpm_var <- data.frame(Latitude = dpm$Latitude[dpm$variable=="Mean"],
                      Longitude = dpm$Longitude[dpm$variable=="Mean"],
                      between_var = between_var,
                      within_var = within_var,
                      rel_within = rel_within,
                      rel_between = rel_between)

# Get Nigerian administrative boundaries for choropleth mapping
nigeria_lgas <- gadm(country = "NGA", level = 2, path = tempdir())
nigeria_lgas_sf <- st_as_sf(nigeria_lgas)
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

# Load bovine data from Excel files
bovine_bct_raw <- read_excel("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
bovine_pcr_raw <- read_excel("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")

# Add test type identifier to each dataset
bovine_bct_raw$Test_Type <- "HCT/BCT"
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
cat("  - BCT tests:", sum(bovine_nigeria_sf$Test_Type == "HCT/BCT"), "\n")
cat("  - PCR tests:", sum(bovine_nigeria_sf$Test_Type == "PCR"), "\n")

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

# Function to aggregate uncertainty data by LGA using spatial join with gap filling
aggregate_uncertainty_by_lga <- function(uncertainty_sf) {
  joined_data <- st_join(uncertainty_sf, nigeria_lgas_sf)
  
  # First, get LGAs with direct data points
  lga_uncertainty <- joined_data %>%
    st_drop_geometry() %>%
    dplyr::filter(!is.na(GID_2)) %>%  # Only keep points that fell within LGAs
    dplyr::group_by(GID_2, NAME_2, NAME_1) %>%
    dplyr::summarise(
      lga_value = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      n_points = n(),
      .groups = 'drop'
    )
  
  # Find LGAs without data (grey areas)
  all_lgas <- nigeria_lgas_sf %>%
    st_drop_geometry() %>%
    dplyr::select(GID_2, NAME_2, NAME_1)
  
  missing_lgas <- all_lgas %>%
    dplyr::anti_join(lga_uncertainty, by = "GID_2")
  
  cat("Found", nrow(missing_lgas), "LGAs without uncertainty data. Filling using spatial interpolation...\n")
  
  if(nrow(missing_lgas) > 0) {
    # For each missing LGA, find closest uncertainty points and interpolate
    for(i in seq_len(nrow(missing_lgas))) {
      lga_id <- missing_lgas$GID_2[i]
      
      # Get centroid of the missing LGA (suppress attribute warning)
      lga_centroid <- nigeria_lgas_sf %>%
        dplyr::filter(GID_2 == lga_id) %>%
        st_centroid(of_largest_polygon = TRUE)
      
      # Find distances to all uncertainty points
      distances <- st_distance(lga_centroid, uncertainty_sf)
      
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
      
      # Add to lga_uncertainty data
      new_row <- data.frame(
        GID_2 = lga_id,
        NAME_2 = missing_lgas$NAME_2[i],
        NAME_1 = missing_lgas$NAME_1[i],
        lga_value = interpolated_value,
        median_value = interpolated_value,
        min_value = interpolated_value,
        max_value = interpolated_value,
        n_points = 0  # Indicate this is interpolated
      )
      
      lga_uncertainty <- rbind(lga_uncertainty, new_row)
    }
  }
  
  return(lga_uncertainty)
}

# Process uncertainty data for choropleth mapping
rel_between_sf <- process_uncertainty_data(dpm_var, "rel_between")
rel_within_sf <- process_uncertainty_data(dpm_var, "rel_within")
between_var_sf <- process_uncertainty_data(dpm_var, "between_var")

# Aggregate by LGA with spatial interpolation
lga_rel_between <- aggregate_uncertainty_by_lga(rel_between_sf)
lga_rel_within <- aggregate_uncertainty_by_lga(rel_within_sf)
lga_between_var <- aggregate_uncertainty_by_lga(between_var_sf)

# Function to create choropleth map for uncertainty measures
create_uncertainty_choropleth <- function(lga_data, title_suffix, legend_title, spec) {
  # Join LGA-level uncertainty back to LGA polygons
  nigeria_lgas_final <- nigeria_lgas_sf %>%
    left_join(lga_data, by = "GID_2")
  
  # Create LGA-level choropleth map
  p <- ggplot(nigeria_lgas_final)
  
  # Add neighboring countries as background context (if available)
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Add Nigeria LGAs with uncertainty data
  p <- p + 
    geom_sf(aes(fill = lga_value), lwd = 0.05, color = "black") +
    geom_sf(data = nigeria_states_sf, fill = NA, color = "black", lwd = 0.5) +
    # Add bovine data points (different colors for BCT and PCR)
    geom_point(data = st_drop_geometry(bovine_nigeria_sf), 
               aes(x = lon, y = lat, size = Number_of_animal_tested, shape = Test_Type), colour="black", fill="white", stroke=0.8,
               alpha = 0.8) +  # Increased from 0.7 to 0.8 for better visibility
    # Set the continuous bar itself explicitly: legend.key.height alone can
    # leave its length dependent on the ggplot2 version's guide defaults.
    uncertainty_fill_scale(spec, name = legend_title,
      guide = ggplot2::guide_colourbar(
        barheight = grid::unit(7.5, "cm"),
        barwidth = grid::unit(0.5, "cm"))) +
    # Add size legend for sample size
    scale_size_continuous(name = "Sample size", 
                         range = c(1, 10),  # Increased from c(1, 4) to make sizes more distinguishable
                         breaks = c(10, 50, 100, 200),
                         labels = c("10", "50", "100", "200+")) +
    # Add color scale for test types
    scale_color_manual(name = "Test type",
                      values = c("HCT/BCT" = "red", "PCR" = "blue"),
                      guide = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
                      scale_shape_manual(name = "Test type",labels = c(
    "HCT/BCT" = "BCT",
    "PCR" = "PCR"
  ),values = c("HCT/BCT" = 21,"PCR" = 22)
)+
    # Add scale bar and north arrow
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                    bar_cols = c("black", "white"), line_width = 1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                          pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                          style = north_arrow_fancy_orienteering(text_size = 8)) +
    theme_void() +
    labs(
      title = paste("Uncertainty analysis by LGA in Nigeria", title_suffix)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.spacing.y = unit(1, "cm"),  # Original spacing between legends
      legend.box.spacing = unit(1, "cm")  # Add spacing between legend boxes
    )
  
  n_clipped <- sum(lga_data$lga_value < spec$limits[1] |
                   lga_data$lga_value > spec$limits[2], na.rm = TRUE)
  if (n_clipped > 0) p <- p + labs(caption = paste(n_clipped,
    "LGA values fall outside the colour range and use its endpoint colours.")) +
    theme(plot.caption = element_text(size = 8, colour = "grey50", hjust = 0.5))
  return(p)
}

# Histograms and maps share the same colour mapping. coord_cartesian crops
# the display only; histogram counts and boxplot summaries use all observations.
create_uncertainty_combined_plot <- function(lga_data, estimate_name, x_label, spec) {
  finite_values <- lga_data$lga_value[is.finite(lga_data$lga_value)]
  stopifnot(length(finite_values) > 0L)
  view_limits <- spec$limits
  axis_breaks <- spec$breaks
  axis_labels <- if (spec$proportion) scales::percent_format(accuracy = 1) else
    scales::label_number()
  p_hist <- ggplot(lga_data, aes(x = lga_value, fill = after_stat(x))) +
    geom_histogram(binwidth = spec$binwidth, boundary = 0, closed = "left",
                   colour = "white", linewidth = 0.2) +
    uncertainty_fill_scale(spec, guide = "none") +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels,
                       expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(xlim = view_limits, expand = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(t = 5, r = 15, b = 0, l = 15)) +
    labs(title = paste("Distribution of", estimate_name, "by LGA in Nigeria"), x = "")
  p_box <- ggplot(lga_data, aes(x = lga_value, y = 1)) +
    geom_boxplot(fill = "grey70", alpha = 0.8, width = 0.8) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels,
                       expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(xlim = view_limits, expand = FALSE) +
    scale_y_continuous(limits = c(0.2, 1.8), expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 10, margin = margin(t = 5)),
          axis.text.x = element_text(size = 9), axis.ticks.x = element_line(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid = element_blank(),
          plot.margin = margin(t = 0, r = 15, b = 5, l = 15)) +
    labs(x = x_label)
  n_outside <- sum(finite_values < view_limits[1] | finite_values > view_limits[2])
  title_text <- paste("Uncertainty Analysis -", estimate_name)
  if (n_outside > 0) title_text <- paste0(title_text, "\n", n_outside,
    " LGA values outside displayed range; retained in summaries")
  title <- grid::textGrob(title_text,
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  gridExtra::arrangeGrob(title, p_hist, p_box, nrow = 3,
                        heights = grid::unit(c(if (n_outside > 0) 1.4 else 0.8, 3.5, 2.5), "cm"))
}

between_spec <- make_uncertainty_spec(relative_between_limits, TRUE, seq(0, 0.25, 0.05))
within_spec <- make_uncertainty_spec(relative_within_limits, TRUE, seq(0.75, 1, 0.05))
variance_limits <- range(lga_between_var$lga_value, na.rm = TRUE)
if (!all(is.finite(variance_limits))) stop("No finite between-dataset variance values.")
if (diff(variance_limits) == 0) variance_limits <- c(0, max(variance_limits[2], 1e-6))
variance_spec <- make_uncertainty_spec(variance_limits, FALSE)

p_rel_between <- create_uncertainty_choropleth(lga_rel_between,
  "\n (relative variability from diagnostic uncertainty)", "Relative\nimportance", between_spec)
p_rel_within <- create_uncertainty_choropleth(lga_rel_within,
  "\n (relative variability from coverage uncertainty)", "Relative\nimportance", within_spec)
p_between_var <- create_uncertainty_choropleth(lga_between_var,
  "\n (Variance from diagnostic uncertainty)", "Variance\n(logit scale)", variance_spec)

p_rel_between_combined <- create_uncertainty_combined_plot(lga_rel_between,
  "relative variability from diagnostic uncertainty",
  "Relative importance of variability from diagnostic uncertainty", between_spec)
p_rel_within_combined <- create_uncertainty_combined_plot(lga_rel_within,
  "relative variability from coverage uncertainty",
  "Relative importance of variability from coverage uncertainty", within_spec)
p_between_var_combined <- create_uncertainty_combined_plot(lga_between_var,
  "variance from diagnostic uncertainty", "Variance from diagnostic uncertainty (logit scale)", variance_spec)

# Save all choropleth maps
dir.create("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA", recursive = TRUE, showWarnings = FALSE)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_rel_between_choropleth.pdf", plot = p_rel_between, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_rel_within_choropleth.pdf", plot = p_rel_within, width = 12, height = 10)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_between_var_choropleth.pdf", plot = p_between_var, width = 12, height = 10)

# Save combined plots (histogram + boxplot for each uncertainty measure)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_rel_between_analysis.pdf", plot = p_rel_between_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_rel_within_analysis.pdf", plot = p_rel_within_combined, width = 10, height = 6)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/lga_between_var_analysis.pdf", plot = p_between_var_combined, width = 10, height = 6)

cat("=== UNCERTAINTY ESTIMATES SUMMARY ===\n\n")

cat("Relative Between-Dataset Variability:\n")
cat("Total LGAs with data:", nrow(lga_rel_between), "\n")
if(nrow(lga_rel_between) > 0) {
  print(summary(lga_rel_between$lga_value))
}

cat("\nRelative Within-Dataset Variability:\n")
cat("Total LGAs with data:", nrow(lga_rel_within), "\n")
if(nrow(lga_rel_within) > 0) {
  print(summary(lga_rel_within$lga_value))
}

cat("\nBetween-Dataset Variance:\n")
cat("Total LGAs with data:", nrow(lga_between_var), "\n")
if(nrow(lga_between_var) > 0) {
  print(summary(lga_between_var$lga_value))
}



# Grid-level diagnostic plots use the same 0-25% display scale.
ggplot(dpm_var) + geom_tile(aes(Latitude, Longitude, fill = rel_between)) +
  coord_fixed(ratio = 1) +
  uncertainty_fill_scale(between_spec,
    name = "Relative importance of\nbetween-dataset variability")
ggsave("Relative_uncertainty_between.png", width = 12, height = 8)

ggplot(dpm_var, aes(x = rel_between)) +
  geom_histogram(aes(fill = after_stat(x)), binwidth = between_spec$binwidth,
                 boundary = 0, closed = "left", colour = "white") +
  uncertainty_fill_scale(between_spec, guide = "none") +
  scale_x_continuous(breaks = between_spec$breaks,
                     labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(xlim = relative_between_limits, expand = FALSE) +
  labs(x = "Relative importance of between-dataset variability", y = "Frequency",
       caption = paste(sum(dpm_var$rel_between > 0.25, na.rm = TRUE),
                       "grid values above 25%; retained in calculations and source data"))
ggsave("Relative_uncertainty_between_histogram.png", width = 8, height = 6)








# ============================================================
# FIGURE 6A SOURCE DATA
#
# Exports the numerical data underlying the Nigerian
# uncertainty maps and associated histogram/boxplot panels.
#
# Sheets:
#   README
#   LGA_values
#   Distribution_values
#   Grid_values
#   Survey_points
#   Plot_settings
#   Summary_statistics
#
# ============================================================

# openxlsx is used only for writing the workbook.
# Explicit package namespaces are used to avoid function masking.

# ------------------------------------------------------------
# 1. Output location
# ------------------------------------------------------------

figure6a_output_dir <- paste0(
  "Code/Prevalence/Bovine BCT and PCR/",
  "Analysis_NGA"
)

dir.create(
  figure6a_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Unique Nigerian LGA lookup
# ------------------------------------------------------------

figure6a_lga_lookup <- nigeria_lgas_sf |>
  sf::st_drop_geometry() |>
  dplyr::transmute(
    GID_2 = as.character(GID_2),
    State = as.character(NAME_1),
    LGA = as.character(NAME_2)
  ) |>
  dplyr::distinct(
    GID_2,
    .keep_all = TRUE
  )

stopifnot(
  nrow(figure6a_lga_lookup) ==
    dplyr::n_distinct(figure6a_lga_lookup$GID_2)
)

# ------------------------------------------------------------
# 3. Standardise each LGA-level uncertainty object
# ------------------------------------------------------------

figure6a_rel_between <- lga_rel_between |>
  dplyr::as_tibble() |>
  dplyr::transmute(
    GID_2 = as.character(GID_2),
    State_from_result = as.character(NAME_1),
    LGA_from_result = as.character(NAME_2),

    Relative_diagnostic_uncertainty =
      as.numeric(lga_value),

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
      "Five-nearest-point inverse-distance interpolation"
    )
  )

figure6a_rel_within <- lga_rel_within |>
  dplyr::as_tibble() |>
  dplyr::transmute(
    GID_2 = as.character(GID_2),

    Relative_coverage_uncertainty =
      as.numeric(lga_value),

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
      "Five-nearest-point inverse-distance interpolation"
    )
  )

figure6a_between_var <- lga_between_var |>
  dplyr::as_tibble() |>
  dplyr::transmute(
    GID_2 = as.character(GID_2),

    Diagnostic_variance =
      as.numeric(lga_value),

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
      "Five-nearest-point inverse-distance interpolation"
    )
  )

# Confirm one row per GID_2 in each result object
stopifnot(
  nrow(figure6a_rel_between) ==
    dplyr::n_distinct(figure6a_rel_between$GID_2),

  nrow(figure6a_rel_within) ==
    dplyr::n_distinct(figure6a_rel_within$GID_2),

  nrow(figure6a_between_var) ==
    dplyr::n_distinct(figure6a_between_var$GID_2)
)

# ------------------------------------------------------------
# 4. Combined LGA-level source-data table
# ------------------------------------------------------------

figure6a_lga_values <- figure6a_lga_lookup |>
  dplyr::left_join(
    figure6a_rel_between,
    by = "GID_2",
    relationship = "one-to-one"
  ) |>
  dplyr::left_join(
    figure6a_rel_within,
    by = "GID_2",
    relationship = "one-to-one"
  ) |>
  dplyr::left_join(
    figure6a_between_var,
    by = "GID_2",
    relationship = "one-to-one"
  ) |>
  dplyr::mutate(
    # Numerical check: the two relative components should sum
    # approximately to one.
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
    State,
    LGA,
    GID_2
  )

stopifnot(
  nrow(figure6a_lga_values) ==
    dplyr::n_distinct(figure6a_lga_values$GID_2)
)

# ------------------------------------------------------------
# 5. Long-format values used by the choropleths,
#    histograms and boxplots
# ------------------------------------------------------------

figure6a_distribution_values <- dplyr::bind_rows(

  figure6a_lga_values |>
    dplyr::transmute(
      GID_2,
      State,
      LGA,

      Uncertainty_measure =
        "Relative diagnostic uncertainty",

      Plot_object =
        "lga_rel_between",

      Source_column =
        "lga_value",

      Plotted_value =
        Relative_diagnostic_uncertainty,

      Number_of_grid_points =
        Number_of_grid_points,

      Interpolated =
        Diagnostic_value_interpolated
    ),

  figure6a_lga_values |>
    dplyr::transmute(
      GID_2,
      State,
      LGA,

      Uncertainty_measure =
        "Relative coverage uncertainty",

      Plot_object =
        "lga_rel_within",

      Source_column =
        "lga_value",

      Plotted_value =
        Relative_coverage_uncertainty,

      Number_of_grid_points =
        Coverage_number_of_grid_points,

      Interpolated =
        Coverage_value_interpolated
    ),

  figure6a_lga_values |>
    dplyr::transmute(
      GID_2,
      State,
      LGA,

      Uncertainty_measure =
        "Diagnostic variance",

      Plot_object =
        "lga_between_var",

      Source_column =
        "lga_value",

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
    State,
    LGA
  )

# ------------------------------------------------------------
# 6. Prediction-grid uncertainty values
# ------------------------------------------------------------

# process_uncertainty_data() swaps Longitude and Latitude before
# plotting. The coordinates below therefore represent the
# corrected plotting coordinates, not merely the original
# column labels in dpm_var.

figure6a_grid_values <- dpm_var |>
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
  exists("bovine_nigeria_sf") &&
  !is.null(bovine_nigeria_sf) &&
  nrow(bovine_nigeria_sf) > 0
) {

  figure6a_survey_points <- bovine_nigeria_sf |>
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

      Test_type =
        as.character(Test_Type)
    ) |>
    dplyr::arrange(
      Test_type,
      Longitude,
      Latitude
    )

} else {

  figure6a_survey_points <- data.frame(
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
# 8. Summary statistics
# ------------------------------------------------------------

make_figure6a_summary <- function(
    values,
    measure_name
) {

  values <- values[
    is.finite(values)
  ]

  data.frame(
    Uncertainty_measure =
      measure_name,

    Number_of_LGAs =
      length(values),

    Mean =
      mean(values, na.rm = TRUE),

    Standard_deviation =
      stats::sd(values, na.rm = TRUE),

    Minimum =
      min(values, na.rm = TRUE),

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
      max(values, na.rm = TRUE),

    stringsAsFactors = FALSE
  )
}

figure6a_summary_statistics <- dplyr::bind_rows(

  make_figure6a_summary(
    figure6a_lga_values$
      Relative_diagnostic_uncertainty,
    "Relative diagnostic uncertainty"
  ),

  make_figure6a_summary(
    figure6a_lga_values$
      Relative_coverage_uncertainty,
    "Relative coverage uncertainty"
  ),

  make_figure6a_summary(
    figure6a_lga_values$
      Diagnostic_variance,
    "Diagnostic variance"
  )
)

# Add interpolation counts to the summary table
figure6a_interpolation_summary <- data.frame(
  Uncertainty_measure = c(
    "Relative diagnostic uncertainty",
    "Relative coverage uncertainty",
    "Diagnostic variance"
  ),

  Directly_aggregated_LGAs = c(
    sum(
      !figure6a_lga_values$
        Diagnostic_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      !figure6a_lga_values$
        Coverage_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      !figure6a_lga_values$
        Variance_value_interpolated,
      na.rm = TRUE
    )
  ),

  Interpolated_LGAs = c(
    sum(
      figure6a_lga_values$
        Diagnostic_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      figure6a_lga_values$
        Coverage_value_interpolated,
      na.rm = TRUE
    ),

    sum(
      figure6a_lga_values$
        Variance_value_interpolated,
      na.rm = TRUE
    )
  ),

  stringsAsFactors = FALSE
)

figure6a_summary_statistics <-
  figure6a_summary_statistics |>
  dplyr::left_join(
    figure6a_interpolation_summary,
    by = "Uncertainty_measure",
    relationship = "one-to-one"
  )

# ------------------------------------------------------------
# 9. Plot settings
# ------------------------------------------------------------

figure6a_plot_settings <- data.frame(
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
    "LGA aggregation",
    "Missing-LGA interpolation",
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
    "GADM level 2 Nigerian Local Government Areas",
    "GID_2",
    as.character(n_datasets),
    as.character(n_units),
    "Logit prevalence scale",
    "Variance across the 1,000 corrected-prevalence model datasets",
    paste(
      "Mean squared within-dataset standard deviation,",
      "estimated from the 2.5th and 97.5th percentiles"
    ),
    "Diagnostic variance plus coverage variance",
    "Diagnostic variance divided by total variance",
    "Coverage variance divided by total variance",
    "Arithmetic mean of prediction-grid values within each LGA",
    paste(
      "LGAs with no prediction-grid points were assigned",
      "an inverse-distance-weighted estimate"
    ),
    "Five nearest prediction-grid points",
    "Inverse distance with 1e-10 added to each distance",
    "Between: 0 to 0.25; within: 0.75 to 1; endpoint colours outside limits",
    paste0(
      format(
        min(
          lga_between_var$lga_value,
          na.rm = TRUE
        ),
        scientific = TRUE,
        digits = 8
      ),
      " to ",
      format(
        max(
          lga_between_var$lga_value,
          na.rm = TRUE
        ),
        scientific = TRUE,
        digits = 8
      )
    ),
    "viridis magma",
    "Relative fractions: 0.01; absolute variance: colour-range width / 25",
    "0",
    "Left-closed",
    "Between: 0 to 0.25; within: 0.75 to 1; display cropped only; full data retained",
    "Between: 0, 0.05, ..., 0.25; within: 0.75, 0.80, ..., 1; variance: scale-specific",
    "0.8",
    "1 to 10",
    "10, 50, 100 and 200",
    "HCT/BCT = 21; PCR = 22",
    "WGS84, EPSG:4326"
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10. README
# ------------------------------------------------------------

figure6a_readme <- data.frame(
  Item = c(
    "Workbook description",
    "Associated figure",
    "LGA_values sheet",
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
      "Source data underlying Figure 6A, including",
      "Nigerian LGA-level uncertainty measures,",
      "prediction-grid uncertainty values, distribution",
      "plot inputs and bovine survey-point overlays."
    ),

    "Figure 6A.",

    paste(
      "One row per unique Nigerian GADM level-2 LGA.",
      "Contains all three mapped uncertainty measures,",
      "within-LGA summaries and interpolation indicators."
    ),

    paste(
      "Long-format LGA values used directly by the",
      "choropleths, histograms and boxplots."
    ),

    paste(
      "Prediction-grid uncertainty decomposition before",
      "aggregation to LGAs."
    ),

    paste(
      "Survey locations plotted over each uncertainty",
      "choropleth."
    ),

    paste(
      "Transformations, interpolation rules, scale limits",
      "and plotting parameters."
    ),

    paste(
      "Descriptive summaries and counts of directly",
      "aggregated and interpolated LGAs."
    ),

    paste(
      "Proportion of total logit-scale variance attributed",
      "to variation across corrected-prevalence datasets."
    ),

    paste(
      "Proportion of total logit-scale variance attributed",
      "to uncertainty within each spatial model projection."
    ),

    paste(
      "Variance across corrected-prevalence model datasets",
      "on the logit prevalence scale."
    ),

    paste(
      "Mean within-dataset variance estimated using the",
      "2.5th and 97.5th percentiles."
    ),

    paste(
      "Indicates whether an LGA value was calculated from",
      "grid points lying directly inside the LGA or by",
      "inverse-distance interpolation."
    ),

    paste(
      "An interpolated LGA has Number_of_grid_points equal",
      "to zero and was assigned a value using the five",
      "nearest prediction-grid points."
    ),

    paste(
      "The input prediction files stored longitude and",
      "latitude in reversed columns. Plot_longitude and",
      "Plot_latitude are the corrected coordinates used",
      "by the plotting functions."
    ),

    paste(
      "Administrative polygon geometry is obtained",
      "separately from GADM and is not embedded in this",
      "workbook."
    ),

    paste(
      "Relative uncertainty values are proportions;",
      "variance values are on the squared logit-prevalence",
      "scale."
    )
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 11. Create workbook
# ------------------------------------------------------------

figure6a_wb <- openxlsx::createWorkbook(
  creator = "AR Kaye",
  title = "Figure 6A source data",
  subject = paste(
    "Nigeria uncertainty decomposition maps,",
    "histograms and boxplots"
  ),
  category = "Research data"
)

openxlsx::addWorksheet(
  figure6a_wb,
  "README"
)

openxlsx::addWorksheet(
  figure6a_wb,
  "LGA_values"
)

openxlsx::addWorksheet(
  figure6a_wb,
  "Distribution_values"
)

openxlsx::addWorksheet(
  figure6a_wb,
  "Grid_values"
)

openxlsx::addWorksheet(
  figure6a_wb,
  "Survey_points"
)

openxlsx::addWorksheet(
  figure6a_wb,
  "Plot_settings"
)

openxlsx::addWorksheet(
  figure6a_wb,
  "Summary_statistics"
)

openxlsx::writeData(
  figure6a_wb,
  "README",
  figure6a_readme
)

openxlsx::writeData(
  figure6a_wb,
  "LGA_values",
  figure6a_lga_values
)

openxlsx::writeData(
  figure6a_wb,
  "Distribution_values",
  figure6a_distribution_values
)

openxlsx::writeData(
  figure6a_wb,
  "Grid_values",
  figure6a_grid_values
)

openxlsx::writeData(
  figure6a_wb,
  "Survey_points",
  figure6a_survey_points
)

openxlsx::writeData(
  figure6a_wb,
  "Plot_settings",
  figure6a_plot_settings
)

openxlsx::writeData(
  figure6a_wb,
  "Summary_statistics",
  figure6a_summary_statistics
)

# ------------------------------------------------------------
# 12. Workbook styles
# ------------------------------------------------------------

figure6a_header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "bottom",
  borderColour = "#FFFFFF"
)

figure6a_decimal_style <- openxlsx::createStyle(
  numFmt = "0.000000"
)

figure6a_variance_style <- openxlsx::createStyle(
  numFmt = "0.0000000000"
)

figure6a_integer_style <- openxlsx::createStyle(
  numFmt = "0"
)

figure6a_percentage_style <- openxlsx::createStyle(
  numFmt = "0.000%"
)

figure6a_wrap_style <- openxlsx::createStyle(
  wrapText = TRUE,
  valign = "top"
)

# ------------------------------------------------------------
# 13. General formatting helper
# ------------------------------------------------------------

figure6a_format_data_sheet <- function(
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
    style = figure6a_header_style,
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

figure6a_format_data_sheet(
  figure6a_wb,
  "LGA_values",
  figure6a_lga_values
)

figure6a_format_data_sheet(
  figure6a_wb,
  "Distribution_values",
  figure6a_distribution_values
)

figure6a_format_data_sheet(
  figure6a_wb,
  "Grid_values",
  figure6a_grid_values
)

figure6a_format_data_sheet(
  figure6a_wb,
  "Survey_points",
  figure6a_survey_points
)

figure6a_format_data_sheet(
  figure6a_wb,
  "Summary_statistics",
  figure6a_summary_statistics
)

# ------------------------------------------------------------
# 14. README formatting
# ------------------------------------------------------------

openxlsx::addStyle(
  figure6a_wb,
  "README",
  figure6a_header_style,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::addStyle(
  figure6a_wb,
  "README",
  figure6a_wrap_style,
  rows = 2:(nrow(figure6a_readme) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::setColWidths(
  figure6a_wb,
  "README",
  cols = 1,
  widths = 32
)

openxlsx::setColWidths(
  figure6a_wb,
  "README",
  cols = 2,
  widths = 90
)

openxlsx::freezePane(
  figure6a_wb,
  "README",
  firstRow = TRUE
)

# ------------------------------------------------------------
# 15. Plot-settings formatting
# ------------------------------------------------------------

openxlsx::addStyle(
  figure6a_wb,
  "Plot_settings",
  figure6a_header_style,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::addStyle(
  figure6a_wb,
  "Plot_settings",
  figure6a_wrap_style,
  rows = 2:(nrow(figure6a_plot_settings) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

openxlsx::setColWidths(
  figure6a_wb,
  "Plot_settings",
  cols = 1,
  widths = 42
)

openxlsx::setColWidths(
  figure6a_wb,
  "Plot_settings",
  cols = 2,
  widths = 90
)

openxlsx::freezePane(
  figure6a_wb,
  "Plot_settings",
  firstRow = TRUE
)

# ------------------------------------------------------------
# 16. Numeric formatting: LGA_values
# ------------------------------------------------------------

figure6a_lga_relative_columns <- match(
  c(
    "Relative_diagnostic_uncertainty",
    "Relative_coverage_uncertainty",
    "Relative_components_sum",
    "Relative_components_difference_from_one"
  ),
  names(figure6a_lga_values)
)

figure6a_lga_relative_columns <-
  figure6a_lga_relative_columns[
    !is.na(figure6a_lga_relative_columns)
  ]

if (
  nrow(figure6a_lga_values) > 0 &&
  length(figure6a_lga_relative_columns) > 0
) {

  openxlsx::addStyle(
    figure6a_wb,
    "LGA_values",
    figure6a_percentage_style,
    rows = 2:(nrow(figure6a_lga_values) + 1),
    cols = figure6a_lga_relative_columns,
    gridExpand = TRUE
  )
}

figure6a_lga_variance_columns <- match(
  c(
    "Diagnostic_variance",
    "Diagnostic_variance_median",
    "Diagnostic_variance_minimum",
    "Diagnostic_variance_maximum"
  ),
  names(figure6a_lga_values)
)

figure6a_lga_variance_columns <-
  figure6a_lga_variance_columns[
    !is.na(figure6a_lga_variance_columns)
  ]

if (
  nrow(figure6a_lga_values) > 0 &&
  length(figure6a_lga_variance_columns) > 0
) {

  openxlsx::addStyle(
    figure6a_wb,
    "LGA_values",
    figure6a_variance_style,
    rows = 2:(nrow(figure6a_lga_values) + 1),
    cols = figure6a_lga_variance_columns,
    gridExpand = TRUE
  )
}

figure6a_lga_integer_columns <- match(
  c(
    "Number_of_grid_points",
    "Coverage_number_of_grid_points",
    "Variance_number_of_grid_points"
  ),
  names(figure6a_lga_values)
)

figure6a_lga_integer_columns <-
  figure6a_lga_integer_columns[
    !is.na(figure6a_lga_integer_columns)
  ]

if (
  nrow(figure6a_lga_values) > 0 &&
  length(figure6a_lga_integer_columns) > 0
) {

  openxlsx::addStyle(
    figure6a_wb,
    "LGA_values",
    figure6a_integer_style,
    rows = 2:(nrow(figure6a_lga_values) + 1),
    cols = figure6a_lga_integer_columns,
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 17. Numeric formatting: Distribution_values
# ------------------------------------------------------------

figure6a_distribution_value_column <- match(
  "Plotted_value",
  names(figure6a_distribution_values)
)

if (
  nrow(figure6a_distribution_values) > 0 &&
  !is.na(figure6a_distribution_value_column)
) {

  openxlsx::addStyle(
    figure6a_wb,
    "Distribution_values",
    figure6a_decimal_style,
    rows = 2:(nrow(figure6a_distribution_values) + 1),
    cols = figure6a_distribution_value_column,
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 18. Numeric formatting: Grid_values
# ------------------------------------------------------------

figure6a_grid_coordinate_columns <- match(
  c(
    "Original_Latitude_column",
    "Original_Longitude_column",
    "Plot_longitude",
    "Plot_latitude"
  ),
  names(figure6a_grid_values)
)

figure6a_grid_coordinate_columns <-
  figure6a_grid_coordinate_columns[
    !is.na(figure6a_grid_coordinate_columns)
  ]

openxlsx::addStyle(
  figure6a_wb,
  "Grid_values",
  figure6a_decimal_style,
  rows = 2:(nrow(figure6a_grid_values) + 1),
  cols = figure6a_grid_coordinate_columns,
  gridExpand = TRUE
)

figure6a_grid_variance_columns <- match(
  c(
    "Diagnostic_variance",
    "Coverage_variance",
    "Total_variance"
  ),
  names(figure6a_grid_values)
)

figure6a_grid_variance_columns <-
  figure6a_grid_variance_columns[
    !is.na(figure6a_grid_variance_columns)
  ]

openxlsx::addStyle(
  figure6a_wb,
  "Grid_values",
  figure6a_variance_style,
  rows = 2:(nrow(figure6a_grid_values) + 1),
  cols = figure6a_grid_variance_columns,
  gridExpand = TRUE
)

figure6a_grid_relative_columns <- match(
  c(
    "Relative_coverage_uncertainty",
    "Relative_diagnostic_uncertainty",
    "Relative_components_sum"
  ),
  names(figure6a_grid_values)
)

figure6a_grid_relative_columns <-
  figure6a_grid_relative_columns[
    !is.na(figure6a_grid_relative_columns)
  ]

openxlsx::addStyle(
  figure6a_wb,
  "Grid_values",
  figure6a_percentage_style,
  rows = 2:(nrow(figure6a_grid_values) + 1),
  cols = figure6a_grid_relative_columns,
  gridExpand = TRUE
)

# ------------------------------------------------------------
# 19. Numeric formatting: Survey_points
# ------------------------------------------------------------

if (nrow(figure6a_survey_points) > 0) {

  openxlsx::addStyle(
    figure6a_wb,
    "Survey_points",
    figure6a_decimal_style,
    rows = 2:(nrow(figure6a_survey_points) + 1),
    cols = c(2, 3),
    gridExpand = TRUE
  )

  openxlsx::addStyle(
    figure6a_wb,
    "Survey_points",
    figure6a_integer_style,
    rows = 2:(nrow(figure6a_survey_points) + 1),
    cols = c(1, 4, 5),
    gridExpand = TRUE
  )

  openxlsx::addStyle(
    figure6a_wb,
    "Survey_points",
    figure6a_percentage_style,
    rows = 2:(nrow(figure6a_survey_points) + 1),
    cols = 6,
    gridExpand = TRUE
  )
}

# ------------------------------------------------------------
# 20. Save workbook
# ------------------------------------------------------------

figure6a_output_file <- file.path(
  figure6a_output_dir,
  "Figure_6A_source_data.xlsx"
)

openxlsx::saveWorkbook(
  figure6a_wb,
  file = figure6a_output_file,
  overwrite = TRUE
)

message(
  "Figure 6A source-data workbook saved to: ",
  figure6a_output_file
)

message(
  "LGA rows exported: ",
  nrow(figure6a_lga_values)
)

message(
  "Distribution values exported: ",
  nrow(figure6a_distribution_values)
)

message(
  "Prediction-grid rows exported: ",
  nrow(figure6a_grid_values)
)

message(
  "Survey points exported: ",
  nrow(figure6a_survey_points)
)
