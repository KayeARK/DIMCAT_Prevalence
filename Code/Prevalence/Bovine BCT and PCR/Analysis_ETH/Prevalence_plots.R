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

##### LOAD DATA FOR COMPUTATION OF NUMBER OF INFECTED ANIMALS

regions <- gadm(country = "ETH", level = 1, path = tempdir())
country_sf <- st_as_sf(regions)


#CATTLE DATA
r_cattle<-raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
# r_buffalo<-raster("Data/Covariates/livestock/buffalo_2015/5_Bf_2015_Da.tif")
# r_goat<-raster("Data/Covariates/livestock/goats_2015/5_Gt_2015_Da.tif")
# r_horse<-raster("Data/Covariates/livestock/horses_2015/5_Ho_2015_Da.tif")
# r_pig<-raster("Data/Covariates/livestock/pigs_2015/5_Pg_2015_Da.tif")
# r_sheep<-raster("Data/Covariates/livestock/sheep_2015/5_Sh_2015_Da.tif")
##############

# Count number of projection files
n_datasets <- length(list.files("Code/Prevalence/Bovine BCT and PCR/Projections_ETH/", pattern = "Projections_model_.*\\.csv"))
cat("Found", n_datasets, "projection files for Ethiopia\n")

# PERFORMANCE OPTIMIZATION: Use terra instead of raster + vectorized operations
library(terra)

for (i in 1:n_datasets){
print(paste("Processing model", i, "of", n_datasets))
dpm <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Projections_ETH/Projections_model_",i,".csv"))

cat("Loaded", nrow(dpm), "prediction points\n")

# OPTIMIZATION 1: Process all three datasets together using data.table-style operations
# Split data by variable but keep in same data structure
dpm_mean <- dpm[dpm$variable == "Mean" & !is.na(dpm$value), ]
dpm_lower <- dpm[dpm$variable == "2.5th percentile" & !is.na(dpm$value), ]
dpm_upper <- dpm[dpm$variable == "97.5th percentile" & !is.na(dpm$value), ]

cat("Data points - Mean:", nrow(dpm_mean), "Lower:", nrow(dpm_lower), "Upper:", nrow(dpm_upper), "\n")

# OPTIMIZATION 2: Check coordinate swapping only once for mean data
lat_in_range <- all(dpm_mean$Latitude >= 3 & dpm_mean$Latitude <= 15, na.rm = TRUE)
lon_in_range <- all(dpm_mean$Longitude >= 33 & dpm_mean$Longitude <= 48, na.rm = TRUE)

# Apply coordinate swapping to all datasets if needed
if (!lat_in_range || !lon_in_range) {
  cat("Swapping coordinates for Ethiopia model", i, "\n")
  
  # Swap all datasets at once
  temp_mean <- dpm_mean$Longitude
  dpm_mean$Longitude <- dpm_mean$Latitude
  dpm_mean$Latitude <- temp_mean
  
  temp_lower <- dpm_lower$Longitude
  dpm_lower$Longitude <- dpm_lower$Latitude
  dpm_lower$Latitude <- temp_lower
  
  temp_upper <- dpm_upper$Longitude
  dpm_upper$Longitude <- dpm_upper$Latitude
  dpm_upper$Latitude <- temp_upper
} else {
  cat("Coordinates already correct for Ethiopia model", i, "\n")
}

# Convert to spatial objects ONCE
dpm_mean_sf <- st_as_sf(dpm_mean, coords = c("Longitude", "Latitude"), crs = 4326)
dpm_lower_sf <- st_as_sf(dpm_lower, coords = c("Longitude", "Latitude"), crs = 4326)
dpm_upper_sf <- st_as_sf(dpm_upper, coords = c("Longitude", "Latitude"), crs = 4326)

# ============================================================================
# UNCERTAINTY FILTERING: Remove pixels with large confidence intervals
# ============================================================================

cat("Applying uncertainty filtering...\n")
cat("Mean data points:", nrow(dpm_mean_sf), "\n")
cat("Lower CI data points:", nrow(dpm_lower_sf), "\n") 
cat("Upper CI data points:", nrow(dpm_upper_sf), "\n")

# Calculate CI width and filter (assuming same coordinate grid)
ci_width <- dpm_upper_sf$value - dpm_lower_sf$value
high_uncertainty_mask <- ci_width > 0.8
n_high_uncertainty <- sum(high_uncertainty_mask, na.rm = TRUE)
n_total <- length(high_uncertainty_mask)

cat("- High uncertainty pixels (CI width > 0.8):", n_high_uncertainty, "of", n_total, 
    "(", round(100 * n_high_uncertainty / n_total, 1), "%)\n")

# Apply filter to all datasets
dpm_mean_sf <- dpm_mean_sf[!high_uncertainty_mask, ]
dpm_lower_sf <- dpm_lower_sf[!high_uncertainty_mask, ]
dpm_upper_sf <- dpm_upper_sf[!high_uncertainty_mask, ]

cat("Retained", nrow(dpm_mean_sf), "pixels for analysis\n")

###### MAJOR PERFORMANCE OPTIMIZATION - VECTORIZED SPATIAL OPERATIONS

cat("Starting optimized spatial processing...\n")

# STEP 1: Convert cattle raster to terra for speed
if (class(r_cattle)[1] != "SpatRaster") {
  r_cattle_terra <- rast(r_cattle)
} else {
  r_cattle_terra <- r_cattle
}

# STEP 2: Batch coordinate transformation (do once for all datasets)
coords_all <- rbind(
  st_coordinates(dpm_mean_sf),
  st_coordinates(dpm_lower_sf), 
  st_coordinates(dpm_upper_sf)
)

# Transform coordinates to raster CRS
coords_all_sf <- st_as_sf(data.frame(coords_all), coords = c("X", "Y"), crs = 4326)
coords_all_transformed <- st_transform(coords_all_sf, crs = st_crs(r_cattle_terra))
coords_all_transformed_xy <- st_coordinates(coords_all_transformed)

# STEP 3: Single batch raster extraction (MAJOR SPEEDUP!)
cat("Extracting cattle density for all", nrow(coords_all_transformed), "points...\n")
cattle_extract_result <- terra::extract(r_cattle_terra, coords_all_transformed_xy, method = "simple")
# Handle the column indexing properly
if (ncol(cattle_extract_result) > 1) {
  cattle_values_all <- cattle_extract_result[,2]
} else {
  cattle_values_all <- cattle_extract_result[,1]
}
cattle_values_all[is.na(cattle_values_all)] <- 0

# Split back into mean, lower, upper
n_points <- nrow(dpm_mean_sf)
cattle_mean <- cattle_values_all[1:n_points]
cattle_lower <- cattle_values_all[(n_points+1):(2*n_points)]
cattle_upper <- cattle_values_all[(2*n_points+1):(3*n_points)]

cat("Extracted cattle density values successfully\n")

# STEP 4: Calculate infected cattle (vectorized)
infected_cattle_mean <- dpm_mean_sf$value * cattle_mean
infected_cattle_lower <- dpm_lower_sf$value * cattle_lower  
infected_cattle_upper <- dpm_upper_sf$value * cattle_upper

# STEP 5: Single batch spatial join with country boundaries  
cat("Performing spatial joins with Ethiopian regions...\n")

# Transform one set of coordinates for region assignment
coords_mean_country <- st_transform(dpm_mean_sf, crs = st_crs(country_sf))
region_assignments <- st_join(coords_mean_country, country_sf, join = st_within)$NAME_1

# Create results data frame (much faster than sf operations)
results_df <- data.frame(
  region = region_assignments,
  mean_value = dpm_mean_sf$value,
  lower_value = dpm_lower_sf$value,  
  upper_value = dpm_upper_sf$value,
  cattle_mean = infected_cattle_mean,
  cattle_lower = infected_cattle_lower,
  cattle_upper = infected_cattle_upper
)

# Remove NA regions
results_df <- results_df[!is.na(results_df$region), ]

cat("Spatial processing complete. Aggregating by region...\n")

# STEP 6: Fast regional aggregation using data.table-style operations
library(dplyr)

# Aggregate by region (much faster than sf operations)
dpm_combined_region <- results_df %>%
  group_by(region) %>%
  summarise(
    mean_value = mean(mean_value, na.rm = TRUE),
    cattle_mean = sum(cattle_mean, na.rm = TRUE),
    lower_value = mean(lower_value, na.rm = TRUE), 
    cattle_lower = sum(cattle_lower, na.rm = TRUE),
    upper_value = mean(upper_value, na.rm = TRUE),
    cattle_upper = sum(cattle_upper, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    # Round livestock numbers
    cattle_mean = round(cattle_mean, 0),
    cattle_lower = round(cattle_lower, 0),
    cattle_upper = round(cattle_upper, 0)
  )

# Add total row
total_row <- data.frame(
  region = "Total",
  mean_value = mean(results_df$mean_value, na.rm = TRUE),
  cattle_mean = sum(results_df$cattle_mean, na.rm = TRUE),
  lower_value = mean(results_df$lower_value, na.rm = TRUE),
  cattle_lower = sum(results_df$cattle_lower, na.rm = TRUE),
  upper_value = mean(results_df$upper_value, na.rm = TRUE),
  cattle_upper = sum(results_df$cattle_upper, na.rm = TRUE)
)

# Round totals
total_row$cattle_mean <- round(total_row$cattle_mean, 0)
total_row$cattle_lower <- round(total_row$cattle_lower, 0) 
total_row$cattle_upper <- round(total_row$cattle_upper, 0)

# Combine results
dpm_combined_region <- rbind(dpm_combined_region, total_row)

cat("Regional aggregation complete!\n")
cat("Processing model", i, "completed successfully\n\n")

#######

write.csv(dpm_combined_region, paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk/Projections_model_", i, ".csv"), row.names = FALSE)



###### RUN CATTLE DATA EXTRACTION - MEAN

# Process Mean
dpm_mean_sf <- st_transform(dpm_mean_sf, crs = crs(r_cattle))
dpm_sp <- as(dpm_mean_sf, "Spatial")
cell_nums <- cellFromXY(r_cattle, coordinates(dpm_sp))
cell_coords <- xyFromCell(r_cattle, cell_nums)
result_mean <- dpm_mean_sf %>%
  mutate(cell_id = cell_nums,
         center_lon = cell_coords[, 1],
         center_lat = cell_coords[, 2])

# for locations in the same cell, take the mean of the values
averaged_result_mean <- result_mean %>%
  group_by(center_lon, center_lat) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# rename center_lon and center_lat to longitude and latitude
dpm_mean_sf <- averaged_result_mean %>%
  rename(longitude = center_lon, latitude = center_lat, 
         value = mean_value)

# Convert back to sf object and join with regions
dpm_mean_sf <- st_as_sf(dpm_mean_sf, coords = c("longitude", "latitude"), crs = crs(r_cattle))
dpm_mean_sf <- st_transform(dpm_mean_sf, crs = st_crs(country_sf))
country_sf <- st_make_valid(country_sf)
dpm_mean_sf <- st_join(dpm_mean_sf, country_sf, join = st_within)

# Select and rename as needed
dpm_mean_sf <- dpm_mean_sf %>%
  dplyr::select(value, NAME_1) %>%
  rename(region = NAME_1)

###### RUN CATTLE DATA EXTRACTION - 2.5TH PERCENTILE

# Process 2.5th percentile
dpm_lower_sf <- st_transform(dpm_lower_sf, crs = crs(r_cattle))
dpm_sp_lower <- as(dpm_lower_sf, "Spatial")
cell_nums_lower <- cellFromXY(r_cattle, coordinates(dpm_sp_lower))
cell_coords_lower <- xyFromCell(r_cattle, cell_nums_lower)
result_lower <- dpm_lower_sf %>%
  mutate(cell_id = cell_nums_lower,
         center_lon = cell_coords_lower[, 1],
         center_lat = cell_coords_lower[, 2])

# for locations in the same cell, take the mean of the values
averaged_result_lower <- result_lower %>%
  group_by(center_lon, center_lat) %>%
  summarise(lower_value = mean(value, na.rm = TRUE), .groups = 'drop')

# rename center_lon and center_lat to longitude and latitude
dpm_lower_sf <- averaged_result_lower %>%
  rename(longitude = center_lon, latitude = center_lat, 
         value = lower_value)

# Convert back to sf object and join with regions
dpm_lower_sf <- st_as_sf(dpm_lower_sf, coords = c("longitude", "latitude"), crs = crs(r_cattle))
dpm_lower_sf <- st_transform(dpm_lower_sf, crs = st_crs(country_sf))
dpm_lower_sf <- st_join(dpm_lower_sf, country_sf, join = st_within)

# Select and rename as needed
dpm_lower_sf <- dpm_lower_sf %>%
  dplyr::select(value, NAME_1) %>%
  rename(region = NAME_1)

###### RUN CATTLE DATA EXTRACTION - 97.5TH PERCENTILE

# Process 97.5th percentile
dpm_upper_sf <- st_transform(dpm_upper_sf, crs = crs(r_cattle))
dpm_sp_upper <- as(dpm_upper_sf, "Spatial")
cell_nums_upper <- cellFromXY(r_cattle, coordinates(dpm_sp_upper))
cell_coords_upper <- xyFromCell(r_cattle, cell_nums_upper)
result_upper <- dpm_upper_sf %>%
  mutate(cell_id = cell_nums_upper,
         center_lon = cell_coords_upper[, 1],
         center_lat = cell_coords_upper[, 2])

# for locations in the same cell, take the mean of the values
averaged_result_upper <- result_upper %>%
  group_by(center_lon, center_lat) %>%
  summarise(upper_value = mean(value, na.rm = TRUE), .groups = 'drop')

# rename center_lon and center_lat to longitude and latitude
dpm_upper_sf <- averaged_result_upper %>%
  rename(longitude = center_lon, latitude = center_lat, 
         value = upper_value)

# Convert back to sf object and join with regions
dpm_upper_sf <- st_as_sf(dpm_upper_sf, coords = c("longitude", "latitude"), crs = crs(r_cattle))
dpm_upper_sf <- st_transform(dpm_upper_sf, crs = st_crs(country_sf))
dpm_upper_sf <- st_join(dpm_upper_sf, country_sf, join = st_within)

# Select and rename as needed
dpm_upper_sf <- dpm_upper_sf %>%
  dplyr::select(value, NAME_1) %>%
  rename(region = NAME_1)

# Cast to points
dpm_mean_sf <- st_cast(dpm_mean_sf, "POINT")
dpm_lower_sf <- st_cast(dpm_lower_sf, "POINT")
dpm_upper_sf <- st_cast(dpm_upper_sf, "POINT")

# OPTIMIZED: Extract cattle data for all datasets at once
# Combine all spatial points for single raster extraction
all_points <- rbind(
  dpm_mean_sf %>% mutate(dataset = "mean"),
  dpm_lower_sf %>% mutate(dataset = "lower"), 
  dpm_upper_sf %>% mutate(dataset = "upper")
)

# Single raster extraction (much faster)
all_points$cattle <- raster::extract(r_cattle, all_points)
all_points$cattle[is.na(all_points$cattle)] <- 0

# Split back into separate datasets
dpm_mean_sf$cattle <- all_points$cattle[all_points$dataset == "mean"]
dpm_lower_sf$cattle <- all_points$cattle[all_points$dataset == "lower"]
dpm_upper_sf$cattle <- all_points$cattle[all_points$dataset == "upper"]
# dpm_upper_sf$buffalo <- raster::extract(r_buffalo, dpm_upper_sf)
# dpm_upper_sf$buffalo[is.na(dpm_upper_sf$buffalo)] <- 0
# dpm_upper_sf$goat <- raster::extract(r_goat, dpm_upper_sf)
# dpm_upper_sf$goat[is.na(dpm_upper_sf$goat)] <- 0
# dpm_upper_sf$horse <- raster::extract(r_horse, dpm_upper_sf)
# dpm_upper_sf$horse[is.na(dpm_upper_sf$horse)] <- 0
# dpm_upper_sf$pig <- raster::extract(r_pig, dpm_upper_sf)
# dpm_upper_sf$pig[is.na(dpm_upper_sf$pig)] <- 0
# dpm_upper_sf$sheep <- raster::extract(r_sheep, dpm_upper_sf)
# dpm_upper_sf$sheep[is.na(dpm_upper_sf$sheep)] <- 0

#set cattle to value*cattle (MEAN)
dpm_mean_sf$cattle <- dpm_mean_sf$value * dpm_mean_sf$cattle
# dpm_mean_sf$buffalo <- dpm_mean_sf$value * dpm_mean_sf$buffalo
# dpm_mean_sf$goat <- dpm_mean_sf$value * dpm_mean_sf$goat
# dpm_mean_sf$horse <- dpm_mean_sf$value * dpm_mean_sf$horse
# dmp_mean_sf$pig <- dpm_mean_sf$value * dpm_mean_sf$pig
# dpm_mean_sf$sheep <- dpm_mean_sf$value * dpm_mean_sf$sheep

#set cattle to value*cattle (2.5TH PERCENTILE)
dpm_lower_sf$cattle <- dpm_lower_sf$value * dpm_lower_sf$cattle
# dpm_lower_sf$buffalo <- dpm_lower_sf$value * dpm_lower_sf$buffalo
# dpm_lower_sf$goat <- dpm_lower_sf$value * dpm_lower_sf$goat
# dpm_lower_sf$horse <- dpm_lower_sf$value * dpm_lower_sf$horse
# dpm_lower_sf$pig <- dpm_lower_sf$value * dpm_lower_sf$pig
# dpm_lower_sf$sheep <- dpm_lower_sf$value * dpm_lower_sf$sheep

#set cattle to value*cattle (97.5TH PERCENTILE)
dpm_upper_sf$cattle <- dpm_upper_sf$value * dpm_upper_sf$cattle
# dpm_upper_sf$buffalo <- dpm_upper_sf$value * dpm_upper_sf$buffalo
# dpm_upper_sf$goat <- dpm_upper_sf$value * dpm_upper_sf$goat
# dpm_upper_sf$horse <- dpm_upper_sf$value * dpm_upper_sf$horse
# dpm_upper_sf$pig <- dpm_upper_sf$value * dpm_upper_sf$pig
# dpm_upper_sf$sheep <- dpm_upper_sf$value * dpm_upper_sf$sheep

#Average by region - MEAN
dpm_mean_region <- dpm_mean_sf %>%
  group_by(region) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            cattle_mean = sum(cattle, na.rm = TRUE))
            # buffalo_mean = sum(buffalo, na.rm = TRUE),
            # goat_mean = sum(goat, na.rm = TRUE),
            # horse_mean = sum(horse, na.rm = TRUE),
            # pig_mean = sum(pig, na.rm = TRUE),
            # sheep_mean = sum(sheep, na.rm = TRUE))

#Average by region - 2.5TH PERCENTILE
dpm_lower_region <- dpm_lower_sf %>%
  group_by(region) %>%
  summarise(lower_value = mean(value, na.rm = TRUE),
            cattle_lower = sum(cattle, na.rm = TRUE))
            # buffalo_lower = sum(buffalo, na.rm = TRUE),
            # goat_lower = sum(goat, na.rm = TRUE),
            # horse_lower = sum(horse, na.rm = TRUE),
            # pig_lower = sum(pig, na.rm = TRUE),
            # sheep_lower = sum(sheep, na.rm = TRUE))

#Average by region - 97.5TH PERCENTILE
dpm_upper_region <- dpm_upper_sf %>%
  group_by(region) %>%
  summarise(upper_value = mean(value, na.rm = TRUE),
            cattle_upper = sum(cattle, na.rm = TRUE))
            # buffalo_upper = sum(buffalo, na.rm = TRUE),
            # goat_upper = sum(goat, na.rm = TRUE),
            # horse_upper = sum(horse, na.rm = TRUE),
            # pig_upper = sum(pig, na.rm = TRUE),
            # sheep_upper = sum(sheep, na.rm = TRUE))

#remove geometry from all region data
dpm_mean_region <- st_drop_geometry(dpm_mean_region)
dpm_lower_region <- st_drop_geometry(dpm_lower_region)
dpm_upper_region <- st_drop_geometry(dpm_upper_region)

#remove NA region from all datasets
dpm_mean_region <- dpm_mean_region[!is.na(dpm_mean_region$region),]
dpm_lower_region <- dpm_lower_region[!is.na(dpm_lower_region$region),]
dpm_upper_region <- dpm_upper_region[!is.na(dpm_upper_region$region),]

# Merge all three datasets by region
dpm_combined_region <- merge(dpm_mean_region, dpm_lower_region, by = "region")
dpm_combined_region <- merge(dpm_combined_region, dpm_upper_region, by = "region")

#round all livestock numbers to 0 decimal places
livestock_cols <- grep("cattle", names(dpm_combined_region), value = TRUE)
dpm_combined_region[livestock_cols] <- lapply(dpm_combined_region[livestock_cols], function(x) round(as.numeric(x), 0))

#add a final row called "Total" with the total number of animals in each category
total_row <- data.frame(
  region = "Total",
  mean_value = mean(dpm_combined_region$mean_value, na.rm = TRUE),
  cattle_mean = sum(dpm_combined_region$cattle_mean, na.rm = TRUE),
  # buffalo_mean = sum(dpm_combined_region$buffalo_mean, na.rm = TRUE),
  # goat_mean = sum(dpm_combined_region$goat_mean, na.rm = TRUE),
  # horse_mean = sum(dpm_combined_region$horse_mean, na.rm = TRUE),
  # pig_mean = sum(dpm_combined_region$pig_mean, na.rm = TRUE),
  # sheep_mean = sum(dpm_combined_region$sheep_mean, na.rm = TRUE),
  lower_value = mean(dpm_combined_region$lower_value, na.rm = TRUE),
  cattle_lower = sum(dpm_combined_region$cattle_lower, na.rm = TRUE),
  # buffalo_lower = sum(dpm_combined_region$buffalo_lower, na.rm = TRUE),
  # goat_lower = sum(dpm_combined_region$goat_lower, na.rm = TRUE),
  # horse_lower = sum(dpm_combined_region$horse_lower, na.rm = TRUE),
  # pig_lower = sum(dpm_combined_region$pig_lower, na.rm = TRUE),
  # sheep_lower = sum(dpm_combined_region$sheep_lower, na.rm = TRUE),
  upper_value = mean(dpm_combined_region$upper_value, na.rm = TRUE),
  cattle_upper = sum(dpm_combined_region$cattle_upper, na.rm = TRUE)
  # buffalo_upper = sum(dpm_combined_region$buffalo_upper, na.rm = TRUE),
  # goat_upper = sum(dpm_combined_region$goat_upper, na.rm = TRUE),
  # horse_upper = sum(dpm_combined_region$horse_upper, na.rm = TRUE),
  # pig_upper = sum(dpm_combined_region$pig_upper, na.rm = TRUE),
  # sheep_upper = sum(dpm_combined_region$sheep_upper, na.rm = TRUE)
)

# Append the total row to the data frame
dpm_combined_region <- rbind(dpm_combined_region, total_row)

#######

write.csv(dpm_combined_region,paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk/Projections_model_",i,".csv"), row.names=FALSE)

}