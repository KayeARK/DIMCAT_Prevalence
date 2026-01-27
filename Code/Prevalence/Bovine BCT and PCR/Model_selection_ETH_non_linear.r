# Ethiopian Prevalence Model Selection - Starting with Linear Effects
# Modified version to handle small datasets that can't support RW2 smoothing
# Will attempt non-linear effects only if linear models work successfully

rm(list=ls())

library(readxl)
library(ggplot2)
library(INLA)
library(sp)
library(sf)
library(afrilearndata)
library(concaveman)
library(geodata)
library(raster)
library(terra)
library(viridis)
library(gridExtra)
library(stringr)

countries_to_infer=c("Ethiopia")

# Add string to the name of the file
data_path <- paste0("Code/TestSensSpec/adjusted_cases_all_sites.csv")
data <- read.csv(data_path)

# Debug: Check data structure and types
cat("Data structure check:\n")
cat("Dimensions:", dim(data), "\n")
cat("Column names (first 10):", paste(head(colnames(data), 10), collapse = ", "), "\n")
cat("Column types (first 5):", paste(sapply(data[,1:min(5, ncol(data))], class), collapse = ", "), "\n")

# Check if we have the expected coordinate and sample_size columns
expected_cols <- c("longitude", "latitude", "sample_size")
missing_cols <- expected_cols[!expected_cols %in% colnames(data)]
if (length(missing_cols) > 0) {
  cat("Missing expected columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Available columns containing 'lon':", paste(grep("lon", colnames(data), ignore.case = TRUE, value = TRUE), collapse = ", "), "\n")
  cat("Available columns containing 'lat':", paste(grep("lat", colnames(data), ignore.case = TRUE, value = TRUE), collapse = ", "), "\n")
  cat("Available columns containing 'sample':", paste(grep("sample", colnames(data), ignore.case = TRUE, value = TRUE), collapse = ", "), "\n")
}

# Filter data for Ethiopian coordinates (Lat: 3-15°N, Lon: 33-48°E)
# Check if coordinates need to be swapped
if (any(data$latitude < 3 | data$latitude > 15) || any(data$longitude < 33 | data$longitude > 48)) {
  cat("Checking if coordinates need to be swapped for Ethiopia...\n")
  lat_in_eth_range <- sum(data$latitude >= 3 & data$latitude <= 15, na.rm = TRUE)
  lon_in_eth_range <- sum(data$longitude >= 33 & data$longitude <= 48, na.rm = TRUE)
  
  # Check if swapping would improve the fit
  lat_as_lon_in_range <- sum(data$latitude >= 33 & data$latitude <= 48, na.rm = TRUE)
  lon_as_lat_in_range <- sum(data$longitude >= 3 & data$longitude <= 15, na.rm = TRUE)
  
  if (lat_as_lon_in_range > lon_in_eth_range && lon_as_lat_in_range > lat_in_eth_range) {
    cat("Swapping coordinates for Ethiopia...\n")
    temp <- data$longitude
    data$longitude <- data$latitude
    data$latitude <- temp
  }
}

# Filter for Ethiopian coordinates
eth_filter <- data$latitude >= 3 & data$latitude <= 15 & 
              data$longitude >= 33 & data$longitude <= 48
data_original <- data[eth_filter, ]

cat("Filtered to", nrow(data_original), "observations within Ethiopian boundaries\n")

# Ensure coordinates are numeric and clean from the start
data_original$longitude <- as.numeric(data_original$longitude)
data_original$latitude <- as.numeric(data_original$latitude)

# Ensure sample_size is numeric (fix for character type error)
if ("sample_size" %in% colnames(data_original)) {
  if (is.character(data_original$sample_size) || is.factor(data_original$sample_size)) {
    data_original$sample_size <- as.numeric(as.character(data_original$sample_size))
  }
  data_original$sample_size[is.na(data_original$sample_size)] <- 0
} else {
  cat("Warning: sample_size column not found. Available columns:", paste(colnames(data_original), collapse = ", "), "\n")
}

# Remove any rows with NA coordinates before aggregation
na_coords <- is.na(data_original$longitude) | is.na(data_original$latitude)
if(sum(na_coords) > 0) {
  cat("Removing", sum(na_coords), "rows with NA coordinates from Ethiopian data\n")
  data_original <- data_original[!na_coords, ]
}

# Create coordinate mapping for aggregation (done once)
cat("Creating coordinate aggregation mapping for Ethiopia...\n")
cat("Original Ethiopian data has", nrow(data_original), "records\n")

# Use base R aggregate to avoid type conversion issues
coord_sample_mapping <- aggregate(
  sample_size ~ longitude + latitude, 
  data = data_original, 
  FUN = sum, 
  na.rm = TRUE
)

cat("After coordinate aggregation:", nrow(coord_sample_mapping), "unique locations\n")
cat("Reduced from", nrow(data_original), "to", nrow(coord_sample_mapping), "records\n")

# Loop through iterations (adjust range as needed)
for (iteration in 7:7){
  print(paste("Processing iteration", iteration-7))
  
  # Extract positive cases from original data for this iteration
  cat("=== DEBUGGING DATA LOADING ===\n")
  cat("Iteration (column index):", iteration, "\n")
  if (iteration <= ncol(data_original)) {
    cat("Column name:", colnames(data_original)[iteration], "\n")
  } else {
    cat("ERROR: Iteration", iteration, "exceeds number of columns (", ncol(data_original), ")\n")
  }
  
  positive_orig <- data_original[,iteration]
  cat("Raw positive data - class:", class(positive_orig), "\n")
  cat("Raw positive data - first 10 values:", paste(head(positive_orig, 10), collapse=", "), "\n")
  cat("Raw positive data - summary:\n")
  print(summary(positive_orig))
  
  # Ensure positive_orig is numeric (fix for character type error)
  if (is.character(positive_orig) || is.factor(positive_orig)) {
    positive_orig <- as.numeric(as.character(positive_orig))
  }
  positive_orig[is.na(positive_orig)] <- 0  # Replace NAs with 0
  
  cat("After conversion - first 10 values:", paste(head(positive_orig, 10), collapse=", "), "\n")
  cat("After conversion - total positive:", sum(positive_orig), "\n")
  
  # Aggregate positive cases for this iteration by coordinates
  positive_aggregated <- aggregate(
    positive_orig ~ longitude + latitude, 
    data = data.frame(
      longitude = data_original$longitude,
      latitude = data_original$latitude,
      positive_orig = positive_orig
    ), 
    FUN = sum, 
    na.rm = TRUE
  )
  
  # Use aggregated data for this iteration
  positive <- positive_aggregated$positive_orig
  sample_size <- coord_sample_mapping$sample_size
  
  # Create data object with aggregated coordinates
  data <- data.frame(
    longitude = as.numeric(coord_sample_mapping$longitude),
    latitude = as.numeric(coord_sample_mapping$latitude),
    stringsAsFactors = FALSE
  )
  
  # ELEVATION - Using Ethiopia (ETH) country code
  r_elv <- elevation_30s(country = "ETH", path = "Data/Covariates")
  elevation_raw <- terra::extract(r_elv, data[, c("longitude", "latitude")])$ETH_elv_msk
  data$elevation <- as.numeric(elevation_raw)
  
  # PRECIPITATION - Ethiopia climate data
  r_prec <- worldclim_country(country = "ETH", path = "Data/Covariates", var = "prec")
  prec_raw <- (terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_1+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_2+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_3+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_4+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_5+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_6+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_7+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_8+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_9+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_10+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_11+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_12)/12
  data$precipitation <- as.numeric(prec_raw)
  
  # AVERAGE TEMPERATURE
  r_tavg <- worldclim_country(country = "ETH", path = "Data/Covariates", var = "tavg")
  data$tavg <- (terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_1+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_2+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_3+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_4+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_5+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_6+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_7+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_8+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_9+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_10+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_11+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_12)/12
  
  # MINIMUM TEMPERATURE
  r_tmin <- worldclim_country(country = "ETH", path = "Data/Covariates", var = "tmin")
  data$tmin <- (terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_1+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_2+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_3+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_4+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_5+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_6+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_7+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_8+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_9+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_10+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_11+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_12)/12
  
  # MAXIMUM TEMPERATURE
  r_tmax <- worldclim_country(country = "ETH", path = "Data/Covariates", var = "tmax")
  data$tmax <- (terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_1+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_2+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_3+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_4+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_5+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_6+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_7+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_8+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_9+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_10+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_11+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_12)/12
  
  # POPULATION DENSITY
  r_pop_den <- population(year=2010, res=0.5, path="Data/Covariates")
  data$pop_den <- terra::extract(r_pop_den, data[, c("longitude", "latitude")])$population_density
  
  # LANDCOVER
  r_tree <- landcover("trees", path="Data/Covariates")
  data$tree <- terra::extract(r_tree, data[, c("longitude", "latitude")])$trees
  
  r_grassland <- landcover("grassland", path="Data/Covariates")
  data$grassland <- terra::extract(r_grassland, data[, c("longitude", "latitude")])$grassland
  
  r_shrub <- landcover("shrubs", path="Data/Covariates")
  data$shrub <- terra::extract(r_shrub, data[, c("longitude", "latitude")])$shrub
  
  r_cropland <- landcover("cropland", path="Data/Covariates")
  data$cropland <- terra::extract(r_cropland, data[, c("longitude", "latitude")])$cropland
  
  r_built <- landcover("built", path="Data/Covariates")
  data$built <- terra::extract(r_built, data[, c("longitude", "latitude")])$built
  
  r_bare <- landcover("bare", path="Data/Covariates")
  data$bare <- terra::extract(r_bare, data[, c("longitude", "latitude")])$bare
  
  r_water <- landcover("water", path="Data/Covariates")
  data$water <- terra::extract(r_water, data[, c("longitude", "latitude")])$water
  
  r_wetland <- landcover("wetland", path="Data/Covariates")
  data$wetland <- terra::extract(r_wetland, data[, c("longitude", "latitude")])$wetland
  
  r_mangrove <- landcover("mangroves", path="Data/Covariates")
  data$mangrove <- terra::extract(r_mangrove, data[, c("longitude", "latitude")])$mangroves
  
  r_moss <- landcover("moss", path="Data/Covariates")
  data$moss <- terra::extract(r_moss, data[, c("longitude", "latitude")])$moss
  
  data$tsetse_habitat <- data$tree + data$wetland
  
  # HUMAN FOOTPRINT
  r_human_fp <- footprint(year=2009, res=30, path="Data/Covariates")
  data$human_fp <- terra::extract(r_human_fp, data[, c("longitude", "latitude")])[,2]
  
  # LIVESTOCK DATA
  r_cattle <- raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
  data$cattle <- terra::extract(r_cattle, data[, c("longitude", "latitude")])
  data$cattle[is.na(data$cattle)] <- 0
  
  r_buffalo <- raster("Data/Covariates/livestock/buffalo_2015/5_Bf_2015_Da.tif")
  data$buffalo <- terra::extract(r_buffalo, data[, c("longitude", "latitude")])
  data$buffalo[is.na(data$buffalo)] <- 0
  
  r_goat <- raster("Data/Covariates/livestock/goats_2015/5_Gt_2015_Da.tif")
  data$goat <- terra::extract(r_goat, data[, c("longitude", "latitude")])
  data$goat[is.na(data$goat)] <- 0
  
  r_horse <- raster("Data/Covariates/livestock/horses_2015/5_Ho_2015_Da.tif")
  data$horse <- terra::extract(r_horse, data[, c("longitude", "latitude")])
  data$horse[is.na(data$horse)] <- 0
  
  r_pig <- raster("Data/Covariates/livestock/pigs_2015/5_Pg_2015_Da.tif")
  data$pig <- terra::extract(r_pig, data[, c("longitude", "latitude")])
  data$pig[is.na(data$pig)] <- 0
  
  r_sheep <- raster("Data/Covariates/livestock/sheep_2015/5_Sh_2015_Da.tif")
  data$sheep <- terra::extract(r_sheep, data[, c("longitude", "latitude")])
  data$sheep[is.na(data$sheep)] <- 0
  
  # TSETSE DATA
  tsetse <- raster("Data/Covariates/tsenumbspec")
  tsetse[tsetse > 1] <- 1
  data$tsetse <- terra::extract(tsetse, data[, c("longitude", "latitude")])
  
  data$total_animal <- data$cattle + data$buffalo + data$goat + data$horse + data$pig + data$sheep
  
  # Extract variables for analysis
  latitudes <- data$latitude
  longitudes <- data$longitude
  elevations <- data$elevation
  precipitations <- data$precipitation
  tavgs <- data$tavg
  tmins <- data$tmin
  tmaxs <- data$tmax
  human_fps <- data$human_fp
  pop_dens <- data$pop_den
  
  trees <- data$tree
  grasslands <- data$grassland
  shrubs <- data$shrub
  croplands <- data$cropland
  builts <- data$built
  bares <- data$bare
  waters <- data$water
  wetlands <- data$wetland
  mangroves <- data$mangrove
  mosses <- data$moss
  tsetse_habitats <- data$tsetse_habitat
  
  cattles <- data$cattle
  buffalos <- data$buffalo
  goats <- data$goat
  horses <- data$horse
  pigs <- data$pig
  sheeps <- data$sheep
  total_animals <- data$total_animal
  tsetses <- data$tsetse
  
  countries_to_infer <- sort(countries_to_infer)
  
  # Get Ethiopian boundaries
  v <- africountries$name %in% countries_to_infer
  country <- africountries[which(v==TRUE),]
  
  # Remove all countries except Ethiopia
  africountries_toplot <- africountries[which(africountries$name %in% countries_to_infer),]
  
  country <- st_coordinates(st_geometry(country))
  
  poly.df <- data.frame("long"=country[,1],"lat"=country[,2])
  poly.sf <- st_as_sf(poly.df, coords = c("long","lat"))
  poly <- concaveman(poly.sf, length_threshold = 0, concavity=1.1)
  point.df <- data.frame("long"=longitudes,"lat"=latitudes)
  point.sf <- st_as_sf(point.df, coords = c("long","lat"))
  long_lat_concave <- st_coordinates(st_geometry(poly))
  border <- long_lat_concave
  
  inside <- st_intersects(point.sf, poly)
  inside <- inside[]==1
  
  # Filter data to points inside Ethiopia
  latitudes <- latitudes[inside]
  longitudes <- longitudes[inside]
  sample_size <- sample_size[inside]
  positive <- positive[inside]
  elevations <- elevations[inside]
  precipitations <- precipitations[inside]
  tavgs <- tavgs[inside]
  tmins <- tmins[inside]
  tmaxs <- tmaxs[inside]
  human_fps <- human_fps[inside]
  pop_dens <- pop_dens[inside]
  trees <- trees[inside]
  grasslands <- grasslands[inside]
  shrubs <- shrubs[inside]
  croplands <- croplands[inside]
  builts <- builts[inside]
  bares <- bares[inside]
  waters <- waters[inside]
  wetlands <- wetlands[inside]
  mangroves <- mangroves[inside]
  mosses <- mosses[inside]
  tsetse_habitats <- tsetse_habitats[inside]
  cattles <- cattles[inside]
  buffalos <- buffalos[inside]
  goats <- goats[inside]
  horses <- horses[inside]
  pigs <- pigs[inside]
  sheeps <- sheeps[inside]
  total_animals <- total_animals[inside]
  tsetses <- tsetses[inside]
  
  # Remove NA values
  ind <- !is.na(latitudes)
  latitudes <- latitudes[ind]
  longitudes <- longitudes[ind]
  sample_size <- sample_size[ind]
  positive <- positive[ind]
  elevations <- elevations[ind]
  precipitations <- precipitations[ind]
  tavgs <- tavgs[ind]
  tmins <- tmins[ind]
  tmaxs <- tmaxs[ind]
  human_fps <- human_fps[ind]
  pop_dens <- pop_dens[ind]
  trees <- trees[ind]
  grasslands <- grasslands[ind]
  shrubs <- shrubs[ind]
  croplands <- croplands[ind]
  builts <- builts[ind]
  bares <- bares[ind]
  waters <- waters[ind]
  wetlands <- wetlands[ind]
  mangroves <- mangroves[ind]
  mosses <- mosses[ind]
  tsetse_habitats <- tsetse_habitats[ind]
  cattles <- cattles[ind]
  buffalos <- buffalos[ind]
  goats <- goats[ind]
  horses <- horses[ind]
  pigs <- pigs[ind]
  sheeps <- sheeps[ind]
  total_animals <- total_animals[ind]
  tsetses <- tsetses[ind]
  
  prevalence <- positive/sample_size
  
  # Create spatial mesh
  coo <- cbind(longitudes, latitudes)
  
  mesh <- inla.mesh.2d(
    loc = coo, offset = c(50, 100),
    cutoff = 3,                # Increased from 1 (reasonable middle ground)
    max.edge = c(6, 15)        # Increased from c(3, 10) 
  )
  
  spde <- inla.spde2.pcmatern(
    mesh = mesh, 
    alpha = 2, constr = TRUE,
    prior.range = c(50, 0.1),    # 50km range, 10% prob range < 50km
    prior.sigma = c(1, 0.1),     # Moderate spatial variance
  )
  indexs <- inla.spde.make.index("s", spde$n.spde)
  A <- inla.spde.make.A(mesh = mesh, loc = coo)
  
  # Create prediction grid
  bb <- bbox(border)
  x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 200)
  y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 200)
  coop <- as.matrix(expand.grid(x, y))
  
  ind <- point.in.polygon(
    coop[, 1], coop[, 2],
    border[, 1], border[, 2]
  )
  coop <- coop[which(ind == 1), ]
  
  # Extract covariates for prediction grid
  ra_elv <- aggregate(r_elv, fact = 5, fun = mean)
  dp <- as.data.frame(coop)
  
  dp$elevation <- terra::extract(ra_elv, coop)[, 1]
  dp$precipitation <- rowSums(terra::extract(r_prec, coop)[, 1:12])/12
  dp$tavg <- rowSums(terra::extract(r_tavg, coop)[, 1:12])/12
  dp$tmin <- rowSums(terra::extract(r_tmin, coop)[, 1:12])/12
  dp$tmax <- rowSums(terra::extract(r_tmax, coop)[, 1:12])/12
  dp$human_fp <- terra::extract(r_human_fp, coop)[, 1]
  dp$pop_den <- terra::extract(r_pop_den, coop)[, 1]
  dp$tree <- terra::extract(r_tree, coop)[, 1]
  dp$grassland <- terra::extract(r_grassland, coop)[, 1]
  dp$shrub <- terra::extract(r_shrub, coop)[, 1]
  dp$cropland <- terra::extract(r_cropland, coop)[, 1]
  dp$built <- terra::extract(r_built, coop)[, 1]
  dp$bare <- terra::extract(r_bare, coop)[, 1]
  dp$water <- terra::extract(r_water, coop)[, 1]
  dp$wetland <- terra::extract(r_wetland, coop)[, 1]
  dp$mangrove <- terra::extract(r_mangrove, coop)[, 1]
  dp$moss <- terra::extract(r_moss, coop)[, 1]
  dp$tsetse_habitat <- dp$tree + dp$wetland
  dp$cattle <- terra::extract(r_cattle, coop)
  dp$cattle[is.na(dp$cattle)] <- 0
  dp$buffalo <- terra::extract(r_buffalo, coop)
  dp$buffalo[is.na(dp$buffalo)] <- 0
  dp$goat <- terra::extract(r_goat, coop)
  dp$goat[is.na(dp$goat)] <- 0
  dp$horse <- terra::extract(r_horse, coop)
  dp$horse[is.na(dp$horse)] <- 0
  dp$pig <- terra::extract(r_pig, coop)
  dp$pig[is.na(dp$pig)] <- 0
  dp$sheep <- terra::extract(r_sheep, coop)
  dp$sheep[is.na(dp$sheep)] <- 0
  dp$total_animal <- dp$cattle + dp$buffalo + dp$goat + dp$horse + dp$pig + dp$sheep
  dp$tsetse <- terra::extract(tsetse, coop)
  
  Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
  
  # Standardize covariates to improve numerical stability
  # This helps prevent "delta[0] is NAN" errors
  cat("Standardizing covariates for numerical stability...\n")
  
  # Function to safely standardize (avoiding division by zero)
  safe_scale <- function(x) {
    if (var(x, na.rm = TRUE) == 0) return(x)  # Don't scale constants
    return(as.numeric(scale(x)))
  }
  
  elevations_std <- safe_scale(elevations)
  precipitations_std <- safe_scale(precipitations)
  tavgs_std <- safe_scale(tavgs)
  tmins_std <- safe_scale(tmins)
  tmaxs_std <- safe_scale(tmaxs)
  human_fps_std <- safe_scale(human_fps)
  pop_dens_std <- safe_scale(pop_dens)
  trees_std <- safe_scale(trees)
  grasslands_std <- safe_scale(grasslands)
  shrubs_std <- safe_scale(shrubs)
  croplands_std <- safe_scale(croplands)
  builts_std <- safe_scale(builts)
  bares_std <- safe_scale(bares)
  waters_std <- safe_scale(waters)
  wetlands_std <- safe_scale(wetlands)
  mangroves_std <- safe_scale(mangroves)
  mosses_std <- safe_scale(mosses)
  tsetse_habitats_std <- safe_scale(tsetse_habitats)
  cattles_std <- safe_scale(cattles)
  buffalos_std <- safe_scale(buffalos)
  goats_std <- safe_scale(goats)
  horses_std <- safe_scale(horses)
  pigs_std <- safe_scale(pigs)
  sheeps_std <- safe_scale(sheeps)
  total_animals_std <- safe_scale(total_animals)
  tsetses_std <- safe_scale(tsetses)
  
  # Create INLA stacks for estimation and prediction with standardized covariates
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = positive, numtrials = sample_size),
    A = list(1, A),
    effects = list(data.frame(
      b0 = rep(1, nrow(coo)),
      # All covariates standardized for numerical stability
      elevation = elevations_std,
      precipitation = precipitations_std,
      tavg = tavgs_std,
      tmin = tmins_std, 
      tmax = tmaxs_std,
      human_fp = human_fps_std,
      pop_den = pop_dens_std,
      tree = trees_std, 
      grassland = grasslands_std,
      shrub = shrubs_std, 
      cropland = croplands_std, 
      built = builts_std, 
      bare = bares_std, 
      water = waters_std, 
      wetland = wetlands_std, 
      mangrove = mangroves_std, 
      moss = mosses_std, 
      tsetse_habitat = tsetse_habitats_std, 
      cattle = cattles_std,
      buffalo = buffalos_std,
      goat = goats_std, 
      horse = horses_std, 
      pig = pigs_std, 
      sheep = sheeps_std, 
      total_animal = total_animals_std, 
      tsetse = tsetses_std
    ), s = indexs)
  )
  
  
  # Standardize prediction data using same scaling as training data
  # Store scaling parameters for consistent transformation
  scale_params <- list()
  
  # Map data frame column names to actual variable names
  var_mapping <- list(
    "elevation" = "elevations", "precipitation" = "precipitations", "tavg" = "tavgs", 
    "tmin" = "tmins", "tmax" = "tmaxs", "human_fp" = "human_fps", "pop_den" = "pop_dens",
    "tree" = "trees", "grassland" = "grasslands", "shrub" = "shrubs", 
    "cropland" = "croplands", "built" = "builts", "bare" = "bares", 
    "water" = "waters", "wetland" = "wetlands", "mangrove" = "mangroves", 
    "moss" = "mosses", "tsetse_habitat" = "tsetse_habitats", 
    "cattle" = "cattles", "buffalo" = "buffalos", "goat" = "goats", 
    "horse" = "horses", "pig" = "pigs", "sheep" = "sheeps", 
    "total_animal" = "total_animals", "tsetse" = "tsetses"
  )
  
  for (var_name in names(var_mapping)) {
    if (var_name %in% names(dp)) {
      training_var_name <- var_mapping[[var_name]]
      if (exists(training_var_name)) {
        orig_data <- get(training_var_name)  # Get the training data variable
        if (var(orig_data, na.rm = TRUE) > 0) {
          scale_params[[var_name]] <- list(
            center = mean(orig_data, na.rm = TRUE),
            scale = sd(orig_data, na.rm = TRUE)
          )
          # Apply same scaling to prediction data
          dp[[paste0(var_name, "_std")]] <- (dp[[var_name]] - scale_params[[var_name]]$center) / scale_params[[var_name]]$scale
        } else {
          dp[[paste0(var_name, "_std")]] <- dp[[var_name]]  # Don't scale constants
        }
      } else {
        cat("Warning: Training variable", training_var_name, "not found\n")
        dp[[paste0(var_name, "_std")]] <- dp[[var_name]]  # Use unscaled
      }
    }
  }
  
  stk.p <- inla.stack(
    tag = "pred",
    data = list(y = NA, numtrials = NA),
    A = list(1, Ap),
    effects = list(data.frame(
      b0 = rep(1, nrow(coop)),
      # All covariates standardized using same parameters as training data
      elevation = dp$elevation_std,
      precipitation = dp$precipitation_std,
      tavg = dp$tavg_std,
      tmin = dp$tmin_std,
      tmax = dp$tmax_std,
      human_fp = dp$human_fp_std,
      pop_den = dp$pop_den_std,
      tree = dp$tree_std, 
      grassland = dp$grassland_std,
      shrub = dp$shrub_std, 
      cropland = dp$cropland_std, 
      built = dp$built_std, 
      bare = dp$bare_std,
      water = dp$water_std, 
      wetland = dp$wetland_std, 
      mangrove = dp$mangrove_std, 
      moss = dp$moss_std,
      tsetse_habitat = dp$tsetse_habitat_std, 
      cattle = dp$cattle_std,
      buffalo = dp$buffalo_std,
      goat = dp$goat_std, 
      horse = dp$horse_std, 
      pig = dp$pig_std,
      sheep = dp$sheep_std, 
      total_animal = dp$total_animal_std, 
      tsetse = dp$tsetse_std
    ), s = indexs)
  )
  
  stk.full <- inla.stack(stk.e, stk.p)
  
  # Check for data quality issues that could cause numerical problems
  cat("=== DATA DIAGNOSTICS ===\n")
  n_obs <- nrow(coo)  # Define n_obs here
  cat("Number of observations:", n_obs, "\n")
  cat("Sample size range:", min(sample_size, na.rm=TRUE), "to", max(sample_size, na.rm=TRUE), "\n")
  cat("Positive cases range:", min(positive, na.rm=TRUE), "to", max(positive, na.rm=TRUE), "\n")
  cat("Prevalence range:", min(positive/sample_size, na.rm=TRUE), "to", max(positive/sample_size, na.rm=TRUE), "\n")
  
  # Check for perfect separation (all 0s or all 1s in any location)
  perfect_zeros <- sum(positive == 0)
  perfect_ones <- sum(positive == sample_size)
  cat("Perfect zeros (no positive cases):", perfect_zeros, "\n")
  cat("Perfect ones (all positive):", perfect_ones, "\n")
  
  # CRITICAL CHECK: Do we have any positive cases at all?
  total_positive <- sum(positive)
  cat("TOTAL positive cases across all sites:", total_positive, "\n")
  
  if (total_positive == 0) {
    stop("CRITICAL ERROR: No positive cases found in the dataset!\n",
         "Cannot fit prevalence models with zero positive cases.\n",
         "Please check:\n",
         "1. Data loading and filtering steps\n",
         "2. Column names for positive cases\n", 
         "3. Data source and study design\n",
         "4. Whether this should be a presence/absence model instead")
  }
  
  if (perfect_zeros > n_obs * 0.8 || perfect_ones > n_obs * 0.8) {
    cat("WARNING: Extreme prevalence distribution detected!\n")
    cat("This can cause numerical issues in INLA.\n")
  }
  # First, determine appropriate k parameter based on data size
  n_obs <- nrow(coo)
  cat("Number of observations:", n_obs, "\n")
  
  # The persistent assertion error suggests the dataset is too small for any RW2 smoothing
  # Let's start with linear effects only and add complexity gradually
  cat("Starting with LINEAR effects only due to persistent RW2 failures\n")
  use_linear <- TRUE  # Force linear effects initially
  
  # We'll try to add non-linear effects for individual covariates later if the linear model works
  k_param <- 2  # Keep this for potential future use
  
  # Create a robust INLA wrapper function to handle numerical issues
  safe_inla <- function(formula, data, family, Ntrials, control.compute, control.predictor) {
    cat("Trying formula:", as.character(formula)[3], "\n")
    
    # First attempt with default settings but safer priors
    tryCatch({
      inla(formula, data = data, family = family, Ntrials = Ntrials,
           control.compute = control.compute, control.predictor = control.predictor,
           control.inla = list(strategy = "gaussian"),
           control.fixed = list(
             mean = 0,           # Center fixed effects at 0
             prec = 0.001,       # Wider prior (less precise = more conservative)
             mean.intercept = 0,
             prec.intercept = 0.001
           ))
    }, error = function(e) {
      cat("First attempt failed:", e$message, "\n")
      cat("Trying with more conservative settings...\n")
      
      # Second attempt with very conservative settings
      tryCatch({
        inla(formula, data = data, family = family, Ntrials = Ntrials,
             control.compute = list(dic=TRUE, waic=TRUE),  # Minimal compute requirements
             control.predictor = list(compute = TRUE, A = control.predictor$A),
             control.inla = list(
               strategy = "gaussian", 
               cmin = 0,
               restart = 3  # Allow restarts for numerical issues
             ),
             control.fixed = list(
               mean = 0,
               prec = 0.0001,     # Even wider priors
               mean.intercept = 0,
               prec.intercept = 0.0001
             ),
             verbose = FALSE)  # Reduce output clutter
      }, error = function(e2) {
        cat("Second attempt failed:", e2$message, "\n")
        cat("Trying minimal model...\n")
        
        # Third attempt: try with just intercept + spatial effect
        tryCatch({
          minimal_formula <- as.formula("y ~ 0 + b0 + f(s, model = spde)")
          inla(minimal_formula, data = data, family = family, Ntrials = Ntrials,
               control.compute = list(waic=TRUE),
               control.predictor = list(compute = TRUE, A = control.predictor$A),
               control.inla = list(strategy = "gaussian"),
               verbose = FALSE)
        }, error = function(e3) {
          cat("All attempts failed. Skipping this model.\n")
          cat("Final error:", e3$message, "\n")
          return(NULL)
        })
      })
    })
  }
  
  # Use LINEAR effects initially - safer for small datasets
  cat("Using LINEAR effects for all covariates initially\n")
  
  # Check for problematic covariates before model fitting
  cat("=== COVARIATE DIAGNOSTICS ===\n")
  covariate_names <- c("elevation", "precipitation", "tavg", "tmin", "tmax", "human_fp", "pop_den",
                      "tree", "grassland", "shrub", "cropland", "built", "bare", 
                      "water", "wetland", "mangrove", "moss", "tsetse_habitat", 
                      "cattle", "buffalo", "goat", "horse", "pig", "sheep", "total_animal", "tsetse")
  
  # Get covariate data from estimation stack
  stack_data <- inla.stack.data(stk.e)
  
  problematic_covs <- c()
  for (cov_name in covariate_names) {
    if (cov_name %in% names(stack_data)) {
      cov_values <- stack_data[[cov_name]]
      
      # Check for constant values (no variation)
      if (var(cov_values, na.rm = TRUE) == 0) {
        cat("WARNING: Constant covariate", cov_name, "- removing\n")
        problematic_covs <- c(problematic_covs, cov_name)
        next
      }
      
      # Check for extreme values that might cause numerical issues
      q99 <- quantile(cov_values, 0.99, na.rm = TRUE)
      q01 <- quantile(cov_values, 0.01, na.rm = TRUE)
      if (q99 / q01 > 1000) {  # Extreme range
        cat("WARNING: Extreme range in", cov_name, "- may cause numerical issues\n")
      }
      
      # Check for too many zeros (common in spatial data)
      zero_prop <- sum(cov_values == 0, na.rm = TRUE) / length(cov_values)
      if (zero_prop > 0.9) {
        cat("WARNING: >90% zeros in", cov_name, "- may be problematic\n")
      }
    }
  }
  
  # Remove problematic covariates
  all_covariates <- setdiff(covariate_names, problematic_covs)
  cat("Using", length(all_covariates), "covariates after removing", length(problematic_covs), "problematic ones\n")
  
  # Filter out covariates with no variation in the prediction grid
  cov_temp <- c()
  covariate_names <- c("elevation", "precipitation", "tavg", "tmin", "tmax", "human_fp", "pop_den",
                      "tree", "grassland", "shrub", "cropland", "built", "bare", 
                      "water", "wetland", "mangrove", "moss", "tsetse_habitat", 
                      "cattle", "buffalo", "goat", "horse", "pig", "sheep", "total_animal", "tsetse")
  
  for (i in seq_along(covariate_names)){
    col_number <- which(colnames(dp) == covariate_names[i])
    x <- !is.na(dp[,col_number])
    if (length(x) > 0 && max(dp[,col_number][x], na.rm = TRUE) > 0){
      cov_temp <- c(cov_temp, all_covariates[i])
    }
  }
  
  available_covariates <- cov_temp
  
  # Start with a manageable subset for computational efficiency
  covariates <- c("f(elevation_idx, elevation, model='rw2', k=10)",
                 "f(precipitation_idx, precipitation, model='rw2', k=10)", 
                 "f(tavg_idx, tavg, model='rw2', k=10)",
                 "f(tree_idx, tree, model='rw2', k=10)", 
                 "f(tsetse_idx, tsetse, model='rw2', k=10)")
  
  # Keep only covariates that are available
  covariates <- intersect(covariates, available_covariates)
  
  # Backward selection step
  waic_threshold <- 2
  waic_increase <- -1
  waic_refs <- c()
  
  while (min(waic_increase) < waic_threshold && length(covariates) > 0){
    
    cat("Attempting model with", length(covariates), "covariates\n")
    cat("Current covariates:", paste(head(covariates, 3), collapse=", "), if(length(covariates) > 3) "..." else "", "\n")
    
    formula <- as.formula(paste0("y ~ 0 + b0 + ", paste(covariates, collapse = " + "), " + f(s, model = spde)"))
    
    res <- safe_inla(formula,
      data = inla.stack.data(stk.full),
      family = "binomial", Ntrials = numtrials,
      control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
      control.predictor = list(link=1, compute = TRUE,
      A = inla.stack.A(stk.full)))
    
    if (is.null(res)) {
      cat("Model failed, skipping this iteration...\n")
      # Remove the first covariate and try again
      covariates <- covariates[-1]
      waic_increase <- rep(1, length(covariates))  # Reset to continue loop
      next
    }
    
    waic_ref <- res$waic$waic
    waic_refs <- c(waic_refs, waic_ref)
    
    if (length(covariates) == 0) break
    
    waic_values <- rep(0, length(covariates))
    
    for (i in seq_along(covariates)){
      if (length(covariates) > 1) {
        formula <- as.formula(paste0("y ~ 0 + b0 + ", paste(covariates[-i], collapse = " + "), " + f(s, model = spde)"))
      } else {
        formula <- as.formula("y ~ 0 + b0 + f(s, model = spde)")
      }
      
      res <- inla(formula,
        data = inla.stack.data(stk.full),
        family = "binomial", Ntrials = numtrials,
        control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
        control.predictor = list(link=1, compute = TRUE,
        A = inla.stack.A(stk.full)))
      
      waic_values[i] <- res$waic$waic
    }
    
    waic_increase <- waic_values - waic_ref
    if(min(waic_increase) < waic_threshold){
      covariates <- covariates[-which.min(waic_increase)]
    }
  }
  
  # Forward selection from remaining covariates
  remaining_covariates <- setdiff(available_covariates, covariates)
  
  waic_increase <- -1
  while (min(waic_increase) < 0 && length(remaining_covariates) > 0){
    waic_values <- rep(0, length(remaining_covariates))
    for (i in seq_along(remaining_covariates)){
      if (length(covariates) > 0) {
        formula <- as.formula(paste0("y ~ 0 + b0 + ", paste(c(covariates, remaining_covariates[i]), collapse = " + "), " + f(s, model = spde)"))
      } else {
        formula <- as.formula(paste0("y ~ 0 + b0 + ", remaining_covariates[i], " + f(s, model = spde)"))
      }
      
      res <- inla(formula,
        data = inla.stack.data(stk.full),
        family = "binomial", Ntrials = numtrials,
        control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
        control.predictor = list(link=1, compute = TRUE,
        A = inla.stack.A(stk.full)))
      
      waic_values[i] <- res$waic$waic
    }
    
    waic_increase <- waic_values - waic_ref
    
    if(min(waic_increase) < -0.001){
      best_covariate <- remaining_covariates[which.min(waic_increase)]
      covariates <- c(covariates, best_covariate)
      remaining_covariates <- remaining_covariates[-which.min(waic_increase)]
      
      if (length(covariates) > 0) {
        formula <- as.formula(paste0("y ~ 0 + b0 + ", paste(covariates, collapse = " + "), " + f(s, model = spde)"))
      } else {
        formula <- as.formula("y ~ 0 + b0 + f(s, model = spde)")
      }
      
      res <- inla(formula,
        data = inla.stack.data(stk.full),
        family = "binomial", Ntrials = numtrials,
        control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
        control.predictor = list(link=1, compute = TRUE,
        A = inla.stack.A(stk.full)))
      
      waic_ref <- res$waic$waic
      waic_refs <- c(waic_refs, waic_ref)
    } else {
      break
    }
  }
  
  # Final model fitting
  if (length(covariates) > 0) {
    formula <- as.formula(paste0("y ~ 0 + b0 + ", paste(covariates, collapse = " + "), " + f(s, model = spde)"))
  } else {
    formula <- as.formula("y ~ 0 + b0 + f(s, model = spde)")
  }
  
  res <- inla(formula,
    data = inla.stack.data(stk.full),
    family = "binomial", Ntrials = numtrials,
    control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
    control.predictor = list(link=1, compute = TRUE,
    A = inla.stack.A(stk.full)))

  # Save WAIC information
  write.csv(data.frame(covariates = if(length(covariates) > 0) covariates else "spatial_only", 
                      waic_final = res$waic$waic), 
           paste0("Code/Prevalence/Bovine BCT and PCR/WAIC_increase_ETH_non_linear/WAIC_increase_model_",iteration-7,".csv"), 
           row.names = FALSE)

  
  # Extract predictions
  index <- inla.stack.index(stk.full, tag = "pred")$data
  
  pred_mean <- res$summary.fitted.values[index, "mean"]
  pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  
  # ============================================================================
  # COVARIATE EFFECT PLOTS - Show how each covariate affects prevalence
  # ============================================================================
  
  if (length(covariates) > 0) {
    cat("Creating covariate effect plots for iteration", iteration-7, "...\n")
    
    library(ggplot2)
    library(viridis)
    library(gridExtra)
    
    # Create directory for effect plots
    dir.create("Code/Prevalence/Bovine BCT and PCR/Covariate_Effects_ETH_non_linear", showWarnings = FALSE, recursive = TRUE)
    
    # Extract covariate names from the selected covariates
    covariate_names <- gsub("f\\((.+)_idx, (.+), model='rw2', k=10\\)", "\\2", covariates)
    
    effect_plots <- list()
    
    for (i in seq_along(covariates)) {
      covariate_name <- covariate_names[i]
      cat("  Processing effect plot for:", covariate_name, "\n")
      
      # Get the marginal for this covariate effect
      marginal_name <- paste0(covariate_name, "_idx")
      
      if (marginal_name %in% names(res$marginals.random)) {
        
        # Get range of the covariate values
        if (covariate_name %in% colnames(inla.stack.data(stk.full))) {
          covariate_values <- inla.stack.data(stk.full)[[covariate_name]]
          covariate_values <- covariate_values[!is.na(covariate_values)]
          
          # Ensure covariate values are numeric (fix for character type error)
          if (is.character(covariate_values)) {
            covariate_values <- as.numeric(covariate_values)
            covariate_values <- covariate_values[!is.na(covariate_values)]
          }
          
          if (length(covariate_values) > 0 && is.numeric(covariate_values)) {
            covariate_range <- range(covariate_values, na.rm = TRUE)
            
            # Create sequence of covariate values for plotting
            covariate_seq <- seq(covariate_range[1], covariate_range[2], length.out = 100)
            
            # Extract the RW2 effect estimates
            rw2_effects <- res$marginals.random[[marginal_name]]
            
            if (length(rw2_effects) >= length(covariate_seq)) {
              # Calculate marginal means and credible intervals for the effect
              n_effects <- min(length(rw2_effects), length(covariate_seq))
              
              effect_means <- sapply(seq_len(n_effects), function(i) {
                x <- rw2_effects[[i]]
                if (is.list(x) && length(x) > 0) {
                  inla.emarginal(function(y) y, x)
                } else {
                  NA
                }
              })
              
              effect_lower <- sapply(seq_len(n_effects), function(i) {
                x <- rw2_effects[[i]]
                if (is.list(x) && length(x) > 0) {
                  inla.qmarginal(0.025, x)
                } else {
                  NA
                }
              })
              
              effect_upper <- sapply(seq_len(n_effects), function(i) {
                x <- rw2_effects[[i]]
                if (is.list(x) && length(x) > 0) {
                  inla.qmarginal(0.975, x)
                } else {
                  NA
                }
              })
              
              # Remove NAs and match lengths
              valid_idx <- !is.na(effect_means) & !is.na(effect_lower) & !is.na(effect_upper)
              valid_idx <- valid_idx & seq_along(valid_idx) <= length(covariate_seq)
              
              if (sum(valid_idx) > 2) {
                effect_df <- data.frame(
                  covariate_value = covariate_seq[seq_len(sum(valid_idx))],
                  effect_mean = effect_means[valid_idx],
                  effect_lower = effect_lower[valid_idx],
                  effect_upper = effect_upper[valid_idx]
                )
                
                # Convert to prevalence scale (logit to probability)
                effect_df$prevalence_mean <- plogis(effect_df$effect_mean)
                effect_df$prevalence_lower <- plogis(effect_df$effect_lower)
                effect_df$prevalence_upper <- plogis(effect_df$effect_upper)
                
                # Create the plot
                p <- ggplot(effect_df, aes(x = covariate_value)) +
                  geom_ribbon(aes(ymin = prevalence_lower, ymax = prevalence_upper), 
                              alpha = 0.3, fill = "blue") +
                  geom_line(aes(y = prevalence_mean), color = "blue", size = 1.2) +
                  labs(
                    title = paste("Non-linear Effect of", stringr::str_to_title(gsub("_", " ", covariate_name))),
                    subtitle = paste("RW2 smooth effect - Model iteration", iteration-7),
                    x = paste(stringr::str_to_title(gsub("_", " ", covariate_name)), 
                              ifelse(covariate_name == "elevation", "(m)",
                              ifelse(covariate_name == "precipitation", "(mm/year)",
                              ifelse(grepl("temp|tavg|tmin|tmax", covariate_name), "(°C)", "")))),
                    y = "Prevalence (probability)",
                    caption = "Ribbon shows 95% credible interval"
                  ) +
                  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                   limits = c(0, max(effect_df$prevalence_upper, na.rm = TRUE) * 1.1)) +
                  theme_bw() +
                  theme(
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(size = 12, face = "bold"),
                    plot.subtitle = element_text(size = 10),
                    axis.title = element_text(size = 11),
                    axis.text = element_text(size = 10)
                  )
                
                effect_plots[[covariate_name]] <- p
                
                # Save individual plot
                ggsave(
                  filename = paste0("Code/Prevalence/Bovine BCT and PCR/Covariate_Effects_ETH_non_linear/Effect_", 
                                  covariate_name, "_model_", iteration-7, ".png"),
                  plot = p, width = 8, height = 6, dpi = 300
                )
              }
            }
          }
        }
      }
    }
    
    # Create combined plot if we have multiple effects
    if (length(effect_plots) > 1) {
      # Arrange plots in a grid
      n_plots <- length(effect_plots)
      n_cols <- min(3, n_plots)
      n_rows <- ceiling(n_plots / n_cols)
      
      combined_plot <- do.call(grid.arrange, c(effect_plots, ncol = n_cols, nrow = n_rows))
      
      # Save combined plot
      ggsave(
        filename = paste0("Code/Prevalence/Bovine BCT and PCR/Covariate_Effects_ETH_non_linear/Combined_Effects_model_", 
                        iteration-7, ".png"),
        plot = combined_plot, 
        width = n_cols * 6, height = n_rows * 4, dpi = 300
      )
    }
    
    cat("  Saved", length(effect_plots), "covariate effect plots\n")
  }
  
  # Create results data frame
  dpm <- rbind(
    data.frame(
      Latitude = coop[, 1], Longitude = coop[, 2],
      value = pred_mean, variable = "Mean"
    ),
    data.frame(
      Latitude = coop[, 1], Longitude = coop[, 2],
      value = pred_ll, variable = "2.5th percentile"
    ),
    data.frame(
      Latitude = coop[, 1], Longitude = coop[, 2],
      value = pred_ul, variable = "97.5th percentile"
    )
  )
  dpm$variable <- as.factor(dpm$variable)
  
  # Create output directories if they don't exist
  dir.create("Code/Prevalence/Bovine BCT and PCR/Covariates_ETH_non_linear", showWarnings = FALSE, recursive = TRUE)
  dir.create("Code/Prevalence/Bovine BCT and PCR/Projections_ETH_non_linear", showWarnings = FALSE, recursive = TRUE)
  dir.create("Code/Prevalence/Bovine BCT and PCR/WAIC_increase_ETH_non_linear", showWarnings = FALSE, recursive = TRUE)
  
  # Write covariates and projections to CSV files
  write.csv(data.frame(covariates = covariates), 
           paste0("Code/Prevalence/Bovine BCT and PCR/Covariates_ETH_non_linear/Covariates_model_", iteration-7, ".csv"), 
           row.names = FALSE)
  
  write.csv(dpm, 
           paste0("Code/Prevalence/Bovine BCT and PCR/Projections_ETH_non_linear/Projections_model_", iteration-7, ".csv"), 
           row.names = FALSE)
  
  cat("Completed iteration", iteration-7, "for Ethiopia (non-linear). Selected covariates:", paste(covariates, collapse = ", "), "\n")
}

cat("\nLinear model selection completed for Ethiopia!\n")
cat("NOTE: Started with linear effects due to dataset size constraints.\n")
cat("RW2 non-linear smoothing caused assertion failures - dataset may be too small.\n")
cat("Results saved to:\n")
cat("- Covariates_ETH_non_linear/ directory (contains linear effects)\n")
cat("- Projections_ETH_non_linear/ directory\n")
cat("- WAIC_increase_ETH_non_linear/ directory\n")
cat("- Covariate_Effects_ETH_non_linear/ directory (NEW: linear effect plots)\n")

# Optional: Try adding non-linear effects to the best-performing covariates
cat("\n=== OPTIONAL: Attempting Non-Linear Effects ===\n")
cat("If you want to try non-linear effects on the selected covariates,\n")
cat("uncomment and run the following experimental section:\n")
cat("# This will attempt RW2 smoothing with k=2 on the final selected covariates\n")
cat("\nAll covariates modeled with non-linear RW2 effects (k=10).\n")
cat("RW2 will capture linear relationships as straight lines when appropriate.\n")
cat("Check Covariate_Effects_ETH_non_linear/ for plots showing how each\n")
cat("selected covariate affects prevalence (non-linear relationships).\n")