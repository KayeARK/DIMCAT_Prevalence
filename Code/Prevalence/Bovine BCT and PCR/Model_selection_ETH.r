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

countries_to_infer=c("Ethiopia")

# Add string to the name of the file
data_path <- paste0("Code/TestSensSpec/adjusted_cases_all_sites.csv")
data <- read.csv(data_path)

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
for (iteration in 207:1007){
  print(paste("Processing iteration", iteration-7))
  
  # Extract positive cases from original data for this iteration
  positive_orig <- data_original[,iteration]
  
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
  data$elevation <- terra::extract(r_elv, data[, c("longitude", "latitude")])$ETH_elv_msk
  
  # PRECIPITATION - Ethiopia climate data
  r_prec <- worldclim_country(country = "ETH", path = "Data/Covariates", var = "prec")
  data$precipitation <- (terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_1+
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
  
  # Create INLA stacks for estimation and prediction
  stk.e <- inla.stack(
    tag = "est",
    data = list(y = positive, numtrials = sample_size),
    A = list(1, A),
    effects = list(data.frame(b0 = rep(1, nrow(coo)), elevation=elevations,
      precipitation=precipitations, tavg=tavgs, tmin=tmins, tmax=tmaxs,
      human_fp=human_fps, pop_den=pop_dens, tree=trees, grassland=grasslands,
      shrub=shrubs, cropland=croplands, built=builts, bare=bares, water=waters, wetland=wetlands,
      mangrove=mangroves, moss=mosses, tsetse_habitat=tsetse_habitats, cattle=cattles, goat=goats,
      horse=horses, pig=pigs, buffalo=buffalos, sheep=sheeps, total_animal=total_animals, tsetse=tsetses), s = indexs)
  )
  
  stk.p <- inla.stack(
    tag = "pred",
    data = list(y = NA, numtrials = NA),
    A = list(1, Ap),
    effects = list(data.frame(b0 = rep(1, nrow(coop)),elevation=dp$elevation,
      precipitation=dp$precipitation, tavg=dp$tavg, tmin=dp$tmin, tmax=dp$tmax,
      human_fp=dp$human_fp, pop_den=dp$pop_den, tree=dp$tree, grassland=dp$grassland,
      shrub=dp$shrub, cropland=dp$cropland, built=dp$built, bare=dp$bare, water=dp$water, 
      wetland=dp$wetland, mangrove=dp$mangrove, moss=dp$moss, tsetse_habitat=dp$tsetse_habitat, cattle=dp$cattle,
      goat=dp$goat, horse=dp$horse, pig=dp$pig, buffalo=dp$buffalo, sheep=dp$sheep,
       total_animal=dp$total_animal, tsetse=dp$tsetse), s = indexs)
  )
  
  stk.full <- inla.stack(stk.e, stk.p)
  
  # Model selection process
  covariates <- c("elevation","precipitation","human_fp","pop_den","tree","grassland",
    "shrub","cropland","built","bare","water","wetland","cattle","tsetse")
  
  # Filter out covariates with no variation
  cov_temp <- c()
  for (i in 1:length(covariates)){
    col_number <- which(colnames(dp) == covariates[i])
    x <- !is.na(dp[,col_number])
    if (max(dp[,col_number][x]) > 0){
      cov_temp <- c(cov_temp, covariates[i])
    }
  }
  
  covariates <- cov_temp
  full_covariates <- covariates
  
  # Backward selection step
  waic_threshold <- 2
  waic_increase <- -1
  waic_refs <- c()
  
  while (min(waic_increase) < waic_threshold && length(covariates) > 0){
    
    formula <- as.formula(paste0("y ~ 0 + b0 +", paste(covariates, collapse = " + "), "+ f(s, model = spde)"))
    
    res <- inla(formula,
      data = inla.stack.data(stk.full),
      family = "binomial", Ntrials = numtrials,
      control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
      control.predictor = list(link=1, compute = TRUE,
      A = inla.stack.A(stk.full)),)
    
    waic_ref <- res$waic$waic
    waic_refs <- c(waic_refs, waic_ref)
    
    waic_values <- rep(0, length(covariates))
    
    for (i in 1:length(covariates)){
      formula <- as.formula(paste0("y ~ 0 + b0 +", paste(covariates[-i], collapse = " + "), "+ f(s, model = spde)"))
      
      res <- inla(formula,
        data = inla.stack.data(stk.full),
        family = "binomial", Ntrials = numtrials,
        control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
        control.predictor = list(link=1, compute = TRUE,
        A = inla.stack.A(stk.full)),)
      
      waic_values[i] <- res$waic$waic
    }
    
    waic_increase <- waic_values - waic_ref
    if(min(waic_increase) < waic_threshold){
      covariates <- covariates[-which.min(waic_increase)]
    }
  }
  
  # Forward selection step
  if (length(covariates) > 0){
    for (i in 1:length(covariates)){
      full_covariates <- full_covariates[-which(full_covariates==covariates[i])]
    }
  }
  
  waic_increase <- -1
  while (min(waic_increase) < 0){
    waic_values <- rep(0, length(full_covariates))
    for (i in 1:length(full_covariates)){
      formula <- as.formula(paste0("y ~ 0 + b0 +", paste(c(covariates,full_covariates[i]), collapse = " + "), "+ f(s, model = spde)"))
      res <- inla(formula,
        data = inla.stack.data(stk.full),
        family = "binomial", Ntrials = numtrials,
        control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
        control.predictor = list(link=1, compute = TRUE,
        A = inla.stack.A(stk.full)),)
      
      waic_values[i] <- res$waic$waic
    }
    
    waic_increase <- waic_values - waic_ref
    
    if((min(waic_increase)) < -0.001){
      covariates <- c(covariates, full_covariates[which.min(waic_increase)])
    }
    full_covariates <- full_covariates[-which.min(waic_increase)]
    
    formula <- as.formula(paste0("y ~ 0 + b0 +", paste(covariates, collapse = " + "), "+ f(s, model = spde)"))
    res <- inla(formula,
      data = inla.stack.data(stk.full),
      family = "binomial", Ntrials = numtrials,
      control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
      control.predictor = list(link=1, compute = TRUE,
      A = inla.stack.A(stk.full)),)
    
    waic_ref <- res$waic$waic
    waic_refs <- c(waic_refs, waic_ref)
  }
  
  # Final removal step
  cov_temp <- c()
  for (i in 1:length(covariates)){
    col_number <- which(colnames(dp) == covariates[i])
    x <- !is.na(dp[,col_number])
    if (max(dp[,col_number][x]) > 0){
      cov_temp <- c(cov_temp, covariates[i])
    }
  }
  
  covariates <- cov_temp
  full_covariates <- covariates
  
  waic_threshold <- 0.01
  waic_increase <- -0.001
  
  while (min(waic_increase) < waic_threshold && length(covariates) > 0){
    
    formula <- as.formula(paste0("y ~ 0 + b0 +", paste(covariates, collapse = " + "), "+ f(s, model = spde)"))
    
    res <- inla(formula,
      data = inla.stack.data(stk.full),
      family = "binomial", Ntrials = numtrials,
      control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
      control.predictor = list(link=1, compute = TRUE,
      A = inla.stack.A(stk.full)),)
    
    waic_ref <- res$waic$waic
    waic_refs <- c(waic_refs, waic_ref)
    
    waic_values <- rep(0, length(covariates))
    
    for (i in 1:length(covariates)){
      formula <- as.formula(paste0("y ~ 0 + b0 +", paste(covariates[-i], collapse = " + "), "+ f(s, model = spde)"))
      
      res <- inla(formula,
        data = inla.stack.data(stk.full),
        family = "binomial", Ntrials = numtrials,
        control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
        control.predictor = list(link=1, compute = TRUE,
        A = inla.stack.A(stk.full)),)
      
      waic_values[i] <- res$waic$waic
    }
    
    waic_increase <- waic_values - waic_ref
    if(min(waic_increase) < waic_threshold){
      covariates <- covariates[-which.min(waic_increase)]
    }
  }

  write.csv(data.frame(covariates,waic_increase), paste0("Code/Prevalence/Bovine BCT and PCR/WAIC_increase_ETH/WAIC_increase_model_",iteration-7,".csv"), row.names = FALSE)

  
  # Extract predictions
  index <- inla.stack.index(stk.full, tag = "pred")$data
  
  pred_mean <- res$summary.fitted.values[index, "mean"]
  pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  
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
  dir.create("Code/Prevalence/Bovine BCT and PCR/Covariates_ETH", showWarnings = FALSE, recursive = TRUE)
  dir.create("Code/Prevalence/Bovine BCT and PCR/Projections_ETH", showWarnings = FALSE, recursive = TRUE)
  
  # Write covariates and projections to CSV files
  write.csv(data.frame(covariates), 
           paste0("Code/Prevalence/Bovine BCT and PCR/Covariates_ETH/Covariates_model_", iteration-7, ".csv"), 
           row.names = FALSE)
  
  write.csv(dpm, 
           paste0("Code/Prevalence/Bovine BCT and PCR/Projections_ETH/Projections_model_", iteration-7, ".csv"), 
           row.names = FALSE)
  
  cat("Completed iteration", iteration-7, "for Ethiopia. Selected covariates:", paste(covariates, collapse = ", "), "\n")
}

cat("\nModel selection completed for Ethiopia!\n")
cat("Results saved to:\n")
cat("- Covariates_ETH/ directory\n")
cat("- Projections_ETH/ directory\n")