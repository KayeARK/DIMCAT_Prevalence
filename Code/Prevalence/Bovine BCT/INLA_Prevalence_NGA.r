library(readxl)
library(ggplot2)
library(INLA)
library(sp)
library(sf)
library(afrilearndata)
library(concaveman)
library(geodata)
library(terra)
library(raster)
library(see)
library(zoom)
library(ggpattern)
library(ggnewscale)
library(viridis)
library(gridExtra)
library(grid)
library(gtable)
library(ggspatial)
library(cowplot)
library(dplyr)
library(tidyr)

countries_to_infer=c("Nigeria")

covariates <- read.csv(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"Covariates.csv"))
#change to array
covariates <- covariates[,1]


#add string to the name of the file

data_path <- paste0("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
data <- read_excel(data_path)

#remove row if longitude or latitude is NA
data <- data[!is.na(data$Longitude) & !is.na(data$Latitude),]

#if number of infections is NA, set to Number_of_animals_tested*TPR
data$Number_of_infections[is.na(data$Number_of_infections)] <- round(data$Number_of_animal_tested[is.na(data$Number_of_infections)]*data$TPR[is.na(data$Number_of_infections)]/100)

positive <- round(data$Number_of_infections)
sample_size <- data$Number_of_animal_tested

positive[positive > sample_size] <- sample_size[positive > sample_size]

#add elevation data


#ELEVATION
r_elv <- elevation_30s(country = "NGA", path = "Data/Covariates")
data$elevation<-terra::extract(r_elv, data[, c("Longitude", "Latitude")])$NGA_elv_msk

#PRECIPITATION
r_prec<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "prec")
data$precipitation<-(terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_1+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_2+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_3+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_4+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_5+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_6+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_7+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_8+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_9+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_10+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_11+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_prec_12)/12


#AVERAGE TEMPERATURE
r_tavg<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tavg")
data$tavg<-(terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_1+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_2+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_3+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_4+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_5+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_6+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_7+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_8+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_9+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_10+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_11+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tavg_12)/12

#MINIMUM TEMPERATURE
r_tmin<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tmin")
data$tmin<-(terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_1+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_2+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_3+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_4+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_5+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_6+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_7+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_8+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_9+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_10+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_11+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmin_12)/12

#MAXIMUM TEMPERATURE
r_tmax<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tmax")
data$tmax<-(terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_1+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_2+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_3+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_4+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_5+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_6+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_7+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_8+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_9+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_10+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_11+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$NGA_wc2.1_30s_tmax_12)/12

#POPULATION DENSITY
r_pop_den <- population(year=2010, res=10, path="Data/Covariates")
data$pop_den <-terra::extract(r_pop_den, data[, c("Longitude", "Latitude")])$population_density

#LANDCOVER
r_tree <- landcover("trees",path="Data/Covariates")
data$tree <- terra::extract(r_tree, data[, c("Longitude", "Latitude")])$trees

r_grassland <- landcover("grassland",path="Data/Covariates")
data$grassland <- terra::extract(r_grassland, data[, c("Longitude", "Latitude")])$grassland

r_shrub <- landcover("shrubs",path="Data/Covariates")
data$shrub <- terra::extract(r_shrub, data[, c("Longitude", "Latitude")])$shrub

r_cropland <- landcover("cropland",path="Data/Covariates")
data$cropland <- terra::extract(r_cropland, data[, c("Longitude", "Latitude")])$cropland

r_built <- landcover("built",path="Data/Covariates")
data$built <- terra::extract(r_built, data[, c("Longitude", "Latitude")])$built

r_bare <- landcover("bare",path="Data/Covariates")
data$bare <- terra::extract(r_bare, data[, c("Longitude", "Latitude")])$bare

r_water <- landcover("water",path="Data/Covariates")
data$water <- terra::extract(r_water, data[, c("Longitude", "Latitude")])$water

r_wetland <- landcover("wetland",path="Data/Covariates")
data$wetland <- terra::extract(r_wetland, data[, c("Longitude", "Latitude")])$wetland

r_mangrove <- landcover("mangroves",path="Data/Covariates")
data$mangrove <- terra::extract(r_mangrove, data[, c("Longitude", "Latitude")])$mangroves

r_moss <- landcover("moss",path="Data/Covariates")
data$moss <- terra::extract(r_moss, data[, c("Longitude", "Latitude")])$moss

data$tsetse_habitat <- data$tree + data$wetland

#HUMAN FOOTPRINT
r_human_fp<-footprint(year=2009, res=30, path="Data/Covariates")
data$human_fp<-terra::extract(r_human_fp, data[, c("Longitude", "Latitude")])[,2]


#CATTLE DATA
r_cattle<-raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
data$cattle<-terra::extract(r_cattle, data[, c("Longitude", "Latitude")])
data$cattle[is.na(data$cattle)] <- 0

r_buffalo<-raster("Data/Covariates/livestock/buffalo_2015/5_Bf_2015_Da.tif")
data$buffalo<-terra::extract(r_buffalo, data[, c("Longitude", "Latitude")])
data$buffalo[is.na(data$buffalo)] <- 0

r_goat<-raster("Data/Covariates/livestock/goats_2015/5_Gt_2015_Da.tif")
data$goat<-terra::extract(r_goat, data[, c("Longitude", "Latitude")])
data$goat[is.na(data$goat)] <- 0

r_horse<-raster("Data/Covariates/livestock/horses_2015/5_Ho_2015_Da.tif")
data$horse<-terra::extract(r_horse, data[, c("Longitude", "Latitude")])
data$horse[is.na(data$horse)] <- 0

r_pig<-raster("Data/Covariates/livestock/pigs_2015/5_Pg_2015_Da.tif")
data$pig<-terra::extract(r_pig, data[, c("Longitude", "Latitude")])
data$pig[is.na(data$pig)] <- 0

r_sheep<-raster("Data/Covariates/livestock/sheep_2015/5_Sh_2015_Da.tif")
data$sheep<-terra::extract(r_sheep, data[, c("Longitude", "Latitude")])
data$sheep[is.na(data$sheep)] <- 0

tsetse<-raster("Data/Covariates/tsenumbspec")
tsetse[tsetse > 1] <- 1
data$tsetse <-terra::extract(tsetse, data[, c("Longitude", "Latitude")])

data$total_animal<-data$cattle+data$buffalo+data$goat+data$horse+data$pig+data$sheep


latitudes<-data$Latitude
longitudes<-data$Longitude
elevations<-data$elevation
precipitations<-data$precipitation
tavgs<-data$tavg
tmins<-data$tmin
tmaxs<-data$tmax
human_fps<-data$human_fp
pop_dens<-data$pop_den

trees<-data$tree
grasslands<-data$grassland
shrubs<-data$shrub
croplands<-data$cropland
builts<-data$built
bares<-data$bare
waters<-data$water
wetlands<-data$wetland
mangroves<-data$mangrove
mosses<-data$moss
tsetse_habitats<-data$tsetse_habitat

cattles<-data$cattle
buffalos<-data$buffalo
goats<-data$goat
horses<-data$horse
pigs<-data$pig
sheeps<-data$sheep
total_animals<-data$total_animal
tsetses<-data$tsetse


countries_to_infer <- sort(countries_to_infer)

#if Africa is in this list then perform the calculations for the whole continent
if("Africa" %in% countries_to_infer){
  border <- st_coordinates(st_geometry(africontinent))

#remove Madagascar
  border <- border[!(border[,1] > 40 & border[,1] < 60 & border[,2] > -30 & border[,2] < -10),]
  border <- border[,1:2]
  #remove Madagascar from africountries
  africountries <- africountries[!(africountries$name=="Madagascar"),]
  country <- africountries
} else{
v=africountries$name %in% countries_to_infer
country <- africountries[which(v==TRUE),]

}

#remove all countries except the ones in countries_to_infer
africountries_toplot <- africountries[which(africountries$name %in% countries_to_infer),]

country <- st_coordinates(st_geometry(country))

poly.df <- data.frame("long"=country[,1],"lat"=country[,2])
poly.sf <- st_as_sf(poly.df, coords = c("long","lat"))
poly <- concaveman(poly.sf, length_threshold = 0,concavity=1.1)
point.df <- data.frame("long"=longitudes,"lat"=latitudes)
point.sf <- st_as_sf(point.df, coords = c("long","lat"))
long_lat_concave <- st_coordinates(st_geometry(poly))
border <- long_lat_concave

inside <- st_intersects(point.sf,poly)
inside <- inside[]==1

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

#remove NA values
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

# Calculate maximum observed prevalence for colorbar scaling
max_observed_prevalence <- max(prevalence, na.rm = TRUE)
cat("Maximum observed prevalence:", round(max_observed_prevalence, 4), "\n")

# put latitude and longitude next to each other
coo <- cbind(longitudes, latitudes)

mesh <- inla.mesh.2d(
  loc = coo, offset = c(50, 100),
  cutoff = 3,                # Increased from 1 (reasonable middle ground)
  max.edge = c(6, 15)        # Increased from c(3, 10) 
)

#plot the mesh object
#save the plot

spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  alpha = 2, constr = TRUE,
  prior.range = c(50, 0.1),    # 50km range, 10% prob range < 50km
  prior.sigma = c(1, 0.1),     # Moderate spatial variance
)

indexs <- inla.spde.make.index("s", spde$n.spde)

A <- inla.spde.make.A(mesh = mesh, loc = coo)

bb <- bbox(border)
x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 200)
y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 200)
coop <- as.matrix(expand.grid(x, y))

ind <- point.in.polygon(
  coop[, 1], coop[, 2],
  border[, 1], border[, 2]
)
coop <- coop[which(ind == 1), ]

ra_elv <- aggregate(r_elv, fact = 5, fun = mean)
dp <- terra::crds(ra_elv)
dp <- as.data.frame(coop)

#coop <- terra::crds(ra)
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

# stack for estimation stk.e
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


# stack for prediction stk.p
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

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

covariates<-covariates[-3]
formula <- as.formula(paste0("y ~ 0 + b0 +",paste(covariates, collapse = " + "),"+ f(s, model = spde)"))

res <- inla(formula,
  data = inla.stack.data(stk.full),
  family = "binomial", Ntrials = numtrials,
  control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
  control.predictor = list(link=1,compute = TRUE,
    A = inla.stack.A(stk.full)
  ),
)

index <- inla.stack.index(stk.full, tag = "pred")$data

pred_mean <- res$summary.fitted.values[index, "mean"]
pred_ll <- res$summary.fitted.values[index, "0.025quant"]
pred_ul <- res$summary.fitted.values[index, "0.975quant"]

# plot

dpm <- rbind(
  data.frame(
    Latitude = coop[, 2], Longitude = coop[, 1],
    value = pred_mean, variable = "Mean"
  ),
  data.frame(
    Latitude = coop[, 2], Longitude = coop[, 1],
    value = pred_ll, variable = "2.5th percentile"
  ),
  data.frame(
    Latitude = coop[, 2], Longitude = coop[, 1],
    value = pred_ul, variable = "97.5th percentile"
  )
)
dpm$variable <- as.factor(dpm$variable)

# Get neighboring countries for geographic context
cat("Loading neighboring countries for geographic context...\n")

# Define neighboring countries based on target country
if(countries_to_infer == "Nigeria") {
  neighboring_countries <- c("BEN", "NER", "TCD", "CMR")  # Benin, Niger, Chad, Cameroon
} else if(countries_to_infer == "Nigeria") {
  neighboring_countries <- c("ERI", "SDN", "SSD", "KEN", "SOM", "DJI", "UGA")  # Eritrea, Sudan, South Sudan, Kenya, Somalia, Djibouti, Uganda
} else {
  neighboring_countries <- c()  # Default: no neighbors for other countries
  cat("No neighboring countries defined for", countries_to_infer, "\n")
}

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

# Get current country boundary as sf object
country_sf <- NULL
tryCatch({
  country_data <- gadm(country = countries_to_infer, level = 0, path = tempdir())
  country_sf <- st_as_sf(country_data)
  # Validate and repair country geometry
  country_sf <- st_make_valid(country_sf)
  cat("Successfully loaded country boundary for", countries_to_infer, "\n")
}, error = function(e) {
  cat("Failed to load country boundary:", e$message, "\n")
})

# Combine neighboring countries if available
neighboring_countries_sf <- NULL
if(length(neighbor_sf_list) > 0) {
  neighboring_countries_full <- do.call(rbind, neighbor_sf_list)
  
  # Validate and repair geometries to avoid intersection errors
  cat("Validating and repairing geometries...\n")
  neighboring_countries_full <- st_make_valid(neighboring_countries_full)
  
  # Create buffer around target country to crop neighbors
  if(!is.null(country_sf)) {
    # Validate country geometry too
    country_sf <- st_make_valid(country_sf)
    country_buffer <- st_buffer(st_union(country_sf), dist = 1.5)  # ~150km buffer
    country_buffer <- st_make_valid(country_buffer)
    
    # Use st_crop with error handling for more robust clipping
    tryCatch({
      neighboring_countries_sf <- st_crop(neighboring_countries_full, country_buffer)
      cat("Successfully loaded and cropped", nrow(neighboring_countries_sf), "neighboring countries\n")
    }, error = function(e) {
      cat("st_crop failed, trying alternative approach:", e$message, "\n")
      # Fallback: use a simpler bbox-based crop
      country_bbox <- st_bbox(country_buffer)
      neighboring_countries_sf <<- st_crop(neighboring_countries_full, country_bbox)
      cat("Successfully loaded and cropped using bbox approach\n")
    })
  }
}

# Get administrative boundaries based on country
admin_sf <- NULL
admin_level <- NULL
admin_name_col <- NULL

if(countries_to_infer == "Nigeria") {
  # Get Nigerian administrative boundaries at LGA level (level 2)
  admin_level <- 2
  admin_name_col <- "NAME_2"
  tryCatch({
    ethiopia_admin <- gadm(country = "NGA", level = admin_level, path = tempdir())
    admin_sf <- st_as_sf(ethiopia_admin)
    admin_sf <- st_make_valid(admin_sf)
    cat("Successfully loaded Nigeria LGA boundaries\n")
  }, error = function(e) {
    cat("Failed to load Nigeria LGA boundaries:", e$message, "\n")
  })
} else if(countries_to_infer == "Nigeria") {
  # Get Nigerian administrative boundaries at LGA level (level 2)
  admin_level <- 2
  admin_name_col <- "NAME_2"
  tryCatch({
    nigeria_admin <- gadm(country = "NGA", level = admin_level, path = tempdir())
    admin_sf <- st_as_sf(nigeria_admin)
    admin_sf <- st_make_valid(admin_sf)
    cat("Successfully loaded Nigeria LGA boundaries\n")
  }, error = function(e) {
    cat("Failed to load Nigeria LGA boundaries:", e$message, "\n")
  })
}

# Function to aggregate prevalence data by administrative units
aggregate_by_admin <- function(projections_data, variable_name) {
  if(is.null(admin_sf)) {
    cat("No administrative boundaries available, returning original data\n")
    return(projections_data)
  }
  
  # Convert projection data to spatial points
  projections_sf <- st_as_sf(projections_data, 
                            coords = c("Longitude", "Latitude"), 
                            crs = 4326)
  
  # Spatial join with administrative boundaries
  joined_data <- st_join(projections_sf, admin_sf)
  
  # Aggregate by administrative unit
  admin_prevalence <- joined_data %>%
    st_drop_geometry() %>%
    dplyr::filter(!is.na(GID_2)) %>%  # Only keep points within admin units
    dplyr::group_by(GID_2, !!sym(admin_name_col), NAME_1) %>%
    dplyr::summarise(
      admin_prevalence = mean(value, na.rm = TRUE),
      median_prevalence = median(value, na.rm = TRUE),
      min_prevalence = min(value, na.rm = TRUE),
      max_prevalence = max(value, na.rm = TRUE),
      n_points = n(),
      .groups = 'drop'
    )
  
  # Handle missing administrative units with spatial interpolation
  all_admin <- admin_sf %>%
    st_drop_geometry() %>%
    dplyr::select(GID_2, !!sym(admin_name_col), NAME_1)
  
  missing_admin <- all_admin %>%
    dplyr::anti_join(admin_prevalence, by = "GID_2")
  
  if(nrow(missing_admin) > 0) {
    cat("Found", nrow(missing_admin), "administrative units without data. Filling using spatial interpolation...\n")
    
    for(i in seq_len(nrow(missing_admin))) {
      admin_id <- missing_admin$GID_2[i]
      
      # Get centroid of the missing admin unit
      admin_centroid <- admin_sf %>%
        dplyr::filter(GID_2 == admin_id) %>%
        st_centroid(of_largest_polygon = TRUE)
      
      # Find distances to all projection points
      distances <- st_distance(admin_centroid, projections_sf)
      distances_numeric <- as.numeric(distances)
      
      # Get indices of 5 closest points
      n_neighbors <- min(5, nrow(projections_sf))
      closest_indices <- order(distances_numeric)[1:n_neighbors]
      
      # Calculate inverse distance weighted average
      closest_points <- projections_sf[closest_indices, ]
      weights <- 1 / (distances_numeric[closest_indices] + 1e-10)
      weights <- weights / sum(weights)
      
      interpolated_prevalence <- sum(st_drop_geometry(closest_points)$value * weights)
      
      # Add to admin_prevalence data
      new_row <- data.frame(
        GID_2 = admin_id,
        admin_prevalence = interpolated_prevalence,
        median_prevalence = interpolated_prevalence,
        min_prevalence = interpolated_prevalence,
        max_prevalence = interpolated_prevalence,
        n_points = 0
      )
      
      # Add the name columns
      new_row[[admin_name_col]] <- missing_admin[[admin_name_col]][i]
      new_row$NAME_1 <- missing_admin$NAME_1[i]
      
      admin_prevalence <- rbind(admin_prevalence, new_row)
    }
  }
  
  return(admin_prevalence)
}

# Aggregate all three datasets by administrative units
dpm_mean <- dpm %>% filter(variable == "Mean")
dpm_lower <- dpm %>% filter(variable == "2.5th percentile") 
dpm_upper <- dpm %>% filter(variable == "97.5th percentile")

admin_prevalence_mean <- aggregate_by_admin(dpm_mean, "Mean")
admin_prevalence_lower <- aggregate_by_admin(dpm_lower, "2.5th percentile")
admin_prevalence_upper <- aggregate_by_admin(dpm_upper, "97.5th percentile")

# Calculate administrative-level uncertainty
admin_uncertainty <- admin_prevalence_lower %>%
  dplyr::select(GID_2, lower = admin_prevalence) %>%
  left_join(admin_prevalence_upper %>% dplyr::select(GID_2, upper = admin_prevalence), by = "GID_2") %>%
  dplyr::mutate(
    ci_width = upper - lower,
    high_uncertainty = ci_width > 0.80
  )

cat("Administrative units with high uncertainty (CI width > 0.80):", 
    sum(admin_uncertainty$high_uncertainty, na.rm = TRUE), "out of", nrow(admin_uncertainty), "\n")

# Function to create administrative-level choropleth map
create_enhanced_map <- function(admin_data, title_suffix, show_uncertainty = FALSE) {
  if(is.null(admin_sf)) {
    cat("No administrative boundaries available, creating pixel-level map\n")
    # Fallback to original pixel-level mapping
    p <- ggplot()
    
    if(!is.null(neighboring_countries_sf)) {
      p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
    }
    
    pixel_data <- dpm %>% filter(variable == gsub("[\\(\\)]", "", title_suffix))
    p <- p + geom_tile(data = pixel_data, aes(Longitude, Latitude, fill = value), alpha = 0.9)
    
    if(!is.null(country_sf)) {
      p <- p + geom_sf(data = country_sf, fill = NA, color = "black", lwd = 0.8)
    }
    
    return(p)
  }
  
  # Join administrative prevalence back to admin polygons
  admin_final <- admin_sf %>%
    left_join(admin_data, by = "GID_2")
  
  # Add uncertainty information if requested
  if(show_uncertainty) {
    admin_final <- admin_final %>%
      left_join(admin_uncertainty, by = "GID_2")
  }
  
  p <- ggplot(admin_final)
  
  # Add neighboring countries as background context
  if(!is.null(neighboring_countries_sf)) {
    p <- p + geom_sf(data = neighboring_countries_sf, fill = "grey95", color = "grey80", lwd = 0.3)
  }
  
  # Add administrative units with prevalence data
  p <- p + 
    geom_sf(aes(fill = admin_prevalence), lwd = 0.1, color = "white") +
    scale_fill_viridis_c(
      name = "AAT\nprevalence",
      option = "viridis",
      direction = 1,
      limits = c(0, max(admin_final$admin_prevalence, na.rm = TRUE)),
      labels = scales::percent_format(accuracy = 1),
      na.value = "grey90"
    )
  
  # Add high uncertainty overlay if requested
  if(show_uncertainty && "high_uncertainty" %in% names(admin_final)) {
    high_uncertainty_admin <- admin_final %>% 
      dplyr::filter(high_uncertainty == TRUE)
    
    if(nrow(high_uncertainty_admin) > 0) {
      p <- p + 
        geom_sf(data = high_uncertainty_admin, fill = "white", alpha = 0.4, 
                color = "red", lwd = 0.3)
    }
  }
  
  # Add country boundary (state/regional level for context)
  if(!is.null(country_sf)) {
    p <- p + geom_sf(data = country_sf, fill = NA, color = "black", lwd = 0.8)
  }
  
  # Add scale bar and north arrow
  p <- p +
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                    bar_cols = c("black", "white"), line_width = 1) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                          pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                          style = north_arrow_fancy_orienteering(text_size = 8)) +
    theme_void() +
    labs(
      title = paste("AAT prevalence by", 
                    ifelse(countries_to_infer == "Nigeria", "LGA", "LGA"), 
                    "in", countries_to_infer, title_suffix)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
    )
  
  # Add uncertainty caption if shown
  if(show_uncertainty && "high_uncertainty" %in% names(admin_final) && 
     any(admin_final$high_uncertainty, na.rm = TRUE)) {
    p <- p + labs(caption = "Areas with white overlay and red border indicate high uncertainty (CI width > 0.80)")
  }
  
  return(p)
}

# Create administrative-level choropleth maps
p_mean <- create_enhanced_map(admin_prevalence_mean, "(mean)", show_uncertainty = FALSE)
p_lower <- create_enhanced_map(admin_prevalence_lower, "(2.5th percentile)", show_uncertainty = TRUE)
p_upper <- create_enhanced_map(admin_prevalence_upper, "(97.5th percentile)", show_uncertainty = TRUE)

# Create combined histogram and summary plots for administrative-level data
create_histogram_plot <- function(admin_data, estimate_name) {
  # Histogram with viridis color scale matching maps
  p_hist <- ggplot(admin_data, aes(x = admin_prevalence, fill = after_stat(x))) +
    geom_histogram(binwidth = 0.04, color = "white", linewidth = 0.2,
                   boundary = 0, closed = "left") +
    scale_fill_viridis_c(option = "viridis", guide = "none", 
                         limits = c(0, max(admin_data$admin_prevalence, na.rm = TRUE))) +
    scale_x_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.25),
                       limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = margin(t = 5, r = 15, b = 0, l = 15),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(
      title = paste("Distribution of", estimate_name, "prevalence by", 
                    ifelse(countries_to_infer == "Nigeria", "LGA", "LGA")),
      x = ""
    )
  
  # Boxplot with matching viridis color
  viridis_color <- viridis::viridis(1, begin = 0.5, end = 0.5, option = "viridis")
  
  p_box <- ggplot(admin_data, aes(x = admin_prevalence, y = 1)) +
    geom_boxplot(fill = viridis_color, alpha = 0.8, width = 0.8) +
    scale_x_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.25),
                       limits = c(0, 1), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0.2, 1.8), expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 5)),
      axis.text.x = element_text(size = 9),
      axis.ticks.x = element_line(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 15, b = 5, l = 15)
    ) +
    labs(
      x = "AAT prevalence"
    )
  
  # Create title and combine plots
  admin_unit <- ifelse(countries_to_infer == "Nigeria", "LGA", "LGA")
  title <- grid::textGrob(paste("AAT prevalence analysis by", admin_unit, "-", estimate_name), 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  combined_plot <- grid.arrange(
    title,
    p_hist,
    p_box,
    nrow = 3,
    heights = unit(c(0.8, 3.5, 2.5), "cm")
  )
  
  return(combined_plot)
}

# Create histogram plots for administrative-level data
p_mean_hist <- create_histogram_plot(admin_prevalence_mean, "mean")
p_lower_hist <- create_histogram_plot(admin_prevalence_lower, "2.5th percentile") 
p_upper_hist <- create_histogram_plot(admin_prevalence_upper, "97.5th percentile")

# Save enhanced maps
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_enhanced_mean.pdf"), 
       plot = p_mean, width = 12, height = 10, dpi = 300)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_enhanced_lower.pdf"), 
       plot = p_lower, width = 12, height = 10, dpi = 300)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_enhanced_upper.pdf"), 
       plot = p_upper, width = 12, height = 10, dpi = 300)

# Save histogram plots
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_hist_mean.pdf"), 
       plot = p_mean_hist, width = 10, height = 6, dpi = 300)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_hist_lower.pdf"), 
       plot = p_lower_hist, width = 10, height = 6, dpi = 300)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_hist_upper.pdf"), 
       plot = p_upper_hist, width = 10, height = 6, dpi = 300)

# Display the enhanced plots
print(p_mean)
print(p_lower)
print(p_upper)

# Display histogram plots
grid.draw(p_mean_hist)
grid.draw(p_lower_hist)
grid.draw(p_upper_hist)

# Print summary statistics for administrative-level estimates
admin_unit <- ifelse(countries_to_infer == "Nigeria", "LGA", "LGA")
cat("=== AAT PREVALENCE ESTIMATES SUMMARY BY", toupper(admin_unit), "===\n\n")

cat("Mean Prevalence:\n")
cat("Total", admin_unit, "units with data:", nrow(admin_prevalence_mean), "\n")
if(nrow(admin_prevalence_mean) > 0) {
  print(summary(admin_prevalence_mean$admin_prevalence))
}

cat("\n2.5th Percentile Prevalence:\n") 
cat("Total", admin_unit, "units with data:", nrow(admin_prevalence_lower), "\n")
if(nrow(admin_prevalence_lower) > 0) {
  print(summary(admin_prevalence_lower$admin_prevalence))
}

cat("\n97.5th Percentile Prevalence:\n")
cat("Total", admin_unit, "units with data:", nrow(admin_prevalence_upper), "\n")
if(nrow(admin_prevalence_upper) > 0) {
  print(summary(admin_prevalence_upper$admin_prevalence))
}

cat("\nUncertainty Analysis:\n")
cat("Total", admin_unit, "units:", nrow(admin_uncertainty), "\n")
cat("High uncertainty", admin_unit, "units (CI width > 0.80):", sum(admin_uncertainty$high_uncertainty, na.rm = TRUE), "\n")
cat("Percentage with high uncertainty:", 
    round(100 * sum(admin_uncertainty$high_uncertainty, na.rm = TRUE) / nrow(admin_uncertainty), 1), "%\n")

# Also save the original combined plot for comparison
ggplot(dpm) + geom_tile(aes(Longitude, Latitude, fill = value)) + 
  geom_polygon(data = data.frame(X = border[,1], Y = border[,2]), aes(x = X, y = Y), color = "black", fill = NA) +
  facet_grid(~factor(variable, levels = c("2.5th percentile","Mean","97.5th percentile"))) +
  coord_fixed(ratio = 1) +
  scale_fill_viridis_c(
    name = "Prevalence",
    option = "viridis",
    limits = c(0, max(dpm$value, na.rm = TRUE))
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) +
  labs(title = paste("AAT prevalence in", countries_to_infer, "- All estimates"))

ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_enhanced_combined.png"), width = 15, height = 6, dpi = 300)


######### EVALUATING THE MODEL ##########

# Get the indices for the estimation data (observations)
index_est <- inla.stack.index(stk.full, tag = "est")$data

# CRITICAL DATA ALIGNMENT VERIFICATION - TRIPLE CHECK
cat("\n=== TRIPLE VERIFICATION: RAW DATA vs MODEL EXTRACTION ===\n")

# Step 1: Track the EXACT original coordinates and prevalences
cat("Step 1: Original raw data before ANY processing\n")
original_coords <- cbind(longitudes, latitudes)  # These are the variables after filtering but before model
original_prevalence <- positive/sample_size  # Raw prevalence calculation
original_sample_size <- sample_size  # Sample sizes

cat("Original data summary:\n")
cat("  N observations:", length(original_prevalence), "\n")
cat("  Coordinate range: Lon [", round(min(original_coords[,1]), 4), ",", round(max(original_coords[,1]), 4), 
    "], Lat [", round(min(original_coords[,2]), 4), ",", round(max(original_coords[,2]), 4), "]\n")

# Step 2: Verify `coo` matches exactly
cat("\nStep 2: Verify coordinate alignment\n")
cat("coo dimensions:", dim(coo), "\n")
cat("original_coords dimensions:", dim(original_coords), "\n")

if(nrow(coo) != nrow(original_coords)) {
  warning("COORDINATE MISMATCH: coo and original_coords have different sizes!")
}

# Check if coordinates are identical
coords_identical <- all.equal(coo, original_coords)
if(isTRUE(coords_identical)) {
  cat("✓ Coordinates are IDENTICAL\n")
} else {
  warning("✗ COORDINATES ARE DIFFERENT!")
  cat("Max coordinate difference:", max(abs(coo - original_coords)), "\n")
}

# Step 3: Verify prevalence alignment  
cat("\nStep 3: Verify prevalence alignment\n")
cat("Current prevalence length:", length(prevalence), "\n")
cat("Original prevalence length:", length(original_prevalence), "\n")

if(length(prevalence) == length(original_prevalence)) {
  prevalence_identical <- all.equal(prevalence, original_prevalence)
  if(isTRUE(prevalence_identical)) {
    cat("✓ Prevalences are IDENTICAL\n")
  } else {
    warning("✗ PREVALENCES ARE DIFFERENT!")
    cat("Max prevalence difference:", max(abs(prevalence - original_prevalence)), "\n")
  }
} else {
  warning("✗ PREVALENCE VECTORS HAVE DIFFERENT LENGTHS!")
}

# Step 4: CRITICAL - Verify INLA extraction matches coordinates
cat("\nStep 4: CRITICAL - Verify model extraction alignment\n")
cat("INLA stack structure:\n")
cat("  Total fitted values:", length(res$summary.fitted.values[,1]), "\n")
cat("  Estimation indices range:", min(index_est), "to", max(index_est), "\n")
cat("  Number of estimation indices:", length(index_est), "\n")
cat("  Number of observations (coo):", nrow(coo), "\n")

if(length(index_est) != nrow(coo)) {
  warning("FATAL ERROR: Estimation indices don't match observation count!")
  cat("This means we're extracting predictions for wrong locations!\n")
}

# Step 5: COORDINATE VERIFICATION PLOTS
cat("\nStep 5: Creating coordinate verification plots...\n")

# Extract the coordinates used to build the INLA mesh/stack
stack_coords_est <- coo  # These should be the coordinates for estimation
  
# For comparison, get prediction coordinates if they exist
index_pred <- inla.stack.index(stk.full, tag = "pred")$data
if(length(index_pred) > 0) {
  # Get prediction grid coordinates (these are different from observation coordinates)
  stack_coords_pred <- coop[1:min(100, nrow(coop)), ]  # Just first 100 for visualization
} else {
  stack_coords_pred <- matrix(NA, 0, 2)
}

# Plot 1: Verify estimation coordinates match original coordinates EXACTLY
p_coord_verify <- ggplot() +
  geom_point(data = data.frame(x = original_coords[,1], y = original_coords[,2]), 
             aes(x = x, y = y), color = "#E74C3C", size = 2, alpha = 0.7, shape = 16) +
  geom_point(data = data.frame(x = stack_coords_est[,1], y = stack_coords_est[,2]), 
             aes(x = x, y = y), color = "#3498DB", size = 1.5, alpha = 0.7, shape = 1) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste("Coordinate Verification -", countries_to_infer),
    x = "Longitude", 
    y = "Latitude",
    subtitle = "Red dots: Original coordinates | Blue circles: Model coordinates",
    caption = "Points should overlap perfectly if alignment is correct"
  )

print(p_coord_verify)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_Coordinate_Verification.png"),
       plot = p_coord_verify, width = 12, height = 8, dpi = 300)

# Plot 2: Direct coordinate comparison (should be perfect diagonal line)
coord_comparison_df <- data.frame(
  orig_lon = original_coords[,1],
  orig_lat = original_coords[,2], 
  model_lon = stack_coords_est[,1],
  model_lat = stack_coords_est[,2]
)

p_coord_diagonal <- ggplot(coord_comparison_df) +
  geom_point(aes(x = orig_lon, y = model_lon), color = "#E74C3C", size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste("Longitude Alignment Check -", countries_to_infer),
    x = "Original Longitude", 
    y = "Model Longitude",
    subtitle = "Points should lie exactly on diagonal line",
    caption = "Perfect alignment = perfect diagonal"
  )

print(p_coord_diagonal)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_Longitude_Alignment.png"),
       plot = p_coord_diagonal, width = 10, height = 10, dpi = 300)

p_coord_diagonal_lat <- ggplot(coord_comparison_df) +
  geom_point(aes(x = orig_lat, y = model_lat), color = "#3498DB", size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste("Latitude Alignment Check -", countries_to_infer),
    x = "Original Latitude", 
    y = "Model Latitude", 
    subtitle = "Points should lie exactly on diagonal line",
    caption = "Perfect alignment = perfect diagonal"
  )

print(p_coord_diagonal_lat)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_Latitude_Alignment.png"),
       plot = p_coord_diagonal_lat, width = 10, height = 10, dpi = 300)

# Quantitative coordinate alignment check
lon_diff <- abs(original_coords[,1] - stack_coords_est[,1])
lat_diff <- abs(original_coords[,2] - stack_coords_est[,2])

cat("Coordinate alignment statistics:\n")
cat("  Max longitude difference:", round(max(lon_diff), 8), "\n")
cat("  Max latitude difference:", round(max(lat_diff), 8), "\n") 
cat("  Mean longitude difference:", round(mean(lon_diff), 8), "\n")
cat("  Mean latitude difference:", round(mean(lat_diff), 8), "\n")

if(max(lon_diff) > 1e-6 || max(lat_diff) > 1e-6) {
  warning("SIGNIFICANT COORDINATE DIFFERENCES DETECTED!")
  cat("This indicates the model predictions are NOT from the same locations as the raw data!\n")
} else {
  cat("✓ Coordinates align within numerical precision\n")
}

# Step 6: POINT-BY-POINT DATA VERIFICATION (MOVED TO AFTER MODES CALCULATION)
# This section was moved to after modes are computed to avoid "object not found" errors

# Step 7: MODEL PERFORMANCE CHECK (MOVED TO AFTER MODES CALCULATION)
# Step 7 content moved to after modes calculation
# cat("\nStep 7: Model performance analysis\n")
# cat("Variable length debugging:\n")
# cat("  Length of prevalence:", length(prevalence), "\n")
# cat("  Length of modes:", length(modes), "\n")
# cat("  Length of original_prevalence:", length(original_prevalence), "\n")
# cat("  Length of original_coords:", nrow(original_coords), "\n")
# Entire Step 7 section commented out - will be restored after modes calculation
# All Step 7 analysis code moved to after modes calculation
# This prevents "object 'modes' not found" errors

# Step 8: CI WIDTH ANALYSIS (MOVED TO AFTER VARIABLE DEFINITION)
# All Step 8 analysis code moved to after lower_ci and upper_ci are calculated
# This prevents "object 'lower_ci' not found" errors

# Extract marginals for the observations directly using the correct approach
# We need to use the INLA built-in functions for proper credible intervals
modes <- numeric(nrow(coo))
means <- numeric(nrow(coo))
lower_ci <- numeric(nrow(coo))
upper_ci <- numeric(nrow(coo))

# Also extract summary statistics for comparison
summary_modes <- numeric(nrow(coo))
summary_means <- numeric(nrow(coo))
summary_lower <- numeric(nrow(coo))
summary_upper <- numeric(nrow(coo))

for(i in 1:nrow(coo)){
  # Get the marginal for this observation using the estimation index
  if(i <= length(index_est)) {
    est_idx <- index_est[i]
    marg <- res$marginals.fitted.values[[est_idx]]
    
    # Also get summary statistics for this same index
    summary_modes[i] <- res$summary.fitted.values[est_idx, "mode"]
    summary_means[i] <- res$summary.fitted.values[est_idx, "mean"]
    summary_lower[i] <- res$summary.fitted.values[est_idx, "0.025quant"]
    summary_upper[i] <- res$summary.fitted.values[est_idx, "0.975quant"]
    
    if(!is.null(marg) && nrow(marg) > 0) {
      # Use INLA's built-in functions for proper calculations
      modes[i] <- inla.mmarginal(marg)  # Mode
      means[i] <- inla.emarginal(function(x) x, marg)  # Mean
      lower_ci[i] <- inla.qmarginal(0.025, marg)  # 2.5th percentile
      upper_ci[i] <- inla.qmarginal(0.975, marg)  # 97.5th percentile
    } else {
      modes[i] <- NA
      means[i] <- NA
      lower_ci[i] <- NA
      upper_ci[i] <- NA
    }
  } else {
    # Index out of bounds
    modes[i] <- NA
    means[i] <- NA
    lower_ci[i] <- NA
    upper_ci[i] <- NA
    summary_modes[i] <- NA
    summary_means[i] <- NA
    summary_lower[i] <- NA
    summary_upper[i] <- NA
  }
}

# Compare marginal vs summary calculations
cat("Calculation comparison (first 3 valid observations):\n")
valid_for_comparison <- which(!is.na(modes) & !is.na(summary_modes))[1:3]
for(i in valid_for_comparison) {
  cat("Obs", i, ":\n")
  cat("  Marginal: mode =", round(modes[i], 4), ", CI = [", round(lower_ci[i], 4), ",", round(upper_ci[i], 4), "]\n")
  cat("  Summary:  mode =", round(summary_modes[i], 4), ", CI = [", round(summary_lower[i], 4), ",", round(summary_upper[i], 4), "]\n")
}



# Debug: Check the credible interval calculations
cat("\n=== DEBUGGING CREDIBLE INTERVALS ===\n")
cat("Number of observations:", nrow(coo), "\n")
cat("Valid modes:", sum(!is.na(modes)), "\n")
cat("Valid lower_ci:", sum(!is.na(lower_ci)), "\n")
cat("Valid upper_ci:", sum(!is.na(upper_ci)), "\n")
cat("CI ranges - Min:", round(min(upper_ci - lower_ci, na.rm=TRUE), 4), 
    "Max:", round(max(upper_ci - lower_ci, na.rm=TRUE), 4), 
    "Mean:", round(mean(upper_ci - lower_ci, na.rm=TRUE), 4), "\n")

# Check a few examples
if(sum(!is.na(modes)) > 0) {
  valid_indices <- which(!is.na(modes))[1:min(5, sum(!is.na(modes)))]
  cat("First few examples:\n")
  for(idx in valid_indices) {
    cat("Obs", idx, ": prevalence =", round(prevalence[idx], 4), 
        ", mode =", round(modes[idx], 4), 
        ", CI = [", round(lower_ci[idx], 4), ",", round(upper_ci[idx], 4), "]\n")
  }
}

# COMPREHENSIVE DATA VALIDATION AND FILTERING
cat("\n=== DATA VALIDATION ===\n")
cat("Original data points:", length(prevalence), "\n")
cat("Sample size distribution:\n")
print(table(cut(sample_size, breaks=c(0,5,10,20,50,100,Inf), labels=c("1-5","6-10","11-20","21-50","51-100",">100"))))

# First filter: basic validity
valid_basic <- !is.na(prevalence) & prevalence >= 0 & !is.na(modes) & !is.na(sample_size) & sample_size > 0
cat("Valid basic observations:", sum(valid_basic), "\n")

# Extract basic valid data
prev_basic <- prevalence[valid_basic]
ss_basic <- sample_size[valid_basic]
modes_basic <- modes[valid_basic]
means_basic <- means[valid_basic]
lower_ci_basic <- lower_ci[valid_basic]
upper_ci_basic <- upper_ci[valid_basic]

# Check for data extraction issues
cat("Prevalence range: [", round(min(prev_basic, na.rm=T), 4), ",", round(max(prev_basic, na.rm=T), 4), "]\n")
cat("Modes range: [", round(min(modes_basic, na.rm=T), 4), ",", round(max(modes_basic, na.rm=T), 4), "]\n")
cat("CI widths - Min:", round(min(upper_ci_basic - lower_ci_basic, na.rm=T), 4), 
    "Max:", round(max(upper_ci_basic - lower_ci_basic, na.rm=T), 4), "\n")

# Check for unrealistic values
unrealistic_prev <- prev_basic < 0 | prev_basic > 1
unrealistic_modes <- modes_basic < 0 | modes_basic > 1
cat("Unrealistic prevalences:", sum(unrealistic_prev, na.rm=T), "\n")
cat("Unrealistic modes:", sum(unrealistic_modes, na.rm=T), "\n")

# COMPREHENSIVE SAMPLE SIZE FILTERING - Test many thresholds including high ones
sample_thresholds <- c(1, 5, 10, 20, 30, 50, 100, 200, 500, 1000)
coverage_results <- data.frame(
  threshold = sample_thresholds,
  n_obs = NA,
  coverage_marginal = NA,
  coverage_summary = NA,
  mean_ci_width = NA
)

cat("\n=== COVERAGE BY SAMPLE SIZE THRESHOLD ===\n")
cat("Threshold\tN_obs\tCov_Marg\tCov_Summ\tMean_CI_Width\n")

for(j in 1:length(sample_thresholds)) {
  min_ss <- sample_thresholds[j]
  ss_filter <- ss_basic >= min_ss
  n_filtered <- sum(ss_filter)
  
  if(n_filtered >= 5) {  # Need at least 5 observations
    # Extract filtered data
    test_prev <- prev_basic[ss_filter]
    test_ss <- ss_basic[ss_filter]
    test_modes <- modes_basic[ss_filter]
    test_lower <- lower_ci_basic[ss_filter]
    test_upper <- upper_ci_basic[ss_filter]
    test_summary_lower <- summary_lower[valid_basic][ss_filter]
    test_summary_upper <- summary_upper[valid_basic][ss_filter]
    
    # Calculate coverage using marginals
    test_within_marg <- (test_prev >= test_lower) & (test_prev <= test_upper)
    test_coverage_marg <- mean(test_within_marg, na.rm = TRUE)
    
    # Calculate coverage using INLA summary
    test_within_summ <- (test_prev >= test_summary_lower) & (test_prev <= test_summary_upper)
    test_coverage_summ <- mean(test_within_summ, na.rm = TRUE)
    
    # Mean CI width
    mean_width <- mean(test_upper - test_lower, na.rm = TRUE)
    
    # Store results
    coverage_results$n_obs[j] <- n_filtered
    coverage_results$coverage_marginal[j] <- test_coverage_marg
    coverage_results$coverage_summary[j] <- test_coverage_summ
    coverage_results$mean_ci_width[j] <- mean_width
    
    cat(sprintf("%d\t\t%d\t\t%.1f%%\t\t%.1f%%\t\t%.4f\n", 
                min_ss, n_filtered, 
                test_coverage_marg*100, test_coverage_summ*100, 
                mean_width))
  } else {
    cat(sprintf("%d\t\t%d\t\tN/A\t\tN/A\t\tN/A\n", min_ss, n_filtered))
  }
}

# Find the best threshold
valid_results <- coverage_results[!is.na(coverage_results$coverage_marginal), ]
if(nrow(valid_results) > 0) {
  # Find threshold that gives coverage closest to 95% with reasonable sample size
  target_coverage <- 0.95
  score <- abs(valid_results$coverage_marginal - target_coverage) + 
           0.1 * abs(valid_results$coverage_summary - target_coverage) +
           0.01 * pmax(0, 50 - valid_results$n_obs)  # Penalty for very small samples
  
  best_idx <- which.min(score)
  best_threshold <- valid_results$threshold[best_idx]
  cat("\nBest sample size threshold:", best_threshold, 
      "(Coverage:", round(valid_results$coverage_marginal[best_idx]*100, 1), "%, N =", 
      valid_results$n_obs[best_idx], ")\n")
} else {
  best_threshold <- 10  # Default fallback
  cat("\nUsing default threshold:", best_threshold, "\n")
}

# OPTION 1: Use optimal threshold (may reduce data significantly)
# ss_threshold <- best_threshold
# OPTION 2: Use minimal filtering to preserve more data points
ss_threshold <- 1  # Keep all valid data points for now

valid_obs_final <- valid_basic & ss_basic >= ss_threshold
cat("\nUsing sample size threshold >=", ss_threshold, "(to preserve data points)\n")
cat("Final valid observations:", sum(valid_obs_final), "\n")

# Apply the final filtering consistently to ALL data
prevalence_valid <- prev_basic[ss_basic >= ss_threshold]
sample_size_valid <- ss_basic[ss_basic >= ss_threshold]
modes_valid <- modes_basic[ss_basic >= ss_threshold]
means_valid <- means_basic[ss_basic >= ss_threshold]
lower_ci_valid <- lower_ci_basic[ss_basic >= ss_threshold]
upper_ci_valid <- upper_ci_basic[ss_basic >= ss_threshold]

cat("Final data verification:\n")
cat("  Prevalence length:", length(prevalence_valid), "\n")
cat("  Sample size length:", length(sample_size_valid), "\n") 
cat("  Modes length:", length(modes_valid), "\n")
cat("  All should match! (", sum(valid_obs_final), ")\n")



# Create sample size categories BEFORE updating the main variables
sample_size_factored <- as.numeric(sample_size_valid)
sample_size_factored[sample_size_valid <= 10] <- 1
sample_size_factored[sample_size_valid > 10 & sample_size_valid <= 20] <- 2
sample_size_factored[sample_size_valid > 20 & sample_size_valid <= 30] <- 3
sample_size_factored[sample_size_valid > 30 & sample_size_valid <= 50] <- 4
sample_size_factored[sample_size_valid > 50 & sample_size_valid <= 100] <- 5
sample_size_factored[sample_size_valid > 100 & sample_size_valid <= 200] <- 6
sample_size_factored[sample_size_valid > 200 & sample_size_valid <= 500] <- 7
sample_size_factored[sample_size_valid > 500 & sample_size_valid <= 1000] <- 8
sample_size_factored[sample_size_valid > 1000] <- 9

# Use the consistently filtered data throughout the analysis
prevalence <- prevalence_valid
sample_size <- sample_size_valid
modes <- modes_valid
means <- means_valid
lower_ci <- lower_ci_valid
upper_ci <- upper_ci_valid

cat("Updated main variables - all have length:", length(prevalence), "\n")
cat("Data preservation summary:\n")
cat("  Original observations:", nrow(coo), "\n")
cat("  After basic filtering:", sum(valid_basic), "\n") 
cat("  Final observations (ss >=", ss_threshold, "):", length(prevalence), "\n")
cat("  Data retention rate:", round(length(prevalence) / nrow(coo) * 100, 1), "%\n")

# Create the main data frame for plotting
dpm <- data.frame(
  prevalence = prevalence, 
  modes = modes, 
  means = means, 
  sample_size = sample_size_factored,
  sample_size_orig = sample_size
)

# Enhanced scatterplot with better aesthetics
p_scatter <- ggplot(dpm, aes(x = prevalence, y = modes)) +
  geom_point(aes(color = factor(sample_size), size = factor(sample_size)), 
             alpha = 0.7, stroke = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", size = 1.2) +
  scale_color_viridis_d(name = "Sample size", option = "plasma",
                        labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  scale_size_manual(values = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                    name = "Sample size",
                    labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Model validation: observed vs predicted prevalence -", countries_to_infer),
    x = "Observed prevalence",
    y = "Predicted prevalence (posterior mode)",
    caption = paste("Correlation:", round(cor(prevalence, modes, use="complete.obs"), 3), 
                   "| Red line: perfect agreement | Blue line: fitted trend")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_scatter)

#save
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"Prevalence_vs_Mode.png"),
       plot = p_scatter, width = 12, height = 8, dpi = 300)

print(cor(prevalence,modes))

#===============================================================================
# ADDITIONAL MODEL VALIDATION PLOTS
#===============================================================================

# 1. RESIDUALS PLOT
residuals <- prevalence - modes
dpm$residuals <- residuals
dpm$fitted <- modes

p_residuals <- ggplot(dpm, aes(x = fitted, y = residuals)) +
  geom_point(aes(color = factor(sample_size), size = factor(sample_size)), 
             alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#E74C3C", size = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "#3498DB", fill = "#3498DB", alpha = 0.2) +
  scale_color_viridis_d(name = "Sample size", option = "plasma",
                        labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  scale_size_manual(values = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                    name = "Sample size",
                    labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Residuals plot -", countries_to_infer),
    x = "Fitted values (predicted prevalence)",
    y = "Residuals (observed - predicted)",
    caption = "Points should be randomly scattered around zero line"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_residuals)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_Residuals_Plot.png"),
       plot = p_residuals, width = 12, height = 8, dpi = 300)

# 2. CREDIBLE INTERVALS COVERAGE
# Use the correctly calculated 95% credible intervals from the posterior marginals
dpm$pred_lower <- lower_ci
dpm$pred_upper <- upper_ci
dpm$within_ci <- (prevalence >= dpm$pred_lower) & (prevalence <= dpm$pred_upper)

coverage_rate <- mean(dpm$within_ci, na.rm = TRUE)

# Debug coverage calculation
cat("\n=== DEBUGGING COVERAGE ===\n")
cat("Observations within CI:", sum(dpm$within_ci, na.rm = TRUE), "\n")
cat("Total valid observations:", sum(!is.na(dpm$within_ci)), "\n")
cat("Coverage rate:", round(coverage_rate * 100, 2), "%\n")

# Check some examples of coverage failures
outside_ci <- which(!dpm$within_ci & !is.na(dpm$within_ci))
if(length(outside_ci) > 0) {
  cat("Examples of observations outside CI:\n")
  n_examples <- min(5, length(outside_ci))
  for(i in 1:n_examples) {
    idx <- outside_ci[i]
    cat("Obs", idx, ": prevalence =", round(prevalence[idx], 4),
        ", CI = [", round(dpm$pred_lower[idx], 4), ",", round(dpm$pred_upper[idx], 4), "], mode =", round(modes[idx], 4), "\n")
  }
}

# ALTERNATIVE ANALYSIS: Use INLA summary statistics directly
cat("\n=== ALTERNATIVE COVERAGE USING INLA SUMMARY ===\n")
alt_modes <- summary_modes[valid_basic][ss_basic >= ss_threshold]
alt_lower <- summary_lower[valid_basic][ss_basic >= ss_threshold]  
alt_upper <- summary_upper[valid_basic][ss_basic >= ss_threshold]

alt_within_ci <- (prevalence >= alt_lower) & (prevalence <= alt_upper)
alt_coverage_rate <- mean(alt_within_ci, na.rm = TRUE)

cat("Alternative coverage rate (INLA summary):", round(alt_coverage_rate * 100, 2), "%\n")

# Compare the two approaches
cat("Coverage comparison:\n")
cat("  Marginal-based:", round(coverage_rate * 100, 2), "%\n")
cat("  Summary-based:", round(alt_coverage_rate * 100, 2), "%\n")

# If summary-based is much better, use that instead
if(!is.na(alt_coverage_rate) && alt_coverage_rate > coverage_rate + 0.1) {
  cat("Using summary-based approach as it shows better coverage.\n")
  dpm$pred_lower <- alt_lower
  dpm$pred_upper <- alt_upper
  dpm$within_ci <- alt_within_ci
  coverage_rate <- alt_coverage_rate
  modes <- alt_modes
  dpm$modes <- alt_modes  # Update the data frame as well
}

# COVERAGE BY SAMPLE SIZE PLOTS
cat("\n=== CREATING COVERAGE PLOTS ===\n")

# Plot coverage vs sample size threshold
valid_coverage_data <- coverage_results[!is.na(coverage_results$coverage_marginal), ]
if(nrow(valid_coverage_data) > 1) {
  
  p_coverage_threshold <- ggplot(valid_coverage_data, aes(x = threshold)) +
    geom_line(aes(y = coverage_marginal * 100, color = "Marginal-based"), size = 1.2) +
    geom_line(aes(y = coverage_summary * 100, color = "Summary-based"), size = 1.2) +
    geom_point(aes(y = coverage_marginal * 100, color = "Marginal-based"), size = 3) +
    geom_point(aes(y = coverage_summary * 100, color = "Summary-based"), size = 3) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "#E74C3C", size = 1) +
    scale_color_viridis_d(name = "Method", option = "plasma", begin = 0.2, end = 0.8) +
    scale_x_log10(breaks = c(1, 5, 10, 20, 50, 100, 200, 500, 1000)) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Coverage Rate by Sample Size Threshold -", countries_to_infer),
      x = "Minimum Sample Size (log scale)",
      y = "Coverage Rate (%)",
      caption = "Target: 95% coverage (red line)"
    ) +
    ylim(0, 100)
  
  print(p_coverage_threshold)
  ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_Coverage_by_SampleSize.png"),
         plot = p_coverage_threshold, width = 12, height = 8, dpi = 300)
  
  # Plot sample sizes vs coverage
  p_samplesize_scatter <- ggplot(valid_coverage_data, aes(x = n_obs, y = coverage_marginal * 100)) +
    geom_point(size = 4, alpha = 0.7, color = "#3498DB") +
    geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", fill = "#E74C3C", alpha = 0.2) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "#E74C3C", size = 1) +
    scale_x_log10() +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11)
    ) +
    labs(
      title = paste("Coverage Rate vs Number of Observations -", countries_to_infer),
      x = "Number of Observations (log scale)",
      y = "Coverage Rate (%)",
      caption = "Target: 95% coverage (red line)"
    ) +
    ylim(0, 100)
  
  print(p_samplesize_scatter)
  ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_Coverage_vs_N_Obs.png"),
         plot = p_samplesize_scatter, width = 12, height = 8, dpi = 300)
}

p_coverage <- ggplot(dpm, aes(x = modes, y = prevalence)) +
  geom_ribbon(aes(ymin = pred_lower, ymax = pred_upper), alpha = 0.3, fill = "#3498DB") +
  geom_point(aes(color = within_ci, size = factor(sample_size)), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", size = 1) +
  scale_color_manual(values = c("FALSE" = "#E74C3C", "TRUE" = "#27AE60"), 
                     name = "Within 95% CI", labels = c("No", "Yes")) +
  scale_size_manual(values = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                    name = "Sample size",
                    labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Credible interval coverage -", countries_to_infer),
    x = "Predicted prevalence",
    y = "Observed prevalence", 
    caption = paste("Coverage rate:", round(coverage_rate * 100, 1), "% (should be ~95%)")
  )

print(p_coverage)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_PI_Coverage.png"),
       plot = p_coverage, width = 12, height = 8, dpi = 300)

# 3. BLAND-ALTMAN PLOT
mean_values <- (prevalence + modes) / 2
diff_values <- prevalence - modes
mean_diff <- mean(diff_values, na.rm = TRUE)
sd_diff <- sd(diff_values, na.rm = TRUE)

dpm$mean_values <- mean_values
dpm$diff_values <- diff_values

p_bland_altman <- ggplot(dpm, aes(x = mean_values, y = diff_values)) +
  geom_point(aes(color = factor(sample_size), size = factor(sample_size)), alpha = 0.7) +
  geom_hline(yintercept = mean_diff, color = "#3498DB", size = 1) +
  geom_hline(yintercept = mean_diff + 1.96 * sd_diff, color = "#E74C3C", linetype = "dashed") +
  geom_hline(yintercept = mean_diff - 1.96 * sd_diff, color = "#E74C3C", linetype = "dashed") +
  scale_color_viridis_d(name = "Sample size", option = "plasma",
                        labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  scale_size_manual(values = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                    name = "Sample size",
                    labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Bland-Altman plot -", countries_to_infer),
    x = "Average of observed and predicted prevalence",
    y = "Difference (observed - predicted)",
    caption = paste("Mean difference:", round(mean_diff, 4), "| 95% limits:", 
                   round(mean_diff - 1.96*sd_diff, 4), "to", round(mean_diff + 1.96*sd_diff, 4))
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_bland_altman)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_Bland_Altman.png"),
       plot = p_bland_altman, width = 12, height = 8, dpi = 300)

# 4. QUANTILE-QUANTILE PLOT
# Create theoretical quantiles vs observed quantiles for model validation
theoretical_quantiles <- qnorm(ppoints(length(prevalence)))
observed_quantiles <- qnorm(rank(prevalence) / (length(prevalence) + 1))
predicted_quantiles <- qnorm(rank(modes) / (length(modes) + 1))

qq_data <- data.frame(
  theoretical = theoretical_quantiles,
  observed = observed_quantiles,
  predicted = predicted_quantiles,
  sample_size = factor(sample_size_factored)
)

p_qq <- ggplot(qq_data, aes(x = observed, y = predicted)) +
  geom_point(aes(color = sample_size, size = sample_size), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", size = 1.2) +
  scale_color_viridis_d(name = "Sample size", option = "plasma",
                        labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  scale_size_manual(values = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                    name = "Sample size",
                    labels = c("1-10", "11-20", "21-30", "31-50", "51-100", "101-200", "201-500", "501-1000", ">1000")) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Q-Q plot: observed vs predicted quantiles -", countries_to_infer),
    x = "Observed quantiles (normal scores)",
    y = "Predicted quantiles (normal scores)", 
    caption = "Points should follow the diagonal line if distributions match"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_qq)
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"_QQ_Plot.png"),
       plot = p_qq, width = 12, height = 8, dpi = 300)

# 5. MODEL PERFORMANCE METRICS
cat("\n=== MODEL VALIDATION METRICS ===\n")
cat("Sample size:", length(prevalence), "observations\n")
cat("Correlation (observed vs predicted):", round(cor(prevalence, modes, use="complete.obs"), 4), "\n")
cat("Mean Absolute Error (MAE):", round(mean(abs(prevalence - modes), na.rm = TRUE), 4), "\n")
cat("Root Mean Square Error (RMSE):", round(sqrt(mean((prevalence - modes)^2, na.rm = TRUE)), 4), "\n")
cat("Mean Bias (observed - predicted):", round(mean(prevalence - modes, na.rm = TRUE), 4), "\n")
cat("95% Credible Interval Coverage:", round(coverage_rate * 100, 1), "% (should be ~95%)\n")

# Calculate R-squared
ss_res <- sum((prevalence - modes)^2, na.rm = TRUE)
ss_tot <- sum((prevalence - mean(prevalence, na.rm = TRUE))^2, na.rm = TRUE)
r_squared <- 1 - (ss_res / ss_tot)
cat("R-squared:", round(r_squared, 4), "\n")

# Concordance correlation coefficient (better for agreement)
mean_obs <- mean(prevalence, na.rm = TRUE)
mean_pred <- mean(modes, na.rm = TRUE)
var_obs <- var(prevalence, na.rm = TRUE)
var_pred <- var(modes, na.rm = TRUE)
cov_obs_pred <- cov(prevalence, modes, use = "complete.obs")
ccc <- (2 * cov_obs_pred) / (var_obs + var_pred + (mean_obs - mean_pred)^2)
cat("Concordance Correlation Coefficient:", round(ccc, 4), "\n")

# 6. SAMPLE SIZE STRATIFIED PLOTS (IMPROVED)
for (s in c(10, 20, 30, 50, 100, 200, 500)){
  modes_s <- modes[sample_size >= s]
  prevalence_s <- prevalence[sample_size >= s]
  n_points <- length(modes_s)
  correlation_s <- if(n_points > 2) round(cor(prevalence_s, modes_s, use="complete.obs"), 3) else NA
  
  dpm_s <- data.frame(prevalence=prevalence_s, modes=modes_s)
  
  p_subset <- ggplot(dpm_s, aes(x = prevalence, y = modes)) +
    geom_point(color = "#3498DB", size = 2.5, alpha = 0.7, stroke = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C", size = 1.2) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste("Model validation: sample size ≥", s, "-", countries_to_infer),
      x = "Observed prevalence", 
      y = "Predicted prevalence (posterior mode)",
      caption = paste("N =", n_points, "| Correlation:", 
                     if(is.na(correlation_s)) "N/A" else correlation_s,
                     "| Red line: perfect agreement")
    ) +
    xlim(0, max(c(prevalence_s, modes_s), na.rm = TRUE) * 1.05) +
    ylim(0, max(c(prevalence_s, modes_s), na.rm = TRUE) * 1.05)
  
  print(p_subset)
  ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"Prevalence_vs_Mode_ss",s,".png"),
         plot = p_subset, width = 10, height = 8, dpi = 300)
}


######### SHADING OUT HIGH UNCERTAINTY #########

#calculate the uncertainty as the difference between the 97.5th and 2.5th percentiles
uncertainty <- pred_ul - pred_ll

#remove the uncertainty values that are greater than 0.9
pred_ul[uncertainty > 0.9] <- NA
pred_ll[uncertainty > 0.9] <- NA
pred_mean[uncertainty > 0.9] <- NA

# plot
dpm <- rbind(
  data.frame(
    Latitude = coop[, 2], Longitude = coop[, 1],
    value = pred_mean, variable = "Mean"
  ),
  data.frame(
    Latitude = coop[, 2], Longitude = coop[, 1],
    value = pred_ll, variable = "2.5th percentile"
  ),
  data.frame(
    Latitude = coop[, 2], Longitude = coop[, 1],
    value = pred_ul, variable = "97.5th percentile"
  )
)

dpm$variable <- as.factor(dpm$variable)

#plot, with hashes where dp$water > 0.5

#water cover is True if dp$water > 0.5
water_cover <- as.numeric(dp$water > 0.5)

#make a dataframe with the water cover and Latitude and Longitude
water_cover <- data.frame(Longitude = coop[, 1], Latitude = coop[, 2], water_cover = water_cover)

#write.csv(water_cover, "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/NGA_water_cover.csv", row.names = FALSE)


water_cover <- do.call(rbind, lapply(c("2.5th percentile", "Mean", "97.5th percentile"), function(v) {
  tmp <- water_cover
  tmp$variable <- factor(v, levels = c("2.5th percentile", "Mean", "97.5th percentile"))
  tmp
}))

# Convert the 'value' column to a factor
water_cover$water_cover <- factor(water_cover$water_cover, levels = c(0, 1), labels = c("High uncertainty", "Water"))


#plot the watercover of Uganda
ggplot(dpm) + geom_tile(aes(Longitude, Latitude, fill = value)) +
  facet_grid(~factor(variable,levels=c("2.5th percentile","Mean","97.5th percentile"))) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Prevalence",
    low = "blue", high = "orange",
    limits = c(0, 1),
    guide = guide_colorbar(order = 2)
  ) +
  new_scale_fill() +  # Add a new fill scale
  geom_tile(data=water_cover, aes(Longitude, Latitude,fill=water_cover,alpha=water_cover),inherit.aes = FALSE)+
  scale_fill_manual(
  name = "",  # Colorbar title
  values = c("High uncertainty" = "grey50", "Water" = "steelblue1"),  # Specific colors for categories
  limits = c("High uncertainty", "Water"),
  guide = guide_legend(order = 1)
)+
guides(alpha = "none")+
  geom_polygon(data =  border, aes(x = X, y = Y), color = "black",fill=NA)+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.box = "vertical")
ggsave(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"Predictions_Shaded.png"),width = 12, height = 8)

#save dpm
write.csv(dpm, paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"Predictions.csv"), row.names = FALSE)
#keep dpm values correspoding to the mean
dpm_mean <- dpm[dpm$variable == "Mean",]
#keep only the values where the mean is not NA
dpm_mean <- dpm_mean[!is.na(dpm_mean$value),]
