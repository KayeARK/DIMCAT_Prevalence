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
library(FNN)

n_iterations <- 20

diagnostics_initialised <- FALSE

countries_to_infer=c("Nigeria")

#add string to the name of the file

data_path <- paste0("Code/TestSensSpec/adjusted_cases_all_sites.csv")
data_original <- read.csv(data_path)

# Ensure coordinates are numeric and clean from the start
data_original$longitude <- as.numeric(data_original$longitude)
data_original$latitude <- as.numeric(data_original$latitude)

# Remove any rows with NA coordinates before aggregation
na_coords <- is.na(data_original$longitude) | is.na(data_original$latitude)
if(sum(na_coords) > 0) {
  cat("Removing", sum(na_coords), "rows with NA coordinates from original data\n")
  data_original <- data_original[!na_coords, ]
}

# Create coordinate mapping for aggregation (done once)
cat("Creating coordinate aggregation mapping...\n")
cat("Original data has", nrow(data_original), "records\n")

# Use base R aggregate to avoid type conversion issues
coord_sample_mapping <- aggregate(
  sample_size ~ longitude + latitude, 
  data = data_original, 
  FUN = sum, 
  na.rm = TRUE
)

cat("After coordinate aggregation:", nrow(coord_sample_mapping), "unique locations\n")
cat("Reduced from", nrow(data_original), "to", nrow(coord_sample_mapping), "records\n")

#i should go from 8 to 1007

for (iteration in (1+7):(1+7+n_iterations-1)){
print(iteration-7)

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
#add elevation data


#ELEVATION
r_elv <- elevation_30s(country = "NGA", path = "Data/Covariates")
data$elevation<-terra::extract(r_elv, data[, c("longitude", "latitude")])$NGA_elv_msk

#PRECIPITATION
r_prec<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "prec")
data$precipitation<-(terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_1+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_2+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_3+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_4+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_5+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_6+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_7+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_8+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_9+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_10+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_11+
terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_12)/12


#AVERAGE TEMPERATURE
r_tavg<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tavg")
data$tavg<-(terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_1+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_2+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_3+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_4+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_5+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_6+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_7+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_8+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_9+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_10+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_11+
terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_12)/12

#MINIMUM TEMPERATURE
r_tmin<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tmin")
data$tmin<-(terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_1+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_2+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_3+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_4+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_5+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_6+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_7+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_8+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_9+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_10+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_11+
terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_12)/12

#MAXIMUM TEMPERATURE
r_tmax<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tmax")
data$tmax<-(terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_1+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_2+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_3+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_4+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_5+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_6+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_7+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_8+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_9+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_10+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_11+
terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_12)/12

#POPULATION DENSITY
r_pop_den <- population(year=2010, res=0.5, path="Data/Covariates")
data$pop_den <-terra::extract(r_pop_den, data[, c("longitude", "latitude")])$population_density

#LANDCOVER
r_tree <- landcover("trees",path="Data/Covariates")
data$tree <- terra::extract(r_tree, data[, c("longitude", "latitude")])$trees

r_grassland <- landcover("grassland",path="Data/Covariates")
data$grassland <- terra::extract(r_grassland, data[, c("longitude", "latitude")])$grassland

r_shrub <- landcover("shrubs",path="Data/Covariates")
data$shrub <- terra::extract(r_shrub, data[, c("longitude", "latitude")])$shrub

r_cropland <- landcover("cropland",path="Data/Covariates")
data$cropland <- terra::extract(r_cropland, data[, c("longitude", "latitude")])$cropland

r_built <- landcover("built",path="Data/Covariates")
data$built <- terra::extract(r_built, data[, c("longitude", "latitude")])$built

r_bare <- landcover("bare",path="Data/Covariates")
data$bare <- terra::extract(r_bare, data[, c("longitude", "latitude")])$bare

r_water <- landcover("water",path="Data/Covariates")
data$water <- terra::extract(r_water, data[, c("longitude", "latitude")])$water

r_wetland <- landcover("wetland",path="Data/Covariates")
data$wetland <- terra::extract(r_wetland, data[, c("longitude", "latitude")])$wetland

r_mangrove <- landcover("mangroves",path="Data/Covariates")
data$mangrove <- terra::extract(r_mangrove, data[, c("longitude", "latitude")])$mangroves

r_moss <- landcover("moss",path="Data/Covariates")
data$moss <- terra::extract(r_moss, data[, c("longitude", "latitude")])$moss

data$tsetse_habitat <- data$tree + data$wetland

#HUMAN FOOTPRINT
r_human_fp<-footprint(year=2009, res=30, path="Data/Covariates")
data$human_fp<-terra::extract(r_human_fp, data[, c("longitude", "latitude")])[,2]


#CATTLE DATA
r_cattle<-raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
data$cattle<-terra::extract(r_cattle, data[, c("longitude", "latitude")])
data$cattle[is.na(data$cattle)] <- 0

r_buffalo<-raster("Data/Covariates/livestock/buffalo_2015/5_Bf_2015_Da.tif")
data$buffalo<-terra::extract(r_buffalo, data[, c("longitude", "latitude")])
data$buffalo[is.na(data$buffalo)] <- 0

r_goat<-raster("Data/Covariates/livestock/goats_2015/5_Gt_2015_Da.tif")
data$goat<-terra::extract(r_goat, data[, c("longitude", "latitude")])
data$goat[is.na(data$goat)] <- 0

r_horse<-raster("Data/Covariates/livestock/horses_2015/5_Ho_2015_Da.tif")
data$horse<-terra::extract(r_horse, data[, c("longitude", "latitude")])
data$horse[is.na(data$horse)] <- 0

r_pig<-raster("Data/Covariates/livestock/pigs_2015/5_Pg_2015_Da.tif")
data$pig<-terra::extract(r_pig, data[, c("longitude", "latitude")])
data$pig[is.na(data$pig)] <- 0

r_sheep<-raster("Data/Covariates/livestock/sheep_2015/5_Sh_2015_Da.tif")
data$sheep<-terra::extract(r_sheep, data[, c("longitude", "latitude")])
data$sheep[is.na(data$sheep)] <- 0

tsetse<-raster("Data/Covariates/tsenumbspec")
tsetse[tsetse > 1] <- 1
data$tsetse <-terra::extract(tsetse, data[, c("longitude", "latitude")])

data$total_animal<-data$cattle+data$buffalo+data$goat+data$horse+data$pig+data$sheep


latitudes<-data$latitude
longitudes<-data$longitude
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

# pdf("Code/AbsencePresence/Bovine BCT/Results/DataPlot_AT_temp.pdf")
# #plot africountries from afrilearndata
# plot(africountries$geometry, border="black", col="white")
# points(longitudes, latitudes, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
# #increase the y limit
# legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
# dev.off()

# put latitude and longitude next to each other
coo <- cbind(longitudes, latitudes)

# Reasonable mesh size - reduces vertices from ~8000 to ~800-1500
mesh <- inla.mesh.2d(
  loc = coo, 
  offset = c(30, 60),        # Slightly reduced offset
  cutoff = 3,                # Increased from 1 (reasonable middle ground)
  max.edge = c(6, 15)        # Increased from c(3, 10) 
)

# Informative priors to prevent spatial overfitting
spde <- inla.spde2.pcmatern(
  mesh = mesh, 
  alpha = 2, 
  prior.range = c(50, 0.1),    # 50km range, 10% prob range < 50km
  prior.sigma = c(1, 0.1),     # Moderate spatial variance
  constr = TRUE
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

cov_file <- paste0("Code/Prevalence/Bovine BCT and PCR/Covariates_NGA/","Covariates_model_", iteration-7, ".csv")
covariates <- read.csv(cov_file)$covariates

formula <- as.formula(paste0("y ~ 0 + b0 +",paste(covariates, collapse = " + "),"+ f(s, model = spde)"))

res <- inla(formula,
  data = inla.stack.data(stk.full),
  family = "binomial", Ntrials = numtrials,
  control.compute=list(return.marginals.predictor=TRUE, dic=TRUE, waic=TRUE),
  control.predictor = list(link=1,compute = TRUE,
  A = inla.stack.A(stk.full)),)


  #### DIAGNOSTICS ####

  ############################################################
################ INITIALISE STORAGE ########################
############################################################

if(!diagnostics_initialised){

  n_pred <- nrow(coop)

  sum_prev_mean <- rep(0, n_pred)
  sum_prev_var  <- rep(0, n_pred)

  sum_eta_mean <- rep(0, n_pred)
  sum_eta_var  <- rep(0, n_pred)

  sum_spatial_mean <- rep(0, n_pred)
  sum_spatial_var  <- rep(0, n_pred)

  sum_nogp <- rep(0, n_pred)

  ##########################################################
  ################ DISTANCE TO DATA ########################
  ##########################################################

  dist_to_data <- get.knnx(
    coo,
    coop,
    k = 1
  )$nn.dist[,1]

  diagnostics_initialised <- TRUE
}

############################################################
################ PREDICTION INDEX ##########################
############################################################

index <- inla.stack.index(
  stk.full,
  tag = "pred"
)$data

############################################################
################ PREVALENCE DIAGNOSTICS ####################
############################################################

pred_mean <- res$summary.fitted.values[
  index,
  "mean"
]

pred_sd <- res$summary.fitted.values[
  index,
  "sd"
]

sum_prev_mean <- sum_prev_mean + pred_mean

sum_prev_var <- sum_prev_var + pred_sd^2

############################################################
################ LATENT LINEAR PREDICTOR ###################
############################################################

eta_mean <- res$summary.linear.predictor[
  index,
  "mean"
]

eta_sd <- res$summary.linear.predictor[
  index,
  "sd"
]

sum_eta_mean <- sum_eta_mean + eta_mean

sum_eta_var <- sum_eta_var + eta_sd^2

############################################################
################ SPATIAL FIELD #############################
############################################################

spatial_mean <- as.vector(
  Ap %*%
  res$summary.random$s$mean
)

spatial_sd <- sqrt(
  as.vector(
    (Ap^2) %*%
    (res$summary.random$s$sd^2)
  )
)

sum_spatial_mean <-
  sum_spatial_mean + spatial_mean

sum_spatial_var <-
  sum_spatial_var + spatial_sd^2

############################################################
################ NO-SPATIAL MODEL ##########################
############################################################

formula_nogp <- as.formula(
  paste0(
    "y ~ 0 + b0 + ",
    paste(covariates, collapse = " + ")
  )
)

res_nogp <- inla(
  formula_nogp,

  data = inla.stack.data(stk.full),

  family = "binomial",

  Ntrials = numtrials,

  control.predictor = list(
    link = 1,
    compute = TRUE,
    A = inla.stack.A(stk.full)
  ),

  verbose = FALSE,

  inla.mode = "compact"
)

pred_nogp <- res_nogp$summary.fitted.values[
  index,
  "mean"
]

sum_nogp <- sum_nogp + pred_nogp

############################################################
################ CLEAN MEMORY ##############################
############################################################

rm(
  res_nogp,
  pred_mean,
  pred_sd,
  eta_mean,
  eta_sd,
  spatial_mean,
  spatial_sd,
  pred_nogp
)

gc()

}

############################################################
################ COMPUTE AVERAGES ##########################
############################################################

avg_prev_mean <- sum_prev_mean / n_iterations

avg_prev_sd <- sqrt(
  sum_prev_var / n_iterations
)

avg_eta_mean <- sum_eta_mean / n_iterations

avg_eta_sd <- sqrt(
  sum_eta_var / n_iterations
)

avg_spatial_mean <- sum_spatial_mean / n_iterations

avg_spatial_sd <- sqrt(
  sum_spatial_var / n_iterations
)

############################################################
################ FIXED EFFECT COMPONENT ####################
############################################################

fixed_effect_mean <- avg_eta_mean - avg_spatial_mean

avg_nogp <- sum_nogp / n_iterations

############################################################
################ SAVE DIAGNOSTICS ##########################
############################################################

diagnostics <- data.frame(

  lon = coop[,1],
  lat = coop[,2],

  ##########################################################
  ################ PREVALENCE ##############################
  ##########################################################

  prevalence_mean = avg_prev_mean,
  prevalence_sd   = avg_prev_sd,

  ##########################################################
  ################ LATENT SCALE ############################
  ##########################################################

  eta_mean = avg_eta_mean,
  eta_sd   = avg_eta_sd,

  ##########################################################
  ################ SPATIAL FIELD ###########################
  ##########################################################

  spatial_mean = avg_spatial_mean,
  spatial_sd   = avg_spatial_sd,

  ##########################################################
  ################ NO GP ###################################
  ##########################################################

  prevalence_no_gp = avg_nogp,

  ##########################################################
  ################ DATA COVERAGE ###########################
  ##########################################################

  dist_to_data = dist_to_data
)

# ############################################################
# ################ SANITY CHECK ##############################
# ############################################################

# cat("\n")
# cat("=====================================\n")
# cat("FIXED EFFECTS VS SPATIAL FIELD\n")
# cat("=====================================\n")

# cat("\nSDs:\n")

# cat(
#   "SD(Xbeta) =",
#   round(sd(fixed_effect_mean, na.rm = TRUE), 3),
#   "\n"
# )

# cat(
#   "SD(spatial field) =",
#   round(sd(avg_spatial_mean, na.rm = TRUE), 3),
#   "\n"
# )

# cat("\nRanges:\n")

# cat(
#   "Range(Xbeta) = [",
#   round(min(fixed_effect_mean, na.rm = TRUE), 3),
#   ", ",
#   round(max(fixed_effect_mean, na.rm = TRUE), 3),
#   "]\n",
#   sep = ""
# )

# cat(
#   "Range(spatial field) = [",
#   round(min(avg_spatial_mean, na.rm = TRUE), 3),
#   ", ",
#   round(max(avg_spatial_mean, na.rm = TRUE), 3),
#   "]\n",
#   sep = ""
# )

# cat("=====================================\n")

# cat(
#   "SD(spatial field) / SD(Xbeta) = ",
#   round(
#     sd(avg_spatial_mean, na.rm = TRUE) /
#     sd(fixed_effect_mean, na.rm = TRUE),
#     2
#   ),
#   "\n"
# )

write.csv(
  diagnostics,
  "Code/Prevalence/Bovine BCT and PCR/Diagnostics NGA/Diagnostics.csv",
  row.names = FALSE
)



