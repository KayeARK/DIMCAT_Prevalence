library(readxl)
library(ggplot2)
library(INLA)
library(sp)
library(sf)
library(afrilearndata)
library(concaveman)
library(geodata)
library(raster)
library(ggregplot)
library(see)

countries_to_infer=c("Ethiopia")
trypanosome<- c("Tc")

covariates <- read.csv(paste0("Code/AbsencePresence/Results/",trypanosome,"/",countries_to_infer,"Covariates.csv"))
#change to array
covariates <- covariates[,1]


#add string to the name of the file

data_path <- paste0("Data/ContAtlas_v2/250318_",trypanosome,"_Presence.xls")
data <- read_excel(data_path)

data <- data[!is.na(data$Trypanosome_detection) | !is.na(data$Number_of_infections),]

data$Trypanosome_detection <- as.factor(data$Trypanosome_detection)
data$Trypanosome_detection <- factor(data$Trypanosome_detection, levels = c("No", "Yes"), labels = c("No", "Yes"))
data$Trypanosome_detection <- as.numeric(data$Trypanosome_detection) - 1
presence <- data$Trypanosome_detection

#add elevation data


#ELEVATION
r_elv <- elevation_30s(country = "ETH", path = "Data/Covariates")
data$elevation<-terra::extract(r_elv, data[, c("Longitude", "Latitude")])$ETH_elv_msk

#PRECIPITATION
r_prec<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "prec")
data$precipitation<-(terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_1+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_2+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_3+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_4+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_5+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_6+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_7+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_8+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_9+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_10+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_11+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_prec_12)/12


#AVERAGE TEMPERATURE
r_tavg<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "tavg")
data$tavg<-(terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_1+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_2+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_3+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_4+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_5+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_6+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_7+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_8+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_9+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_10+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_11+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tavg_12)/12

#MINIMUM TEMPERATURE
r_tmin<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "tmin")
data$tmin<-(terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_1+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_2+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_3+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_4+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_5+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_6+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_7+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_8+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_9+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_10+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_11+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmin_12)/12

#MAXIMUM TEMPERATURE
r_tmax<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "tmax")
data$tmax<-(terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_1+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_2+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_3+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_4+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_5+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_6+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_7+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_8+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_9+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_10+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_11+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$ETH_wc2.1_30s_tmax_12)/12

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
presence <- presence[inside]
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

#remove NA values
ind <- !is.na(latitudes)
latitudes <- latitudes[ind]
longitudes <- longitudes[ind]
presence <- presence[ind]
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

pdf("Code/AbsencePresence/Results/Images/DataPlot_AT_temp.pdf")
#plot africountries from afrilearndata
plot(africountries$geometry, border="black", col="white")
points(longitudes, latitudes, pch=19, col=ifelse(presence == 1, "red", "black"), cex=0.5)
#increase the y limit
legend("bottomleft", legend=c("Present","Absent"), fill=rev(c("black","red")))
dev.off()

# put latitude and longitude next to each other
coo <- cbind(longitudes, latitudes)

mesh <- inla.mesh.2d(
  loc = coo, offset = c(50, 100),
  cutoff = 1, max.edge = c(30, 60)
)

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

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
dp$elevation <- extract(ra_elv, coop)[, 1]
dp$precipitation <- rowSums(extract(r_prec, coop)[, 1:12])/12
dp$tavg <- rowSums(extract(r_tavg, coop)[, 1:12])/12
dp$tmin <- rowSums(extract(r_tmin, coop)[, 1:12])/12
dp$tmax <- rowSums(extract(r_tmax, coop)[, 1:12])/12
dp$human_fp <- extract(r_human_fp, coop)[, 1]
dp$pop_den <- extract(r_pop_den, coop)[, 1]
dp$tree <- extract(r_tree, coop)[, 1]
dp$grassland <- extract(r_grassland, coop)[, 1]
dp$shrub <- extract(r_shrub, coop)[, 1]
dp$cropland <- extract(r_cropland, coop)[, 1]
dp$built <- extract(r_built, coop)[, 1]
dp$bare <- extract(r_bare, coop)[, 1]
dp$water <- extract(r_water, coop)[, 1]
dp$wetland <- extract(r_wetland, coop)[, 1]
dp$mangrove <- extract(r_mangrove, coop)[, 1]
dp$moss <- extract(r_moss, coop)[, 1]
dp$tsetse_habitat <- dp$tree + dp$wetland
dp$cattle <- extract(r_cattle, coop)
dp$cattle[is.na(dp$cattle)] <- 0
dp$buffalo <- extract(r_buffalo, coop)
dp$buffalo[is.na(dp$buffalo)] <- 0
dp$goat <- extract(r_goat, coop)
dp$goat[is.na(dp$goat)] <- 0
dp$horse <- extract(r_horse, coop)
dp$horse[is.na(dp$horse)] <- 0
dp$pig <- extract(r_pig, coop)
dp$pig[is.na(dp$pig)] <- 0
dp$sheep <- extract(r_sheep, coop)
dp$sheep[is.na(dp$sheep)] <- 0
dp$total_animal <- dp$cattle + dp$buffalo + dp$goat + dp$horse + dp$pig + dp$sheep

Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = presence),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo)), elevation=elevations,
  precipitation=precipitations, tavg=tavgs, tmin=tmins, tmax=tmaxs,
  human_fp=human_fps, pop_den=pop_dens, tree=trees, grassland=grasslands,
  shrub=shrubs, cropland=croplands, built=builts, bare=bares, water=waters, wetland=wetlands,
  mangrove=mangroves, moss=mosses, tsetse_habitat=tsetse_habitats, cattle=cattles, goat=goats,
  horse=horses, pig=pigs, buffalo=buffalos, sheep=sheeps, total_animal=total_animals), s = indexs)
)


# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop)),elevation=dp$elevation,
  precipitation=dp$precipitation, tavg=dp$tavg, tmin=dp$tmin, tmax=dp$tmax,
  human_fp=dp$human_fp, pop_den=dp$pop_den, tree=dp$tree, grassland=dp$grassland,
  shrub=dp$shrub, cropland=dp$cropland, built=dp$built, bare=dp$bare, water=dp$water, 
  wetland=dp$wetland, mangrove=dp$mangrove, moss=dp$moss, tsetse_habitat=dp$tsetse_habitat, cattle=dp$cattle,
  goat=dp$goat, horse=dp$horse, pig=dp$pig, buffalo=dp$buffalo, sheep=dp$sheep,
   total_animal=dp$total_animal), s = indexs)
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

formula <- as.formula(paste0("y ~ 0 + b0 +",paste(covariates, collapse = " + "),"+ f(s, model = spde)"))

res <- inla(formula,
  data = inla.stack.data(stk.full),
  family = "binomial",
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

ggplot(dpm) + geom_tile(aes(Latitude, Longitude, fill = value)) + geom_polygon(data =  border, aes(x = X, y = Y), color = "black",fill=NA)+
  facet_grid(~factor(variable,levels=c("2.5th percentile","Mean","97.5th percentile"))) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Prevalence",
    low = "blue", high = "orange",
    limits = c(0, 1)
  ) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) 
ggsave("Code/AbsencePresence/Results/Images/Predictions_SpecificCountry.png",width = 12, height = 8)

#####
log_tau <- -1.22  # from INLA
log_kappa <- -2.31

tau <- exp(log_tau)
kappa <- exp(log_kappa)
nu <- 1  # Fixed in SPDE model

sigma2 <- 1 / (4 * pi * kappa^2 * tau^2)   # marginal variance
range <- sqrt(8 * nu) / kappa             # practical range

matern_cov <- function(h, sigma2, kappa, nu = 1) {
  d <- h * kappa
  const <- sigma2 * (2^(1 - nu)) / gamma(nu)
  cov <- const * (d^nu) * besselK(d, nu)
  cov[is.nan(cov)] <- sigma2  # handle h = 0
  return(cov)
}

h <- seq(0, 5 * range, length.out = 200)  # distances
cov_values <- matern_cov(h, sigma2, kappa)

#plot and save h and cov_value in ggplot
cov_df <- data.frame(h = h, cov = cov_values)
ggplot(cov_df, aes(x = h, y = cov)) +
  geom_line() +
  labs(x = "Distance (km)", y = "Covariance") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())
    ggsave("Code/AbsencePresence/Results/Images/Covariance.png",width = 12, height = 8)


closest_point_index <- numeric(nrow(coo))
for(i in 1:nrow(coo)){
  closest_point_index[i] <- which.min(rowSums((coop-t(matrix(coo[i,],2,nrow(coop))))^2))
}

marginals <-  array(NA,c(43,2,nrow(coo)))
for(i in 1:nrow(coo)){
  marg <- res$marginals.fitted.values[index][[closest_point_index[i]]]
  marginals[,,i] <- marg
}


x0=marginals[,,presence==0]
x1=marginals[,,presence==1]

modes0<-numeric(length(x0[1,1,]))
for (i in (1:(length(x0[1,1,])))){
  modes0[i]<-x0[which.max(x0[,2,i]),1,i]
}

modes1<-numeric(length(x1[1,1,]))
for (i in (1:(length(x1[1,1,])))){
  modes1[i]<-x1[which.max(x1[,2,i]),1,i]
}


#plot violin plot of modes
modes<-c(modes0,modes1)
modes<-as.data.frame(modes)
modes$presence<-c(rep(0,length(modes0)),rep(1,length(modes1)))

ggplot(modes, aes(x = factor(presence), y = modes, fill=factor(presence))) +
  geom_violin() +
  labs(x = "True presence", y = "Mode of the prevalence posterior") +
  theme_bw(base_size=18) +
  scale_fill_okabeito()+
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none")
ggsave("Code/AbsencePresence/Results/Images/Modes.png",width = 12, height = 8)


