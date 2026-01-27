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

r_elv <- elevation_30s(country = "MWI", path = "Data/Covariates")

#make a grid of longitude and latitude values between extent of r_elv
lon <- seq(ext(r_elv)[1], ext(r_elv)[2], by = 0.1)
lat <- seq(ext(r_elv)[3], ext(r_elv)[4], by = 0.1)

#make a data frame of all combinations of lon and lat
data <- expand.grid(Longitude = lon, Latitude = lat)


data$elevation<-terra::extract(r_elv, data[, c("Longitude", "Latitude")])$MWI_elv_msk
#ELEVATION
r_elv <- elevation_30s(country = "MWI", path = "Data/Covariates")
data$elevation<-terra::extract(r_elv, data[, c("Longitude", "Latitude")])$MWI_elv_msk

#PRECIPITATION
r_prec<- worldclim_country(country = "MWI",path = "Data/Covariates", var = "prec")
data$precipitation<-(terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_1+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_2+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_3+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_4+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_5+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_6+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_7+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_8+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_9+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_10+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_11+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_prec_12)/12


#AVERAGE TEMPERATURE
r_tavg<- worldclim_country(country = "MWI",path = "Data/Covariates", var = "tavg")
data$tavg<-(terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_1+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_2+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_3+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_4+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_5+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_6+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_7+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_8+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_9+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_10+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_11+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tavg_12)/12

#MINIMUM TEMPERATURE
r_tmin<- worldclim_country(country = "MWI",path = "Data/Covariates", var = "tmin")
data$tmin<-(terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_1+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_2+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_3+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_4+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_5+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_6+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_7+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_8+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_9+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_10+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_11+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmin_12)/12

#MAXIMUM TEMPERATURE
r_tmax<- worldclim_country(country = "MWI",path = "Data/Covariates", var = "tmax")
data$tmax<-(terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_1+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_2+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_3+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_4+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_5+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_6+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_7+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_8+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_9+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_10+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_11+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$MWI_wc2.1_30s_tmax_12)/12

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

#remove rows with NA values
data<-na.omit(data)

# var1<-data$elevation
# var2<-data$precipitation
# correlation<-cor(var1, var2)

# Create a scatter plot
# ggplot(data, aes(x = var1, y = var2)) +
#   geom_point() +
#   labs(title = paste("Correlation between Elevation and Precipitation in Ethiopia: ", round(correlation, 2)),
#        x = "Elevation (m)",
#        y = "Precipitation (mm)")
# ggsave("Code/AbsencePresence/Correlation.png", width = 8, height = 6)

#compute correlation between all variables
correlation_matrix <- cor(data[, c("elevation", "precipitation", "tavg", "tmin", "tmax", "pop_den", "tree", "grassland", "shrub", "cropland", "built", "bare", "water", "wetland", "human_fp", "cattle", "goat", "horse", "pig", "sheep")], use = "pairwise.complete.obs")

# Create a heatmap of the correlation matrix
library(ggcorrplot)
ggcorrplot(correlation_matrix, 
           type = "lower", 
           lab = TRUE, 
           title = "Correlation matrix heatmap",
           colors = c("blue", "white", "red"),
           tl.cex = 10,
           tl.col = "black",
           tl.srt = 45,
           outline.color = "black") +
  theme_bw()
ggsave("Code/AbsencePresence/Results/Images/Correlations/Correlation_heatmap_MWI.png", width = 10, height = 8)

