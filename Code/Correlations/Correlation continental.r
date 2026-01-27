library(readxl)
library(ggplot2)
library(INLA)
library(sp)
library(sf)
library(afrilearndata)
library(concaveman)
library(geodata)
library(raster)

r_elv <- elevation_global(0.5,path = "Data/Covariates")

#make a grid of longitude and latitude values between extent of r_elv
lon <- seq(ext(r_elv)[1], ext(r_elv)[2], by = 0.1)
lat <- seq(ext(r_elv)[3], ext(r_elv)[4], by = 0.1)

#make a data frame of all combinations of lon and lat
data <- expand.grid(Longitude = lon, Latitude = lat)


data$elevation<-terra::extract(r_elv, data[, c("Longitude", "Latitude")])$wc2.1_30s_elev

#PRECIPITATION
r_prec<- worldclim_global(0.5,path = "Data/Covariates", var = "prec")
data$precipitation<-(terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_01+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_02+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_03+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_04+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_05+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_06+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_07+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_08+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_09+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_10+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_11+
terra::extract(r_prec, data[, c("Longitude", "Latitude")])$wc2.1_30s_prec_12)/12


#AVERAGE TEMPERATURE
r_tavg<- worldclim_global(0.5,path = "Data/Covariates", var = "tavg")
data$tavg<-(terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_01+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_02+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_03+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_04+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_05+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_06+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_07+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_08+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_09+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_10+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_11+
terra::extract(r_tavg, data[, c("Longitude", "Latitude")])$wc2.1_30s_tavg_12)/12

#MINIMUM TEMPERATURE
r_tmin<- worldclim_global(0.5,path = "Data/Covariates", var = "tmin")
data$tmin<-(terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_01+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_02+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_03+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_04+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_05+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_06+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_07+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_08+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_09+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_10+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_11+
terra::extract(r_tmin, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmin_12)/12

#MAXIMUM TEMPERATURE
r_tmax<- worldclim_global(0.5,path = "Data/Covariates", var = "tmax")
data$tmax<-(terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_01+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_02+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_03+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_04+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_05+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_06+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_07+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_08+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_09+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_10+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_11+
terra::extract(r_tmax, data[, c("Longitude", "Latitude")])$wc2.1_30s_tmax_12)/12

#POPULATION DENSITY
r_pop_den <- population(year=2010, res=0.5, path="Data/Covariates")
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

r_tsetse<-raster("Data/Covariates/tsenumbspec")
tsetse[tsetse > 1] <- 1
data$tsetse <-terra::extract(tsetse, data[, c("Longitude", "Latitude")])

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
correlation_matrix <- cor(data[, c("elevation", "precipitation", "tavg", "tmin", "tmax", "pop_den", "tree", "grassland", "shrub", "cropland", "built", "bare", "water", "wetland","mangrove", "human_fp", "cattle", "goat", "horse", "pig", "sheep","tsetse")], use = "pairwise.complete.obs")

# Function to rename covariates for better readability
rename_covariates <- function(names) {
  name_mapping <- c(
    "pop_den" = "Population density",
    "tmax" = "Maximum temperature", 
    "tmin" = "Minimum temperature",
    "tavg" = "Average temperature",
    "human_fp" = "Human footprint",
    "elevation" = "Elevation",
    "precipitation" = "Precipitation",
    "tree" = "Tree cover",
    "grassland" = "Grassland",
    "shrub" = "Shrubland", 
    "cropland" = "Cropland",
    "built" = "Built area",
    "bare" = "Bare ground",
    "water" = "Water",
    "wetland" = "Wetland",
    "mangrove" = "Mangrove",
    "cattle" = "Cattle density",
    "goat" = "Goat density",
    "horse" = "Horse density", 
    "pig" = "Pig density",
    "sheep" = "Sheep density",
    "tsetse" = "Tsetse presence"
  )
  
  for(old_name in names(name_mapping)) {
    names <- gsub(old_name, name_mapping[old_name], names, fixed = TRUE)
  }
  return(names)
}

# Rename correlation matrix row and column names
rownames(correlation_matrix) <- rename_covariates(rownames(correlation_matrix))
colnames(correlation_matrix) <- rename_covariates(colnames(correlation_matrix))

# Create an enhanced heatmap of the correlation matrix
library(ggcorrplot)
library(viridis)

p <- ggcorrplot(correlation_matrix, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 2.5,
           title = "Continental correlation matrix heatmap",
           colors = c("#2166AC", "white", "#B2182B"),  # Better diverging color scheme
           outline.color = "white") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold", margin = margin(b = 20)),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                       midpoint = 0, limit = c(-1, 1), 
                       name = "Correlation\ncoefficient",
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"))

print(p)
ggsave("Code/Correlations/Correlation_heatmap_Continental.png", plot = p, 
       width = 14, height = 12, dpi = 300, bg = "white")

