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

states <- gadm(country = "NGA", level = 1, path = tempdir())
country_sf <- st_as_sf(states)


#CATTLE DATA
r_cattle<-raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
# r_buffalo<-raster("Data/Covariates/livestock/buffalo_2015/5_Bf_2015_Da.tif")
# r_goat<-raster("Data/Covariates/livestock/goats_2015/5_Gt_2015_Da.tif")
# r_horse<-raster("Data/Covariates/livestock/horses_2015/5_Ho_2015_Da.tif")
# r_pig<-raster("Data/Covariates/livestock/pigs_2015/5_Pg_2015_Da.tif")
# r_sheep<-raster("Data/Covariates/livestock/sheep_2015/5_Sh_2015_Da.tif")


##############

# Count number of projection files
n_datasets <- length(list.files("Code/Prevalence/Bovine BCT and PCR/Projections_NGA/", pattern = "Projections_model_.*\\.csv"))
cat("Found", n_datasets, "projection files for Nigeria\n")

for (i in 1:n_datasets){
print(i)
dpm <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Projections_NGA/Projections_model_",i,".csv"))

# Prepare for cattle data extraction - Mean
dpm_mean <- dpm[dpm$variable == "Mean",]
#keep only the values where the mean is not NA
dpm_mean <- dpm_mean[!is.na(dpm_mean$value),]

# Prepare for cattle data extraction - 2.5th percentile
dpm_lower <- dpm[dpm$variable == "2.5th percentile",]
dpm_lower <- dpm_lower[!is.na(dpm_lower$value),]

# Prepare for cattle data extraction - 97.5th percentile
dpm_upper <- dpm[dpm$variable == "97.5th percentile",]
dpm_upper <- dpm_upper[!is.na(dpm_upper$value),]

# Swap longitude and latitude (they are reversed in the data) - Mean
temp <- dpm_mean$Longitude
dpm_mean$Longitude <- dpm_mean$Latitude
dpm_mean$Latitude <- temp

# Swap longitude and latitude - 2.5th percentile
temp <- dpm_lower$Longitude
dpm_lower$Longitude <- dpm_lower$Latitude
dpm_lower$Latitude <- temp

# Swap longitude and latitude - 97.5th percentile
temp <- dpm_upper$Longitude
dpm_upper$Longitude <- dpm_upper$Latitude
dpm_upper$Latitude <- temp

# Convert to spatial objects
dpm_mean_sf <- st_as_sf(dpm_mean, coords = c("Longitude", "Latitude"), crs = 4326)
dpm_lower_sf <- st_as_sf(dpm_lower, coords = c("Longitude", "Latitude"), crs = 4326)
dpm_upper_sf <- st_as_sf(dpm_upper, coords = c("Longitude", "Latitude"), crs = 4326)

###### RUN CATTLE DATA EXTRACTION - MEAN

dpm_mean_sf <- st_transform(dpm_mean_sf, crs = crs(r_cattle))
dpm_sp <- as(dpm_mean_sf, "Spatial")
cell_nums <- cellFromXY(r_cattle, coordinates(dpm_sp))
cell_coords <- xyFromCell(r_cattle, cell_nums)
result <- dpm_mean_sf %>%
  mutate(cell_id = cell_nums,
         center_lon = cell_coords[, 1],
         center_lat = cell_coords[, 2])

# for locations in the same cell, take the mean of the values
averaged_result <- result %>%
  group_by(center_lon, center_lat) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# rename center_lon and center_lat to longitude and latitude
dpm_mean_sf <- averaged_result %>%
  rename(longitude = center_lon, latitude = center_lat, 
         value = mean_value)

###### RUN CATTLE DATA EXTRACTION - 2.5TH PERCENTILE

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
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# rename center_lon and center_lat to longitude and latitude
dpm_lower_sf <- averaged_result_lower %>%
  rename(longitude = center_lon, latitude = center_lat, 
         value = mean_value)

###### RUN CATTLE DATA EXTRACTION - 97.5TH PERCENTILE

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
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# rename center_lon and center_lat to longitude and latitude
dpm_upper_sf <- averaged_result_upper %>%
  rename(longitude = center_lon, latitude = center_lat, 
         value = mean_value)

# Convert back to sf objects (coordinates are in cattle raster CRS)
dpm_mean_sf <- st_as_sf(dpm_mean_sf, coords = c("longitude", "latitude"), crs = crs(r_cattle))
dpm_lower_sf <- st_as_sf(dpm_lower_sf, coords = c("longitude", "latitude"), crs = crs(r_cattle))
dpm_upper_sf <- st_as_sf(dpm_upper_sf, coords = c("longitude", "latitude"), crs = crs(r_cattle))

dpm_mean_sf <- st_transform(dpm_mean_sf, crs = st_crs(country_sf))
dpm_lower_sf <- st_transform(dpm_lower_sf, crs = st_crs(country_sf))
dpm_upper_sf <- st_transform(dpm_upper_sf, crs = st_crs(country_sf))

# Join the data with the state polygons
country_sf <- st_make_valid(country_sf)
dpm_mean_sf <- st_join(dpm_mean_sf, country_sf, join = st_within)
dpm_lower_sf <- st_join(dpm_lower_sf, country_sf, join = st_within)
dpm_upper_sf <- st_join(dpm_upper_sf, country_sf, join = st_within)

# Now select and rename as needed
dpm_mean_sf <- dpm_mean_sf %>%
  dplyr::select(value, NAME_1) %>%
  rename(state = NAME_1)

dpm_lower_sf <- dpm_lower_sf %>%
  dplyr::select(value, NAME_1) %>%
  rename(state = NAME_1)

dpm_upper_sf <- dpm_upper_sf %>%
  dplyr::select(value, NAME_1) %>%
  rename(state = NAME_1)

dpm_mean_sf <- st_cast(dpm_mean_sf, "POINT")
dpm_lower_sf <- st_cast(dpm_lower_sf, "POINT")
dpm_upper_sf <- st_cast(dpm_upper_sf, "POINT")

#extract cattle data and add to dpm_mean_sf (MEAN)
dpm_mean_sf$cattle <- raster::extract(r_cattle, dpm_mean_sf)
dpm_mean_sf$cattle[is.na(dpm_mean_sf$cattle)] <- 0
# dpm_mean_sf$buffalo <- raster::extract(r_buffalo, dpm_mean_sf)
# dpm_mean_sf$buffalo[is.na(dpm_mean_sf$buffalo)] <- 0
# dpm_mean_sf$goat <- raster::extract(r_goat, dpm_mean_sf)
# dpm_mean_sf$goat[is.na(dpm_mean_sf$goat)] <- 0
# dpm_mean_sf$horse <- raster::extract(r_horse, dpm_mean_sf)
# dpm_mean_sf$horse[is.na(dpm_mean_sf$horse)] <- 0
# dpm_mean_sf$pig <- raster::extract(r_pig, dpm_mean_sf)
# dpm_mean_sf$pig[is.na(dpm_mean_sf$pig)] <- 0
# dpm_mean_sf$sheep <- raster::extract(r_sheep, dpm_mean_sf)
# dpm_mean_sf$sheep[is.na(dpm_mean_sf$sheep)] <- 0

#extract cattle data and add to dpm_lower_sf (2.5TH PERCENTILE)
dpm_lower_sf$cattle <- raster::extract(r_cattle, dpm_lower_sf)
dpm_lower_sf$cattle[is.na(dpm_lower_sf$cattle)] <- 0
# dpm_lower_sf$buffalo <- raster::extract(r_buffalo, dpm_lower_sf)
# dpm_lower_sf$buffalo[is.na(dpm_lower_sf$buffalo)] <- 0
# dpm_lower_sf$goat <- raster::extract(r_goat, dpm_lower_sf)
# dpm_lower_sf$goat[is.na(dpm_lower_sf$goat)] <- 0
# dpm_lower_sf$horse <- raster::extract(r_horse, dpm_lower_sf)
# dpm_lower_sf$horse[is.na(dpm_lower_sf$horse)] <- 0
# dpm_lower_sf$pig <- raster::extract(r_pig, dpm_lower_sf)
# dpm_lower_sf$pig[is.na(dpm_lower_sf$pig)] <- 0
# dpm_lower_sf$sheep <- raster::extract(r_sheep, dpm_lower_sf)
# dpm_lower_sf$sheep[is.na(dpm_lower_sf$sheep)] <- 0

#extract cattle data and add to dpm_upper_sf (97.5TH PERCENTILE)
dpm_upper_sf$cattle <- raster::extract(r_cattle, dpm_upper_sf)
dpm_upper_sf$cattle[is.na(dpm_upper_sf$cattle)] <- 0
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
# dpm_mean_sf$pig <- dpm_mean_sf$value * dpm_mean_sf$pig
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

#Average by state - MEAN
dpm_mean_state <- dpm_mean_sf %>%
  group_by(state) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            cattle_mean = sum(cattle, na.rm = TRUE))
            # buffalo_mean = sum(buffalo, na.rm = TRUE),
            # goat_mean = sum(goat, na.rm = TRUE),
            # horse_mean = sum(horse, na.rm = TRUE),
            # pig_mean = sum(pig, na.rm = TRUE),
            # sheep_mean = sum(sheep, na.rm = TRUE))

#Average by state - 2.5TH PERCENTILE
dpm_lower_state <- dpm_lower_sf %>%
  group_by(state) %>%
  summarise(lower_value = mean(value, na.rm = TRUE),
            cattle_lower = sum(cattle, na.rm = TRUE))
            # buffalo_lower = sum(buffalo, na.rm = TRUE),
            # goat_lower = sum(goat, na.rm = TRUE),
            # horse_lower = sum(horse, na.rm = TRUE),
            # pig_lower = sum(pig, na.rm = TRUE),
            # sheep_lower = sum(sheep, na.rm = TRUE))

#Average by state - 97.5TH PERCENTILE
dpm_upper_state <- dpm_upper_sf %>%
  group_by(state) %>%
  summarise(upper_value = mean(value, na.rm = TRUE),
            cattle_upper = sum(cattle, na.rm = TRUE))
            # buffalo_upper = sum(buffalo, na.rm = TRUE),
            # goat_upper = sum(goat, na.rm = TRUE),
            # horse_upper = sum(horse, na.rm = TRUE),
            # pig_upper = sum(pig, na.rm = TRUE),
            # sheep_upper = sum(sheep, na.rm = TRUE))

#remove geometry from all state data
dpm_mean_state <- st_drop_geometry(dpm_mean_state)
dpm_lower_state <- st_drop_geometry(dpm_lower_state)
dpm_upper_state <- st_drop_geometry(dpm_upper_state)

#remove NA state from all datasets
dpm_mean_state <- dpm_mean_state[!is.na(dpm_mean_state$state),]
dpm_lower_state <- dpm_lower_state[!is.na(dpm_lower_state$state),]
dpm_upper_state <- dpm_upper_state[!is.na(dpm_upper_state$state),]

# Merge all three datasets by state
dpm_combined_state <- merge(dpm_mean_state, dpm_lower_state, by = "state")
dpm_combined_state <- merge(dpm_combined_state, dpm_upper_state, by = "state")

#round all livestock numbers to 0 decimal places
livestock_cols <- grep("cattle", names(dpm_combined_state), value = TRUE)
dpm_combined_state[livestock_cols] <- lapply(dpm_combined_state[livestock_cols], function(x) round(as.numeric(x), 0))

#add a final row called "Total" with the total number of animals in each category
total_row <- data.frame(
  state = "Total",
  mean_value = mean(dpm_combined_state$mean_value, na.rm = TRUE),
  cattle_mean = sum(dpm_combined_state$cattle_mean, na.rm = TRUE),
  # buffalo_mean = sum(dpm_combined_state$buffalo_mean, na.rm = TRUE),
  # goat_mean = sum(dpm_combined_state$goat_mean, na.rm = TRUE),
  # horse_mean = sum(dpm_combined_state$horse_mean, na.rm = TRUE),
  # pig_mean = sum(dpm_combined_state$pig_mean, na.rm = TRUE),
  # sheep_mean = sum(dpm_combined_state$sheep_mean, na.rm = TRUE),
  lower_value = mean(dpm_combined_state$lower_value, na.rm = TRUE),
  cattle_lower = sum(dpm_combined_state$cattle_lower, na.rm = TRUE),
  # buffalo_lower = sum(dpm_combined_state$buffalo_lower, na.rm = TRUE),
  # goat_lower = sum(dpm_combined_state$goat_lower, na.rm = TRUE),
  # horse_lower = sum(dpm_combined_state$horse_lower, na.rm = TRUE),
  # pig_lower = sum(dpm_combined_state$pig_lower, na.rm = TRUE),
  # sheep_lower = sum(dpm_combined_state$sheep_lower, na.rm = TRUE),
  upper_value = mean(dpm_combined_state$upper_value, na.rm = TRUE),
  cattle_upper = sum(dpm_combined_state$cattle_upper, na.rm = TRUE)
  # buffalo_upper = sum(dpm_combined_state$buffalo_upper, na.rm = TRUE),
  # goat_upper = sum(dpm_combined_state$goat_upper, na.rm = TRUE),
  # horse_upper = sum(dpm_combined_state$horse_upper, na.rm = TRUE),
  # pig_upper = sum(dpm_combined_state$pig_upper, na.rm = TRUE),
  # sheep_upper = sum(dpm_combined_state$sheep_upper, na.rm = TRUE)
)

# Append the total row to the data frame
dpm_combined_state <- rbind(dpm_combined_state, total_row)

#######


# sym <- check_symmetry(dpm_mean, dpm_lower, dpm_upper)
# print(sym$frac_flagged)  # % of units with noticeable skew


# ggplot(dpm) + geom_tile(aes(Latitude, Longitude, fill = value)) +
#   facet_grid(~factor(variable,levels=c("2.5th percentile","Mean","97.5th percentile"))) +
#   coord_fixed(ratio = 1) +
#   scale_fill_gradient(
#     name = "Prevalence",
#     low = "blue", high = "orange",
#     limits = c(0, 1),
#     guide = guide_colorbar(order = 2)
#   ) +
#   new_scale_fill() +  # Add a new fill scale
#   #geom_tile(data=water_cover, aes(Longitude, Latitude,fill=water_cover,alpha=water_cover),inherit.aes = FALSE)+
#   scale_fill_manual(
#   name = "",  # Colorbar title
#   values = c("High uncertainty" = "grey50", "Water" = "steelblue1"),  # Specific colors for categories
#   limits = c("High uncertainty", "Water"),
#   guide = guide_legend(order = 1)
# )+
# guides(alpha = "none")+
#   #geom_polygon(data =  border, aes(x = X, y = Y), color = "black",fill=NA)+
#   theme_bw()+
#   theme(axis.line = element_line(colour = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     legend.box = "vertical")
# ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/Prevalence_plots/Projections_model_",i,".png"),width = 12, height = 8)

write.csv(dpm_combined_state,paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA_all_PCR/Cattle_at_risk/Projections_model_",i,".csv"), row.names=FALSE)


}

