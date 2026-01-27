library(dplyr)
library(terra)
library(sf)
library(readxl)
library(raster)
library(geodata)

#clear workspace
rm(list = ls())

countries_to_infer=c("Nigeria")
dpm <- read.csv(paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"Predictions.csv"))

# Select relevant columns
dpm_mean <- dpm[dpm$variable == "Mean",]
#keep only the values where the mean is not NA
dpm_mean <- dpm_mean[!is.na(dpm_mean$value),]

# Convert to a spatial object
dpm_mean_sf <- st_as_sf(dpm_mean, coords = c("Longitude", "Latitude"), crs = 4326)


# Load state data
states <- gadm(country = "NGA", level = 1, path = tempdir())
country_sf <- st_as_sf(states)


#CATTLE DATA
r_cattle<-raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
r_buffalo<-raster("Data/Covariates/livestock/buffalo_2015/5_Bf_2015_Da.tif")
r_goat<-raster("Data/Covariates/livestock/goats_2015/5_Gt_2015_Da.tif")
r_horse<-raster("Data/Covariates/livestock/horses_2015/5_Ho_2015_Da.tif")
r_pig<-raster("Data/Covariates/livestock/pigs_2015/5_Pg_2015_Da.tif")
r_sheep<-raster("Data/Covariates/livestock/sheep_2015/5_Sh_2015_Da.tif")

dpm_mean_sf <- st_transform(dpm_mean_sf, crs = crs(r_cattle))
dpm_sp <- as(dpm_mean_sf, "Spatial")
cell_nums <- cellFromXY(r_cattle, coordinates(dpm_sp))
cell_coords <- xyFromCell(r_cattle, cell_nums)
result <- dpm_mean_sf %>%
  dplyr::mutate(cell_id = cell_nums,
         center_lon = cell_coords[, 1],
         center_lat = cell_coords[, 2])

# for locations in the same cell, take the mean of the values
averaged_result <- result %>%
  dplyr::group_by(center_lon, center_lat) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

# Convert back to sf object using the center coordinates
dpm_mean_sf <- st_as_sf(averaged_result, coords = c("center_lon", "center_lat"), crs = crs(r_cattle))
dpm_mean_sf <- st_transform(dpm_mean_sf, crs = st_crs(country_sf))
# Join the data with the state polygons
country_sf <- st_make_valid(country_sf)
dpm_mean_sf <- st_join(dpm_mean_sf, country_sf, join = st_within)

# Ensure we have point geometries before extracting coordinates
dpm_mean_sf <- st_cast(dpm_mean_sf, "POINT")

# Extract coordinates and add as columns
coords <- st_coordinates(dpm_mean_sf)
dpm_mean_sf$longitude <- coords[, 1]
dpm_mean_sf$latitude <- coords[, 2]

# Now select and rename as needed
dpm_mean_sf <- dpm_mean_sf %>%
  dplyr::select(longitude, latitude, mean_value, NAME_1) %>%
  dplyr::rename(state = NAME_1, value = mean_value)
#extract cattle data and add to dpm_mean_sf
dpm_mean_sf$cattle_total <- raster::extract(r_cattle, dpm_mean_sf)
dpm_mean_sf$cattle_total[is.na(dpm_mean_sf$cattle_total)] <- 0 # set NA values to 0
dpm_mean_sf$buffalo_total <- raster::extract(r_buffalo, dpm_mean_sf)
dpm_mean_sf$buffalo_total[is.na(dpm_mean_sf$buffalo_total)] <- 0 # set NA values to 0
dpm_mean_sf$goat_total <- raster::extract(r_goat, dpm_mean_sf)
dpm_mean_sf$goat_total[is.na(dpm_mean_sf$goat_total)] <- 0 # set NA values to 0
dpm_mean_sf$horse_total <- raster::extract(r_horse, dpm_mean_sf)
dpm_mean_sf$horse_total[is.na(dpm_mean_sf$horse_total)] <- 0 # set NA values to 0
dpm_mean_sf$pig_total <- raster::extract(r_pig, dpm_mean_sf)
dpm_mean_sf$pig_total[is.na(dpm_mean_sf$pig_total)] <- 0 # set NA values to 0
dpm_mean_sf$sheep_total <- raster::extract(r_sheep, dpm_mean_sf)
dpm_mean_sf$sheep_total[is.na(dpm_mean_sf$sheep_total)] <- 0 # set NA values to 0

#set cattle to value*cattle
dpm_mean_sf$cattle <- dpm_mean_sf$value * dpm_mean_sf$cattle_total
dpm_mean_sf$buffalo <- dpm_mean_sf$value * dpm_mean_sf$buffalo_total
dpm_mean_sf$goat <- dpm_mean_sf$value * dpm_mean_sf$goat_total
dpm_mean_sf$horse <- dpm_mean_sf$value * dpm_mean_sf$horse_total
dpm_mean_sf$pig <- dpm_mean_sf$value * dpm_mean_sf$pig_total
dpm_mean_sf$sheep <- dpm_mean_sf$value * dpm_mean_sf$sheep_total


#Average by state
dpm_mean_state <- dpm_mean_sf %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
            longitude = mean(longitude, na.rm = TRUE),
            latitude = mean(latitude, na.rm = TRUE),
            cattle = sum(cattle, na.rm = TRUE),
            cattle_total = sum(cattle_total, na.rm = TRUE),
            buffalo = sum(buffalo, na.rm = TRUE),
            buffalo_total = sum(buffalo_total, na.rm = TRUE),
            goat = sum(goat, na.rm = TRUE),
            goat_total = sum(goat_total, na.rm = TRUE),
            horse = sum(horse, na.rm = TRUE),
            horse_total = sum(horse_total, na.rm = TRUE),
            pig = sum(pig, na.rm = TRUE),
            pig_total = sum(pig_total, na.rm = TRUE),
            sheep = sum(sheep, na.rm = TRUE),
            sheep_total = sum(sheep_total, na.rm = TRUE))

#remove geometry from dpm_mean_state
dpm_mean_state <- st_drop_geometry(dpm_mean_state)
#remove NA state
dpm_mean_state <- dpm_mean_state[!is.na(dpm_mean_state$state),]






# Load the AAT data
data_path <- paste0("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
data <- read_excel(data_path)

#remove row if longitude or latitude is NA
data <- data[!is.na(data$Longitude) & !is.na(data$Latitude),]

#if number of infections is NA, set to Number_of_animals_tested*TPR
data$Number_of_infections[is.na(data$Number_of_infections)] <- round(data$Number_of_animal_tested[is.na(data$Number_of_infections)]*data$TPR[is.na(data$Number_of_infections)]/100)

positive <- round(data$Number_of_infections)
sample_size <- data$Number_of_animal_tested

data$prevalence <- positive/sample_size
data$number_of_animals_tested <- data$Number_of_animal_tested

#work out which state each point in data is in
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)
data_sf <- st_transform(data_sf, crs = st_crs(country_sf))
# Join the data with the state polygons
data_sf <- st_join(data_sf, country_sf, join = st_within)
# Select relevant columns
data_sf$longitude <- st_coordinates(data_sf)[, 1]
data_sf$latitude <- st_coordinates(data_sf)[, 2]

# Now select and rename as needed
data_sf <- data_sf %>%
  dplyr::select(longitude, latitude, prevalence, number_of_animals_tested,NAME_1) %>%
  dplyr::rename(state = NAME_1)

#Average by state and include number of data points
data_mean_state <- data_sf %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(mean_prevalence = mean(prevalence, na.rm = TRUE),
            longitude = mean(longitude, na.rm = TRUE),
            latitude = mean(latitude, na.rm = TRUE),
            number_of_studies = n(),
            number_of_animals_tested = sum(number_of_animals_tested, na.rm = TRUE))

#add mean_prevalence from data_mean_state to dpm_mean_state according to state
dpm_mean_state <- dpm_mean_state %>%
  dplyr::left_join(data_mean_state, by = "state") %>%
  dplyr::rename(mean_value_data = mean_prevalence)

#remove longitude.y, latitude.y
dpm_mean_state <- dpm_mean_state %>%
  dplyr::select(-longitude.y, -latitude.y)

#change all NA values in n to 0
dpm_mean_state$number_of_studies[is.na(dpm_mean_state$number_of_studies)] <- 0

#move mean_value_data to after mean_value
dpm_mean_state <- dpm_mean_state %>%
  dplyr::select(state,mean_value_data, mean_value, cattle, cattle_total, buffalo, buffalo_total, goat, goat_total, horse, horse_total, pig, pig_total, sheep, sheep_total, number_of_animals_tested,number_of_studies)

#add a final row called "Total" with the total number of animals in each category
total_row <- data.frame(
  state = "Total",
  mean_value_data = mean(dpm_mean_state$mean_value_data, na.rm = TRUE),
  mean_value = mean(dpm_mean_state$mean_value, na.rm = TRUE),
  cattle = sum(dpm_mean_state$cattle, na.rm = TRUE),
  cattle_total = sum(dpm_mean_state$cattle_total, na.rm = TRUE),
  buffalo = sum(dpm_mean_state$buffalo, na.rm = TRUE),
  buffalo_total = sum(dpm_mean_state$buffalo_total, na.rm = TRUE),
  goat = sum(dpm_mean_state$goat, na.rm = TRUE),
  goat_total = sum(dpm_mean_state$goat_total, na.rm = TRUE),
  horse = sum(dpm_mean_state$horse, na.rm = TRUE),
  horse_total = sum(dpm_mean_state$horse_total, na.rm = TRUE),
  pig = sum(dpm_mean_state$pig, na.rm = TRUE),
  pig_total = sum(dpm_mean_state$pig_total, na.rm = TRUE),
  sheep = sum(dpm_mean_state$sheep, na.rm = TRUE),
  sheep_total = sum(dpm_mean_state$sheep_total, na.rm = TRUE),
  number_of_studies = sum(dpm_mean_state$number_of_studies, na.rm = TRUE),
  number_of_animals_tested = sum(dpm_mean_state$number_of_animals_tested, na.rm = TRUE)
)

# Append the total row to the data frame
dpm_mean_state <- rbind(dpm_mean_state, total_row)

#round mean_value_data to 3 significant figures
dpm_mean_state$mean_value_data <- round(as.numeric(dpm_mean_state$mean_value_data), 3)
dpm_mean_state$mean_value <- round(as.numeric(dpm_mean_state$mean_value), 3)

#round cattle, buffalo, goat, horse, pig, sheep to 0 decimal places
dpm_mean_state$cattle <- round(as.numeric(dpm_mean_state$cattle), 0)
dpm_mean_state$cattle_total <- round(as.numeric(dpm_mean_state$cattle_total), 0)
dpm_mean_state$buffalo <- round(as.numeric(dpm_mean_state$buffalo), 0)
dpm_mean_state$buffalo_total <- round(as.numeric(dpm_mean_state$buffalo_total), 0)
dpm_mean_state$goat <- round(as.numeric(dpm_mean_state$goat), 0)
dpm_mean_state$goat_total <- round(as.numeric(dpm_mean_state$goat_total), 0)
dpm_mean_state$horse <- round(as.numeric(dpm_mean_state$horse), 0)
dpm_mean_state$horse_total <- round(as.numeric(dpm_mean_state$horse_total), 0)
dpm_mean_state$pig <- round(as.numeric(dpm_mean_state$pig), 0)
dpm_mean_state$pig_total <- round(as.numeric(dpm_mean_state$pig_total), 0)
dpm_mean_state$sheep <- round(as.numeric(dpm_mean_state$sheep), 0)
dpm_mean_state$sheep_total <- round(as.numeric(dpm_mean_state$sheep_total), 0)

#change NA values in mean_value_data to "No data"
dpm_mean_state$mean_value_data[is.na(dpm_mean_state$mean_value_data)] <- "No data"
dpm_mean_state$number_of_animals_tested[is.na(dpm_mean_state$number_of_animals_tested)] <- 0

#write to csv
write.csv(dpm_mean_state,paste0("Code/Prevalence/Bovine BCT/Results/",countries_to_infer,"/",countries_to_infer,"PredictionsByState.csv"), row.names=FALSE)
