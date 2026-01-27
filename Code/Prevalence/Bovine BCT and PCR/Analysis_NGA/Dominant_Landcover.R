rm(list=ls())

library(terra)
library(geodata)
library(ggplot2)
library(viridis)
library(dplyr)
library(sf)
library(raster)

# Set country for analysis
country_code <- "NGA"
country_name <- "Nigeria"

cat("Loading land cover data for", country_name, "...\n")

# Load all land cover types using the same approach as Model_selection_NGA.r
r_tree <- landcover("trees", path="Data/Covariates")
r_grassland <- landcover("grassland", path="Data/Covariates")
r_shrub <- landcover("shrubs", path="Data/Covariates")
r_cropland <- landcover("cropland", path="Data/Covariates")
r_built <- landcover("built", path="Data/Covariates")
r_bare <- landcover("bare", path="Data/Covariates")
r_water <- landcover("water", path="Data/Covariates")
r_wetland <- landcover("wetland", path="Data/Covariates")
r_mangrove <- landcover("mangroves", path="Data/Covariates")
r_moss <- landcover("moss", path="Data/Covariates")

# Get country boundary for cropping
country_bounds <- geodata::gadm(country = country_code, level = 0, path = "Data/Covariates")

# Crop all rasters to country boundary
cat("Cropping rasters to", country_name, "boundary...\n")
r_tree_crop <- terra::crop(r_tree, country_bounds)
r_grassland_crop <- terra::crop(r_grassland, country_bounds)
r_shrub_crop <- terra::crop(r_shrub, country_bounds)
r_cropland_crop <- terra::crop(r_cropland, country_bounds)
r_built_crop <- terra::crop(r_built, country_bounds)
r_bare_crop <- terra::crop(r_bare, country_bounds)
r_water_crop <- terra::crop(r_water, country_bounds)
r_wetland_crop <- terra::crop(r_wetland, country_bounds)
r_mangrove_crop <- terra::crop(r_mangrove, country_bounds)
r_moss_crop <- terra::crop(r_moss, country_bounds)

# Mask to country boundaries (set pixels outside country to NA)
r_tree_mask <- terra::mask(r_tree_crop, country_bounds)
r_grassland_mask <- terra::mask(r_grassland_crop, country_bounds)
r_shrub_mask <- terra::mask(r_shrub_crop, country_bounds)
r_cropland_mask <- terra::mask(r_cropland_crop, country_bounds)
r_built_mask <- terra::mask(r_built_crop, country_bounds)
r_bare_mask <- terra::mask(r_bare_crop, country_bounds)
r_water_mask <- terra::mask(r_water_crop, country_bounds)
r_wetland_mask <- terra::mask(r_wetland_crop, country_bounds)
r_mangrove_mask <- terra::mask(r_mangrove_crop, country_bounds)
r_moss_mask <- terra::mask(r_moss_crop, country_bounds)

# Stack all land cover rasters
landcover_stack <- c(r_tree_mask, r_grassland_mask, r_shrub_mask, r_cropland_mask, 
                     r_built_mask, r_bare_mask, r_water_mask, r_wetland_mask, 
                     r_mangrove_mask, r_moss_mask)

# Set names for the stack
names(landcover_stack) <- c("trees", "grassland", "shrubs", "cropland", "built", 
                           "bare", "water", "wetland", "mangroves", "moss")

cat("Calculating dominant land cover type for each pixel...\n")

# Find the dominant land cover type (maximum value) for each pixel
dominant_landcover <- terra::which.max(landcover_stack)

# Create a data frame for plotting
dominant_df <- terra::as.data.frame(dominant_landcover, xy = TRUE)
colnames(dominant_df) <- c("longitude", "latitude", "dominant_type")

# Remove NA values
dominant_df <- dominant_df[!is.na(dominant_df$dominant_type), ]

# Convert numeric codes to land cover names
landcover_names <- c("Trees", "Grassland", "Shrubs", "Cropland", "Built", 
                     "Bare", "Water", "Wetland", "Mangroves", "Moss")
dominant_df$landcover <- factor(landcover_names[dominant_df$dominant_type], 
                               levels = landcover_names)

cat("Creating visualization...\n")

# Define colors for each land cover type
landcover_colors <- c(
  "Trees" = "#228B22",        # Forest Green
  "Grassland" = "#9ACD32",    # Yellow Green
  "Shrubs" = "#8FBC8F",       # Dark Sea Green
  "Cropland" = "#FFD700",     # Gold
  "Built" = "#696969",        # Dim Gray
  "Bare" = "#D2B48C",         # Tan
  "Water" = "#4169E1",        # Royal Blue
  "Wetland" = "#00CED1",      # Dark Turquoise
  "Mangroves" = "#008B8B",    # Dark Cyan
  "Moss" = "#32CD32"          # Lime Green
)

# Create the heatmap
p <- ggplot(dominant_df, aes(x = longitude, y = latitude, fill = landcover)) +
  geom_raster() +
  scale_fill_manual(values = landcover_colors, name = "Dominant\nland cover") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, hjust = 0.5),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = paste("Dominant land cover type -", country_name),
    x = "Longitude",
    y = "Latitude",
    caption = "Each pixel colored by the land cover type with highest proportion"
  )

# Print summary statistics
cat("\nSummary of dominant land cover types in", country_name, ":\n")
summary_table <- table(dominant_df$landcover)
summary_df <- data.frame(
  LandCover = names(summary_table),
  PixelCount = as.numeric(summary_table),
  Percentage = round(as.numeric(summary_table) / sum(summary_table) * 100, 2)
)
print(summary_df)

# Display the plot
print(p)


ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/Dominant_Landcover_", country_code, ".pdf"), 
       plot = p, width = 12, height = 8)

cat("\nPlot saved as Dominant_Landcover_", country_code, ".png and .pdf in Analysis_NGA folder\n")

#===============================================================================
# ADDITIONAL COVARIATE VISUALIZATIONS
#===============================================================================

cat("\n=== Loading and visualizing other covariates ===\n")

# ELEVATION
cat("Processing elevation data...\n")
r_elv <- elevation_30s(country = country_code, path = "Data/Covariates")
r_elv_crop <- terra::crop(r_elv, country_bounds)
r_elv_mask <- terra::mask(r_elv_crop, country_bounds)

elv_df <- terra::as.data.frame(r_elv_mask, xy = TRUE)
colnames(elv_df) <- c("longitude", "latitude", "elevation")
elv_df <- elv_df[!is.na(elv_df$elevation), ]

p_elv <- ggplot(elv_df, aes(x = longitude, y = latitude, fill = elevation)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Elevation\n(m)") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = paste("Elevation -", country_name), x = "Longitude", y = "Latitude")

ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/Elevation_", country_code, ".pdf"), plot = p_elv, width = 12, height = 8)

# PRECIPITATION
cat("Processing precipitation data...\n")
r_prec <- worldclim_country(country = country_code, path = "Data/Covariates", var = "prec")
# Calculate annual average precipitation
r_prec_annual <- sum(r_prec) / 12
r_prec_crop <- terra::crop(r_prec_annual, country_bounds)
r_prec_mask <- terra::mask(r_prec_crop, country_bounds)

prec_df <- terra::as.data.frame(r_prec_mask, xy = TRUE)
colnames(prec_df) <- c("longitude", "latitude", "precipitation")
prec_df <- prec_df[!is.na(prec_df$precipitation), ]

p_prec <- ggplot(prec_df, aes(x = longitude, y = latitude, fill = precipitation)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Annual\nprecipitation\n(mm)", option = "B") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = paste("Annual precipitation -", country_name), x = "Longitude", y = "Latitude")

ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/Precipitation_", country_code, ".pdf"), plot = p_prec, width = 12, height = 8)

# TEMPERATURE (Average)
cat("Processing temperature data...\n")
r_tavg <- worldclim_country(country = country_code, path = "Data/Covariates", var = "tavg")
r_tavg_annual <- sum(r_tavg) / 12
r_tavg_crop <- terra::crop(r_tavg_annual, country_bounds)
r_tavg_mask <- terra::mask(r_tavg_crop, country_bounds)

tavg_df <- terra::as.data.frame(r_tavg_mask, xy = TRUE)
colnames(tavg_df) <- c("longitude", "latitude", "temperature")
tavg_df <- tavg_df[!is.na(tavg_df$temperature), ]

p_tavg <- ggplot(tavg_df, aes(x = longitude, y = latitude, fill = temperature)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Average\ntemperature\n(°C)", option = "C") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = paste("Average temperature -", country_name), x = "Longitude", y = "Latitude")

ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/Temperature_", country_code, ".pdf"), plot = p_tavg, width = 12, height = 8)

# POPULATION DENSITY
cat("Processing population density data...\n")
r_pop_den <- population(year = 2010, res = 0.5, path = "Data/Covariates")
r_pop_crop <- terra::crop(r_pop_den, country_bounds)
r_pop_mask <- terra::mask(r_pop_crop, country_bounds)

pop_df <- terra::as.data.frame(r_pop_mask, xy = TRUE)
colnames(pop_df) <- c("longitude", "latitude", "population_density")
pop_df <- pop_df[!is.na(pop_df$population_density), ]

# Add small constant to avoid log(0) and calculate log population density
pop_df$log_pop_density <- log10(pop_df$population_density + 1)

p_pop <- ggplot(pop_df, aes(x = longitude, y = latitude, fill = log_pop_density)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Log of\npopulation\ndensity", na.value = "transparent", option = "E") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  ) +
  labs(title = paste("Log of population density -", country_name), x = "Longitude", y = "Latitude")

# Save with high DPI for sharp rasters
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/PopulationDensity_", country_code, ".png"), 
       plot = p_pop, width = 12, height = 8, dpi = 600, type = "cairo-png")
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/PopulationDensity_", country_code, ".tiff"), 
       plot = p_pop, width = 12, height = 8, dpi = 600, compression = "none")
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/PopulationDensity_", country_code, ".pdf"), 
       plot = p_pop, width = 12, height = 8, device = "pdf", dpi = 600, 
       useDingbats = FALSE, compress = FALSE)

# HUMAN FOOTPRINT
cat("Processing human footprint data...\n")
r_human_fp <- footprint(year = 2009, res = 30, path = "Data/Covariates")
r_human_fp_crop <- terra::crop(r_human_fp, country_bounds)
r_human_fp_mask <- terra::mask(r_human_fp_crop, country_bounds)

hfp_df <- terra::as.data.frame(r_human_fp_mask, xy = TRUE)
colnames(hfp_df) <- c("longitude", "latitude", "human_footprint")
hfp_df <- hfp_df[!is.na(hfp_df$human_footprint), ]

p_hfp <- ggplot(hfp_df, aes(x = longitude, y = latitude, fill = human_footprint)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Human\nfootprint\nindex", option = "D") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = paste("Human footprint -", country_name), x = "Longitude", y = "Latitude")

ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/HumanFootprint_", country_code, ".pdf"), plot = p_hfp, width = 12, height = 8)

# CATTLE DENSITY
cat("Processing cattle density data...\n")
r_cattle <- raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
r_cattle <- terra::rast(r_cattle)  # Convert to SpatRaster
r_cattle_crop <- terra::crop(r_cattle, country_bounds)
r_cattle_mask <- terra::mask(r_cattle_crop, country_bounds)

cattle_df <- terra::as.data.frame(r_cattle_mask, xy = TRUE)
colnames(cattle_df) <- c("longitude", "latitude", "cattle_density")
cattle_df <- cattle_df[!is.na(cattle_df$cattle_density), ]

# Convert country bounds to sf for plotting
country_sf <- sf::st_as_sf(country_bounds)

p_cattle <- ggplot(cattle_df, aes(x = longitude, y = latitude, fill = cattle_density)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cattle\ndensity", trans = "sqrt", na.value = "transparent") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  ) +
  labs(title = paste("Cattle density -", country_name), x = "Longitude", y = "Latitude")

# Save with high DPI for sharp rasters - multiple formats for best quality
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/CattleDensity_", country_code, ".png"), 
       plot = p_cattle, width = 12, height = 8, dpi = 600, type = "cairo-png")
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/CattleDensity_", country_code, ".tiff"), 
       plot = p_cattle, width = 12, height = 8, dpi = 600, compression = "none")
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/CattleDensity_", country_code, ".pdf"), 
       plot = p_cattle, width = 12, height = 8, device = "pdf", dpi = 600, 
       useDingbats = FALSE, compress = FALSE)

# TSETSE PRESENCE
cat("Processing tsetse data...\n")
tsetse <- raster("Data/Covariates/tsenumbspec")
tsetse[tsetse > 1] <- 1
tsetse <- terra::rast(tsetse)  # Convert to SpatRaster
tsetse_crop <- terra::crop(tsetse, country_bounds)
tsetse_mask <- terra::mask(tsetse_crop, country_bounds)

tsetse_df <- terra::as.data.frame(tsetse_mask, xy = TRUE)
colnames(tsetse_df) <- c("longitude", "latitude", "tsetse_presence")
tsetse_df <- tsetse_df[!is.na(tsetse_df$tsetse_presence), ]
tsetse_df$tsetse_presence <- factor(tsetse_df$tsetse_presence, levels = c(0, 1), labels = c("Absent", "Present"))

p_tsetse <- ggplot(tsetse_df, aes(x = longitude, y = latitude, fill = tsetse_presence)) +
  geom_tile() +
  scale_fill_manual(values = c("Absent" = "#FEE5D9", "Present" = "#A50F15"), name = "Tsetse\nPresence", na.value = "transparent") +
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  ) +
  labs(title = paste("Tsetse fly presence -", country_name), x = "Longitude", y = "Latitude")

# Save with high DPI for sharp rasters - multiple formats for best quality
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/TsetsePresence_", country_code, ".png"), 
       plot = p_tsetse, width = 12, height = 8, dpi = 600, type = "cairo-png")
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/TsetsePresence_", country_code, ".tiff"), 
       plot = p_tsetse, width = 12, height = 8, dpi = 600, compression = "none")
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_maps/TsetsePresence_", country_code, ".pdf"), 
       plot = p_tsetse, width = 12, height = 8, device = "pdf", dpi = 600, 
       useDingbats = FALSE, compress = FALSE)

cat("\n=== All covariate visualizations completed for", country_name, "===\n")
cat("Files saved to Analysis_NGA/ folder\n")

