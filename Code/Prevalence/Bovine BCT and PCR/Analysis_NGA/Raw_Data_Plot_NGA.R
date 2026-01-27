library(ggplot2)
library(dplyr)
library(readr)
library(sf)
library(geodata)
library(viridis)
library(gridExtra)
library(grid)
library(gtable)
library(readxl)
library(ggspatial)
library(cowplot)

# Get Nigerian administrative boundaries
cat("Loading Nigerian administrative boundaries...\n")
nigeria_states <- gadm(country = "NGA", level = 1, path = tempdir()) 
nigeria_states_sf <- st_as_sf(nigeria_states)

# Get neighboring countries for geographic context
cat("Loading neighboring countries for geographic context...\n")
neighboring_countries <- c("BEN", "NER", "TCD", "CMR")  # Benin, Niger, Chad, Cameroon
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

# Combine all neighboring countries into one sf object
if(length(neighbor_sf_list) > 0) {
  neighboring_countries_full <- do.call(rbind, neighbor_sf_list)
  
  # Validate and repair geometries to avoid intersection errors
  cat("Validating and repairing geometries...\n")
  neighboring_countries_full <- st_make_valid(neighboring_countries_full)
  nigeria_states_valid <- st_make_valid(nigeria_states_sf)
  
  # Create a buffer around Nigeria to crop neighboring countries to border regions only
  nigeria_buffer <- st_buffer(st_union(nigeria_states_valid), dist = 1.5)  # ~150km buffer
  nigeria_buffer <- st_make_valid(nigeria_buffer)
  
  # Use st_crop instead of st_intersection for more robust clipping
  tryCatch({
    neighboring_countries <- st_crop(neighboring_countries_full, nigeria_buffer)
    cat("Successfully loaded and cropped", nrow(neighboring_countries), "neighboring country border regions\n")
  }, error = function(e) {
    cat("st_crop failed, trying alternative approach:", e$message, "\n")
    # Fallback: use a simpler bbox-based crop
    nigeria_bbox <- st_bbox(nigeria_buffer)
    neighboring_countries <<- st_crop(neighboring_countries_full, nigeria_bbox)
    cat("Successfully loaded and cropped using bbox approach\n")
  })
  
  neighboring_borders <- neighboring_countries  # Keep for backwards compatibility
} else {
  neighboring_countries <- NULL
  neighboring_borders <- NULL
}

cat("Geographic boundaries prepared successfully!\n")

# Load bovine data from Excel files
cat("Loading AAT bovine data...\n")
bovine_bct_raw <- readxl::read_excel("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
bovine_pcr_raw <- readxl::read_excel("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")

cat("BCT dataset columns:", ncol(bovine_bct_raw), "\n")
cat("PCR dataset columns:", ncol(bovine_pcr_raw), "\n")

# Check column names
cat("BCT columns:", paste(names(bovine_bct_raw), collapse = ", "), "\n")
cat("PCR columns:", paste(names(bovine_pcr_raw), collapse = ", "), "\n")

# Standardize the datasets to have common columns
# Handle BCT dataset
bovine_bct_clean <- bovine_bct_raw %>%
  dplyr::mutate(Test_Type = "BCT/HCT")

# Add prevalence rate column if TPR exists, otherwise use NA
if("TPR" %in% names(bovine_bct_clean)) {
  bovine_bct_clean$Prevalence_Rate <- bovine_bct_clean$TPR
} else {
  bovine_bct_clean$Prevalence_Rate <- NA
}

# Select standard columns
bovine_bct_clean <- bovine_bct_clean %>%
  dplyr::select(
    Longitude, Latitude, Number_of_animal_tested, 
    Number_of_infections, Prevalence_Rate, Test_Type
  )

# Handle PCR dataset
bovine_pcr_clean <- bovine_pcr_raw %>%
  dplyr::mutate(Test_Type = "PCR")

# Add prevalence rate column if T_ATPR exists, otherwise use NA  
if("T_ATPR" %in% names(bovine_pcr_clean)) {
  bovine_pcr_clean$Prevalence_Rate <- bovine_pcr_clean$T_ATPR
} else {
  bovine_pcr_clean$Prevalence_Rate <- NA
}

# Add infections column if it doesn't exist
if(!"Number_of_infections" %in% names(bovine_pcr_clean)) {
  bovine_pcr_clean$Number_of_infections <- NA
}

# Select standard columns
bovine_pcr_clean <- bovine_pcr_clean %>%
  dplyr::select(
    Longitude, Latitude, Number_of_animal_tested, 
    Number_of_infections, Prevalence_Rate, Test_Type
  )

# Combine both datasets using bind_rows for safer joining
bovine_data_raw <- dplyr::bind_rows(bovine_bct_clean, bovine_pcr_clean)

cat("Combined dataset has", nrow(bovine_data_raw), "rows\n")
cat("Test type distribution:\n")
print(table(bovine_data_raw$Test_Type, useNA = "ifany"))

# Process bovine data for mapping
cat("Processing bovine data for Nigeria...\n")
bovine_clean <- bovine_data_raw %>%
  # Remove rows with missing coordinates
  dplyr::filter(!is.na(Longitude), !is.na(Latitude)) %>%
  # Remove rows with missing or invalid sample sizes
  dplyr::filter(!is.na(Number_of_animal_tested), Number_of_animal_tested > 0) %>%
  # Calculate prevalence
  dplyr::mutate(
    # Handle missing infections - calculate from prevalence rate if available
    Number_of_infections = case_when(
      !is.na(Number_of_infections) ~ Number_of_infections,
      !is.na(Prevalence_Rate) ~ round(Number_of_animal_tested * Prevalence_Rate / 100),
      TRUE ~ 0
    ),
    # Handle missing prevalence rate
    Prevalence_Rate = case_when(
      !is.na(Prevalence_Rate) ~ Prevalence_Rate,
      !is.na(Number_of_infections) & Number_of_animal_tested > 0 ~ 
        (Number_of_infections / Number_of_animal_tested) * 100,
      TRUE ~ 0
    ),
    # Ensure infections don't exceed sample size
    Number_of_infections = pmin(Number_of_infections, Number_of_animal_tested, na.rm = TRUE),
    # Calculate prevalence as proportion
    Prevalence = ifelse(Number_of_animal_tested > 0, 
                       Number_of_infections / Number_of_animal_tested, 0),
    # Create sample size categories for point sizing (wider ranges for smaller dots)
    Sample_Size_Category = case_when(
      Number_of_animal_tested <= 25 ~ "1-25",
      Number_of_animal_tested <= 50 ~ "26-50", 
      Number_of_animal_tested <= 100 ~ "51-100",
      Number_of_animal_tested <= 200 ~ "101-200",
      TRUE ~ ">200"
    ),
    Sample_Size_Category = factor(Sample_Size_Category, 
                                 levels = c("1-25", "26-50", "51-100", "101-200", ">200"))
  )

# Convert to spatial points
bovine_sf <- st_as_sf(bovine_clean, 
                      coords = c("Longitude", "Latitude"), 
                      crs = 4326)

# Filter for points within Nigeria using spatial intersection
bovine_nigeria <- st_filter(bovine_sf, nigeria_states_sf)

# Add back coordinate columns for plotting
coords <- st_coordinates(bovine_nigeria)
bovine_nigeria$lon <- coords[, 1]
bovine_nigeria$lat <- coords[, 2]

cat("Processed", nrow(bovine_nigeria), "observations in Nigeria\n")
cat("Test types available:", unique(bovine_nigeria$Test_Type), "\n")
cat("Sample size distribution:\n")
print(table(bovine_nigeria$Sample_Size_Category, bovine_nigeria$Test_Type))

# Calculate maximum prevalence for scale capping
max_prevalence <- max(bovine_nigeria$Prevalence, na.rm = TRUE)
cat("Maximum prevalence:", round(max_prevalence * 100, 1), "%\n")

# Create the map
cat("Creating raw data visualization...\n")

# Define point sizes for sample size categories (smaller overall)
size_values <- c("1-25" = 1.0, "26-50" = 1.5, "51-100" = 2.5, 
                 "101-200" = 3.5, ">200" = 5.0)

# Define shapes for test types (circle for BCT/HCT, triangle for PCR)
shape_values <- c("BCT/HCT" = 16, "PCR" = 17)

# Get Nigeria bounds for plot limits
nigeria_bbox <- st_bbox(nigeria_states_sf)
lon_range <- nigeria_bbox[c("xmin", "xmax")]
lat_range <- nigeria_bbox[c("ymin", "ymax")]

# Add buffer to bounds (matching LGA_choropleth approach)
lon_buffer <- diff(lon_range) * 0.1  # Moderate buffer to match choropleth plots
lat_buffer <- diff(lat_range) * 0.1

p <- ggplot() +
  # Add neighboring countries (matching LGA_choropleth style)
  {if(!is.null(neighboring_countries) && nrow(neighboring_countries) > 0) {
    geom_sf(data = neighboring_countries, fill = "grey95", color = "grey80", lwd = 0.3)
  }} +
  
  # Add Nigeria state boundaries with single clean border
  geom_sf(data = nigeria_states_sf, fill = "white", color = "gray50", size = 0.3) +
  
  # Add bovine data points
  geom_point(data = as.data.frame(bovine_nigeria), 
             aes(x = lon, y = lat, 
                 color = Prevalence,
                 size = Sample_Size_Category,
                 shape = Test_Type),
             alpha = 0.8, stroke = 0.5) +
  
  # Use viridis green for prevalence color scale (capped at maximum)
  scale_color_viridis_c(name = "Prevalence", 
                        option = "viridis",
                        limits = c(0, max_prevalence),
                        labels = scales::percent_format(accuracy = 1),
                        guide = guide_colorbar(title.position = "top",
                                             title.hjust = 0.5,
                                             barwidth = unit(0.7, "cm"),
                                             barheight = unit(4, "cm"),
                                             order = 1)) +
  
  # Set point sizes
  scale_size_manual(name = "Sample size", 
                    values = size_values,
                    guide = guide_legend(title.position = "top",
                                       title.hjust = 0.5,
                                       override.aes = list(color = "black", alpha = 1),
                                       order = 2)) +
  
  # Set point shapes for test types  
  scale_shape_manual(name = "Test type",
                     values = shape_values,
                     guide = guide_legend(title.position = "top",
                                        title.hjust = 0.5,
                                        override.aes = list(size = 4, alpha = 1),
                                        order = 3)) +
  
  # Set coordinate system and limits
  coord_sf(xlim = lon_range + c(-lon_buffer, lon_buffer),
           ylim = lat_range + c(-lat_buffer, lat_buffer),
           expand = FALSE) +
  
  # Add scale bar and north arrow (matching CI plots)
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8, 
                  bar_cols = c("black", "white"), line_width = 1) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                        pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
                        style = north_arrow_fancy_orienteering(text_size = 8)) +
  
  # Theme and labels
  theme_void() +
  labs(
    title = "Raw AAT survey data - Nigeria",
    subtitle = paste0(nrow(bovine_nigeria), " observations | ", 
                     "BCT/HCT: ", sum(bovine_nigeria$Test_Type == "BCT/HCT"),
                     " | PCR: ", sum(bovine_nigeria$Test_Type == "PCR")),
    caption = "Size indicates sample size, colour indicates prevalence rate, shape indicates test type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 10, color = "grey60"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    legend.margin = margin(l = 20)
  )

print(p)

# Save the map plot
output_file <- "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Raw_AAT_Data_Nigeria.png"
ggsave(output_file, plot = p, width = 14, height = 10, dpi = 300, bg = "white")
cat("Map plot saved to:", output_file, "\n")

# Create prevalence histogram
cat("Creating prevalence histogram...\n")

hist_data <- bovine_nigeria %>%
  st_drop_geometry() %>%
  dplyr::filter(!is.na(Prevalence))

h <- ggplot(hist_data, aes(x = Prevalence * 100)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white", size = 0.3) +
  facet_wrap(~Test_Type, scales = "free_y", ncol = 1) +
  scale_x_continuous(name = "Prevalence (%)", 
                     limits = c(0, max_prevalence * 100),
                     expand = c(0.02, 0)) +
  scale_y_continuous(name = "Number of observations", 
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(
    title = "Distribution of AAT Prevalence - Nigeria"
  )

print(h)

# Save histogram
hist_file <- "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/AAT_Prevalence_Histogram_Nigeria.png"
ggsave(hist_file, plot = h, width = 10, height = 8, dpi = 300, bg = "white")
cat("Histogram saved to:", hist_file, "\n")

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total observations in Nigeria:", nrow(bovine_nigeria), "\n")
cat("By test type:\n")
test_summary <- bovine_nigeria %>%
  st_drop_geometry() %>%
  group_by(Test_Type) %>%
  summarise(
    n_obs = n(),
    mean_prevalence = mean(Prevalence, na.rm = TRUE),
    median_sample_size = median(Number_of_animal_tested, na.rm = TRUE),
    total_animals_tested = sum(Number_of_animal_tested, na.rm = TRUE),
    total_positive = sum(Number_of_infections, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    overall_prevalence = total_positive / total_animals_tested
  )

print(test_summary)