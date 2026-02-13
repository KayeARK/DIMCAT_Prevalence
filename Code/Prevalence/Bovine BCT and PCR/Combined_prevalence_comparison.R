rm(list=ls())

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggridges)
library(viridis)
library(ggtext)
library(sf)
library(geodata)
library(cowplot)  # For combining plots
library(ggrepel)  # For repelling text labels

# Function to load prevalence data for a country
load_country_data <- function(country_code, country_name) {
  cat("Loading data for", country_name, "...\n")
  
  # Initialize arrays
  i <- 1
  data <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/Cattle_at_risk/Projections_model_", i, ".csv"))
  locations <- data[,1]
  n_locations <- length(locations)
  n_models <- 1000
  
  prevalence <- array(NA, dim = c(n_locations, n_models))
  n_valid_models <- 0
  
  # Load all model files
  for(i in 1:n_models) {
    file_path <- paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_", country_code, "/Cattle_at_risk/Projections_model_", i, ".csv")
    if(file.exists(file_path)) {
      data <- read.csv(file_path)
      prevalence[, i] <- data[, 2]  # Assuming prevalence is in second column
      n_valid_models <- n_valid_models + 1
    }
  }
  
  cat("Loaded", n_valid_models, "model files for", country_name, "\n")
  
  # Convert to tibble
  prevalence_tibble <- as_tibble(prevalence[, 1:n_valid_models])
  colnames(prevalence_tibble) <- paste0("Sim_", 1:n_valid_models)
  prevalence_tibble <- prevalence_tibble %>%
    mutate(Location = locations,
           Country = country_name) %>%
    pivot_longer(cols = starts_with("Sim_"), names_to = "Simulation", values_to = "Prevalence")
  
  return(prevalence_tibble)
}

# Load data for both countries
nigeria_data <- load_country_data("NGA", "Nigeria")
ethiopia_data <- load_country_data("ETH", "Ethiopia")

# Create small reference maps showing regional boundaries
create_reference_map <- function(country_code, country_name, prevalence_data, admin_level = 1) {
  tryCatch({
    # Load country administrative boundaries
    country_admin <- gadm(country = country_code, level = admin_level, path = tempdir())
    country_sf <- st_as_sf(country_admin)
    
    # Apply the same name standardization as used in the main data
    country_sf <- country_sf %>%
      mutate(
        standardized_name = case_when(
          get(paste0("NAME_", admin_level)) == "Benshangul-Gumaz" ~ "Benishangul-Gumuz",
          get(paste0("NAME_", admin_level)) == "Gambela Peoples" ~ "Gambela",
          get(paste0("NAME_", admin_level)) == "Harari People" ~ "Harari",
          get(paste0("NAME_", admin_level)) == "Southern Nations, Nationalities" ~ "SNNPR",
          TRUE ~ get(paste0("NAME_", admin_level))
        )
      )
    
    # Extract prevalence data for this country
    country_prevalence <- prevalence_data %>%
      filter(Country == country_name) %>%
      group_by(Location) %>%
      summarise(mean_prev = mean(Prevalence, na.rm = TRUE), .groups = 'drop')
    
    # Join prevalence data with spatial data
    country_sf <- country_sf %>%
      left_join(country_prevalence, by = c("standardized_name" = "Location"))
    
    # Create simple reference map with repelled labels
    # Extract centroids for text positioning
    country_centroids <- country_sf %>%
      mutate(
        lon = st_coordinates(st_centroid(.))[,1],
        lat = st_coordinates(st_centroid(.))[,2]
      ) %>%
      st_drop_geometry() %>%
      # Manual adjustment for Rivers label
      mutate(
        lon = ifelse(standardized_name == "Rivers", lon - 0.1, lon),  # Move left
        lat = ifelse(standardized_name == "Rivers", lat + 0.2, lat),  # Move up
        lon = ifelse(standardized_name == "Oromia", lon + 1, lon)  # Move left
      )
    
    ref_map <- ggplot(country_sf) +
      geom_sf(aes(fill = mean_prev), color = "white", size = 0.3) +
      scale_fill_viridis_c(name = "Prevalence", option = "viridis", trans = "sqrt", 
                          na.value = "lightgray", limits = c(0, 1)) +
      geom_text_repel(data = country_centroids,
                      aes(x = lon, y = lat, label = standardized_name), 
                      size = 1.6, 
                      color = "black",      # Black text with white outline for visibility
                      fontface = "bold",    # Bold for better readability
                      bg.color = "white",   # White background behind text
                      bg.r = 0.1,          # Radius of background
                      box.padding = 0.05,
                      point.padding = 0.05,
                      force = 0,
                      force_pull = 20,      # Strong pull to keep labels near original position
                      max.overlaps = Inf,
                      segment.color = NA,  # Hide connecting lines
                      max.iter = 3000,     # Fewer iterations for minimal movement
                      seed = 42) +         # Reproducible positioning
      theme_void() +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5, "pt"),
        legend.position = "none"  # Hide legend since main plot has it
      ) +
      labs(title = "")
    
    return(ref_map)
  }, error = function(e) {
    cat("Could not create reference map for", country_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Create reference maps will be created after data processing

zero_prev_region <- tibble(
  Location = "",
  Country = "Nigeria",
  Simulation = paste0("Sim_", 1:200),
  Prevalence = -1
)

one_prev_region <- tibble(
  Location = "",
  Country = "Ethiopia",
  Simulation = paste0("Sim_", 1:200),
  Prevalence = 2
)

# Combine the data
combined_data <- bind_rows(nigeria_data, ethiopia_data, zero_prev_region, one_prev_region)

combined_data <- combined_data %>%
  mutate(Location = case_when(
    Location == "Benshangul-Gumaz" ~ "Benishangul-Gumuz",
    Location == "Gambela Peoples" ~ "Gambela",
    Location == "Harari People" ~ "Harari",
    Location == "Southern Nations, Nationalities" ~ "SNNPR",
    TRUE ~ Location
  ))

# Calculate mean prevalence for ordering locations within each country
combined_data <- combined_data %>%
  group_by(Country, Location) %>%
  mutate(mean_prev = mean(Prevalence, na.rm = TRUE)) %>%
  ungroup() %>%
  # Order locations by mean prevalence within each country
  group_by(Country) %>%
  mutate(Location = fct_reorder(Location, mean_prev)) %>%
  ungroup()

# Create unique location identifiers to avoid duplicates
combined_data <- combined_data %>%
  mutate(Location_Unique = paste(Location, Country, sep = " (") %>% paste0(")"))

# Create reference maps with prevalence-based coloring
cat("Creating reference maps with prevalence coloring...\n")
nigeria_ref_map <- create_reference_map("NGA", "Nigeria", combined_data, 1)  # State level
ethiopia_ref_map <- create_reference_map("ETH", "Ethiopia", combined_data, 1)  # Regional level

# Get ordered locations by country
nigeria_locations <- combined_data %>% 
  filter(Country == "Nigeria") %>% 
  arrange(mean_prev) %>% 
  pull(Location_Unique) %>% 
  unique()

ethiopia_locations <- combined_data %>% 
  filter(Country == "Ethiopia") %>% 
  arrange(mean_prev) %>% 
  pull(Location_Unique) %>% 
  unique()

# Add gap by creating spacer labels
gap_size <- 2  # Number of empty rows between countries
spacer_labels <- paste("SPACER", 1:gap_size)

# Create the full y-axis order: Ethiopia (bottom), gap, Nigeria (top)
y_axis_order <- c(ethiopia_locations, spacer_labels, nigeria_locations)

# Update the combined data with the new y-axis ordering
combined_data <- combined_data %>%
  mutate(Location_Display = factor(Location_Unique, levels = y_axis_order))

# Create stacked ridge plot with gap
combined_ridge_plot <- ggplot(combined_data, aes(x = Prevalence, y = Location_Display, fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_viridis(name = "Prevalence", option = "viridis", trans = "sqrt") +
  scale_y_discrete(
    breaks = y_axis_order[!grepl("SPACER", y_axis_order)],  # Remove spacer from axis labels
    labels = function(x) gsub(" \\([^)]+\\)$", "", x)  # Remove country suffix from display labels
  ) +
  # Add country labels as annotations
  annotate("text", x = 0.625, y = length(ethiopia_locations)/2 + 0.5, 
           label = "ETHIOPIA", angle = 90, size = 4, fontface = "bold", 
           hjust = 0.5, color = "gray30") +
  annotate("text", x = 0.625, y = length(ethiopia_locations) + gap_size + length(nigeria_locations)/2 + 0.5, 
           label = "NIGERIA", angle = 90, size = 4, fontface = "bold", 
           hjust = 0.5, color = "gray30") +
  labs(title = 'AAT Prevalence Comparison: Nigeria vs Ethiopia',
       subtitle = 'Distribution of prevalence estimates across administrative regions',
       x = 'Prevalence',
       y = '') +
  xlim(0, 1) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    #panel.grid.major.y = element_blank(),  # Remove horizontal grid lines for cleaner look
    panel.grid.minor.y = element_blank()
  )

# Save the combined plot with reference maps
if(!is.null(nigeria_ref_map) && !is.null(ethiopia_ref_map)) {
  # Calculate positions to align with country labels
  # Nigeria label is at y = length(ethiopia_locations) + gap_size + length(nigeria_locations)/2 + 0.5
  # Ethiopia label is at y = length(ethiopia_locations)/2 + 0.5
  
  # Convert y-axis positions to plot coordinates (0-1 scale)
  total_regions <- length(ethiopia_locations) + gap_size + length(nigeria_locations)
  
  # Position for Nigeria inset (aligned with Nigeria label)
  nigeria_y_pos <- (length(ethiopia_locations) + gap_size + length(nigeria_locations)/2 + 0.5) / total_regions
  nigeria_y_centered <- nigeria_y_pos - 0.2  # Center the map around the label position
  
  # Position for Ethiopia inset (aligned with Ethiopia label)  
  ethiopia_y_pos <- (length(ethiopia_locations)/2 + 0.5) / total_regions
  ethiopia_y_centered <- ethiopia_y_pos - 0.1  # Center the map around the label position
  
  # Create final combined plot with reference maps as insets
  final_plot <- ggdraw(combined_ridge_plot) +
    # Add Nigeria reference map (aligned with Nigeria label)
    draw_plot(nigeria_ref_map, x = 0.67, y = nigeria_y_centered, width = 0.35, height = 0.35) +
    # Add Ethiopia reference map (aligned with Ethiopia label)
    draw_plot(ethiopia_ref_map, x = 0.67, y = ethiopia_y_centered, width = 0.35, height = 0.35)

  ggsave("Code/Prevalence/Bovine BCT and PCR/Combined_prevalence_comparison.pdf", 
         plot = final_plot, width = 12, height = 10, dpi = 300)
} else {
  # Fallback: save just the ridge plot if reference maps failed
  ggsave("Code/Prevalence/Bovine BCT and PCR/Combined_prevalence_comparison.pdf", 
         plot = combined_ridge_plot, width = 12, height = 10, dpi = 300)
}
