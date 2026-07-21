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
      # Identify overlapping locations in southern Nigeria (lat < 7.5°N approximately)
      mutate(
        is_southern_nigeria = country_name == "Nigeria" & lat < 7.5,
        # Define which labels are likely to overlap - expanded list
        needs_external_label = standardized_name %in% c("Rivers", "Cross River", "Akwa Ibom", "Bayelsa", 
                                                        "Kwara", "Osun", "Anambra", "Ebonyi", "Harari", 
                                                        "Ondo", "Abia"),
        # Manual adjustment for overlapping areas - move labels outside and show lines
        lon_adjusted = case_when(
          # Original southern Nigeria overlaps
          standardized_name == "Rivers" ~ lon - 0.3,  # Move far left
          standardized_name == "Cross River" ~ lon + 1.2,  # Move right
          standardized_name == "Akwa Ibom" ~ lon + 0.8,  # Move slightly right
          standardized_name == "Bayelsa" ~ lon - 1.0,  # Move left
          # Additional Nigeria overlaps
          standardized_name == "Kwara" ~ lon - 1.2,  # Move left
          standardized_name == "Osun" ~ lon - 1.5,  # Move far left
          standardized_name == "Anambra" ~ lon - 1.3,  # Move right
          standardized_name == "Ebonyi" ~ lon + 1.3,  # Move far right
          standardized_name == "Ondo" ~ lon - 1.1,  # Move left
          standardized_name == "Abia" ~ lon + 0.9,  # Move right
          # Ethiopia overlap
          standardized_name == "Harari" ~ lon - 0,  # Move left
          standardized_name == "Oromia" ~ lon + 1,  # Move Oromia (Ethiopia) left
          TRUE ~ lon
        ),
        lat_adjusted = case_when(
          # Original southern Nigeria overlaps
          standardized_name == "Rivers" ~ lat - 0.9,  # Move up
          standardized_name == "Cross River" ~ lat + 0.5,  # Move up
          standardized_name == "Akwa Ibom" ~ lat - 0.7,  # Move down
          standardized_name == "Bayelsa" ~ lat - 0.5,  # Move down
          # Additional Nigeria overlaps
          standardized_name == "Kwara" ~ lat + 0.8,  # Move up
          standardized_name == "Osun" ~ lat ,  # Move down
          standardized_name == "Anambra" ~ lat ,  # Move up
          standardized_name == "Ebonyi" ~ lat + 0.6,  # Move up
          standardized_name == "Ondo" ~ lat - 0.8,  # Move down
          standardized_name == "Abia" ~ lat ,  # Move down
          # Ethiopia overlap
          standardized_name == "Harari" ~ lat - 0.5,  # Move up
          TRUE ~ lat
        )
      )
    
    ref_map <- ggplot(country_sf) +
      geom_sf(aes(fill = mean_prev), color = "white", size = 0.3) +
      scale_fill_viridis_c(name = "Prevalence", option = "viridis", trans = "sqrt", 
                          na.value = "lightgray", limits = c(0, 1)) +
      # Add connecting lines for externally positioned labels
      geom_segment(data = country_centroids %>% filter(needs_external_label),
                   aes(x = lon, y = lat, xend = lon_adjusted, yend = lat_adjusted),
                   color = "gray50", size = 0.3, alpha = 0.6) +
      # Add small points at original centroids for labels that were moved outside
      geom_point(data = country_centroids %>% filter(needs_external_label),
                 aes(x = lon, y = lat), size = 1, color = "gray30", alpha = 0.7) +
      # Add text labels at adjusted positions
      geom_text_repel(data = country_centroids,
                      aes(x = lon_adjusted, y = lat_adjusted, label = standardized_name), 
                      size = 2.6, 
                      color = "black",      # Black text with white outline for visibility
                      fontface = "bold",    # Bold for better readability
                      bg.color = "white",   # White background behind text
                      bg.r = 0.1,          # Radius of background
                      box.padding = 0.05,
                      point.padding = 0.05,
                      force = 0,           # Slight force to avoid overlaps
                      force_pull = 20,      # Reduced pull to allow more movement
                      max.overlaps = Inf,
                      segment.color = NA,  # Hide ggrepel's own segments since we draw our own
                      max.iter = 1000,     # Allow more iterations for better positioning
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







# ============================================================
# FIGURE 4 SOURCE DATA
#
# Source data for:
#   1. Administrative-area prevalence ridge distributions
#   2. Nigeria state-level reference choropleth
#   3. Ethiopia region-level reference choropleth
#
# Everything is written to one Excel workbook.
# ============================================================

library(openxlsx)
library(dplyr)
library(tidyr)
library(forcats)

output_dir <- "Code/Prevalence/Bovine BCT and PCR"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Extract genuine model outputs only
# ------------------------------------------------------------

# Exclude the two artificial plotting rows:
# Nigeria prevalence = -1
# Ethiopia prevalence = 2
#
# Also exclude blank administrative-area names.

figure4_prevalence_draws <- combined_data %>%
  mutate(
    Location = as.character(Location),
    Country = as.character(Country),
    Simulation = as.character(Simulation)
  ) %>%
  filter(
    !is.na(Location),
    Location != "",
    !is.na(Prevalence),
    Prevalence >= 0,
    Prevalence <= 1
  ) %>%
  transmute(
    Country = Country,
    Administrative_area = Location,
    Model_draw = Simulation,
    Prevalence = Prevalence
  ) %>%
  arrange(
    Country,
    Administrative_area,
    readr::parse_number(Model_draw)
  )

# ------------------------------------------------------------
# 2. Verify that rows are uniquely identified
# ------------------------------------------------------------

duplicate_draws <- figure4_prevalence_draws %>%
  count(
    Country,
    Administrative_area,
    Model_draw,
    name = "number_of_rows"
  ) %>%
  filter(number_of_rows > 1)

if (nrow(duplicate_draws) > 0) {
  stop(
    paste(
      "Duplicate country-area-model combinations were found.",
      "Check the source projection files before exporting."
    )
  )
}

# ------------------------------------------------------------
# 3. Mean prevalence values used for the choropleths
# ------------------------------------------------------------

figure4_map_values <- figure4_prevalence_draws %>%
  group_by(
    Country,
    Administrative_area
  ) %>%
  summarise(
    Mean_prevalence = mean(Prevalence, na.rm = TRUE),
    Number_of_model_draws = sum(!is.na(Prevalence)),
    .groups = "drop"
  ) %>%
  arrange(
    Country,
    Administrative_area
  )

# ------------------------------------------------------------
# 4. Plot ordering used for the ridge distributions
# ------------------------------------------------------------

# Reconstruct the ordering directly from the genuine data.
# Locations are ordered by mean prevalence within each country.

figure4_plot_order <- figure4_map_values %>%
  group_by(Country) %>%
  arrange(
    Mean_prevalence,
    Administrative_area,
    .by_group = TRUE
  ) %>%
  mutate(
    Order_within_country = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    Country_plot_block = case_when(
      Country == "Ethiopia" ~ 1L,
      Country == "Nigeria" ~ 2L,
      TRUE ~ 3L
    )
  ) %>%
  arrange(
    Country_plot_block,
    Order_within_country
  ) %>%
  mutate(
    Overall_plot_order = row_number()
  ) %>%
  select(
    Country,
    Administrative_area,
    Mean_prevalence,
    Order_within_country,
    Overall_plot_order
  )

# ------------------------------------------------------------
# 5. Summary statistics for reconstruction checks
# ------------------------------------------------------------

summary_parameter <- function(x) {
  tibble(
    Mean = mean(x, na.rm = TRUE),
    Standard_deviation = sd(x, na.rm = TRUE),
    Minimum = min(x, na.rm = TRUE),
    Lower_2.5_percentile = unname(quantile(x, 0.025, na.rm = TRUE)),
    Median = median(x, na.rm = TRUE),
    Upper_97.5_percentile = unname(quantile(x, 0.975, na.rm = TRUE)),
    Maximum = max(x, na.rm = TRUE)
  )
}

figure4_country_summary <- figure4_prevalence_draws %>%
  group_by(Country) %>%
  group_modify(
    ~ summary_parameter(.x$Prevalence)
  ) %>%
  ungroup() %>%
  left_join(
    figure4_map_values %>%
      group_by(Country) %>%
      summarise(
        Number_of_administrative_areas = n(),
        Minimum_draws_per_area = min(Number_of_model_draws),
        Maximum_draws_per_area = max(Number_of_model_draws),
        .groups = "drop"
      ),
    by = "Country"
  ) %>%
  select(
    Country,
    Number_of_administrative_areas,
    Minimum_draws_per_area,
    Maximum_draws_per_area,
    Mean,
    Standard_deviation,
    Minimum,
    Lower_2.5_percentile,
    Median,
    Upper_97.5_percentile,
    Maximum
  )

# ------------------------------------------------------------
# 6. Check Excel row limit
# ------------------------------------------------------------

excel_max_data_rows <- 1048575L

if (nrow(figure4_prevalence_draws) > excel_max_data_rows) {
  stop(
    paste0(
      "The Prevalence_draws table contains ",
      format(nrow(figure4_prevalence_draws), big.mark = ","),
      " rows, which exceeds the Excel worksheet limit of ",
      format(excel_max_data_rows, big.mark = ","),
      " data rows."
    )
  )
}

# ------------------------------------------------------------
# 7. README
# ------------------------------------------------------------

figure4_readme <- data.frame(
  Item = c(
    "Workbook description",
    "Associated figure",
    "Prevalence_draws sheet",
    "Country",
    "Administrative_area",
    "Model_draw",
    "Prevalence",
    "Map_values sheet",
    "Mean_prevalence",
    "Number_of_model_draws",
    "Plot_order sheet",
    "Order_within_country",
    "Overall_plot_order",
    "Country_summary sheet",
    "Administrative levels",
    "Artificial plotting rows",
    "Geographic boundaries",
    "Units",
    "Missing values"
  ),
  Description = c(
    paste(
      "Source data underlying the comparison of modelled AAT",
      "prevalence distributions and administrative-area",
      "reference maps for Nigeria and Ethiopia."
    ),
    "Figure 4.",
    paste(
      "One row per country, administrative area and retained",
      "model draw. These prevalence values are the observations",
      "used to estimate the ridge density for each area."
    ),
    "Country represented in the figure.",
    paste(
      "Administrative-area name used to label the ridge",
      "distribution and join prevalence values to the map."
    ),
    paste(
      "Identifier of the model realisation from which the",
      "prevalence estimate was obtained."
    ),
    paste(
      "Modelled prevalence estimate for the corresponding",
      "administrative area and model draw."
    ),
    paste(
      "One row per mapped administrative area. Contains the",
      "mean prevalence value used to colour the reference maps."
    ),
    paste(
      "Arithmetic mean of all retained prevalence model draws",
      "for the administrative area."
    ),
    paste(
      "Number of retained model prevalence estimates contributing",
      "to the administrative-area mean."
    ),
    paste(
      "Administrative-area ordering used to construct the vertical",
      "axis of the ridge plot."
    ),
    paste(
      "Rank of the administrative area by mean prevalence within",
      "its country, from lowest to highest."
    ),
    paste(
      "Position in the complete plot ordering, with Ethiopia",
      "followed by Nigeria. The visual gap between countries is",
      "a plotting feature and is not represented as data."
    ),
    paste(
      "Country-level summary statistics supplied to facilitate",
      "verification of the source-data reconstruction."
    ),
    paste(
      "Nigeria is represented at GADM administrative level 1",
      "(states). Ethiopia is represented at GADM administrative",
      "level 1 (regions)."
    ),
    paste(
      "Artificial prevalence values of -1 and 2 used in the",
      "plotting script have been excluded because they are not",
      "model outputs and fall outside the plotted prevalence range."
    ),
    paste(
      "Administrative boundaries are obtained separately from",
      "GADM using the geodata R package and are not duplicated",
      "in this workbook."
    ),
    "All prevalence values are proportions ranging from 0 to 1.",
    "Blank cells represent unavailable values."
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 8. Create one Excel workbook
# ------------------------------------------------------------

wb <- createWorkbook(
  creator = "AR Kaye",
  title = "Figure 4 source data",
  subject = paste(
    "Nigeria and Ethiopia administrative-area prevalence",
    "distributions and reference-map values"
  ),
  category = "Research data"
)

addWorksheet(wb, "README")
addWorksheet(wb, "Prevalence_draws")
addWorksheet(wb, "Map_values")
addWorksheet(wb, "Plot_order")
addWorksheet(wb, "Country_summary")

writeData(
  wb,
  sheet = "README",
  x = figure4_readme
)

writeData(
  wb,
  sheet = "Prevalence_draws",
  x = figure4_prevalence_draws
)

writeData(
  wb,
  sheet = "Map_values",
  x = figure4_map_values
)

writeData(
  wb,
  sheet = "Plot_order",
  x = figure4_plot_order
)

writeData(
  wb,
  sheet = "Country_summary",
  x = figure4_country_summary
)

# ------------------------------------------------------------
# 9. Workbook styles
# ------------------------------------------------------------

header_style <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#1F4E78",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "bottom",
  borderColour = "#FFFFFF"
)

prevalence_style <- createStyle(
  numFmt = "0.000000"
)

integer_style <- createStyle(
  numFmt = "0"
)

wrap_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

# ------------------------------------------------------------
# 10. Format README
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "README",
  style = header_style,
  rows = 1,
  cols = 1:ncol(figure4_readme),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "README",
  style = wrap_style,
  rows = 2:(nrow(figure4_readme) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "README",
  cols = 1,
  widths = 28
)

setColWidths(
  wb,
  sheet = "README",
  cols = 2,
  widths = 85
)

setRowHeights(
  wb,
  sheet = "README",
  rows = 2:(nrow(figure4_readme) + 1),
  heights = 32
)

freezePane(
  wb,
  sheet = "README",
  firstRow = TRUE
)

# ------------------------------------------------------------
# 11. Format prevalence-draw sheet
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Prevalence_draws",
  style = header_style,
  rows = 1,
  cols = 1:ncol(figure4_prevalence_draws),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Prevalence_draws",
  style = prevalence_style,
  rows = 2:(nrow(figure4_prevalence_draws) + 1),
  cols = 4,
  gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "Prevalence_draws",
  cols = 1,
  widths = 14
)

setColWidths(
  wb,
  sheet = "Prevalence_draws",
  cols = 2,
  widths = 30
)

setColWidths(
  wb,
  sheet = "Prevalence_draws",
  cols = 3,
  widths = 16
)

setColWidths(
  wb,
  sheet = "Prevalence_draws",
  cols = 4,
  widths = 16
)

freezePane(
  wb,
  sheet = "Prevalence_draws",
  firstRow = TRUE
)

addFilter(
  wb,
  sheet = "Prevalence_draws",
  row = 1,
  cols = 1:ncol(figure4_prevalence_draws)
)

# ------------------------------------------------------------
# 12. Format map-values sheet
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Map_values",
  style = header_style,
  rows = 1,
  cols = 1:ncol(figure4_map_values),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Map_values",
  style = prevalence_style,
  rows = 2:(nrow(figure4_map_values) + 1),
  cols = 3,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Map_values",
  style = integer_style,
  rows = 2:(nrow(figure4_map_values) + 1),
  cols = 4,
  gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "Map_values",
  cols = 1:ncol(figure4_map_values),
  widths = "auto"
)

freezePane(
  wb,
  sheet = "Map_values",
  firstRow = TRUE
)

addFilter(
  wb,
  sheet = "Map_values",
  row = 1,
  cols = 1:ncol(figure4_map_values)
)

# ------------------------------------------------------------
# 13. Format plot-order sheet
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Plot_order",
  style = header_style,
  rows = 1,
  cols = 1:ncol(figure4_plot_order),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Plot_order",
  style = prevalence_style,
  rows = 2:(nrow(figure4_plot_order) + 1),
  cols = 3,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Plot_order",
  style = integer_style,
  rows = 2:(nrow(figure4_plot_order) + 1),
  cols = 4:5,
  gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "Plot_order",
  cols = 1:ncol(figure4_plot_order),
  widths = "auto"
)

freezePane(
  wb,
  sheet = "Plot_order",
  firstRow = TRUE
)

addFilter(
  wb,
  sheet = "Plot_order",
  row = 1,
  cols = 1:ncol(figure4_plot_order)
)

# ------------------------------------------------------------
# 14. Format summary sheet
# ------------------------------------------------------------

addStyle(
  wb,
  sheet = "Country_summary",
  style = header_style,
  rows = 1,
  cols = 1:ncol(figure4_country_summary),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Country_summary",
  style = integer_style,
  rows = 2:(nrow(figure4_country_summary) + 1),
  cols = 2:4,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Country_summary",
  style = prevalence_style,
  rows = 2:(nrow(figure4_country_summary) + 1),
  cols = 5:ncol(figure4_country_summary),
  gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "Country_summary",
  cols = 1:ncol(figure4_country_summary),
  widths = "auto"
)

freezePane(
  wb,
  sheet = "Country_summary",
  firstRow = TRUE
)

# ------------------------------------------------------------
# 15. Save workbook
# ------------------------------------------------------------

figure4_output_file <- file.path(
  output_dir,
  "Figure_4_source_data.xlsx"
)

saveWorkbook(
  wb,
  file = figure4_output_file,
  overwrite = TRUE
)

message(
  "Figure 4 source-data workbook saved to: ",
  figure4_output_file
)

message(
  "Prevalence draws exported: ",
  format(nrow(figure4_prevalence_draws), big.mark = ",")
)

message(
  "Administrative areas exported: ",
  format(nrow(figure4_map_values), big.mark = ",")
)