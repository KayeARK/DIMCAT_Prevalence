rm(list=ls())

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggridges)
library(viridis)
library(ggtext)

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
  scale_fill_viridis(name = "Prevalence", option = "viridis") +
  scale_y_discrete(
    breaks = y_axis_order[!grepl("SPACER", y_axis_order)],  # Remove spacer from axis labels
    labels = function(x) gsub(" \\([^)]+\\)$", "", x)  # Remove country suffix from display labels
  ) +
  # Add country labels as annotations
  annotate("text", x = 0.95, y = length(ethiopia_locations)/2 + 0.5, 
           label = "ETHIOPIA", angle = 90, size = 4, fontface = "bold", 
           hjust = 0.5, color = "gray30") +
  annotate("text", x = 0.95, y = length(ethiopia_locations) + gap_size + length(nigeria_locations)/2 + 0.5, 
           label = "NIGERIA", angle = 90, size = 4, fontface = "bold", 
           hjust = 0.5, color = "gray30") +
  labs(title = 'AAT Prevalence Comparison: Nigeria vs Ethiopia',
       subtitle = 'Distribution of prevalence estimates across administrative regions',
       x = 'Prevalence',
       y = '') +
  xlim(0, 1) +
  theme_minimal() +
  theme(
    legend.position = "right",
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

# Save the combined plot
ggsave("Code/Prevalence/Bovine BCT and PCR/Combined_prevalence_comparison.pdf", 
       plot = combined_ridge_plot, width = 12, height = 10, dpi = 300)
