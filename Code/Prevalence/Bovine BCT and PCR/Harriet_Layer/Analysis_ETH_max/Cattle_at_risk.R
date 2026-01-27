rm(list=ls())

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggridges)
library(viridis)
library(ggtext)

i<-1
data <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk/Projections_model_",i,".csv"))

locations <-  data[,1]

#set up an empty array based on actual data dimensions
n_locations <- length(locations)
n_models <- 200  # Adjust based on actual number of models

cattle_at_risk <- array(NA, dim = c(n_locations, n_models))
prevalence <- array(NA, dim = c(n_locations, n_models))

# Load data from all model projections
for (i in 1:n_models){
  
  file_path <- paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk/Projections_model_",i,".csv")
  
  # Check if file exists before trying to read
  if(file.exists(file_path)) {
    data <- read.csv(file_path)
    cattle_at_risk[,i] <- data[,3]  # Assuming column 3 is cattle at risk
    prevalence[,i] <- data[,2]      # Assuming column 2 is prevalence
  } else {
    cat("Warning: File", file_path, "not found\n")
  }
}

# Remove columns with all NAs (missing files)
valid_cols <- !apply(is.na(prevalence), 2, all)
prevalence <- prevalence[, valid_cols]
cattle_at_risk <- cattle_at_risk[, valid_cols]
n_valid_models <- sum(valid_cols)

cat("Loaded data from", n_valid_models, "model files\n")
cat("Number of locations:", n_locations, "\n")

#change prevalence to a tibble for visualization
prevalence_tibble <- as_tibble(prevalence)
colnames(prevalence_tibble) <- paste0("Sim_", 1:n_valid_models)
prevalence_tibble <- prevalence_tibble %>%
  mutate(Location = locations) %>%
  pivot_longer(cols = starts_with("Sim_"), names_to = "Simulation", values_to = "Prevalence")

# Calculate mean prevalence for ordering locations
prevalence_tibble <- prevalence_tibble %>%
  group_by(Location) %>%
  mutate(mean_prev = mean(Prevalence, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Location = fct_reorder(Location, mean_prev))

# Create ridge plot for prevalence across Ethiopian regions/zones
ridge_plot <- ggplot(prevalence_tibble, aes(x = Prevalence, y = Location, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Prevalence", option = "C") +
  labs(title = 'Prevalence across regions in Ethiopia') +
  xlim(0, 1) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 14),
    plot.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11.5)
  )

# Save the ridge plot
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Prevalence_plot.pdf", 
       plot = ridge_plot, width=8, height=6)

# Create individual histograms for cattle at risk by location
for (i in 1:n_locations){
  df <- data.frame(cattle_at_risk = cattle_at_risk[i,])
  # Remove NA values
  df <- df[!is.na(df$cattle_at_risk), , drop=FALSE]
  
  if(nrow(df) > 0) {
    p <- ggplot(df, aes(x=cattle_at_risk)) + 
      geom_histogram(bins=20, fill="#D55E00", color="black", alpha=0.7) + 
      ggtitle(paste0("Cattle at risk - ", locations[i])) + 
      xlab("Number of cattle at risk") + 
      ylab("Frequency") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_plots/Histogram_cattle_",
                  gsub(" ", "_", locations[i]), ".png"), 
           plot=p, width=8, height=6)
  }
}

# Create individual histograms for prevalence by location  
for (i in 1:n_locations){
  df <- data.frame(prevalence = prevalence[i,])
  # Remove NA values
  df <- df[!is.na(df$prevalence), , drop=FALSE]
  
  if(nrow(df) > 0) {
    p <- ggplot(df, aes(x=prevalence)) + 
      geom_histogram(bins=20, fill="#009E73", color="black", alpha=0.7) + 
      ggtitle(paste0("Prevalence distribution - ", locations[i])) + 
      xlab("Prevalence") + 
      ylab("Frequency") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk_plots/Histogram_prevalence_",
                  gsub(" ", "_", locations[i]), ".png"), 
           plot=p, width=8, height=6)
  }
}

# Summary statistics
cat("\n=== ETHIOPIA CATTLE AT RISK ANALYSIS ===\n")
cat("Summary statistics across all locations and models:\n")

# Calculate summary stats
prevalence_summary <- prevalence_tibble %>%
  group_by(Location) %>%
  summarise(
    mean_prevalence = mean(Prevalence, na.rm = TRUE),
    median_prevalence = median(Prevalence, na.rm = TRUE),
    q25_prevalence = quantile(Prevalence, 0.25, na.rm = TRUE),
    q75_prevalence = quantile(Prevalence, 0.75, na.rm = TRUE),
    min_prevalence = min(Prevalence, na.rm = TRUE),
    max_prevalence = max(Prevalence, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_prevalence))

print(prevalence_summary)

# Calculate cattle at risk summary
cattle_summary <- data.frame(
  Location = locations,
  mean_cattle_at_risk = apply(cattle_at_risk, 1, mean, na.rm = TRUE),
  median_cattle_at_risk = apply(cattle_at_risk, 1, median, na.rm = TRUE),
  q25_cattle_at_risk = apply(cattle_at_risk, 1, function(x) quantile(x, 0.25, na.rm = TRUE)),
  q75_cattle_at_risk = apply(cattle_at_risk, 1, function(x) quantile(x, 0.75, na.rm = TRUE)),
  min_cattle_at_risk = apply(cattle_at_risk, 1, min, na.rm = TRUE),
  max_cattle_at_risk = apply(cattle_at_risk, 1, max, na.rm = TRUE)
) %>%
  arrange(desc(mean_cattle_at_risk))

cat("\nCattle at risk summary:\n")
print(cattle_summary)

# Save summary tables
write.csv(prevalence_summary, "Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/ethiopia_prevalence_summary.csv", row.names = FALSE)
write.csv(cattle_summary, "Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/ethiopia_cattle_at_risk_summary.csv", row.names = FALSE)

cat("\nAnalysis complete!\n")
cat("Generated files:\n")
cat("- Prevalence_plot.pdf (ridge plot)\n") 
cat("- Individual histograms in Cattle_at_risk_plots/\n")
cat("- ethiopia_prevalence_summary.csv\n")
cat("- ethiopia_cattle_at_risk_summary.csv\n")