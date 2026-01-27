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
data <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk/Projections_model_",i,".csv"))

locations <-  data[,1]

#set up an empty array of size 34*1000
cattle_at_risk <- array(NA, dim = c(length(locations),1000))
prevalence <- array(NA, dim = c(length(locations),1000))

for (i in 1:1000){



data <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk/Projections_model_",i,".csv"))
cattle_at_risk[,i] <- data[,3]
prevalence[,i] <- data[,2]

}

#change prevalence to a tibble (just like the lincoln weather data)
prevalence_tibble <- as_tibble(prevalence)
colnames(prevalence_tibble) <- paste0("Sim_", 1:1000)
prevalence_tibble <- prevalence_tibble %>%
  mutate(Location = locations) %>%
  pivot_longer(cols = starts_with("Sim_"), names_to = "Simulation", values_to = "Prevalence")

prevalence_tibble <- prevalence_tibble %>%
  group_by(Location) %>%
  mutate(mean_prev = mean(Prevalence, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Location = fct_reorder(Location, mean_prev))

# Plot
ggplot(prevalence_tibble, aes(x = `Prevalence`, y = `Location`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Prevalence", option = "C") +
  labs(title = 'Prevalence across states in Nigeria') +
  xlim(0, 1) +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    )
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Prevalence_plot.pdf", width=8, height=6)

for (i in 1:38){
df <- data.frame(cattle_at_risk=cattle_at_risk[i,])
p <- ggplot(df, aes(x=cattle_at_risk)) + geom_histogram(bins=20, fill="#D55E00", color="black") + ggtitle(paste0("Histogram of cattle at risk for location ",locations[i])) + xlab("Cattle at risk") + ylab("Count")
ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_plots/Histogram_location_",locations[i],".png"), plot=p)
}

for (i in 1:38){
  df <- data.frame(prevalence=prevalence[i,])
  p <- ggplot(df, aes(x=prevalence)) + geom_histogram(bins=20, fill="#009E73", color="black") + ggtitle(paste0("Histogram of prevalence for location ",locations[i])) + xlab("Prevalence") + ylab("Count")
  ggsave(paste0("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_plots/Histogram_prevalence_location_",locations[i],".png"), plot=p)
}
