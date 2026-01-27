
rm(list=ls())

covariates_all<-c("elevation","precipitation","tavg","tmax","human_fp","pop_den","tree","grassland",
"shrub","cropland","built","bare","water","wetland","mangrove","cattle","tsetse")

covariates_present <-  rep(0, length(covariates_all))

for (i in 1:200){
data <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Covariates_NGA/Covariates_model_",i,".csv"))
covariates <- data[,1]

#compare covariates to covariates_all, if there, add a 1 to covariates_present
for (j in 1:length(covariates_all)){
  if (covariates_all[j] %in% covariates){
    covariates_present[j] <- covariates_present[j] + 1
  }
}
}
#make empty array of length 5

print(covariates_present)

#plot using ggplot
library(ggplot2)
library(reshape2)
df <- data.frame(covariates_all, covariates_present)
colnames(df) <- c("Covariate", "Count")
df$Covariate <- factor(df$Covariate, levels = df$Covariate[order(df$Count, decreasing = TRUE)])
p <- ggplot(df, aes(x=Covariate, y=Count)) + 
  geom_bar(stat="identity", fill="steelblue")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title="Covariate Selection Frequency in Bovine BCT and PCR Models",
       x="Covariate",
       y="Number of Times Selected (out of 100 models)") +
  ylim(0, 200)
ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_Selection_Frequency.png", plot = p, width = 10, height = 6)

# ====== ALTERNATIVE VISUALIZATIONS FOR COVARIATE IMPORTANCE ======

# 1. Relative Importance (as percentage)
library(dplyr)
library(tidyr)
library(viridis)

df$Percentage <- (df$Count / 200) * 100
df$Importance_Category <- cut(df$Percentage, 
                             breaks = c(0, 25, 50, 75, 100), 
                             labels = c("Low (0-25%)", "Medium (25-50%)", "High (50-75%)", "Very High (75-100%)"),
                             include.lowest = TRUE)

# Horizontal bar chart with color coding
p2 <- ggplot(df, aes(x = reorder(Covariate, Percentage), y = Percentage, fill = Importance_Category)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis_d(name = "Importance Level", option = "plasma") +
  labs(title = "Covariate Importance in Model Selection",
       subtitle = "Percentage of models (out of 200) containing each covariate",
       x = "Covariate",
       y = "Selection Percentage (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11))

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_Importance_Percentage.png", plot = p2, width = 10, height = 8)

# 2. Heatmap showing covariate co-occurrence patterns
# Create a matrix to track which covariates appear together
covariate_matrix <- matrix(0, nrow = length(covariates_all), ncol = length(covariates_all))
rownames(covariate_matrix) <- covariates_all
colnames(covariate_matrix) <- covariates_all

# Fill the co-occurrence matrix
for (i in 1:200) {
  data <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Covariates_NGA/Covariates_model_",i,".csv"))
  covariates <- data[,1]
  
  # For each pair of covariates that appear together, increment the counter
  for (j in 1:length(covariates)) {
    for (k in 1:length(covariates)) {
      if (j != k) {  # Don't count self-pairs
        covar1_idx <- which(covariates_all == covariates[j])
        covar2_idx <- which(covariates_all == covariates[k])
        if (length(covar1_idx) > 0 && length(covar2_idx) > 0) {
          covariate_matrix[covar1_idx, covar2_idx] <- covariate_matrix[covar1_idx, covar2_idx] + 1
        }
      }
    }
  }
}

# Convert to data frame for ggplot
library(reshape2)
heatmap_data <- melt(covariate_matrix)
colnames(heatmap_data) <- c("Covariate1", "Covariate2", "Co_occurrence")

p3 <- ggplot(heatmap_data, aes(x = Covariate1, y = Covariate2, fill = Co_occurrence)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Co-occurrence\nCount", option = "inferno") +
  labs(title = "Covariate Co-occurrence Heatmap",
       subtitle = "How often pairs of covariates appear together in models",
       x = "Covariate", y = "Covariate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11))

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_Cooccurrence_Heatmap.png", plot = p3, width = 12, height = 10)

# 3. Dot plot with confidence intervals (showing variability in selection)
# Calculate some basic statistics about selection patterns
set.seed(42)  # For reproducible "uncertainty" estimates

df$Lower_CI <- pmax(0, df$Count - sqrt(df$Count) * 1.96)
df$Upper_CI <- pmin(200, df$Count + sqrt(df$Count) * 1.96)

p4 <- ggplot(df, aes(x = reorder(Covariate, Count), y = Count)) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI, color = Importance_Category), size = 0.8) +
  coord_flip() +
  scale_color_viridis_d(name = "Importance Level", option = "plasma") +
  labs(title = "Covariate Selection with Uncertainty",
       subtitle = "Frequency of selection with approximate confidence intervals",
       x = "Covariate",
       y = "Number of Times Selected (out of 200 models)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11))

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_Selection_Uncertainty.png", plot = p4, width = 10, height = 8)

# 4. Circular/Radial plot showing covariate importance
library(ggplot2)
# Prepare data for circular plot
df$angle <- seq(0, 2*pi, length.out = nrow(df) + 1)[1:nrow(df)]
df$x <- cos(df$angle) * df$Percentage
df$y <- sin(df$angle) * df$Percentage

p5 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(size = Percentage, color = Importance_Category), alpha = 0.7) +
  geom_text(aes(label = Covariate), hjust = 0.5, vjust = -1.5, size = 3) +
  scale_size_continuous(name = "Selection %", range = c(2, 10)) +
  scale_color_viridis_d(name = "Importance Level", option = "plasma") +
  labs(title = "Covariate Importance - Radial View",
       subtitle = "Distance from center indicates selection frequency") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom") +
  coord_equal()

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_Radial_Plot.png", plot = p5, width = 10, height = 10)

# 5. Model complexity analysis - how many covariates per model?
model_complexity <- numeric(200)
for (i in 1:200) {
  data <- read.csv(paste0("Code/Prevalence/Bovine BCT and PCR/Covariates_NGA/Covariates_model_",i,".csv"))
  model_complexity[i] <- nrow(data)
}

complexity_df <- data.frame(
  Model = 1:200,
  Number_of_Covariates = model_complexity
)

p6 <- ggplot(complexity_df, aes(x = Number_of_Covariates)) +
  geom_histogram(bins = 15, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(aes(xintercept = mean(Number_of_Covariates)), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Model Complexity",
       subtitle = paste("Average number of covariates per model:", round(mean(model_complexity), 1)),
       x = "Number of Covariates per Model",
       y = "Frequency (Number of Models)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11))

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Model_Complexity_Distribution.png", plot = p6, width = 10, height = 6)

# 6. Covariate categories analysis (grouping related variables)
df$Category <- case_when(
  df$Covariate %in% c("elevation") ~ "Topography",
  df$Covariate %in% c("precipitation", "tavg", "tmax") ~ "Climate",
  df$Covariate %in% c("human_fp", "pop_den") ~ "Human factors",
  df$Covariate %in% c("tree", "grassland", "shrub", "cropland", "built", "bare", "water", "wetland", "mangrove") ~ "Land use/cover",
  df$Covariate %in% c("cattle") ~ "Livestock",
  df$Covariate %in% c("tsetse") ~ "Vector",
  TRUE ~ "Other"
)

category_summary <- df %>%
  group_by(Category) %>%
  summarise(
    Mean_Selection_Rate = mean(Percentage),
    Total_Selections = sum(Count),
    Number_of_Variables = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Selection_Rate))

p7 <- ggplot(category_summary, aes(x = reorder(Category, Mean_Selection_Rate), y = Mean_Selection_Rate)) +
  geom_col(aes(fill = Number_of_Variables), alpha = 0.8) +
  coord_flip() +
  scale_fill_viridis_c(name = "Number of\nVariables", option = "viridis") +
  labs(title = "Covariate Category Importance",
       subtitle = "Average selection rate by covariate category",
       x = "Covariate Category",
       y = "Average Selection Rate (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11))

ggsave("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Covariate_Category_Importance.png", plot = p7, width = 10, height = 6)

# Print summary statistics
cat("\n=== COVARIATE IMPORTANCE ANALYSIS SUMMARY ===\n")
cat("Total models analyzed:", 200, "\n")
cat("Average number of covariates per model:", round(mean(model_complexity), 1), "\n")
cat("Range of model complexity:", min(model_complexity), "to", max(model_complexity), "covariates\n")
cat("\nTop 5 most important covariates:\n")
top_covariates <- df[order(df$Count, decreasing = TRUE)[1:5], c("Covariate", "Count", "Percentage")]
print(top_covariates)
cat("\nCovariate category summary:\n")
print(category_summary)

cat("\nVisualization files created:\n")
cat("1. Covariate_Selection_Frequency.png - Original frequency bar chart\n")
cat("2. Covariate_Importance_Percentage.png - Horizontal percentage chart with importance levels\n")
cat("3. Covariate_Cooccurrence_Heatmap.png - Shows which covariates appear together\n")
cat("4. Covariate_Selection_Uncertainty.png - Dot plot with confidence intervals\n")
cat("5. Covariate_Radial_Plot.png - Circular visualization of importance\n")
cat("6. Model_Complexity_Distribution.png - Distribution of number of covariates per model\n")
cat("7. Covariate_Category_Importance.png - Importance by covariate category\n")
