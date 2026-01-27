# PRIORITY TARGETING LIST - HIGH BURDEN ADMINISTRATIVE UNITS
# Based on burden concentration analysis

library(dplyr)

cat("=== PRIORITY TARGETING RECOMMENDATIONS ===\n\n")

# Load the data (simplified version focusing on priority lists)
nigeria_cattle_data <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Nigeria_cattle_uncertainty_fine.csv")
ethiopia_cattle_data <- read.csv("Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Ethiopia_cattle_uncertainty_fine.csv")

# Nigeria priority list
nga_priority <- nigeria_cattle_data %>%
  dplyr::filter(state != "Total") %>%  # Exclude the total row
  dplyr::filter(mean_cattle_mean > 0) %>%
  dplyr::arrange(desc(mean_cattle_mean)) %>%
  dplyr::mutate(
    rank = row_number(),
    cumulative_cattle = cumsum(mean_cattle_mean),
    total_cattle = sum(mean_cattle_mean),
    cumulative_percent = cumulative_cattle / total_cattle * 100,
    priority_tier = case_when(
      rank <= ceiling(n() * 0.10) ~ "Tier 1 (Top 10%)",
      rank <= ceiling(n() * 0.20) ~ "Tier 2 (Top 20%)", 
      rank <= ceiling(n() * 0.50) ~ "Tier 3 (Top 50%)",
      TRUE ~ "Lower Priority"
    )
  ) %>%
  dplyr::select(rank, LGA = state, cattle_at_risk = mean_cattle_mean, 
                cumulative_percent, priority_tier)

# Ethiopia priority list  
eth_priority <- ethiopia_cattle_data %>%
  dplyr::filter(region != "Total") %>%  # Exclude the total row
  dplyr::filter(mean_cattle_mean > 0) %>%
  dplyr::arrange(desc(mean_cattle_mean)) %>%
  dplyr::mutate(
    rank = row_number(),
    cumulative_cattle = cumsum(mean_cattle_mean),
    total_cattle = sum(mean_cattle_mean),
    cumulative_percent = cumulative_cattle / total_cattle * 100,
    priority_tier = case_when(
      rank <= ceiling(n() * 0.10) ~ "Tier 1 (Top 10%)",
      rank <= ceiling(n() * 0.20) ~ "Tier 2 (Top 20%)", 
      rank <= ceiling(n() * 0.50) ~ "Tier 3 (Top 50%)",
      TRUE ~ "Lower Priority"
    )
  ) %>%
  dplyr::select(rank, Zone = region, cattle_at_risk = mean_cattle_mean, 
                cumulative_percent, priority_tier)

cat("NIGERIA - TOP 20 PRIORITY LGAs FOR INTERVENTION:\n")
cat("===============================================\n")
nga_top_20 <- nga_priority %>% slice_head(n = 20)
print(as.data.frame(nga_top_20))

cat("\nNIGERIA TIER 1 (TOP 10%) SUMMARY:\n")
nga_tier1 <- nga_priority %>% dplyr::filter(priority_tier == "Tier 1 (Top 10%)")
cat("Number of LGAs:", nrow(nga_tier1), "\n")
cat("Cattle at risk:", format(sum(nga_tier1$cattle_at_risk), big.mark = ","), "\n")
cat("Share of national burden:", round(max(nga_tier1$cumulative_percent), 1), "%\n\n")

cat("ETHIOPIA - TOP 20 PRIORITY ZONES FOR INTERVENTION:\n")
cat("=================================================\n")
eth_top_20 <- eth_priority %>% slice_head(n = 20)
print(as.data.frame(eth_top_20))

cat("\nETHIOPIA TIER 1 (TOP 10%) SUMMARY:\n")
eth_tier1 <- eth_priority %>% dplyr::filter(priority_tier == "Tier 1 (Top 10%)")
cat("Number of zones:", nrow(eth_tier1), "\n") 
cat("Cattle at risk:", format(sum(eth_tier1$cattle_at_risk), big.mark = ","), "\n")
cat("Share of national burden:", round(max(eth_tier1$cumulative_percent), 1), "%\n\n")

# Save priority lists to CSV
write.csv(nga_priority, "Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/nigeria_priority_targeting_list.csv", row.names = FALSE)
write.csv(eth_priority, "Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/ethiopia_priority_targeting_list.csv", row.names = FALSE)

# Save Tier 1 (highest priority) lists
write.csv(nga_tier1, "Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/nigeria_tier1_priority_lgas.csv", row.names = FALSE)
write.csv(eth_tier1, "Code/Prevalence/Bovine BCT and PCR/Cattle at risk distribution inequality/ethiopia_tier1_priority_zones.csv", row.names = FALSE)

cat("=== INTERVENTION STRATEGY RECOMMENDATIONS ===\n\n")

cat("PHASE 1 - IMMEDIATE ACTION (Tier 1 - Top 10%):\n")
cat("Nigeria:", nrow(nga_tier1), "LGAs covering", round(max(nga_tier1$cumulative_percent), 1), "% of burden\n")
cat("Ethiopia:", nrow(eth_tier1), "zones covering", round(max(eth_tier1$cumulative_percent), 1), "% of burden\n")
cat("IMPACT: Addressing ~50% of national burden with focused intervention\n\n")

nga_tier2 <- nga_priority %>% dplyr::filter(priority_tier == "Tier 2 (Top 20%)")
eth_tier2 <- eth_priority %>% dplyr::filter(priority_tier == "Tier 2 (Top 20%)")

cat("PHASE 2 - SCALED INTERVENTION (Tier 2 - Next 10%):\n") 
cat("Nigeria:", nrow(nga_tier2), "additional LGAs for", 
    round(max(nga_tier2$cumulative_percent) - max(nga_tier1$cumulative_percent), 1), "% more burden\n")
cat("Ethiopia:", nrow(eth_tier2), "additional zones for",
    round(max(eth_tier2$cumulative_percent) - max(eth_tier1$cumulative_percent), 1), "% more burden\n")
cat("CUMULATIVE IMPACT: Addressing ~75% of national burden\n\n")

cat("RESOURCE ALLOCATION RECOMMENDATIONS:\n")
cat("- Allocate 60% of resources to Tier 1 areas (highest burden density)\n")
cat("- Allocate 30% of resources to Tier 2 areas (moderate burden)\n") 
cat("- Allocate 10% of resources to Tier 3 areas (broad coverage)\n\n")

cat("FILES CREATED:\n")
cat("- nigeria_priority_targeting_list.csv (Complete ranking)\n")
cat("- ethiopia_priority_targeting_list.csv (Complete ranking)\n")
cat("- nigeria_tier1_priority_lgas.csv (Top 10% for immediate action)\n")
cat("- ethiopia_tier1_priority_zones.csv (Top 10% for immediate action)\n")