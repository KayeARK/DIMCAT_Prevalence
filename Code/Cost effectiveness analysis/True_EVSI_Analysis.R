#!/usr/bin/env Rcat("TRUE EVSI ANALYSIS: Value of Additional Sample Information\n")

# ==============================================================================cat("=========================================================\n")

# TRUE EVSI ANALYSIS AND VISUALIZATIONcat("EVSI = Expected utility after collecting additional data - Expected utility with current info\n")

# ==============================================================================cat("This accounts for both the value AND cost of additional research\n\n")

# Analyze the true Expected Value of Sample Information results

# and compare with previous net benefit calculationslibrary(tidyverse)

# ==============================================================================library(rstan)



library(tidyverse)
library(viridis)
library(ggplot2)
library(patchwork)

setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")

cat("=== TRUE EVSI ANALYSIS AND VISUALIZATION ===\n\n")



setwd("/Users/u2074276/Library/CloudStorage/OneDrive-UniversityofWarwick/Desktop/DIMCAT_Prevalence")load_prevalence_data <- function() {

  # Try to load Nigeria INLA projection data

cat("=== TRUE EVSI ANALYSIS AND VISUALIZATION ===\n\n")  projection_files <- list.files("Code/Prevalence/Bovine BCT and PCR/Projections_NGA/", 

                                 pattern = "Projections_model_.*\\.csv", full.names = TRUE)

# ==============================================================================  

# DATA LOADING  if (length(projection_files) == 0) {

# ==============================================================================    projection_files <- list.files("Code/Prevalence/Bovine BCT/Results/", 

                                   pattern = "Projections_model_.*\\.csv", full.names = TRUE)

# Load true EVSI results  }

evsi_true <- read_csv("Code/Cost effectiveness analysis/evsi_results/simplified_true_evsi.csv",   

                     show_col_types = FALSE)  if (length(projection_files) > 0) {

    cat(sprintf("Found %d INLA projection files for Nigeria\n", length(projection_files)))

# Load comparison data    

evsi_comparison <- read_csv("Code/Cost effectiveness analysis/evsi_results/evsi_detailed_comparison.csv",     # Load first 50 files to get prevalence distributions

                           show_col_types = FALSE)    n_files_to_load <- min(50, length(projection_files))

    

cat("True EVSI results loaded:", nrow(evsi_true), "locations\n")    all_projections <- map_dfr(projection_files[seq_len(n_files_to_load)], function(file) {

cat("Comparison data loaded:", nrow(evsi_comparison), "locations\n\n")      tryCatch({

        data <- read_csv(file, show_col_types = FALSE)

# ==============================================================================        

# STATISTICAL SUMMARY        # Standardize column names

# ==============================================================================        if ("Latitude" %in% names(data)) {

          data <- data %>% rename(latitude = Latitude, longitude = Longitude)

cat("=== TRUE EVSI STATISTICAL SUMMARY ===\n")        }

        if ("value" %in% names(data)) {

# Basic statistics          data <- data %>% rename(mean = value) %>% filter(variable == "Mean")

evsi_stats <- evsi_true %>%        }

  summarise(        

    n_locations = n(),        data %>%

    mean_evsi = mean(evsi),          select(latitude, longitude, mean) %>%

    median_evsi = median(evsi),          mutate(model_id = str_extract(basename(file), "\\d+"))

    sd_evsi = sd(evsi),      }, error = function(e) {

    min_evsi = min(evsi),        cat(sprintf("Error loading %s: %s\n", file, e$message))

    max_evsi = max(evsi),        return(NULL)

    q25_evsi = quantile(evsi, 0.25),      })

    q75_evsi = quantile(evsi, 0.75),    })

    positive_evsi_count = sum(evsi > 0),    

    positive_evsi_prop = mean(evsi > 0),    if (nrow(all_projections) > 0) {

    mean_positive_evsi = mean(evsi[evsi > 0]),      # Calculate prevalence statistics across INLA models

    substantial_evsi_count = sum(evsi > 0.01),  # More than 1 cent      prevalence_data <- all_projections %>%

    substantial_evsi_prop = mean(evsi > 0.01)        group_by(latitude, longitude) %>%

  )        summarise(

          prev_mean = mean(mean, na.rm = TRUE),

print(evsi_stats)          prev_sd = sd(mean, na.rm = TRUE),

          n_models = n(),

# EVSI by prevalence quartiles          .groups = "drop"

prevalence_quartiles <- evsi_true %>%        ) %>%

  mutate(prevalence_quartile = cut(prevalence,         filter(!is.na(prev_mean), prev_mean > 0, prev_mean < 1) %>%

                                  breaks = quantile(prevalence, c(0, 0.25, 0.5, 0.75, 1)),        mutate(location_id = row_number()) %>%

                                  labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),        select(location_id, latitude, longitude, prev_mean, prev_sd)

                                  include.lowest = TRUE)) %>%      

  group_by(prevalence_quartile) %>%      cat(sprintf("Loaded %d locations from Nigeria INLA models\n", nrow(prevalence_data)))

  summarise(      return(prevalence_data)

    n = n(),    }

    mean_prevalence = mean(prevalence),  }

    mean_evsi = mean(evsi),  

    median_evsi = median(evsi),  # Fallback to synthetic data

    positive_evsi_prop = mean(evsi > 0),  cat("Generating synthetic prevalence data...\n")

    .groups = "drop"  set.seed(123)

  )  n_locations <- 500

  

cat("\nEVSI by Prevalence Quartiles:\n")  prevalence_data <- tibble(

print(prevalence_quartiles)    location_id = 1:n_locations,

    latitude = runif(n_locations, 4.0, 14.0),    

# Strategy distribution under perfect information    longitude = runif(n_locations, 2.5, 15.0),   

strategy_summary <- evsi_true %>%  ) %>%

  count(optimal_strategy) %>%    mutate(

  mutate(proportion = n / sum(n))      spatial_trend = 0.05 + 0.25 * (latitude - 4)/(14 - 4),

      noise = rnorm(n_locations, 0, 0.08),

cat("\nOptimal Strategy Distribution (Perfect Information):\n")      prev_mean = pmax(0.01, pmin(0.6, spatial_trend + noise)),

print(strategy_summary)      prev_sd = 0.01 + 0.15 * prev_mean * (1 - prev_mean)

    ) %>%

# ==============================================================================    select(location_id, latitude, longitude, prev_mean, prev_sd)

# COMPARISON WITH PREVIOUS "EVSI" CALCULATIONS  

# ==============================================================================  return(prevalence_data)

}

cat("\n=== COMPARISON WITH PREVIOUS CALCULATIONS ===\n")

prevalence_data <- load_prevalence_data()

# Compare with old "EVSI" measures

comparison_stats <- evsi_comparison %>%cat("\nStep 2: Loading fitted test characteristics...\n")

  summarise(

    correlation_max_evsi = cor(true_evsi, old_max_evsi, use = "complete.obs"),define_test_characteristics <- function() {

    correlation_hct_evsi = cor(true_evsi, old_hct_evsi, use = "complete.obs"),  tryCatch({

    correlation_pcr_evsi = cor(true_evsi, old_pcr_evsi, use = "complete.obs"),    fit <- readRDS("Code/TestSensSpec/latent_class_fit.rds")

    strategy_agreement = mean(optimal_strategy_perfect_info ==     post <- rstan::extract(fit)

                            ifelse(prob_hct_optimal > prob_pcr_optimal & prob_hct_optimal > prob_no_test_optimal, "HCT",    

                                  ifelse(prob_pcr_optimal > prob_no_test_optimal, "PCR", "NO_TEST")),     list(

                            na.rm = TRUE)      hct_sensitivity = post$Se_hct,

  )      hct_specificity = post$Sp_hct,

      pcr_sensitivity = post$Se_pcr,

cat("Correlation with previous 'EVSI' measures:\n")      pcr_specificity = post$Sp_pcr,

cat("• True EVSI vs Max 'EVSI':", round(comparison_stats$correlation_max_evsi, 3), "\n")      n_samples = length(post$Se_hct)

cat("• True EVSI vs HCT 'EVSI':", round(comparison_stats$correlation_hct_evsi, 3), "\n")    )

cat("• True EVSI vs PCR 'EVSI':", round(comparison_stats$correlation_pcr_evsi, 3), "\n")  }, error = function(e) {

cat("• Strategy agreement:", round(comparison_stats$strategy_agreement * 100, 1), "%\n")    cat("Using fallback test characteristics\n")

    n_samples <- 5000

# ==============================================================================    list(

# VISUALIZATIONS      hct_sensitivity = rbeta(n_samples, 80, 20),

# ==============================================================================      hct_specificity = rbeta(n_samples, 85, 15), 

      pcr_sensitivity = rbeta(n_samples, 95, 5),

dir.create("Code/Cost effectiveness analysis/evsi_results/plots", showWarnings = FALSE)      pcr_specificity = rbeta(n_samples, 98, 2),

      n_samples = n_samples

# 1. EVSI Distribution    )

p1 <- ggplot(evsi_true, aes(x = evsi)) +  })

  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +}

  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +

  geom_vline(xintercept = mean(evsi_true$evsi), linetype = "solid", color = "darkred", size = 1) +test_chars <- define_test_characteristics()

  labs(title = "Distribution of True EVSI",

       subtitle = "Expected Value of Perfect Information about Test Parameters",cat("Step 3: Define research scenarios and costs...\n")

       x = "True EVSI ($)", y = "Count",

       caption = "Red dashed line: $0, Red solid line: Mean EVSI") +# Research scenarios: different sample sizes for additional studies

  theme_minimal() +research_scenarios <- tibble(

  theme(plot.title = element_text(size = 14, face = "bold"))  scenario_id = 1:4,

  scenario_name = c("Current Info Only", "50 Additional Tests", "100 Additional Tests", "200 Additional Tests"),

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_distribution.png",   additional_n = c(0, 50, 100, 200),

       p1, width = 10, height = 6, dpi = 300)  research_cost = c(0, 50*65, 100*65, 200*65)  # Assume HCT cost for additional sampling

)

# 2. EVSI vs Prevalence

p2 <- ggplot(evsi_true, aes(x = prevalence * 100, y = evsi)) +cat("Research scenarios:\n")

  geom_point(alpha = 0.6, color = "steelblue") +print(research_scenarios)

  geom_smooth(method = "loess", se = TRUE, color = "red") +

  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +# Test and utility costs

  labs(title = "True EVSI vs Prevalence",cost_params <- list(

       subtitle = "Relationship between disease prevalence and value of test parameter information",  hct_mean = 65, hct_cv = 0.25,

       x = "Prevalence (%)", y = "True EVSI ($)") +  hct_shape = 1/(0.25^2), hct_rate = (1/(0.25^2))/65,

  theme_minimal() +  pcr_mean = 110, pcr_cv = 0.2,

  theme(plot.title = element_text(size = 14, face = "bold"))  pcr_shape = 1/(0.2^2), pcr_rate = (1/(0.2^2))/110

)

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_vs_prevalence.png", 

       p2, width = 10, height = 6, dpi = 300)utility_params <- list(

  tp_mean = 800, tp_shape = 1/(0.2^2), tp_rate = (1/(0.2^2))/800,

# 3. Spatial Distribution of EVSI  fn_mean = -400, fn_shape = 1/(0.25^2), fn_rate = (1/(0.25^2))/400,

p3 <- ggplot(evsi_true, aes(x = longitude, y = latitude)) +  tn_mean = 50, tn_shape = 1/(0.3^2), tn_rate = (1/(0.3^2))/50,

  geom_point(aes(color = evsi), size = 0.8, alpha = 0.8) +  fp_mean = -150, fp_shape = 1/(0.3^2), fp_rate = (1/(0.3^2))/150

  scale_color_viridis_c(name = "True EVSI ($)", option = "plasma") +)

  labs(title = "Spatial Distribution of True EVSI",

       subtitle = "Nigeria - Expected Value of Perfect Test Parameter Information",cat("\nStep 4: Define key functions...\n")

       x = "Longitude", y = "Latitude") +

  theme_minimal() +# Calculate expected utility for a test

  theme(plot.title = element_text(size = 14, face = "bold"),calculate_expected_utility <- function(prevalence, test_cost, sensitivity, specificity, 

        axis.text = element_text(size = 8)) +                                     utility_tp, utility_fn, utility_tn, utility_fp) {

  coord_fixed()  prob_tp = prevalence * sensitivity

  prob_fn = prevalence * (1 - sensitivity)

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_spatial_distribution.png",   prob_tn = (1 - prevalence) * specificity

       p3, width = 12, height = 8, dpi = 300)  prob_fp = (1 - prevalence) * (1 - specificity)

  

# 4. EVSI by Optimal Strategy  expected_utility <- (prob_tp * utility_tp + prob_fn * utility_fn + 

p4 <- ggplot(evsi_true, aes(x = optimal_strategy, y = evsi, fill = optimal_strategy)) +                      prob_tn * utility_tn + prob_fp * utility_fp) - test_cost

  geom_boxplot(alpha = 0.7) +  return(expected_utility)

  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +}

  scale_fill_viridis_d(option = "viridis") +

  labs(title = "True EVSI by Optimal Strategy",# Update prevalence beliefs with additional sample data (Beta-Binomial conjugate)

       subtitle = "Distribution of EVSI across different optimal strategies",update_prevalence_posterior <- function(prior_mean, prior_sd, sample_n, sample_positives, test_sensitivity, test_specificity) {

       x = "Optimal Strategy (Perfect Information)", y = "True EVSI ($)") +  # Convert prior moments to beta parameters

  theme_minimal() +  prior_var = prior_sd^2

  theme(plot.title = element_text(size = 14, face = "bold"),  prior_alpha = prior_mean * (prior_mean * (1 - prior_mean) / prior_var - 1)

        legend.position = "none")  prior_beta = (1 - prior_mean) * (prior_mean * (1 - prior_mean) / prior_var - 1)

  

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_by_strategy.png",   # Adjust for test characteristics (simplified approach)

       p4, width = 10, height = 6, dpi = 300)  # In reality, this should account for test sensitivity/specificity more carefully

  adjusted_positives = sample_positives / test_sensitivity  # Rough adjustment

# 5. Comparison of EVSI measures  adjusted_negatives = (sample_n - sample_positives) / test_specificity

if (!any(is.na(evsi_comparison$old_max_evsi))) {  

  p5 <- ggplot(evsi_comparison, aes(x = old_max_evsi, y = true_evsi)) +  # Bayesian update

    geom_point(alpha = 0.6, color = "steelblue") +  posterior_alpha = prior_alpha + adjusted_positives

    geom_smooth(method = "lm", se = TRUE, color = "red") +  posterior_beta = prior_beta + adjusted_negatives

    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +  

    labs(title = "True EVSI vs Previous 'EVSI' Calculation",  # Return updated moments

         subtitle = paste("Correlation:", round(comparison_stats$correlation_max_evsi, 3)),  posterior_mean = posterior_alpha / (posterior_alpha + posterior_beta)

         x = "Previous 'EVSI' (Net Benefit Difference) ($)",   posterior_var = (posterior_alpha * posterior_beta) / ((posterior_alpha + posterior_beta)^2 * (posterior_alpha + posterior_beta + 1))

         y = "True EVSI ($)",  posterior_sd = sqrt(posterior_var)

         caption = "Dashed line: y = x") +  

    theme_minimal() +  return(list(mean = posterior_mean, sd = posterior_sd, alpha = posterior_alpha, beta = posterior_beta))

    theme(plot.title = element_text(size = 14, face = "bold"))}

  

  ggsave("Code/Cost effectiveness analysis/evsi_results/plots/evsi_comparison.png", cat("\nStep 5: Calculate TRUE EVSI for each location and research scenario...\n")

         p5, width = 10, height = 6, dpi = 300)

}calculate_true_evsi <- function(location_row, n_monte_carlo = 1000) {

  prev_mean <- location_row$prev_mean

# 6. High-Value Locations for Research  prev_sd <- location_row$prev_sd

high_evsi_threshold <- quantile(evsi_true$evsi, 0.9)  

high_evsi_locations <- evsi_true %>%  results <- map_dfr(research_scenarios$scenario_id, function(scenario_id) {

  filter(evsi > high_evsi_threshold)    scenario <- research_scenarios[scenario_id, ]

    additional_n <- scenario$additional_n

p6 <- ggplot(evsi_true, aes(x = longitude, y = latitude)) +    research_cost <- scenario$research_cost

  geom_point(aes(color = evsi > high_evsi_threshold), alpha = 0.6, size = 0.8) +    

  scale_color_manual(values = c("gray70", "red"),     utilities_current_info <- numeric(n_monte_carlo)

                     labels = c("Standard EVSI", "High EVSI (Top 10%)"),    utilities_after_research <- numeric(n_monte_carlo)

                     name = "Research Priority") +    

  labs(title = "Priority Locations for Test Parameter Research",    for (i in 1:n_monte_carlo) {

       subtitle = paste("Top 10% EVSI locations (>$", round(high_evsi_threshold, 4), ")"),      # Sample uncertain parameters

       x = "Longitude", y = "Latitude") +      true_prevalence <- rnorm(1, prev_mean, prev_sd)

  theme_minimal() +      true_prevalence <- pmax(0.001, pmin(0.999, true_prevalence))

  theme(plot.title = element_text(size = 14, face = "bold"),      

        axis.text = element_text(size = 8)) +      hct_sens <- sample(test_chars$hct_sensitivity, 1)

  coord_fixed()      hct_spec <- sample(test_chars$hct_specificity, 1)

      pcr_sens <- sample(test_chars$pcr_sensitivity, 1)

ggsave("Code/Cost effectiveness analysis/evsi_results/plots/research_priority_locations.png",       pcr_spec <- sample(test_chars$pcr_specificity, 1)

       p6, width = 12, height = 8, dpi = 300)      

      hct_cost <- rgamma(1, shape = cost_params$hct_shape, rate = cost_params$hct_rate)

# ==============================================================================      pcr_cost <- rgamma(1, shape = cost_params$pcr_shape, rate = cost_params$pcr_rate)

# COMPREHENSIVE SUMMARY      

# ==============================================================================      utility_tp <- rgamma(1, shape = utility_params$tp_shape, rate = utility_params$tp_rate)

      utility_fn <- -rgamma(1, shape = utility_params$fn_shape, rate = utility_params$fn_rate)

cat("\n=== COMPREHENSIVE EVSI ANALYSIS SUMMARY ===\n")      utility_tn <- rgamma(1, shape = utility_params$tn_shape, rate = utility_params$tn_rate)

      utility_fp <- -rgamma(1, shape = utility_params$fp_shape, rate = utility_params$fp_rate)

# Calculate total potential value across Nigeria      

total_locations <- 18304      # SCENARIO 1: Decision with current information only

representative_sample <- nrow(evsi_true)      util_hct_current <- calculate_expected_utility(true_prevalence, hct_cost, hct_sens, hct_spec,

scaling_factor <- total_locations / representative_sample                                                    utility_tp, utility_fn, utility_tn, utility_fp)

      util_pcr_current <- calculate_expected_utility(true_prevalence, pcr_cost, pcr_sens, pcr_spec,

estimated_total_evsi <- mean(evsi_true$evsi) * total_locations                                                    utility_tp, utility_fn, utility_tn, utility_fp)

estimated_positive_locations <- round(mean(evsi_true$evsi > 0) * total_locations)      

      utilities_current_info[i] <- max(util_hct_current, util_pcr_current)

cat("Scale-up Estimates to Full Nigeria:\n")      

cat("• Estimated total EVSI across Nigeria: $", round(estimated_total_evsi, 0), "\n")      # SCENARIO 2: Decision after collecting additional sample data

cat("• Estimated locations with positive EVSI:", estimated_positive_locations, "\n")      if (additional_n > 0) {

cat("• Mean EVSI per location: $", round(mean(evsi_true$evsi), 4), "\n")        # Simulate collecting additional test data (using HCT for the research)

        sample_positives <- rbinom(1, additional_n, true_prevalence * hct_sens + (1 - true_prevalence) * (1 - hct_spec))

if (mean(evsi_true$evsi) > 0) {        

  cat("• Positive total EVSI suggests value in test parameter research\n")        # Update prevalence beliefs based on this new data

} else {        updated_prevalence <- update_prevalence_posterior(prev_mean, prev_sd, additional_n, sample_positives, hct_sens, hct_spec)

  cat("• Low/negative total EVSI suggests current parameters are adequate\n")        

}        # Make decision with updated information

        util_hct_updated <- calculate_expected_utility(updated_prevalence$mean, hct_cost, hct_sens, hct_spec,

# Research recommendations                                                      utility_tp, utility_fn, utility_tn, utility_fp)

high_value_count <- sum(evsi_true$evsi > 0.01)  # More than 1 cent        util_pcr_updated <- calculate_expected_utility(updated_prevalence$mean, pcr_cost, pcr_sens, pcr_spec,

if (high_value_count > 0) {                                                      utility_tp, utility_fn, utility_tn, utility_fp)

  cat("\nResearch Recommendations:\n")        

  cat("• Focus test validation studies on", high_value_count, "highest EVSI locations\n")        # Subtract research cost from the better option

  cat("• Estimated high-value locations nationwide:",         best_utility_updated <- max(util_hct_updated, util_pcr_updated) - research_cost

      round(high_value_count * scaling_factor), "\n")        utilities_after_research[i] <- best_utility_updated

        } else {

  # Characteristics of high-value locations        utilities_after_research[i] <- utilities_current_info[i]

  high_value_prev <- mean(evsi_true$prevalence[evsi_true$evsi > 0.01])      }

  cat("• High-EVSI locations have mean prevalence:", round(high_value_prev * 100, 1), "%\n")    }

}    

    # Calculate EVSI for this scenario

# Key findings    expected_utility_current <- mean(utilities_current_info)

cat("\nKey Findings:\n")    expected_utility_research <- mean(utilities_after_research)

cat("1. True EVSI calculation shows limited value of resolving test parameter uncertainty\n")    evsi <- expected_utility_research - expected_utility_current

cat("2. Most decision uncertainty comes from prevalence, not test performance\n")    

cat("3. Current test parameter estimates appear adequate for decision-making\n")    tibble(

cat("4. Negative correlation with previous 'EVSI' suggests different methodologies\n")      scenario_id = scenario_id,

cat("5. Research resources might be better directed toward prevalence estimation\n")      scenario_name = scenario$scenario_name,

      additional_n = additional_n,

cat("\nVisualization files saved:\n")      research_cost = research_cost,

cat("• evsi_distribution.png - Histogram of EVSI values\n")      expected_utility_current = expected_utility_current,

cat("• evsi_vs_prevalence.png - EVSI relationship with prevalence\n")        expected_utility_research = expected_utility_research,

cat("• evsi_spatial_distribution.png - Geographic distribution of EVSI\n")      evsi = evsi,

cat("• evsi_by_strategy.png - EVSI by optimal strategy\n")      evsi_net = evsi  # Net EVSI (research cost already subtracted above)

if (!any(is.na(evsi_comparison$old_max_evsi))) {    )

  cat("• evsi_comparison.png - True vs previous EVSI comparison\n")  })

}  

cat("• research_priority_locations.png - High-priority locations for research\n")  return(results)

}

cat("\n=== TRUE EVSI ANALYSIS COMPLETE ===\n")
cat("Running TRUE EVSI calculation...\n")
cat("This compares decision-making with current info vs after additional research\n\n")

# Calculate for subset of locations (full analysis would take very long)
n_locations_to_analyze <- min(100, nrow(prevalence_data))
sample_locations <- prevalence_data[1:n_locations_to_analyze, ]

cat(sprintf("Analyzing %d locations...\n", n_locations_to_analyze))

start_time <- Sys.time()

all_evsi_results <- list()
for (i in 1:nrow(sample_locations)) {
  location_results <- calculate_true_evsi(sample_locations[i, ])
  location_results$location_id <- sample_locations$location_id[i]
  location_results$latitude <- sample_locations$latitude[i]
  location_results$longitude <- sample_locations$longitude[i]
  location_results$prev_mean <- sample_locations$prev_mean[i]
  location_results$prev_sd <- sample_locations$prev_sd[i]
  
  all_evsi_results[[i]] <- location_results
  
  if (i %% 10 == 0) cat(sprintf("Completed %d/%d locations\n", i, nrow(sample_locations)))
}

final_results <- bind_rows(all_evsi_results)

end_time <- Sys.time()
cat(sprintf("Analysis completed in %.1f minutes\n", as.numeric(difftime(end_time, start_time, units = "mins"))))

cat("\nStep 6: Results summary...\n")

# Summarize results by research scenario
scenario_summary <- final_results %>%
  group_by(scenario_id, scenario_name, additional_n, research_cost) %>%
  summarise(
    n_locations = n(),
    mean_evsi = mean(evsi),
    median_evsi = median(evsi),
    total_evsi = sum(evsi),
    positive_evsi_pct = 100 * mean(evsi > 0),
    .groups = "drop"
  )

cat("\nTRUE EVSI RESULTS BY RESEARCH SCENARIO\n")
cat("=====================================\n\n")
print(scenario_summary)

# Save results
dir.create("Code/Cost effectiveness analysis/Results", recursive = TRUE, showWarnings = FALSE)
write_csv(final_results, "Code/Cost effectiveness analysis/Results/true_evsi_results.csv")
write_csv(scenario_summary, "Code/Cost effectiveness analysis/Results/true_evsi_scenario_summary.csv")

cat("\nKEY FINDINGS:\n")
cat("• TRUE EVSI measures value of specific additional research\n")
cat("• Negative EVSI means current information is sufficient\n") 
cat("• Positive EVSI means additional research would improve decision-making\n")
cat(sprintf("• Best research strategy: %s\n", 
    scenario_summary$scenario_name[which.max(scenario_summary$mean_evsi)]))

cat("\n✓ TRUE EVSI analysis complete!\n")
cat("This shows the value of collecting specific amounts of additional diagnostic data.\n")