# Fully Bayesian prevalence adjustment using full posterior of Se & Sp
# - Inputs: y (positives), n (sample size)
# - Se_draws, Sp_draws: numeric vectors (posterior draws from Stan), length K
# - Returns: grid posterior for p, posterior samples of p and T=n*p, plots

library(ggplot2)
library(stats)
library(readxl)

# theta(p, Se, Sp) = p*Se + (1-p)*(1-Sp)
theta_fun <- function(p, Se, Sp) p*Se + (1 - p)*(1 - Sp)

# Build per-site grid posteriors conditional on each Se/Sp draw,
# then sample joint draws with shared Se/Sp.
adjust_prevalence_multi_site_grid <- function(
  y, n,
  Se_samples_HCT, Sp_samples_HCT,
  Se_samples_PCR, Sp_samples_PCR,
  test_types,
  prior_a = 1, prior_b = 1,
  grid_n = 1000,          # 1000 is usually plenty; 3000 is heavy
  n_draws = 1000,
  eps = 1e-12
) {
  I <- length(y)
  stopifnot(length(n) == I, length(test_types) == I)
  stopifnot(length(Se_samples_HCT) == length(Sp_samples_HCT))
  stopifnot(length(Se_samples_PCR) == length(Sp_samples_PCR))
  K <- length(Se_samples_HCT)
  stopifnot(length(Se_samples_PCR) == K)

  # common grid + prior
  p_grid <- seq(eps, 1 - eps, length.out = grid_n)
  log_prior <- dbeta(p_grid, prior_a, prior_b, log = TRUE)

  # preallocate outputs
  p_samples <- matrix(NA_real_, nrow = n_draws, ncol = I)
  T_samples <- matrix(NA_real_, nrow = n_draws, ncol = I)
  colnames(p_samples) <- paste0("site_", seq_len(I))
  colnames(T_samples) <- paste0("site_", seq_len(I))

  # choose which (Se,Sp) draw each posterior draw t will use
  # (this enforces cross-site coupling within each draw)
  k_idx <- sample.int(K, size = n_draws, replace = (n_draws > K))

  # process each unique k once
  uniq_k <- sort(unique(k_idx))
  count=1
  for (kk in uniq_k) {
    print(count)
    count=count+1
    # indices of posterior draws that use this k
    t_sel <- which(k_idx == kk)

    # pull Se/Sp for both tests
    Se_HCT_k <- Se_samples_HCT[kk]; Sp_HCT_k <- Sp_samples_HCT[kk]
    Se_PCR_k <- Se_samples_PCR[kk]; Sp_PCR_k <- Sp_samples_PCR[kk]

    # theta on grid by test type
    theta_HCT <- pmin(pmax(p_grid * Se_HCT_k + (1 - p_grid) * (1 - Sp_HCT_k), eps), 1 - eps)
    theta_PCR <- pmin(pmax(p_grid * Se_PCR_k + (1 - p_grid) * (1 - Sp_PCR_k), eps), 1 - eps)

    # for each site, build discrete posterior on the fly, then sample for all t in t_sel
    for (i in seq_len(I)) {
      theta <- if (test_types[i] == "HCT") theta_HCT else if (test_types[i] == "PCR") theta_PCR else stop("Unknown test type")
      # log posterior over grid
      loglik <- dbinom(y[i], size = n[i], prob = theta, log = TRUE)
      logpost <- log_prior + loglik

      # normalize stably (log-sum-exp)
      m <- max(logpost)
      w <- exp(logpost - m)
      prob <- w / sum(w)   # sums to 1

      # sample p for all draws t that use this k (vectorized via CDF)
      # build CDF once per site
      cdf <- cumsum(prob)
      u <- runif(length(t_sel))
      # integer index via CDF
      j <- findInterval(u, cdf) + 1L
      # guard against last-bin boundary
      j[j < 1L] <- 1L; j[j > grid_n] <- grid_n

      p_draws <- p_grid[j]
      p_samples[t_sel, i] <- p_draws

      # posterior predictive T = Binomial(n_i, p_i)
      T_samples[t_sel, i] <- rbinom(length(t_sel), size = n[i], prob = p_draws)
    }
    # discard per-k matrices immediately — nothing retained
  }

  # summaries
  summarize <- function(x) c(
    mean = mean(x), median = median(x),
    q025 = quantile(x, 0.025), q975 = quantile(x, 0.975)
  )
  p_summary <- apply(p_samples, 2, summarize)
  T_summary <- apply(T_samples, 2, summarize)

  list(
    p_grid = p_grid,
    p_samples = p_samples,
    T_samples = T_samples,
    p_summary = p_summary,
    T_summary = T_summary,
    k_idx = k_idx            # which (Se,Sp) draw each row used (for reproducibility)
  )
}




# -----------------------------
# EXAMPLE USAGE
# -----------------------------
# Replace these with your real posterior draws:
# Se_draws <- posterior$Se_pcr   # numeric vector length K
# Sp_draws <- posterior$Sp_pcr

fit_HCT<-readRDS("Code/TestSensSpec/latent_class_fit.rds")
post <- rstan::extract(fit_HCT)
Se_HCT_empirical <- post$Se_hct
Sp_HCT_empirical <- post$Sp_hct
Se_PCR_empirical <- post$Se_pcr
Sp_PCR_empirical <- post$Sp_pcr




data_path <- paste0("Data/ContAtlas_v3/Bovine data/AAT_PR_bovine_BCT-HCT_data_table.xls")
data <- read_excel(data_path)
#remove row if longitude or latitude is NA
data <- data[!is.na(data$Longitude) & !is.na(data$Latitude),]
#if number of infections is NA, set to Number_of_animals_tested*TPR
data$Number_of_infections[is.na(data$Number_of_infections)] <- round(data$Number_of_animal_tested[is.na(data$Number_of_infections)]*data$TPR[is.na(data$Number_of_infections)]/100)
positive_HCT <- round(data$Number_of_infections)
sample_size_HCT <- data$Number_of_animal_tested
positive_HCT[positive_HCT > sample_size_HCT] <- sample_size_HCT[positive_HCT > sample_size_HCT]
longitude_HCT <- data$Longitude
latitude_HCT <- data$Latitude
##repeat "HCT" for length(positive_HCT) times
label_HCT <- rep("HCT", length(positive_HCT))


data_path <- paste0("Data/ContAtlas_v3/Bovine data/AT_PREV_bovine_PCR_Table.xls")
data <- read_excel(data_path)
#remove row if longitude or latitude is NA
data <- data[!is.na(data$Longitude) & !is.na(data$Latitude),]
#if number of infections is NA, set to Number_of_animals_tested*TPR
data$Number_of_infections[is.na(data$Number_of_infections)] <- round(data$Number_of_animal_tested[is.na(data$Number_of_infections)]*data$T_ATPR[is.na(data$Number_of_infections)]/100)
positive_PCR <- round(data$Number_of_infections)
sample_size_PCR <- data$Number_of_animal_tested
longitude_PCR <- data$Longitude
latitude_PCR <- data$Latitude
label_PCR <- rep("PCR", length(positive_PCR))


#combine the two datasets
number_of_positives <- c(positive_HCT, positive_PCR)
number_of_tests <- c(sample_size_HCT, sample_size_PCR)
test_types <- c(label_HCT, label_PCR)
longitudes <- c(longitude_HCT, longitude_PCR)
latitudes <- c(latitude_HCT, latitude_PCR)




# number_of_tests <- c(2065,1417)
# number_of_positives <- c(236,185)
# test_types <- c("HCT","HCT")  # specify test type for each site

res <- adjust_prevalence_multi_site_grid(y = number_of_positives, n = number_of_tests,
                               Se_samples_HCT = Se_HCT_empirical,
                               Sp_samples_HCT = Sp_HCT_empirical,
                               Se_samples_PCR = Se_PCR_empirical,
                               Sp_samples_PCR = Sp_PCR_empirical,
                               test_types = test_types)



adjusted_cases<-res$T_samples
adjusted_prevalence <- res$p_samples

# #plot joint distribution of adjusted_cases[1] and adjusted_cases[2]
# p1 <- ggplot(data = as.data.frame(adjusted_prevalence), aes(x = site_1, y = site_2)) +
#   geom_point(alpha = 0.1) +
#   labs(title = "",
#        x = "Niamina East adjusted prevalence",
#        y = "Katangini adjusted prevalence") +
#   theme_minimal()
# ggsave("Code/TestSensSpec/joint_adjusted.pdf", width=7, height=4)
# # ggsave("posterior_T.png", res$T_plot, width=7, height=4)

#save adjusted_prevalence as a csv, including longitudes and latitudes and sample sizes adjusting for the
# fact that adjusted_prevalence is a matrix with n_draws rows and I columns
adjusted_prevalence_df <- as.data.frame(adjusted_prevalence)
#transpose the dataframe
adjusted_prevalence_df <- as.data.frame(t(adjusted_prevalence_df))
adjusted_prevalence_df$longitude <- longitudes
adjusted_prevalence_df$latitude <- latitudes
adjusted_prevalence_df$sample_size <- number_of_tests
adjusted_prevalence_df$positive_tests <- number_of_positives
adjusted_prevalence_df$original_prevalence <- number_of_positives/number_of_tests
adjusted_prevalence_df$test_type <- test_types
#put longitudes, latitudes, sample sizes, number of positives, original prevalence and test types as the first columns
adjusted_prevalence_df <- adjusted_prevalence_df[, c(ncol(adjusted_prevalence_df)-5, ncol(adjusted_prevalence_df)-4, ncol(adjusted_prevalence_df)-3, ncol(adjusted_prevalence_df)-2, ncol(adjusted_prevalence_df)-1, ncol(adjusted_prevalence_df), 1:(ncol(adjusted_prevalence_df)-6))]
write.csv(adjusted_prevalence_df, "Code/TestSensSpec/adjusted_prevalence_all_sites.csv", row.names = TRUE)


adjusted_cases_df <- as.data.frame(adjusted_cases)
#transpose the dataframe
adjusted_cases_df <- as.data.frame(t(adjusted_cases_df))
adjusted_cases_df$longitude <- longitudes
adjusted_cases_df$latitude <- latitudes
adjusted_cases_df$sample_size <- number_of_tests
adjusted_cases_df$positive_tests <- number_of_positives
adjusted_cases_df$original_prevalence <- number_of_positives/number_of_tests
adjusted_cases_df$test_type <- test_types
#put longitudes, latitudes, sample sizes, number of positives, original prevalence and test types as the first columns
adjusted_cases_df <- adjusted_cases_df[, c(ncol(adjusted_cases_df)-5, ncol(adjusted_cases_df)-4, ncol(adjusted_cases_df)-3, ncol(adjusted_cases_df)-2, ncol(adjusted_cases_df)-1, ncol(adjusted_cases_df), 1:(ncol(adjusted_cases_df)-6))]
write.csv(adjusted_cases_df, "Code/TestSensSpec/adjusted_cases_all_sites.csv", row.names = TRUE)
