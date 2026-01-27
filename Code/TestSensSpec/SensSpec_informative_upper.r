library(readxl)
library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(coda)

#define logit function
logit <- function(p) {
  log(p / (1 - p))
}

data <- read_excel("Data/To_infer_sens_spec_v2.xlsx")



#change all entries in data$Diagnostic_CAT_2 which are PCR multi-species to PCR
data$Parasit_Diagnostic_CAT_2[data$Parasit_Diagnostic_CAT_2 == "BCT"] <- "HCT"
data$Parasit_Diagnostic_CAT_2[data$Parasit_Diagnostic_CAT_2 == "HCT * Animal Inoculation"] <- "HCT"
data$Parasit_Diagnostic_CAT_2[data$Parasit_Diagnostic_CAT_2 == "BCT * Blood smears * Mice Inoculation"] <- "HCT"
data$Parasit_Diagnostic_CAT_2[data$Parasit_Diagnostic_CAT_2 == "BCT + blood smears"] <- "BCT"


data$Molec_Diagnostic_CAT_2[data$Molec_Diagnostic_CAT_2 == "PCR multi-species"] <- "PCR"
data$Molec_Diagnostic_CAT_2[data$Molec_Diagnostic_CAT_2 == "PCR species-specific"] <- "PCR"
data$Molec_Diagnostic_CAT_2[data$Molec_Diagnostic_CAT_2 == "PCR multi-species + PCR species-specific"] <- "PCR"

#remove rows where SPECIES_AN_en is not "Cattle"
data <- data[data$SPECIES_AN_en == "Cattle", ]

#extract unique SURVEY_ID  
ids <- unique(data$LOCATION_ID)

HCTs <- data$Parasit_infections
PCRs <- data$Molec_infections
sample_size_arr <- data$SAMPLE_SIZE




# latent_class_analysis.R
# REQUIREMENTS:
# install.packages(c("rstan","bayesplot","ggplot2","coda","dplyr"))
# rstan setup: see https://mc-stan.org/rstan/install.html


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ---------------------------
# Replace these with your data or ensure they are in the environment:
# HCTs         : vector (counts or proportions) of HCT positives per site
# PCRs         : vector (counts or proportions) of PCR positives per site
# sample_size_arr : vector of sample sizes per site (integers)
# ---------------------------

# Example placeholder (comment out when using your own arrays)
#HCTs <- c(5, 20, 3, 0, 12)
#PCRs <- c(6, 18, 2, 0, 10)
#sample_size_arr <- c(100, 200, 80, 50, 150)

#generate synthetic HCTs and PCRs
# synthetic_HCT_sensitivity <- 0.25
# synthetic_PCR_sensitivity <- 0.75
# synthetic_HCT_specificity <- 0.60
# synthetic_PCR_specificity <- 0.10
# synthetic_prevalences <- runif(100, 0.01, 0.5) # example prevalences for 5 sites
# synthetic_sample_sizes <- sample(50:200, 100, replace = TRUE) # random sample sizes between 50 and 200
# set.seed(123) # for reproducibility
# HCT_prob_pos <- synthetic_HCT_sensitivity * synthetic_prevalences + (1 - synthetic_HCT_specificity) * (1 - synthetic_prevalences)
# PCR_prob_pos <- synthetic_PCR_sensitivity * synthetic_prevalences + (1 - synthetic_PCR_specificity) * (1 - synthetic_prevalences)
# HCTs <- rbinom(length(synthetic_sample_sizes), synthetic_sample_sizes, HCT_prob_pos)
# PCRs <- rbinom(length(synthetic_sample_sizes), synthetic_sample_sizes, PCR_prob_pos)
# sample_size_arr <- synthetic_sample_sizes



# Basic checks
if(!(exists("HCTs") && exists("PCRs") && exists("sample_size_arr"))){
  stop("Please provide HCTs, PCRs, and sample_size_arr in the environment before running.")
}

I <- length(sample_size_arr)
if(length(HCTs) != I || length(PCRs) != I){
  stop("Lengths of HCTs, PCRs, and sample_size_arr must be equal.")
}

# If HCTs/PCRs look like proportions (all between 0 and 1), convert to counts:
is_proportion_vec <- function(x) all(x >= 0 & x <= 1) && any(x > 0 & x < 1)
if(is_proportion_vec(HCTs) || is_proportion_vec(PCRs)){
  if(is_proportion_vec(HCTs)){
    HCT_counts <- round(HCTs * sample_size_arr)
    message("Converted HCTs (proportions) -> counts by rounding.")
  } else {
    HCT_counts <- as.integer(HCTs)
  }
  if(is_proportion_vec(PCRs)){
    PCR_counts <- round(PCRs * sample_size_arr)
    message("Converted PCRs (proportions) -> counts by rounding.")
  } else {
    PCR_counts <- as.integer(PCRs)
  }
} else {
  HCT_counts <- as.integer(HCTs)
  PCR_counts <- as.integer(PCRs)
}

# Final sanity
if(any(HCT_counts < 0) || any(PCR_counts < 0) || any(sample_size_arr <= 0)){
  stop("Counts must be non-negative and sample sizes positive.")
}
if(any(HCT_counts > sample_size_arr) || any(PCR_counts > sample_size_arr)){
  warning("Some counts exceed sample size at a site - check your inputs.")
}

# ---------------------------
# Stan model (latent-class, marginals only)
# ---------------------------
# stan_model_code <- "
# data {
#   int<lower=1> I;
#   int<lower=0> n[I];
#   int<lower=0> y_pcr[I];
#   int<lower=0> y_hct[I];
# }
# parameters {
#   real<lower=0,upper=1> Se_pcr;
#   real<lower=0,upper=1> Sp_pcr;
#   real<lower=0,upper=1> Se_hct;
#   real<lower=0,upper=1> Sp_hct;

#   vector[I] z;
#   real mu;
#   real<lower=0> sigma;
# }
# transformed parameters {
#   vector[I] logit_p;
#   vector[I] p;
#   vector[I] theta_pcr;
#   vector[I] theta_hct;

#   logit_p = mu + sigma * z;   // non-centered
#   p = inv_logit(logit_p);

#   for (i in 1:I) p[i] = inv_logit(logit_p[i]);

#   for (i in 1:I) {
#     theta_pcr[i] = p[i]*Se_pcr + (1 - p[i])*(1 - Sp_pcr);
#     theta_hct[i] = p[i]*Se_hct + (1 - p[i])*(1 - Sp_hct);
#   }
# }
# model {
#   // Priors (adjust as needed)
#   Se_pcr ~ beta(24,5);
#   Sp_pcr ~ uniform(0.75,1);
#   Se_hct ~ beta(5,24);
#   Sp_hct ~ beta(24,5);

#   mu ~ normal(-2, 2);
#   sigma ~ normal(0, 1);
#   logit_p ~ normal(mu, sigma);

#   // Likelihood
#   y_pcr ~ binomial(n, theta_pcr);
#   y_hct ~ binomial(n, theta_hct);
# }
# generated quantities {
#   real rel_Se = Se_hct / Se_pcr;
#   real rel_Sp = Sp_hct / Sp_pcr;

#   // Posterior predictive draws
#   int y_pcr_rep[I];
#   int y_hct_rep[I];
#   for (i in 1:I) {
#     y_pcr_rep[i] = binomial_rng(n[i], theta_pcr[i]);
#     y_hct_rep[i] = binomial_rng(n[i], theta_hct[i]);
#   }
# }
# "

stan_model_code <- "
data {
  int<lower=1> I;
  int<lower=0> n[I];
  int<lower=0> y_pcr[I];
  int<lower=0> y_hct[I];
}
parameters {
  // Se/Sp on logit scale (unconstrained reparam)
  real mu_logit_Se_pcr;
  real mu_logit_Sp_pcr;
  real mu_logit_Se_hct;
  real mu_logit_Sp_hct;

  // hierarchical prevalence (non-centered)
  real mu;                 // mean logit prevalence
  real<lower=0> sigma;     // sd of logit prevalence (problematic param)
  vector[I] z;             // standard normal auxiliary

}
transformed parameters {
  vector[I] logit_p;
  vector[I] p;
  vector[I] theta_pcr;
  vector[I] theta_hct;

  // map logits back to probabilities
  real<lower=0,upper=1> Se_pcr = inv_logit(mu_logit_Se_pcr);
  real<lower=0,upper=1> Sp_pcr = inv_logit(mu_logit_Sp_pcr);
  real<lower=0,upper=1> Se_hct = inv_logit(mu_logit_Se_hct);
  real<lower=0,upper=1> Sp_hct = inv_logit(mu_logit_Sp_hct);

  // non-centered parameterization
  logit_p = mu + sigma * z;
  p = inv_logit(logit_p);

  for (i in 1:I) {
    theta_pcr[i] = p[i]*Se_pcr + (1 - p[i])*(1 - Sp_pcr);
    theta_hct[i] = p[i]*Se_hct + (1 - p[i])*(1 - Sp_hct);
  }
}
model {
  // ---------- Priors ----------
  // Priors for Se/Sp on logit scale (tune to literature)
  mu_logit_Se_pcr ~ normal(2, 1);  // centers Se_pcr ~ 0.5
  mu_logit_Sp_pcr ~ normal(2, 1);  // centers Sp_pcr ~ 0.5
  mu_logit_Se_hct ~ normal(2, 1); // centers Se_hct ~ 0.5
  mu_logit_Sp_hct ~ normal(2, 1); // centers Sp_hct ~ 0.5

  // Prior for hierarchical prevalence mean & SD
  mu ~ normal(logit(0.1), 1);         // centers global prevalence at 0.1
  sigma ~ normal(0, 1);     // standard half-normal prior on sigma

  // standard normal for non-centered z
  z ~ normal(0, 1);

  // ---------- Likelihood ----------
  y_pcr ~ binomial(n, theta_pcr);
  y_hct ~ binomial(n, theta_hct);
}
generated quantities {
  real rel_Se = Se_hct / Se_pcr;
  real rel_Sp = Sp_hct / Sp_pcr;

  int y_pcr_rep[I];
  int y_hct_rep[I];
  for (i in 1:I) {
    y_pcr_rep[i] = binomial_rng(n[i], theta_pcr[i]);
    y_hct_rep[i] = binomial_rng(n[i], theta_hct[i]);
  }
}
"



# Write and compile Stan model
stan_file <- tempfile(fileext = ".stan")
writeLines(stan_model_code, con = stan_file)
sm <- stan_model(stan_file)

# Data for Stan
stan_data <- list(
  I = I,
  n = as.integer(sample_size_arr),
  y_pcr = as.integer(PCR_counts),
  y_hct = as.integer(HCT_counts)
)

init_fun <- function(chain_id = 1) {
  list(
    mu_logit_Se_pcr = runif(1,0,5),
    mu_logit_Sp_pcr = runif(1,1.1,5),
    mu_logit_Se_hct = runif(1,0,5),
    mu_logit_Sp_hct = runif(1,1.1,5),
    mu = rnorm(1,logit(0.1),1),
    sigma = runif(1,0,1),
    z = rnorm(I, 0, 1)
  )
}

# Fit model
fit <- sampling(sm,
                data = stan_data,
                init = init_fun,
                chains = 5,
                iter = 10000,
                warmup = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                seed = 1234)

# Print core summaries
print(fit, pars = c("Se_pcr", "Sp_pcr", "Se_hct", "Sp_hct", "rel_Se", "rel_Sp", "mu","sigma"), probs = c(0.025,0.5,0.975))

# Extract posterior arrays
post <- rstan::extract(fit)
# Compute Pr(Se_hct > Se_pcr) and Pr(Sp_hct > Sp_pcr)
prob_Se_higher <- mean(post$Se_hct > post$Se_pcr)
prob_Sp_higher <- mean(post$Sp_hct > post$Sp_pcr)

cat(sprintf("\\nPr(Se_hct > Se_pcr) = %.3f\\n", prob_Se_higher))
cat(sprintf("Pr(Sp_hct > Sp_pcr) = %.3f\\n", prob_Sp_higher))

# Posterior summaries table:
summ_df <- data.frame(
  param = c("Se_pcr","Se_hct","Sp_pcr","Sp_hct","rel_Se","rel_Sp"),
  mean = c(mean(post$Se_pcr), mean(post$Se_hct), mean(post$Sp_pcr), mean(post$Sp_hct), mean(post$rel_Se), mean(post$rel_Sp)),
  median = c(median(post$Se_pcr), median(post$Se_hct), median(post$Sp_pcr), median(post$Sp_hct), median(post$rel_Se), median(post$rel_Sp)),
  ci2.5 = c(quantile(post$Se_pcr, .025), quantile(post$Se_hct, .025), quantile(post$Sp_pcr, .025), quantile(post$Sp_hct, .025), quantile(post$rel_Se, .025), quantile(post$rel_Sp, .025)),
  ci97.5 = c(quantile(post$Se_pcr, .975), quantile(post$Se_hct, .975), quantile(post$Sp_pcr, .975), quantile(post$Sp_hct, .975), quantile(post$rel_Se, .975), quantile(post$rel_Sp, .975))
)
print(summ_df)

# Simple PPC: compare observed vs posterior predictive mean per site
y_pcr_rep <- post$y_pcr_rep   # dims: iterations x I
y_hct_rep <- post$y_hct_rep

ppc_summary <- data.frame(
  site = 1:I,
  n = sample_size_arr,
  obs_pcr = PCR_counts,
  pred_pcr_mean = apply(y_pcr_rep, 2, mean),
  obs_hct = HCT_counts,
  pred_hct_mean = apply(y_hct_rep, 2, mean)
)

print(ppc_summary)

# Plot posterior densities for sensitivities/specificities and relative metrics
df_plot <- data.frame(
  Se_pcr = post$Se_pcr,
  Se_hct = post$Se_hct,
  Sp_pcr = post$Sp_pcr,
  Sp_hct = post$Sp_hct,
  rel_Se = post$rel_Se,
  rel_Sp = post$rel_Sp
)

# Save fit for later
saveRDS(fit, file = "Code/TestSensSpec/latent_class_fit_upper.rds")
message("Fit saved to latent_class_fit_upper.rds")

