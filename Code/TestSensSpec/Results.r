library(ggplot2)
library(tidyr)
library(dplyr)
library(rstan)
library(reshape2)
library(posterior)
library(bayesplot)

fit<-readRDS("Code/TestSensSpec/latent_class_fit.rds")

post <- rstan::extract(fit)

df_plot <- data.frame(
  Se_pcr = post$Se_pcr,
  Se_hct = post$Se_hct,
  Sp_pcr = post$Sp_pcr,
  Sp_hct = post$Sp_hct,
  rel_Se = post$rel_Se,
  rel_Sp = post$rel_Sp
)

df_post <- melt(df_plot[, c("Se_pcr", "Se_hct")])

# Build analytical prior curve
eps <- 1e-6
x_vals <- seq(eps, 1 - eps, length.out = 1000)
y_vals <- dnorm(qlogis(x_vals), mean = 0, sd = 1) / (x_vals * (1 - x_vals))

# Convert analytical curve into "fake" data points for plotting
# by sampling proportionally to density (so it looks like a KDE)
n_prior_points <- nrow(df_post) / 2  # match posterior size roughly
prior_samples <- sample(x_vals, size = 1000000, replace = TRUE, prob = y_vals)
df_prior <- data.frame(value = prior_samples, variable = "Prior")

# Combine posterior and prior into one data frame
df_all <- rbind(df_post, df_prior)

# Plot all as densities with fill
p1 <- ggplot(df_all, aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.35) +
  scale_color_manual(
    name = "",
    values = c("Se_pcr" = "#D55E00", "Se_hct" = "#0072B2", "Prior" = "darkgreen"),
    labels = c("Se_pcr" = "PCR",
               "Se_hct" = "BCT/HCT",
               "Prior" = "Prior")
  ) +
  scale_fill_manual(
    name = "",
    values = c("Se_pcr" = "#D55E00", "Se_hct" = "#0072B2", "Prior" = "darkgreen"),
    labels = c("Se_pcr" = "PCR",
               "Se_hct" = "BCT/HCT",
               "Prior" = "Prior")
  ) +
  labs(title = "",
       x = "Sensitivity", y = "Density") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 16),       # axis titles
    axis.text = element_text(size = 14),     # legend title
    legend.text = element_text(size = 14)
  )
ggsave("Code/TestSensSpec/sensitivity_density.pdf", plot = p1, width = 8, height = 6)

eps <- 1e-6
x_vals <- seq(0.93, 1 - eps, length.out = 1000)
y_vals <- dnorm(qlogis(x_vals), mean = 0, sd = 1) / (x_vals * (1 - x_vals))

# Convert analytical curve into "fake" data points for plotting
# by sampling proportionally to density (so it looks like a KDE)
n_prior_points <- nrow(df_post) / 2  # match posterior size roughly
df_post <- melt(df_plot[, c("Sp_pcr", "Sp_hct")])

p1 <- ggplot(df_post, aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.35) +
  scale_color_manual(
    name = "",
    values = c("Sp_pcr" = "#D55E00", "Sp_hct" = "#0072B2"),
    labels = c("Sp_pcr" = "PCR",
               "Sp_hct" = "BCT/HCT")
  ) +
  scale_fill_manual(
    name = "",
    values = c("Sp_pcr" = "#D55E00", "Sp_hct" = "#0072B2"),
    labels = c("Sp_pcr" = "PCR",
               "Sp_hct" = "BCT/HCT")
  ) +
  labs(title = "",
       x = "Specificity", y = "Density") +
  coord_cartesian(xlim = c(0.98, 1))+
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 16),       # axis titles
    axis.text = element_text(size = 14),     # legend title
    legend.text = element_text(size = 14)
  )


print(p1)
ggsave("Code/TestSensSpec/specificity_density.pdf", plot = p1, width = 8, height = 6)



draws <- as_draws_df(as.array(fit))

variables(draws) <- gsub("^Se_pcr$", "PCR sensitivity", variables(draws))
variables(draws) <- gsub("^Sp_pcr$", "PCR specificity", variables(draws))
variables(draws) <- gsub("^Se_hct$", "BCT/HCT sensitivity", variables(draws))
variables(draws) <- gsub("^Sp_hct$", "BCT/HCT specificity", variables(draws))

mcmc_trace(draws, pars = c("PCR sensitivity", "PCR specificity",
                           "BCT/HCT sensitivity", "BCT/HCT specificity"))
ggsave("Code/TestSensSpec/mcmc_trace.pdf", width = 8, height = 6)




print(fit, pars = c("Se_pcr","Sp_pcr","Se_hct","Sp_hct","mu","sigma"), 
      digits_summary = 10)

  
#THIS CREATES PLOTS WITH A TRANSPARENT BACKGROUND FOR USE IN A POSTER

# Create identical plots with transparent background and custom font color
font_color <- "#EEE5D4"

# Recreate sensitivity plot with transparent background and custom font color
df_post_sens <- melt(df_plot[, c("Se_pcr", "Se_hct")])
df_all_sens <- rbind(df_post_sens, df_prior)

p1_transparent <- ggplot(df_all_sens, aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.35) +
  scale_color_manual(
    name = "",
    values = c("Se_pcr" = "#FFAA00", "Se_hct" = "#00FFFF", "Prior" = "#00FF00"),
    labels = c("Se_pcr" = "PCR",
               "Se_hct" = "BCT/HCT",
               "Prior" = "Prior")
  ) +
  scale_fill_manual(
    name = "",
    values = c("Se_pcr" = "#FFAA00", "Se_hct" = "#00FFFF", "Prior" = "#00FF00"),
    labels = c("Se_pcr" = "PCR",
               "Se_hct" = "BCT/HCT",
               "Prior" = "Prior")
  ) +
  labs(title = "",
       x = "Sensitivity", y = "Density") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = font_color),
    axis.title = element_text(size = 16, color = font_color),
    axis.text = element_text(size = 14, color = font_color),
    legend.text = element_text(size = 14, color = font_color),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )
ggsave("Code/TestSensSpec/sensitivity_density_transparent.pdf", plot = p1_transparent, 
       width = 8, height = 6, bg = "transparent")

# Recreate specificity plot with transparent background and custom font color
df_post_spec <- melt(df_plot[, c("Sp_pcr", "Sp_hct")])

p2_transparent <- ggplot(df_post_spec, aes(x = value, color = variable, fill = variable)) +
  geom_density(alpha = 0.35) +
  scale_color_manual(
    name = "",
    values = c("Sp_pcr" = "#FFAA00", "Sp_hct" = "#00FFFF"),
    labels = c("Sp_pcr" = "PCR",
               "Sp_hct" = "BCT/HCT")
  ) +
  scale_fill_manual(
    name = "",
    values = c("Sp_pcr" = "#FFAA00", "Sp_hct" = "#00FFFF"),
    labels = c("Sp_pcr" = "PCR",
               "Sp_hct" = "BCT/HCT")
  ) +
  labs(title = "",
       x = "Specificity", y = "Density") +
  coord_cartesian(xlim = c(0.98, 1)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = font_color),
    axis.title = element_text(size = 16, color = font_color),
    axis.text = element_text(size = 14, color = font_color),
    legend.text = element_text(size = 14, color = font_color),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )

print(p2_transparent)
ggsave("Code/TestSensSpec/specificity_density_transparent.pdf", plot = p2_transparent, 
       width = 8, height = 6, bg = "transparent")



