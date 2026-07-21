library(ggplot2)
library(tidyr)
library(dplyr)
library(rstan)
library(reshape2)
library(posterior)
library(bayesplot)
library(openxlsx)

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



library(openxlsx)
library(dplyr)

# ------------------------------------------------------------
# Create source-data tables
# ------------------------------------------------------------

posterior_draws <- data.frame(
  posterior_draw = seq_along(post$Se_pcr),
  PCR_sensitivity = post$Se_pcr,
  HCT_BCT_sensitivity = post$Se_hct,
  PCR_specificity = post$Sp_pcr,
  HCT_BCT_specificity = post$Sp_hct
)

# Analytical LogitNormal(0,1) prior curve
eps <- 1e-6

prior_curve <- data.frame(
  parameter_value = seq(eps, 1 - eps, length.out = 10000)
) |>
  mutate(
    prior_density =
      dnorm(qlogis(parameter_value), mean = 0, sd = 1) /
      (parameter_value * (1 - parameter_value)),
    prior_distribution = "LogitNormal(0,1)"
  )

# Posterior summary function
summarise_parameter <- function(x, parameter_name) {
  data.frame(
    parameter = parameter_name,
    posterior_mean = mean(x),
    posterior_median = median(x),
    posterior_sd = sd(x),
    lower_95_CrI = unname(quantile(x, 0.025)),
    upper_95_CrI = unname(quantile(x, 0.975)),
    number_of_draws = length(x)
  )
}

posterior_summary <- bind_rows(
  summarise_parameter(post$Se_pcr, "PCR sensitivity"),
  summarise_parameter(post$Se_hct, "HCT/BCT sensitivity"),
  summarise_parameter(post$Sp_pcr, "PCR specificity"),
  summarise_parameter(post$Sp_hct, "HCT/BCT specificity")
)

readme <- data.frame(
  item = c(
    "Workbook description",
    "Posterior_draws",
    "Prior_curve",
    "Posterior_summary",
    "Prior specification",
    "Credible intervals"
  ),
  description = c(
    "Source data underlying the diagnostic sensitivity and specificity plots.",
    paste(
      "Retained posterior samples for PCR and HCT/BCT sensitivity",
      "and specificity. Each row is one posterior draw."
    ),
    paste(
      "Analytical probability-density curve for the",
      "LogitNormal(0,1) prior used for sensitivity and specificity."
    ),
    paste(
      "Posterior means, medians, standard deviations and",
      "equal-tailed 95% credible intervals."
    ),
    "logit(parameter) ~ Normal(0,1).",
    "2.5th and 97.5th percentiles of the posterior draws."
  )
)

# ------------------------------------------------------------
# Create one workbook containing all outputs
# ------------------------------------------------------------

wb <- createWorkbook(
  creator = "AR Kaye",
  title = "Diagnostic test performance source data",
  subject = "Source data for sensitivity and specificity figure"
)

addWorksheet(wb, "README")
addWorksheet(wb, "Posterior_draws")
addWorksheet(wb, "Prior_curve")
addWorksheet(wb, "Posterior_summary")

writeData(wb, "README", readme)
writeData(wb, "Posterior_draws", posterior_draws)
writeData(wb, "Prior_curve", prior_curve)
writeData(wb, "Posterior_summary", posterior_summary)

# ------------------------------------------------------------
# Formatting
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

decimal_style <- createStyle(numFmt = "0.000000")
integer_style <- createStyle(numFmt = "0")
wrap_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

# README
addStyle(
  wb,
  sheet = "README",
  style = header_style,
  rows = 1,
  cols = 1:ncol(readme),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "README",
  style = wrap_style,
  rows = 2:(nrow(readme) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

setColWidths(wb, "README", cols = 1, widths = 24)
setColWidths(wb, "README", cols = 2, widths = 80)
freezePane(wb, "README", firstRow = TRUE)

# Posterior draws
addStyle(
  wb,
  sheet = "Posterior_draws",
  style = header_style,
  rows = 1,
  cols = 1:ncol(posterior_draws),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Posterior_draws",
  style = integer_style,
  rows = 2:(nrow(posterior_draws) + 1),
  cols = 1,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Posterior_draws",
  style = decimal_style,
  rows = 2:(nrow(posterior_draws) + 1),
  cols = 2:5,
  gridExpand = TRUE
)

setColWidths(
  wb,
  "Posterior_draws",
  cols = 1:ncol(posterior_draws),
  widths = "auto"
)

freezePane(wb, "Posterior_draws", firstRow = TRUE)
addFilter(
  wb,
  "Posterior_draws",
  row = 1,
  cols = 1:ncol(posterior_draws)
)

# Prior curve
addStyle(
  wb,
  sheet = "Prior_curve",
  style = header_style,
  rows = 1,
  cols = 1:ncol(prior_curve),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Prior_curve",
  style = decimal_style,
  rows = 2:(nrow(prior_curve) + 1),
  cols = 1:2,
  gridExpand = TRUE
)

setColWidths(
  wb,
  "Prior_curve",
  cols = 1:ncol(prior_curve),
  widths = "auto"
)

freezePane(wb, "Prior_curve", firstRow = TRUE)
addFilter(
  wb,
  "Prior_curve",
  row = 1,
  cols = 1:ncol(prior_curve)
)

# Posterior summaries
addStyle(
  wb,
  sheet = "Posterior_summary",
  style = header_style,
  rows = 1,
  cols = 1:ncol(posterior_summary),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Posterior_summary",
  style = decimal_style,
  rows = 2:(nrow(posterior_summary) + 1),
  cols = 2:6,
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "Posterior_summary",
  style = integer_style,
  rows = 2:(nrow(posterior_summary) + 1),
  cols = 7,
  gridExpand = TRUE
)

setColWidths(
  wb,
  "Posterior_summary",
  cols = 1:ncol(posterior_summary),
  widths = "auto"
)

freezePane(wb, "Posterior_summary", firstRow = TRUE)
addFilter(
  wb,
  "Posterior_summary",
  row = 1,
  cols = 1:ncol(posterior_summary)
)

# ------------------------------------------------------------
# Save the single Excel workbook
# ------------------------------------------------------------

output_file <- "Code/TestSensSpec/Figure_1_source_data.xlsx"

saveWorkbook(
  wb,
  file = output_file,
  overwrite = TRUE
)

message("Saved source-data workbook to: ", output_file)
