library(ggplot2)

data <- read.csv("Code/AbsencePresence/Results/Burkina FasoWAICs.csv")

waic_refs <- data$x

waic_df <- data.frame(waic = waic_refs, covariates = 1:length(waic_refs))
ggplot(waic_df, aes(x = covariates, y = waic)) +
  geom_line() +
  geom_point() +
  labs(x = "Iteration", y = "WAIC") +
  theme_bw(base_size = 18) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())
ggsave("Code/AbsencePresence/Results/Images/WAICs.png",width = 12, height = 8)
