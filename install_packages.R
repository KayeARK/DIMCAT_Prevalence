# Install remaining packages for DIMCAT Prevalence project

# Install devtools if needed
if (!require(devtools, quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org/")
}

# Install afrilearndata from GitHub
library(devtools)
cat("Installing afrilearndata from GitHub...\n")
devtools::install_github("afrimapr/afrilearndata")

cat("Installation complete!\n")