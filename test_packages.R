# Test key packages for DIMCAT Prevalence project
packages_to_test <- c("ggplot2", "dplyr", "tidyr", "forcats", "ggridges", "viridis", "ggtext", "readxl", "INLA", "sp", "sf", "terra", "raster", "geodata", "afrilearndata", "concaveman", "rstan", "bayesplot", "coda")

cat("Testing package installations:\n")
for (pkg in packages_to_test) {
  result <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    "OK"
  }, error = function(e) {
    "FAILED"
  })
  cat(sprintf("%-15s: %s\n", pkg, result))
}
cat("\nPackage installation summary complete!\n")