############################################################
################ NGA DATA HISTOGRAMS #######################
############################################################

rm(list = ls())

library(ggplot2)
library(dplyr)

############################################################
################ LOAD DATA #################################
############################################################

data_path <- "Code/TestSensSpec/adjusted_prevalence_all_sites.csv"

data_adjusted <- read.csv(data_path)

############################################################
################ STORE ADJUSTED PREVALENCE #################
############################################################

adjusted_prev <- rowMeans(
  data_adjusted[,8:1007],
  na.rm = TRUE
)

############################################################
################ LOAD NIGERIA BORDER #######################
############################################################

library(sf)
library(afrilearndata)

nga <- africountries[
  africountries$name == "Nigeria",
]

############################################################
################ CONVERT TO SF #############################
############################################################

data_sf <- st_as_sf(
  data_adjusted,
  coords = c("longitude", "latitude"),   # CHANGE IF NEEDED
  crs = 4326
)

############################################################
################ FILTER TO NIGERIA #########################
############################################################

data_sf <- st_intersection(
  data_sf,
  nga
)

############################################################
################ KEEP MATCHING ROWS ########################
############################################################

keep_rows <- as.numeric(
  rownames(data_sf)
)

############################################################
################ CREATE HISTOGRAM DATA #####################
############################################################

raw_data <- data.frame(
  prevalence = data_adjusted[keep_rows, 6]
)

adjusted_data <- data.frame(
  prevalence = adjusted_prev[keep_rows]
)
############################################################
################ LOAD GEOSTATISTICAL DATA ##################
############################################################


# Load and average all projection files in Projections_NGA directory
projections_dir <- "Code/Prevalence/Bovine BCT and PCR/Projections_NGA/"
projection_files <- list.files(projections_dir, pattern = "Projections_model_.*\\.csv", full.names = TRUE)

cat("Found", length(projection_files), "projection files to average\n")

# Function to load and process a single projection file
load_projection_file <- function(file_path) {
  tryCatch({
    data <- read.csv(file_path)
    # Add model identifier
    data$model_id <- basename(file_path)
    return(data)
  }, error = function(e) {
    cat("Error loading", file_path, ":", e$message, "\n")
    return(NULL)
  })
}

# Load all projection files
cat("Loading all projection files...\n")
all_projections <- do.call(rbind, lapply(projection_files, load_projection_file))

# Remove any NULL entries (failed loads)
all_projections <- all_projections[!is.null(all_projections), ]

cat("Total rows loaded:", nrow(all_projections), "\n")
cat("Unique models:", length(unique(all_projections$model_id)), "\n")

# Average across all models for each coordinate and variable combination
cat("Computing ensemble averages...\n")
ensemble_projections <- all_projections %>%
  dplyr::group_by(Latitude, Longitude, variable) %>%
  dplyr::summarise(
    value = mean(value, na.rm = TRUE),
    n_models = n(),
    .groups = 'drop'
  )

cat("Ensemble projections created with", nrow(ensemble_projections), "unique coordinate-variable combinations\n")
cat("Variables available:", unique(ensemble_projections$variable), "\n")

# Prepare data for mean, 2.5th and 97.5th percentiles from ensemble
projections_mean <- ensemble_projections %>%
  filter(variable == "Mean") %>%
  dplyr::select(Longitude, Latitude, prevalence = value)


geo_data <- data.frame(
  prevalence = projections_mean[,3]
)

############################################################
################ HISTOGRAM FUNCTION ########################
############################################################

plot_histogram <- function(
  data,
  fill_colour,
  title
){

  ##########################################################
  ################ COMPUTE MEAN ############################
  ##########################################################

  mean_value <- mean(
    data$prevalence,
    na.rm = TRUE
  )

  ##########################################################
  ################ CREATE PLOT #############################
  ##########################################################

  p <- ggplot(
    data,
    aes(x = prevalence)
  ) +

    geom_histogram(
      bins = 30,
      colour = "black",
      fill = fill_colour
    ) +

    ########################################################
    ################ MEAN LINE #############################
    ########################################################

    geom_vline(
      xintercept = mean_value,
      linewidth = 1,
      linetype = "dashed",
      colour = "red"
    ) +

    ########################################################
    ################ MEAN LABEL ############################
    ########################################################

    annotate(
      "text",
      x = mean_value,
      y = Inf,
      label = paste0(
        "Mean = ",
        round(mean_value, 3)
      ),
      vjust = 2,
      hjust = -0.1,
      colour = "red",
      fontface = "bold",
      size = 5
    ) +

    theme_bw() +

    labs(
      title = title,
      x = "Prevalence",
      y = "Count"
    ) +

    theme(

      plot.title = element_text(
        face = "bold",
        size = 16
      )
    )

  return(p)
}


############################################################
################ RAW DATA HISTOGRAM ########################
############################################################

p_raw <- plot_histogram(
  data = raw_data,
  fill_colour = "steelblue",
  title = "Raw NGA Data"
)

print(p_raw)

############################################################
################ ADJUSTED DATA HISTOGRAM ###################
############################################################

p_adjusted <- plot_histogram(
  data = adjusted_data,
  fill_colour = "darkorange",
  title = "Adjusted NGA Data"
)

print(p_adjusted)

############################################################
################ GEOSTATISTICAL HISTOGRAM ##################
############################################################

p_geo <- plot_histogram(
  data = geo_data,
  fill_colour = "forestgreen",
  title = "Geostatistical NGA Data"
)

print(p_geo)

############################################################
################ SAVE FIGURES ##############################
############################################################

ggsave(
  "Code/Prevalence/Bovine BCT and PCR/Diagnostics NGA/Figure_Raw_Histogram.png",
  p_raw,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "Code/Prevalence/Bovine BCT and PCR/Diagnostics NGA/Figure_Adjusted_Histogram.png",
  p_adjusted,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "Code/Prevalence/Bovine BCT and PCR/Diagnostics NGA/Figure_Geostatistical_Histogram.png",
  p_geo,
  width = 7,
  height = 5,
  dpi = 300
)