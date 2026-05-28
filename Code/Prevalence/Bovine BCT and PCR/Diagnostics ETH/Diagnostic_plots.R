############################################################
################ REVIEWER 4 DIAGNOSTIC PLOTS ###############
############################################################

rm(list = ls())

library(ggplot2)
library(sf)
library(viridis)
library(afrilearndata)
library(ggspatial)
library(cowplot)
library(geodata)
library(dplyr)
library(units)

############################################################
################ LOAD DIAGNOSTICS ##########################
############################################################

diagnostics <- read.csv(
  "Code/Prevalence/Bovine BCT and PCR/Diagnostics ETH/Diagnostics.csv"
)

diagnostics$dist_to_data_km <- diagnostics$dist_to_data * 111



############################################################
################ LOAD STATE BOUNDARIES #####################
############################################################

nigeria_states <- gadm(
  country = "ETH",
  level = 1,
  path = tempdir()
)

nigeria_states_sf <- st_as_sf(
  nigeria_states
)

############################################################
################ LOAD LGA BOUNDARIES #######################
############################################################

nigeria_lgas <- gadm(
  country = "ETH",
  level = 3,
  path = tempdir()
)

nigeria_lgas_sf <- st_as_sf(
  nigeria_lgas
)

############################################################
################ CONVERT DIAGNOSTICS TO SF #################
############################################################

diagnostics_sf <- st_as_sf(
  diagnostics,
  coords = c("lon", "lat"),
  crs = 4326
)

############################################################
################ AGGREGATE TO LGA ##########################
############################################################

aggregate_diagnostic <- function(variable){

  ##########################################################
  ################ SPATIAL JOIN ############################
  ##########################################################

  joined <- st_join(
    diagnostics_sf,
    nigeria_lgas_sf
  )

  ##########################################################
  ################ DIRECT AGGREGATION ######################
  ##########################################################

  aggregated <- joined %>%

    st_drop_geometry() %>%

    filter(!is.na(GID_3)) %>%

    group_by(
      GID_3,
      NAME_2,
      NAME_1
    ) %>%

    summarise(
      value = mean(
        .data[[variable]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  ##########################################################
  ################ IDENTIFY MISSING LGAS ###################
  ##########################################################

  all_lgas <- nigeria_lgas_sf %>%

    st_drop_geometry() %>%

    select(
      GID_3,
      NAME_2,
      NAME_1
    )

  missing_lgas <- all_lgas %>%

    anti_join(
      aggregated,
      by = "GID_3"
    )

  ##########################################################
  ################ INTERPOLATE MISSING LGAS ################
  ##########################################################

  if(nrow(missing_lgas) > 0){

    cat(
      "Interpolating",
      nrow(missing_lgas),
      "missing LGAs for",
      variable,
      "\n"
    )

    for(i in seq_len(nrow(missing_lgas))){

      ######################################################
      ################ LGA CENTROID ########################
      ######################################################

      lga_id <- missing_lgas$GID_3[i]

      lga_centroid <- nigeria_lgas_sf %>%

        filter(
          GID_3 == lga_id
        ) %>%

        st_centroid(
          of_largest_polygon = TRUE
        )

      ######################################################
      ################ DISTANCES ###########################
      ######################################################

      distances <- st_distance(
        lga_centroid,
        diagnostics_sf
      )

      distances_numeric <- as.numeric(
        distances
      )

      ######################################################
      ################ CLOSEST POINTS ######################
      ######################################################

      n_neighbors <- min(
        5,
        nrow(diagnostics_sf)
      )

      closest_indices <- order(
        distances_numeric
      )[1:n_neighbors]

      ######################################################
      ################ IDW WEIGHTS #########################
      ######################################################

      weights <- 1 / (
        distances_numeric[closest_indices] + 1e-10
      )

      weights <- weights / sum(weights)

      ######################################################
      ################ INTERPOLATED VALUE ##################
      ######################################################

      interpolated_value <- sum(
        diagnostics_sf[[variable]][closest_indices] *
          weights,
        na.rm = TRUE
      )

      ######################################################
      ################ APPEND ##############################
      ######################################################

      new_row <- data.frame(

        GID_3 = lga_id,

        NAME_2 = missing_lgas$NAME_2[i],

        NAME_1 = missing_lgas$NAME_1[i],

        value = interpolated_value
      )

      aggregated <- rbind(
        aggregated,
        new_row
      )
    }
  }

  ##########################################################
  ################ JOIN BACK TO POLYGONS ##################
  ##########################################################

  map_data <- nigeria_lgas_sf %>%

    left_join(
      aggregated,
      by = c(
        "GID_3",
        "NAME_2",
        "NAME_1"
      )
    )

  return(map_data)
}

############################################################
################ CHOROPLETH FUNCTION #######################
############################################################

plot_map <- function(
  variable,
  title,
  filename,
  limits = NULL,
  viridis_option = "viridis"
){

  ##########################################################
  ################ AGGREGATE DATA ##########################
  ##########################################################

  map_data <- aggregate_diagnostic(
    variable
  )

  ##########################################################
  ################ INITIALISE PLOT #########################
  ##########################################################

  p <- ggplot()



  ##########################################################
  ################ CHOROPLETH ##############################
  ##########################################################

  p <- p +

    geom_sf(
      data = map_data,
      aes(fill = value),
      colour = "white",
      linewidth = 0.05
    )

  ##########################################################
  ################ STATE BORDERS ###########################
  ##########################################################

  p <- p +

    geom_sf(
      data = nigeria_states_sf,
      fill = NA,
      colour = "black",
      linewidth = 0.4
    )

 

  ##########################################################
  ################ COLOUR SCALE ############################
  ##########################################################

  if(is.null(limits)){

    p <- p +

      scale_fill_viridis_c(
        option = viridis_option,
        na.value = "grey90"
      )

  } else {

    p <- p +

      scale_fill_viridis_c(
        option = viridis_option,
        limits = limits,
        na.value = "grey90"
      )
  }

  ##########################################################
  ################ FORMATTING ##############################
  ##########################################################

  p <- p +

    coord_sf(
      expand = FALSE
    ) +

    labs(
      title = title,
      fill = NULL
    ) +

    theme_void() +

    theme(

      plot.title = element_text(
        size = 24,
        face = "bold",
        hjust = 0.5
      ),

      legend.position = "right",

      legend.key.height = unit(
        1.5,
        "cm"
      ),

      legend.key.width = unit(
        0.5,
        "cm"
      )
    )

  ##########################################################
  ################ SAVE ####################################
  ##########################################################

  ggsave(
    filename,
    p,
    width = 12,
    height = 10,
    dpi = 300
  )

  print(p)
}

############################################################
################ PREVALENCE UNCERTAINTY ####################
############################################################

plot_map(
  variable = "prevalence_sd",
  title = "Posterior standard deviation of prevalence (linear scale)",
  filename = "Code/Prevalence/Bovine BCT and PCR/Diagnostics ETH/Figure_PrevalenceSD.png"
)

############################################################
################ SPATIAL FIELD #############################
############################################################

plot_map(
  variable = "spatial_mean",
  title = "Posterior mean spatial field",
  filename = "Code/Prevalence/Bovine BCT and PCR/Diagnostics ETH/Figure_SpatialMean.png"
)

############################################################
################ NO GP PREVALENCE ##########################
############################################################

plot_map(
  variable = "prevalence_no_gp",
  title = "Average predicted prevalence without spatial field",
  filename = "Code/Prevalence/Bovine BCT and PCR/Diagnostics ETH/Figure_NoGP.png",
  limits = c(0,1)
)

############################################################
################ DISTANCE TO DATA ##########################
############################################################

plot_map(
  variable = "dist_to_data_km",
  title = "Average distance to nearest survey location (km)",
  filename = "Code/Prevalence/Bovine BCT and PCR/Diagnostics ETH/Figure_DistanceToData.png"
)