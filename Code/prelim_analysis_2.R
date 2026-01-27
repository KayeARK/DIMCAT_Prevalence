# prelim_analysis.R: loads in the Continental Atlas and produces
# initial visualisations of the raw data.

# Load in data manipulation and visualisation libraries.
library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(sf)

# Define the name of the continental atlas.
cont_atlas_path <- "Data"
cont_atlas_filename <- "240708_FAO_Continental_Atlas_AT_1990_2020.xlsx"

# Load in the raw data.
cont_atlas_loc <- read_excel(
  path = paste(cont_atlas_path, cont_atlas_filename, sep = "/"),
  sheet = 3,
  col_names = TRUE,
  col_types = "guess"
)
cont_atlas_epi <- read_excel(
  path = paste(cont_atlas_path, cont_atlas_filename, sep = "/"),
  sheet = 4,
  col_names = TRUE,
  col_types = "guess"
)

# Load in Africa country boundaries.
#africa <- read_sf(
#  "../data/shapefiles/afr_g2014_2013_0.shp"
#)

# Join the locations of all surveys.
cont_atlas_epi <- cont_atlas_epi %>%
  left_join(cont_atlas_loc, by = c("LOCATION_ID", "SOURCE_ID"))

# Ammend all starting years of the study.
# Filter out any years that are ambigious.
cont_atlas <- cont_atlas_epi %>%
  mutate(YEAR_ST = gsub("\\?", "", YEAR_ST)) %>%
  filter(grepl("^[0-9]{4}$", YEAR_ST)) %>%
  mutate(YEAR_ST = as.numeric(YEAR_ST))

# Visualise the location of all surveys.
gg_surveys <- cont_atlas %>%
  ggplot() +
  geom_sf(
    data = africa
  ) +
  geom_point(
    mapping = aes(
      x = LONG,
      y = LAT,
      fill = YEAR_ST,
      size = SAMPLE_SIZE
    ),
    colour = "grey50",
    shape = 21
  ) +
  xlab("") +
  ylab("") +
  scale_fill_fermenter(
    "Year of study",
    n.breaks = 8,
    palette = "Accent",
    show.limits = TRUE
  ) +
  scale_radius(
    "Number of animals\nsampled in the study",
    range = c(1.5, 6)
  ) +
  coord_sf(
    xlim = c(-20, 55),
    ylim = c(-40, 40),
    expand = FALSE
  ) +
  theme(
    legend.key.height = unit(1, "cm")
  ) +
  guides(
    fill = guide_colorsteps(
      frame.colour = "black",
      ticks = TRUE,
      ticks.colour = "black"
    )
  ) +
  labs(
    title = paste(
      "Location, timing and sample size of all studies",
      "in the Continental Atlas"
    ),
    caption = paste(
      "The figure shows the location, timing and sample size of",
      nrow(cont_atlas) -
        sum(is.na(cont_atlas$LONG)) -
        sum(is.na(cont_atlas$SAMPLE_SIZE)),
      "studies.\n",
      "Note that all studies shown cover a range",
      "of different animal hosts and parasites.\n",
      "Multi-year studies or studies with ambiguous year information",
      "were removed, corresponding to the removal of",
      nrow(cont_atlas_epi) - nrow(cont_atlas),
      "studies.\n",
      "Not shown in this figure:",
      sum(is.na(cont_atlas$LONG)),
      "studies did not have any corresponding location information.\n",
      sum(is.na(cont_atlas$SAMPLE_SIZE)),
      "studies did not have any corresponding sample size information.\n"
    )
  )

# Save the figure.
ggsave(
  filename = "all_surveys.png",
  dpi = 300,
  width = 9,
  height = 9,
  path = "figures",
)

gg_cattle <- cont_atlas %>%
  filter(grepl("cattle", SPECIES_AN, ignore.case = TRUE)) %>%
  {
    ggplot(.) +
      geom_sf(
        data = africa
      ) +
      geom_point(
        mapping = aes(
          x = LONG,
          y = LAT,
          fill = YEAR_ST,
          size = SAMPLE_SIZE
        ),
        colour = "grey50",
        shape = 21
      ) +
      xlab("") +
      ylab("") +
      scale_fill_fermenter(
        "Year of study",
        n.breaks = 8,
        palette = "Accent",
        show.limits = TRUE
      ) +
      scale_radius(
        "Number of animals\nsampled in the study",
        range = c(1.5, 6)
      ) +
      coord_sf(
        xlim = c(-20, 55),
        ylim = c(-40, 40),
        expand = FALSE
      ) +
      theme(
        legend.key.height = unit(1, "cm")
      ) +
      guides(
        fill = guide_colorsteps(
          frame.colour = "black",
          ticks = TRUE,
          ticks.colour = "black"
        )
      ) +
      labs(
        title = paste(
          "Location, timing and sample size of all studies",
          "in the Continental Atlas that were conducted in cattle"
        ),
        caption = paste(
          "The figure shows the location, timing and sample size of",
          nrow(.) -
            sum(is.na(.$LONG)) -
            sum(is.na(.$SAMPLE_SIZE)),
          "studies which feature cattle.\n",
          "Note that all studies shown cover a range",
          "of different parasites.\nSome studies may be",
          "grouped with other animal hosts (e.g. sheep and goats).\n",
          "Multi-year studies or studies with ambiguous year information",
          "were removed,corresponding to the removal of",
          sum(grepl("cattle", cont_atlas_epi$SPECIES_AN, ignore.case = TRUE)) - 
            nrow(.),
          "studies.\n",
          "Not shown in this figure:",
          sum(is.na(.$LONG)),
          "studies did not have any corresponding location information.\n",
          sum(is.na(.$SAMPLE_SIZE)),
          "studies did not have any corresponding sample size information.\n"
        )
      )
  }

# Save the figure.
ggsave(
  filename = "cattle_surveys.png",
  dpi = 300,
  width = 9,
  height = 9,
  path = "figures",
)

gg_camels <- cont_atlas %>%
  filter(grepl("camels", SPECIES_AN, ignore.case = TRUE)) %>%
  {
    ggplot(.) +
      geom_sf(
        data = africa
      ) +
      geom_point(
        mapping = aes(
          x = LONG,
          y = LAT,
          fill = YEAR_ST,
          size = SAMPLE_SIZE
        ),
        colour = "grey50",
        shape = 21
      ) +
      xlab("") +
      ylab("") +
      scale_fill_fermenter(
        "Year of study",
        n.breaks = 8,
        palette = "Accent",
        show.limits = TRUE
      ) +
      scale_radius(
        "Number of animals\nsampled in the study",
        range = c(1.5, 6)
      ) +
      coord_sf(
        xlim = c(-20, 55),
        ylim = c(-40, 40),
        expand = FALSE
      ) +
      theme(
        legend.key.height = unit(1, "cm")
      ) +
      guides(
        fill = guide_colorsteps(
          frame.colour = "black",
          ticks = TRUE,
          ticks.colour = "black"
        )
      ) +
      labs(
        title = paste(
          "Location, timing and sample size of all studies",
          "in the Continental Atlas that were conducted in camels"
        ),
        caption = paste(
          "The figure shows the location, timing and sample size of",
          nrow(.) -
            sum(is.na(.$LONG)) -
            sum(is.na(.$SAMPLE_SIZE)),
          "studies which feature camels.\n",
          "Note that all studies shown cover a range",
          "of different parasites.\nSome studies may be",
          "grouped with other animal hosts (e.g. horses and donkeys).\n",
          "Multi-year studies or studies with ambiguous year information",
          "were removed,corresponding to the removal of",
          sum(grepl("camels", cont_atlas_epi$SPECIES_AN, ignore.case = TRUE)) - 
            nrow(.),
          "studies.\n",
          "Not shown in this figure:",
          sum(is.na(.$LONG)),
          "studies did not have any corresponding location information.\n",
          sum(is.na(.$SAMPLE_SIZE)),
          "studies did not have any corresponding sample size information.\n"
        )
      )
  }

# Save the figure.
ggsave(
  filename = "camel_surveys.png",
  dpi = 300,
  width = 9,
  height = 9,
  path = "figures",
)

gg_pigs <- cont_atlas %>%
  filter(grepl("pigs", SPECIES_AN, ignore.case = TRUE)) %>%
  {
    ggplot(.) +
      geom_sf(
        data = africa
      ) +
      geom_point(
        mapping = aes(
          x = LONG,
          y = LAT,
          fill = YEAR_ST,
          size = SAMPLE_SIZE
        ),
        colour = "grey50",
        shape = 21
      ) +
      xlab("") +
      ylab("") +
      scale_fill_fermenter(
        "Year of study",
        n.breaks = 8,
        palette = "Accent",
        show.limits = TRUE
      ) +
      scale_radius(
        "Number of animals\nsampled in the study",
        range = c(1.5, 6)
      ) +
      coord_sf(
        xlim = c(-20, 55),
        ylim = c(-40, 40),
        expand = FALSE
      ) +
      theme(
        legend.key.height = unit(1, "cm")
      ) +
      guides(
        fill = guide_colorsteps(
          frame.colour = "black",
          ticks = TRUE,
          ticks.colour = "black"
        )
      ) +
      labs(
        title = paste(
          "Location, timing and sample size of all studies",
          "in the Continental Atlas that were conducted in pigs"
        ),
        caption = paste(
          "The figure shows the location, timing and sample size of",
          nrow(.) -
            sum(is.na(.$LONG)) -
            sum(is.na(.$SAMPLE_SIZE)),
          "studies which feature pigs.\n",
          "Note that all studies shown cover a range",
          "of different parasites.\nSome studies may be",
          "grouped with other animal hosts (e.g. goats and sheep).\n",
          "Multi-year studies or studies with ambiguous year information",
          "were removed,corresponding to the removal of",
          sum(grepl("pigs", cont_atlas_epi$SPECIES_AN, ignore.case = TRUE)) - 
            nrow(.),
          "studies.\n",
          "Not shown in this figure:",
          sum(is.na(.$LONG)),
          "studies did not have any corresponding location information.\n",
          sum(is.na(.$SAMPLE_SIZE)),
          "studies did not have any corresponding sample size information.\n"
        )
      )
  }

# Save the figure.
ggsave(
  filename = "pigs_surveys.png",
  dpi = 300,
  width = 9,
  height = 9,
  path = "figures",
)