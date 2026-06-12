rm(list = ls())

library(sf)
library(raster)
library(terra)
library(dplyr)
library(readxl)
library(geodata)

cat("Loading Nigerian state boundaries...\n")

# State boundaries
states <- gadm(country = "NGA", level = 1, path = tempdir())
states_sf <- st_as_sf(states)

# ------------------------------------------------------------------
# Load cattle raster
# ------------------------------------------------------------------

cattle_path <- "Data/Covariates/livestock/cattle_2020/5_Ct_2020_Da.tif"

cat("Loading cattle raster...\n")
r_cattle <- raster(cattle_path)

# Same conversion used in cattle-at-risk script
if (grepl("2020", cattle_path)) {

  res_deg <- res(r_cattle)

  res_km <- res_deg * 111.32

  cell_area_km2 <- res_km[1] * res_km[2]

  conversion_factor <- cell_area_km2

  cat("2020 density raster detected\n")
  cat("Cell area:", round(cell_area_km2, 2), "km²\n")

} else {

  conversion_factor <- 1

  cat("2015 absolute cattle raster detected\n")
}

# Convert density to cattle per cell
r_cattle_total <- r_cattle * conversion_factor

# ------------------------------------------------------------------
# Convert raster cells to points
# ------------------------------------------------------------------

cat("Extracting raster cells...\n")

cattle_points <- rasterToPoints(r_cattle_total)

cattle_df <- data.frame(
  x = cattle_points[,1],
  y = cattle_points[,2],
  cattle = cattle_points[,3]
)

cattle_df <- cattle_df[!is.na(cattle_df$cattle), ]

cat("Valid cattle cells:", nrow(cattle_df), "\n")

# Convert to sf
cattle_sf <- st_as_sf(
  cattle_df,
  coords = c("x", "y"),
  crs = crs(r_cattle)
)

# Match CRS
cattle_sf <- st_transform(
  cattle_sf,
  st_crs(states_sf)
)

# ------------------------------------------------------------------
# Assign raster cells to states
# ------------------------------------------------------------------

cat("Assigning cattle cells to states...\n")

join_result <- st_join(
  cattle_sf,
  states_sf,
  join = st_within
)

# Remove cells outside Nigeria
join_result <- join_result[!is.na(join_result$NAME_1), ]

# ------------------------------------------------------------------
# Aggregate raster cattle by state
# ------------------------------------------------------------------

cat("Calculating raster-derived cattle totals...\n")

raster_totals <- join_result %>%
  st_drop_geometry() %>%
  group_by(NAME_1) %>%
  summarise(
    raster_cattle = sum(cattle, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(state = NAME_1)

# ------------------------------------------------------------------
# Load official cattle numbers
# ------------------------------------------------------------------

cat("Loading official cattle counts...\n")

official_totals <- read_excel(
  "Data/Cattle_numbers_State.xlsx"
)

official_totals <- official_totals %>%
  dplyr::select(
    state = NAME,
    official_cattle = N0
  ) %>%
  dplyr::mutate(
    state = dplyr::case_when(
      state == "Nassarawa" ~ "Nasarawa",
      state == "Abuja" ~ "Federal Capital Territory",
      TRUE ~ state
    )
  )
# ------------------------------------------------------------------
# Calculate correction factors
# ------------------------------------------------------------------

correction_table <- raster_totals %>%
  left_join(
    official_totals,
    by = "state"
  ) %>%
  mutate(
    correction_factor =
      official_cattle / raster_cattle
  ) %>%
  arrange(desc(correction_factor))

# ------------------------------------------------------------------
# Add national total
# ------------------------------------------------------------------

national_total <- data.frame(
  state = "Total",
  raster_cattle = sum(correction_table$raster_cattle, na.rm = TRUE),
  official_cattle = sum(correction_table$official_cattle, na.rm = TRUE),
  correction_factor =
    sum(correction_table$official_cattle, na.rm = TRUE) /
    sum(correction_table$raster_cattle, na.rm = TRUE)
)

correction_table <- bind_rows(
  correction_table,
  national_total
)

# ------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------

write.csv(
  correction_table,
  "Code/Prevalence/Bovine BCT and PCR/Analysis_NGA/Cattle_at_risk_correction_factors.csv",
  row.names = FALSE
)

cat("\nCorrection factors saved\n")

print(correction_table)

cat("\nSummary:\n")
cat("National raster cattle:",
    round(sum(correction_table$raster_cattle[correction_table$state != 'Total'])),
    "\n")

cat("National official cattle:",
    round(sum(correction_table$official_cattle[correction_table$state != 'Total'])),
    "\n")

cat("National correction factor:",
    round(
      sum(correction_table$official_cattle[correction_table$state != 'Total']) /
      sum(correction_table$raster_cattle[correction_table$state != 'Total']),
      3
    ),
    "\n")