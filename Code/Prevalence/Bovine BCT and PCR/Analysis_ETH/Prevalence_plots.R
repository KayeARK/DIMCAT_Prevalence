# Faster version of the supplied ETH cattle-at-risk post-processing script.
# Preserves the FINAL CSV calculation in the original (which wrote each file twice).
# Preserves coordinate swapping, CI-width > 0.8 exclusion, raster-cell averaging,
# multipoint/region membership, missing cattle -> 0, rounding and Total-row rules.
# No model fitting, additional checkpoint files or changes to raster coverage.
# Run from the same project working directory as the original.

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(raster)
})

# The supplied attachment is ETH, not NGA.
projection_dir <- "Code/Prevalence/Bovine BCT and PCR/Projections_ETH/"
output_dir <- "Code/Prevalence/Bovine BCT and PCR/Analysis_ETH/Cattle_at_risk/"
cattle_file <- "Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif"

# Cache only geometry-dependent results in memory. Values are recalculated for
# every file. A changed retained subset of a cattle cell receives its own entry.
make_spatial_lookup <- function(r_cattle, country_sf) {
  country_sf <- sf::st_make_valid(country_sf)
  point_keys <- character()
  point_geometry <- sf::st_sfc(crs = sf::st_crs(raster::crs(r_cattle)))
  point_cells <- numeric()
  spatial_cache <- new.env(hash = TRUE, parent = emptyenv())
  cache_entries <- 0L
  coord_key <- function(x, y) paste(sprintf("%.17g", x), sprintf("%.17g", y), sep = ":")

  function(d) {
    if (!nrow(d)) return(data.frame(value=numeric(), region=character(), cattle=numeric()))
    if (any(!is.finite(d$Longitude)) || any(!is.finite(d$Latitude)))
      stop("Non-finite coordinates: the original spatial conversion also requires valid coordinates.")
    keys <- coord_key(d$Longitude, d$Latitude)
    new_rows <- which(!duplicated(keys) & !keys %in% point_keys)
    if (length(new_rows)) {
      pts <- sf::st_as_sf(d[new_rows, c("Longitude", "Latitude"), drop=FALSE],
                         coords=c("Longitude", "Latitude"), crs=4326)
      pts <- sf::st_transform(pts, crs=raster::crs(r_cattle))
      point_cells <<- c(point_cells, raster::cellFromXY(r_cattle, sf::st_coordinates(pts)))
      point_geometry <<- c(point_geometry, sf::st_geometry(pts))
      point_keys <<- c(point_keys, keys[new_rows])
    }
    ids <- match(keys, point_keys)
    cells <- point_cells[ids]
    centers <- raster::xyFromCell(r_cattle, cells)
    # Same grouping and means as the original, including duplicate input rows.
    groups <- data.frame(center_lon=centers[,1], center_lat=centers[,2],
                         point_id=ids, value=d$value) %>%
      group_by(center_lon, center_lat) %>%
      summarise(value=mean(value, na.rm=TRUE),
                point_ids=list(sort(unique(point_id))), .groups="drop")
    group_keys <- vapply(groups$point_ids, paste, collapse=",", FUN.VALUE=character(1))
    missing <- which(!vapply(group_keys, exists, logical(1), envir=spatial_cache, inherits=FALSE))
    if (length(missing)) {
      message("  Preparing ", length(missing), " new retained cell geometries...")
      # Same sf union -> transform -> within join -> POINT cast as the original.
      # Joining the MULTIPOINT (before casting) preserves its boundary behaviour.
      lens <- lengths(groups$point_ids[missing])
      idx <- unlist(groups$point_ids[missing], use.names=FALSE)
      to_group <- sf::st_sf(cache_key=rep(group_keys[missing], lens),
                           geometry=point_geometry[idx])
      grouped <- to_group %>% group_by(cache_key) %>% summarise(.groups="drop")
      grouped <- sf::st_transform(grouped, crs=sf::st_crs(country_sf))
      joined <- sf::st_join(grouped, country_sf[, "NAME_1", drop=FALSE], join=sf::st_within)
      expanded <- suppressWarnings(sf::st_cast(joined, "POINT"))
      expanded$cattle <- as.numeric(raster::extract(r_cattle, expanded))
      expanded$cattle[is.na(expanded$cattle)] <- 0
      tables <- split(sf::st_drop_geometry(expanded)[, c("NAME_1", "cattle")], expanded$cache_key)
      for (key in group_keys[missing]) {
        z <- tables[[key]]
        if (is.null(z)) stop("A retained geometry produced no points: ", key)
        names(z)[1] <- "region"
        assign(key, z, envir=spatial_cache)
      }
      cache_entries <<- cache_entries + length(missing)
    }
    pieces <- mget(group_keys, envir=spatial_cache, inherits=FALSE)
    counts <- vapply(pieces, nrow, integer(1))
    out <- do.call(rbind, unname(pieces))
    out$value <- rep(groups$value, counts)
    out$cattle <- out$value * out$cattle
    rownames(out) <- NULL
    # Bound memory across unusually variable masks; never write disk checkpoints.
    if (cache_entries > 50000L) {
      rm(list=ls(spatial_cache, all.names=TRUE), envir=spatial_cache)
      cache_entries <<- 0L
      message("  Cleared the in-memory geometry cache to bound memory use.")
    }
    out[, c("value", "region", "cattle")]
  }
}

prepare_draw <- function(dpm, label="") {
  needed <- c("variable", "value", "Longitude", "Latitude")
  if (!all(needed %in% names(dpm))) stop("Missing columns in ", label, ": ", paste(setdiff(needed,names(dpm)),collapse=", "))
  slices <- lapply(c("Mean", "2.5th percentile", "97.5th percentile"), function(v)
    dpm[!is.na(dpm$variable) & dpm$variable == v & !is.na(dpm$value), needed, drop=FALSE])
  # The original subtracted/filtered these by row position. Fail clearly if that
  # assumption is not true instead of recycling or silently pairing locations.
  if (length(unique(vapply(slices,nrow,integer(1)))) != 1L)
    stop(label, ": mean/lower/upper have different non-missing lengths; cannot preserve the original row-paired filter.")
  for (j in 2:3) if (!identical(as.matrix(slices[[1]][,c("Longitude","Latitude")]),
                               as.matrix(slices[[j]][,c("Longitude","Latitude")]))) {
    # Matrix row names can differ merely because variables occupy different rows.
    if (!isTRUE(all.equal(unname(as.matrix(slices[[1]][,c("Longitude","Latitude")])),
                          unname(as.matrix(slices[[j]][,c("Longitude","Latitude")])))))
      stop(label, ": mean/lower/upper coordinates are not aligned. No result was written.")
  }
  dm <- slices[[1]]
  lat_in_range <- all(dm$Latitude >= 3 & dm$Latitude <= 15, na.rm=TRUE)
  lon_in_range <- all(dm$Longitude >= 33 & dm$Longitude <= 48, na.rm=TRUE)
  if (!lat_in_range || !lon_in_range) {
    message("  Swapping coordinates for Ethiopia ", label)
    slices <- lapply(slices, function(x) {tmp<-x$Longitude;x$Longitude<-x$Latitude;x$Latitude<-tmp;x})
  }
  high_uncertainty <- (slices[[3]]$value - slices[[2]]$value) > 0.8
  if (anyNA(high_uncertainty)) stop(label, ": undefined interval widths.")
  message("  Retained ", sum(!high_uncertainty), " of ", length(high_uncertainty), " points (CI width <= 0.8).")
  lapply(slices, function(x) x[!high_uncertainty,,drop=FALSE])
}

summarise_draw <- function(slices, lookup) {
  summaries <- lapply(seq_along(slices), function(j) {
    x <- lookup(slices[[j]]) %>% group_by(region) %>%
      summarise(value=mean(value,na.rm=TRUE), cattle=sum(cattle,na.rm=TRUE), .groups="drop")
    x <- as.data.frame(x[!is.na(x$region),,drop=FALSE])
    names(x)[2:3] <- list(c("mean_value","cattle_mean"),
                        c("lower_value","cattle_lower"),
                        c("upper_value","cattle_upper"))[[j]]
    x
  })
  result <- merge(summaries[[1]], summaries[[2]], by="region")
  result <- merge(result, summaries[[3]], by="region")
  livestock_cols <- grep("cattle",names(result),value=TRUE)
  result[livestock_cols] <- lapply(result[livestock_cols], function(x) round(as.numeric(x),0))
  # Preserve the FINAL original output: mean of region means and sum of rounded
  # region cattle totals (the overwritten preliminary block used different rules).
  total <- data.frame(region="Total", mean_value=mean(result$mean_value,na.rm=TRUE),
                      cattle_mean=sum(result$cattle_mean,na.rm=TRUE),
                      lower_value=mean(result$lower_value,na.rm=TRUE),
                      cattle_lower=sum(result$cattle_lower,na.rm=TRUE),
                      upper_value=mean(result$upper_value,na.rm=TRUE),
                      cattle_upper=sum(result$cattle_upper,na.rm=TRUE))
  rbind(result,total)
}

run_prevalence_plots <- function() {
  files <- list.files(projection_dir, pattern="^Projections_model_[0-9]+\\.csv$", full.names=TRUE)
  if (!length(files)) stop("No projection files found in ", projection_dir)
  ids <- as.numeric(sub("^Projections_model_([0-9]+)\\.csv$", "\\1", basename(files)))
  files <- files[order(ids)]
  if (!file.exists(cattle_file)) stop("Cattle raster not found: ", cattle_file)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  message("Found ", length(files), " projection files for Ethiopia.")
  message("Loading Ethiopian regions and cattle raster once...")
  regions <- geodata::gadm(country="ETH", level=1, path=tempdir())
  country_sf <- sf::st_as_sf(regions)
  r_cattle <- raster::raster(cattle_file)
  lookup <- make_spatial_lookup(r_cattle, country_sf)
  started <- proc.time()[[3]]
  for (i in seq_along(files)) {
    tick <- proc.time()[[3]]
    message("Processing ", basename(files[i]), " (", i, "/", length(files), ")")
    # Base read.csv retains the original parsing and numeric precision.
    dpm <- read.csv(files[i])
    result <- summarise_draw(prepare_draw(dpm, basename(files[i])), lookup)
    write.csv(result, file.path(output_dir, basename(files[i])), row.names=FALSE)
    elapsed <- proc.time()[[3]]-tick
    message(sprintf("  Saved in %.2f s; total elapsed %.1f min.", elapsed, (proc.time()[[3]]-started)/60))
  }
  message("Completed. Results: ", normalizePath(output_dir, winslash="/"))
  invisible(NULL)
}

run_prevalence_plots()
