# ETH: randomly sample complete covariate sets and fit new adjusted-count draws.
# Put this script alongside Results_ETH/ and Projections_ETH/.
# Projections_ETH at that level is understood to contain EXPECTED-prevalence CSVs.
# Run with source("AAT_ETH_random_covariate_projections.R").
# No covariate search. The model engine below is adapted from your attached ETH script.
# Retains quadratic groups, location IID, priors, mesh settings and swapped CSV labels.
# Existing CSVs are never overwritten by different results.
#
# First run: build a frozen empirical catalogue and assignment plan, then fit ONE new draw.
# Change max_new_draws to 20L, 100L or Inf to process more. Resourcing resumes automatically.
# draw_ids uses MODEL DRAW numbers; input CSV column = draw + 7.
# Existing expected CSVs in either input folder identify draws already produced.
# Published results: Results_ETH/Expected_prevalence_ETH/Projections_model_<draw>.csv
# Checkpoints/diagnostics: Results_ETH/Random_covariate_run_ETH/.
# Copying new generated covariate files back into the original donor catalogue is unnecessary.
#
# Random reassignment approximates full per-draw covariate selection; frequencies are not
# Bayesian model probabilities. Validate this shortcut on a subset before manuscript use.

# Additional helpers for the random-covariate continuation workflow.
random_script_directory <- function() {
  ofiles <- lapply(sys.frames(), function(x) x$ofile)
  ofiles <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x), ofiles)
  if (length(ofiles)) return(dirname(normalizePath(tail(ofiles, 1)[[1]], mustWork = TRUE)))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)))
  getwd()
}


random_settings <- list(
  base_dir = random_script_directory(), # script folder containing Results_ETH and Projections_ETH
  project_root = "",                   # auto-find ancestor containing Code/TestSensSpec; override if needed
  data_path = "",                      # blank uses the original Code/TestSensSpec input
  preparation_run_dir = "",            # optional original run containing Checkpoints_ETH; otherwise auto-detect
  total_draws = 1000L,                 # model draws 1:1000 correspond to input columns 8:1007
  draw_ids = 842:1000,                     # NULL = unfinished draws from the frozen plan; e.g. 101:120
  max_new_draws = Inf,                  # FIRST RUN: one new draw. Set 20L, 100L or Inf afterwards
  selection_seed = 20261006L,
  retain_fits = FALSE,                  # one fit checkpoint per draw for later joint aggregation
  retain_joint_samples = FALSE         # large; fits allow resampling if needed
)

random_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f))
  saveRDS(x, f, compress = FALSE, version = 2)
  unname(tools::md5sum(f))
}

random_child_dir <- function(parent, name) {
  ds <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  hits <- ds[tolower(basename(ds)) == tolower(name)]
  if (length(hits) > 1L) stop("Ambiguous directory names: ", paste(hits, collapse = ", "))
  if (length(hits)) hits[1] else file.path(parent, name)
}

random_find_project <- function(base, override = "") {
  if (nzchar(override)) {
    root <- normalizePath(path.expand(override), mustWork = TRUE)
    return(root)
  }
  p <- normalizePath(base, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(p, "Code/TestSensSpec/adjusted_cases_all_sites.csv"))) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  stop("Cannot find the DIMCAT project root above the script. Set random_settings$project_root.")
}

random_file_id <- function(path, stem) {
  pattern <- paste0("^", stem, "([0-9]+)([[:space:]]*\\([0-9]+\\))?\\.csv$")
  if (!all(grepl(pattern, basename(path)))) stop("Unexpected numbered CSV filename.")
  as.integer(sub(pattern, "\\1", basename(path)))
}

random_validate_projection <- function(path) {
  z <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e)
    stop("Cannot read projection ", path, ": ", conditionMessage(e)))
  if (!all(c("Latitude", "Longitude", "value", "variable") %in% names(z)) || !nrow(z))
    stop("Incomplete projection CSV: ", path)
  vars <- c("Mean", "2.5th percentile", "97.5th percentile")
  if (!setequal(unique(z$variable), vars) || !is.numeric(z$value) ||
      any(!is.na(z$value) & (!is.finite(z$value) | z$value < 0 | z$value > 1)) ||
      !any(is.finite(z$value[z$variable == "Mean"]))) stop("Invalid projection values: ", path)
  if (!is.numeric(z$Latitude) || !is.numeric(z$Longitude) ||
      any(!is.finite(z$Latitude)) || any(!is.finite(z$Longitude))) stop("Invalid projection coordinates: ", path)
  keys <- paste(z$Latitude, z$Longitude, sep = "|")
  if (anyDuplicated(paste(keys, z$variable, sep = "|"))) stop("Duplicate grid rows: ", path)
  grids <- lapply(vars, function(v) sort(keys[z$variable == v]))
  if (!identical(grids[[1]], grids[[2]]) || !identical(grids[[1]], grids[[3]]))
    stop("Truncated or inconsistent projection grid: ", path)
  invisible(TRUE)
}

random_projection_inventory <- function(directories) {
  files <- unique(unlist(lapply(directories[dir.exists(directories)], function(d)
    list.files(d, pattern = "^Projections_model_[0-9]+([[:space:]]*\\([0-9]+\\))?\\.csv$",
               full.names = TRUE))))
  if (!length(files)) return(data.frame(draw = integer(), path = character(), md5 = character()))
  ids <- random_file_id(files, "Projections_model_")
  hashes <- unname(tools::md5sum(files))
  if (anyNA(hashes)) stop("Cannot fingerprint a projection CSV. Check OneDrive downloads.")
  for (id in unique(ids)) {
    at <- which(ids == id)
    if (length(unique(hashes[at])) != 1L) stop("Conflicting copies of expected-prevalence draw ", id,
      ". Keep the intended version before continuing: ", paste(files[at], collapse = "; "))
    random_validate_projection(files[at[1]])
  }
  ans <- data.frame(draw = ids, path = files, md5 = hashes)
  ans <- ans[!duplicated(ans$draw), , drop = FALSE]
  ans[order(ans$draw), , drop = FALSE]
}

random_read_catalogue <- function(directory, completed, candidates) {
  files <- list.files(directory,
    pattern = "^Covariates_model_[0-9]+([[:space:]]*\\([0-9]+\\))?\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No Covariates_model_*.csv files in ", directory)
  ids <- random_file_id(files, "Covariates_model_")
  sets <- lapply(files, function(f) {
    z <- read.csv(f, stringsAsFactors = FALSE)
    if (!"covariates" %in% names(z)) stop("Missing covariates column: ", f)
    cc <- trimws(as.character(z$covariates))
    if (anyNA(cc) || any(!nzchar(cc)) || anyDuplicated(cc) || !all(cc %in% candidates))
      stop("Invalid/unknown covariates in ", f, ". Use ETH group names, not quadratic coefficient names.")
    candidates[candidates %in% cc]  # canonical order; squared terms are reconstructed by the engine
  })
  labels <- vapply(sets, paste, collapse = " + ", FUN.VALUE = character(1))
  for (id in unique(ids)) {
    if (length(unique(labels[ids == id])) > 1L) stop("Different covariate sets share donor draw ", id,
      ". Resolve the conflicting files rather than giving one draw extra weight.")
  }
  keep <- !duplicated(ids) & ids %in% completed$draw
  omitted <- sort(unique(ids[!ids %in% completed$draw]))
  if (length(omitted)) message("Excluding covariate files without a complete expected projection: ",
                               paste(omitted, collapse = ", "))
  if (!any(keep)) stop("No covariate sets have matching completed expected-prevalence files.")
  ord <- which(keep)[order(ids[keep])]
  list(table = data.frame(donor_draw = ids[ord], source_file = normalizePath(files[ord]),
    source_md5 = unname(tools::md5sum(files[ord])), covariates = labels[ord]), sets = sets[ord])
}

random_make_plan <- function(catalogue, targets, seed) {
  if (!length(catalogue$sets)) stop("Empty donor catalogue.")
  # Per-target seed makes assignments independent of session length/resume order.
  donor_row <- vapply(targets, function(d) {
    set.seed(as.integer((as.double(seed) + 1000003 * d) %% 2147483646 + 1))
    sample.int(length(catalogue$sets), 1L)
  }, integer(1))
  data.frame(draw = as.integer(targets), input_column = as.integer(targets + 7L),
    donor_row = donor_row, donor_draw = catalogue$table$donor_draw[donor_row],
    covariates = catalogue$table$covariates[donor_row])
}

random_publish <- function(source, destination) {
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(destination)) {
    if (identical(unname(tools::md5sum(source)), unname(tools::md5sum(destination)))) return(invisible(destination))
    stop("Refusing to overwrite a different existing projection: ", destination)
  }
  tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(destination))
  on.exit(unlink(tmp))
  if (!file.copy(source, tmp)) stop("Cannot copy projection for publication.")
  if (!file.rename(tmp, destination)) stop("Cannot publish projection: ", destination)
  invisible(destination)
}

random_preparation_source <- function(cfg, signature) {
  siblings <- list.dirs(dirname(cfg$random_results_dir), recursive = FALSE, full.names = TRUE)
  siblings <- siblings[grepl("^Spatial_joint_ETH", basename(siblings))]
  roots <- unique(c(cfg$random_results_dir, list.dirs(cfg$random_results_dir, recursive = TRUE), siblings))
  if (nzchar(cfg$random_preparation_override))
    roots <- normalizePath(path.expand(cfg$random_preparation_override), mustWork = TRUE)
  roots <- roots[!startsWith(roots, cfg$random_output_root)]
  required <- c("input_md5", "country", "input_filter", "elevation_policy", "spatial_crs",
    "cutoff_m", "max_edge_m", "offset_m", "range_prior_m", "sigma_prior", "package_versions",
    "covariate_degree", "coefficient_prior_sd", "intercept_prior_sd", "use_location_iid", "iid_prior")
  ok <- vapply(roots, function(root) {
    fs <- file.path(root, c("run_manifest.rds", "Checkpoints_ETH/spatial_preparation.rds"))
    if (!all(file.exists(fs))) return(FALSE)
    old <- readRDS(fs[1])
    all(vapply(required, function(nm) identical(old[[nm]], signature[[nm]]), logical(1)))
  }, logical(1))
  roots <- roots[ok]
  if (!length(roots)) {
    if (nzchar(cfg$random_preparation_override)) stop("preparation_run_dir is incomplete or incompatible with the current ETH data/model.")
    return("")
  }
  stamps <- file.info(file.path(roots, "Checkpoints_ETH/spatial_preparation.rds"))$mtime
  roots[order(stamps, decreasing = TRUE)][1]
}

run_random_eth <- function(cfg, random_cfg) {
  base <- normalizePath(path.expand(random_cfg$base_dir), mustWork = TRUE)
  results <- random_child_dir(base, "Results_ETH")
  external <- random_child_dir(base, "Projections_ETH")
  if (!dir.exists(results)) stop("Results_ETH must sit next to this script: ", results)
  cov_dir <- random_child_dir(results, "Covariates_ETH")
  expected <- random_child_dir(results, "Expected_prevalence_ETH")
  output <- file.path(results, "Random_covariate_run_ETH")
  dir.create(output, showWarnings = FALSE)
  lock <- file.path(results, ".random_projection_lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Another random-projection session is running. ",
    "After a killed session, remove .random_projection_lock only once it has stopped.")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  project <- random_find_project(base, random_cfg$project_root)
  oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE); setwd(project)
  if (nzchar(random_cfg$data_path)) cfg$data_path <- path.expand(random_cfg$data_path)
  cfg$data_path <- normalizePath(cfg$data_path, mustWork = TRUE)
  header <- names(read.csv(cfg$data_path, nrows = 0, check.names = FALSE))
  if (length(header) < 8L) stop("No adjusted-count columns found.")
  candidates <- c("elevation", "precipitation", "human_fp", "pop_den", "tree", "grassland",
    "shrub", "cropland", "built", "bare", "water", "wetland", "cattle", "tsetse")
  snapshot_file <- file.path(output, "Random_assignment_snapshot.rds")
  context <- list(input_md5 = unname(tools::md5sum(cfg$data_path)),
    total_draws = random_cfg$total_draws, seed = random_cfg$selection_seed,
    posterior_samples = cfg$posterior_samples, integration_tolerance = cfg$integration_tolerance,
    posterior_seed = cfg$seed)
  if (length(context$total_draws) != 1L || !is.finite(context$total_draws) ||
      context$total_draws < 1 || context$total_draws != as.integer(context$total_draws) ||
      context$total_draws > length(header) - 7L) stop("total_draws exceeds available count columns or is invalid.")
  if (length(context$seed) != 1L || !is.finite(context$seed) || context$seed < 1 ||
      context$seed != as.integer(context$seed)) stop("Invalid selection_seed.")
  if (file.exists(snapshot_file)) {
    snapshot <- readRDS(snapshot_file)
    if (!identical(snapshot$context, context)) stop("Input/seed/total_draws differ from the saved assignment plan. ",
      "Use the original settings to resume this run.")
  } else {
    completed <- random_projection_inventory(c(external, expected))
    catalogue <- random_read_catalogue(cov_dir, completed, candidates)
    targets <- setdiff(seq_len(context$total_draws), completed$draw)
    snapshot <- list(context = context, catalogue = catalogue,
      existing = completed, plan = random_make_plan(catalogue, targets, context$seed))
    atomic_save_rds(snapshot, snapshot_file)
    atomic_write_csv(catalogue$table, file.path(output, "Covariate_catalogue.csv"))
    atomic_write_csv(snapshot$plan, file.path(output, "Random_covariate_assignments.csv"))
    atomic_write_csv(completed, file.path(output, "Existing_expected_projections.csv"))
  }
  requested <- if (is.null(random_cfg$draw_ids)) snapshot$plan$draw else random_cfg$draw_ids
  if (!is.numeric(requested) || anyNA(requested) || any(!is.finite(requested)) ||
      any(requested != as.integer(requested)) || anyDuplicated(requested) ||
      !all(requested %in% snapshot$plan$draw)) stop("draw_ids must be new MODEL DRAW numbers from Random_covariate_assignments.csv.")
  requested <- as.integer(requested)
  complete <- vapply(requested, function(d) {
    marker <- file.path(output, "Checkpoints_ETH", sprintf("draw_%04d", d), "completed.rds")
    if (!file.exists(marker)) return(FALSE)
    z <- readRDS(marker)
    all(file.exists(z$files))
  }, logical(1))
  requested <- requested[!complete]
  limit <- random_cfg$max_new_draws
  if (length(limit) != 1L || is.na(limit) || limit <= 0 ||
      (is.finite(limit) && limit != as.integer(limit))) stop("max_new_draws must be a positive integer or Inf.")
  if (is.finite(limit)) requested <- head(requested, as.integer(limit))
  if (!length(requested)) {
    message("No unfinished requested draws remain. Existing outputs are in ", expected)
    return(invisible(output))
  }
  # Check collisions before any model fitting, including files created since the snapshot.
  for (d in requested) for (folder in c(external, expected)) {
    dest <- file.path(folder, paste0("Projections_model_", d, ".csv"))
    own <- file.path(output, "Expected_prevalence_ETH", basename(dest))
    if (file.exists(dest) && (!file.exists(own) ||
        !identical(unname(tools::md5sum(dest)), unname(tools::md5sum(own)))))
      stop("A new conflicting projection appeared after the plan was saved: ", dest)
  }
  cfg$iteration_columns <- requested + 7L
  cfg$random_snapshot <- snapshot
  cfg$random_catalogue_hash <- random_hash(snapshot$catalogue)
  cfg$random_results_dir <- normalizePath(results)
  cfg$random_output_root <- normalizePath(output)
  cfg$random_expected_dir <- expected
  cfg$random_preparation_override <- random_cfg$preparation_run_dir
  cfg$candidate_source_dir <- ""
  cfg$resume_dir <- if (file.exists(file.path(output, "run_manifest.rds"))) output else ""
  cfg$output_parent <- results
  cfg$save_fit_draws <- if (random_cfg$retain_fits) snapshot$plan$draw else integer(0)
  cfg$save_joint_draws <- if (random_cfg$retain_joint_samples) snapshot$plan$draw else integer(0)
  message("Frozen catalogue: ", nrow(snapshot$catalogue$table), " original completed draws; ",
    length(unique(snapshot$catalogue$table$covariates)), " distinct covariate sets.")
  message("This session will attempt model draws: ", paste(requested, collapse = ", "))
  message("Expected-prevalence CSVs will be published in: ", expected)
  run_aat_quick(cfg)
}
settings <- list(
  data_path = "Code/TestSensSpec/adjusted_cases_all_sites.csv",
  iteration_columns = 8L,                 # same draw 1 default; later 8:12, 8:1007
  mesh_profile = "pilot20km",             # pilot20km | check10km | reference6km
  output_parent = "Code/Prevalence/Bovine BCT and PCR",
  resume_dir = "",                        # NEW quadratic run; never resume old linear scores
  candidate_source_dir = "",              # fresh ETH search; never point at Nigerian runs
  reuse_preparation_dir = "",             # fresh ETH rasters/mesh; resume_dir reuses these
  reference_output_dir = "",              # optional existing output comparison
  coefficient_prior_sd = 2.5,
  covariate_degree = 2L,                 # linear + square; 1L permits a separate linear run
  waic_min_improvement = 2,              # add a group only if WAIC decreases MORE than this
  waic_removal_tolerance = 2,            # remove a group if WAIC increases LESS than this
  intercept_prior_sd = 5,
  use_location_iid = TRUE,
  iid_prior = c(1, 0.1),                  # unchanged P(sigma_IID > 1) = 0.1
  support_tolerance = 1e-10,
  workers = 1L,                          # capped at physical cores minus one
  inla_mode = NULL,                      # inherit current R session; no forced switch
  max_session_minutes = 30000,
  fit_timeout_seconds = 30000,
  sampling_timeout_seconds = 30000,
  prediction_timeout_seconds = 30000,
  control_inla = list(tolerance = 1e-5, tolerance.step = 1e-8),
  waic_recheck_tolerance = 0.05,
  posterior_samples = 4096L,              # per SELECTED model, not count draws
  prediction_chunk_size = 128L,           # limits temporary matrix memory
  integration_tolerance = 1e-4,           # probability units (0.01 percentage point)
  mc_screen_pp = 3,                       # split-half difference screen, in pp
  seed = 20260928L,
  diagnostic_draws = 1L,
  make_diagnostic_plots = TRUE,
  save_fit_draws = 1L,
  save_joint_draws = 1L                   # large; retain draw 1 for audit/reanalysis
)

atomic_save_rds <- function(object, path) {
  tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Could not commit checkpoint: ", path)
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  write.csv(object, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not commit CSV: ", path)
  invisible(path)
}

# Retain the original Ethiopia input orientation check, before aggregation.
# This is independent of the legacy swapped CSV output column labels.
filter_ethiopia_inputs <- function(data) {
  data$longitude <- suppressWarnings(as.numeric(as.character(data$longitude)))
  data$latitude <- suppressWarnings(as.numeric(as.character(data$latitude)))
  keep <- is.finite(data$longitude) & is.finite(data$latitude)
  if (any(!keep)) message("Removing ", sum(!keep), " rows with missing/non-finite coordinates.")
  data <- data[keep, , drop = FALSE]
  lat_ok <- data$latitude >= 3 & data$latitude <= 15
  lon_ok <- data$longitude >= 33 & data$longitude <= 48
  lat_as_lon <- data$latitude >= 33 & data$latitude <= 48
  lon_as_lat <- data$longitude >= 3 & data$longitude <= 15
  if (sum(lat_as_lon) > sum(lon_ok) && sum(lon_as_lat) > sum(lat_ok)) {
    message("Applying the original Ethiopia input coordinate-swap heuristic.")
    tmp <- data$longitude
    data$longitude <- data$latitude
    data$latitude <- tmp
  }
  keep <- data$latitude >= 3 & data$latitude <= 15 &
          data$longitude >= 33 & data$longitude <= 48
  message("Ethiopia bounding-box filter: retained ", sum(keep), " of ", nrow(data),
          " finite-coordinate records. Country-domain masking follows during preparation.")
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) stop("No observations remain in the Ethiopia input bounding box (longitude 33-48, latitude 3-15). Check the input CSV and coordinate orientation.")
  data
}

# Extract the sole elevation layer by position, not a mask-dependent layer name.
# Explicit point CRS also avoids assumptions about the raster's native CRS.
extract_elevation_values <- function(r, lonlat) {
  if (!inherits(r, "SpatRaster") || terra::nlyr(r) != 1L)
    stop("Elevation download did not return a single-layer SpatRaster. Check Data/Covariates/elevation.")
  raster_crs <- terra::crs(r)
  if (!nzchar(raster_crs)) stop("Elevation raster has no CRS.")
  lonlat <- as.matrix(lonlat)
  if (ncol(lonlat) != 2L || any(!is.finite(lonlat))) stop("Invalid elevation extraction coordinates.")
  pts <- terra::vect(lonlat, crs = "OGC:CRS84")
  pts <- terra::project(pts, raster_crs)
  z <- terra::extract(r, pts, ID = FALSE, method = "simple")
  if (nrow(z) != nrow(lonlat) || ncol(z) != 1L)
    stop("Elevation extraction changed row count or returned unexpected columns.")
  z <- as.numeric(z[[1L]])
  z[!is.finite(z)] <- NA_real_
  z
}

check_elevation_sites <- function(values, lonlat, sample_size, mapping_rows, diagnostics_dir) {
  if (length(values) != nrow(lonlat) || length(sample_size) != nrow(lonlat) ||
      length(mapping_rows) != nrow(lonlat)) stop("Elevation audit alignment failure.")
  missing <- !is.finite(values)
  audit <- data.frame(mapping_row = mapping_rows, longitude_deg = lonlat[, 1],
    latitude_deg = lonlat[, 2], sample_size = sample_size, elevation_m = values,
    missing_elevation = missing, source = "ETH_unmasked_30s")
  path <- file.path(diagnostics_dir, "Elevation_observation_extraction.csv")
  atomic_write_csv(audit, path)
  if (any(missing)) stop(sum(missing), " fitted locations still have missing elevation in the unmasked raster. ",
    "No mesh or model has been fitted. See ", path,
    "; inspect these coordinates and the downloaded raster before proceeding. No sites were dropped or imputed.")
  invisible(audit)
}

coordinate_key <- function(longitude, latitude) {
  paste(sprintf("%.17g", longitude), sprintf("%.17g", latitude), sep = "|")
}

aggregate_draw_counts <- function(data_original, iteration, mapping) {
  y <- data_original[[iteration]]
  if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0) ||
      any(y > data_original$sample_size) || any(abs(y - round(y)) > 1e-8)) {
    stop("Invalid binomial counts in input column ", iteration)
  }
  group <- match(coordinate_key(data_original$longitude, data_original$latitude),
                 coordinate_key(mapping$longitude, mapping$latitude))
  if (anyNA(group)) stop("An input row could not be matched to its coordinates.")
  totals <- rowsum(matrix(y, ncol = 1), group = group, reorder = TRUE)
  aligned <- totals[match(as.character(seq_len(nrow(mapping))), rownames(totals)), 1]
  if (anyNA(aligned)) stop("Aggregated response has missing locations.")
  as.numeric(aligned)
}

# Pure search routine; score() handles fitting/checkpoint reuse for one draw.
# Two-point practical parsimony rule; this is NOT a WAIC significance test.
select_per_draw <- function(candidates, score, note = function(...) NULL,
    score_many = function(sets) vapply(sets, score, numeric(1)),
    backward_tolerance = 2, addition_min_improvement = 2, final_backward_tolerance = 2) {
  covs <- as.character(candidates)
  backward <- function(covs, threshold, phase) {
    current <- score(covs)
    repeat {
      if (!length(covs)) return(list(covariates = covs, waic = current, delta = numeric(0)))
      values <- score_many(lapply(seq_along(covs), function(j) covs[-j]))
      delta <- values - current
      j <- which.min(delta)
      if (delta[j] >= threshold) {
        return(list(covariates = covs, waic = current, delta = delta))
      }
      note(phase, "remove", covs[j], delta[j])
      covs <- covs[-j]
      current <- values[j]
    }
  }
  first <- backward(covs, backward_tolerance, paste0("backward_", backward_tolerance))
  covs <- first$covariates
  current <- first$waic
  remaining <- setdiff(candidates, covs)
  while (length(remaining)) {
    values <- score_many(lapply(seq_along(remaining), function(j) c(covs, remaining[j])))
    delta <- values - current
    j <- which.min(delta)
    if (delta[j] >= -addition_min_improvement) break
    note("forward", "add", remaining[j], delta[j])
    covs <- c(covs, remaining[j])
    remaining <- remaining[-j]
    current <- values[j]
  }
  backward(covs, final_backward_tolerance, paste0("backward_", final_backward_tolerance))
}


# Compare original-format output files without interpreting swapped column names.
compare_run_outputs <- function(current_dir, reference_dir, draw, destination) {
  pred_name <- paste0("Projections_model_", draw, ".csv")
  cov_name <- paste0("Covariates_model_", draw, ".csv")
  files <- c(file.path(current_dir, "Projections_ETH", pred_name),
             file.path(reference_dir, "Projections_ETH", pred_name),
             file.path(current_dir, "Covariates_ETH", cov_name),
             file.path(reference_dir, "Covariates_ETH", cov_name))
  if (!all(file.exists(files))) return(invisible(FALSE))
  a <- read.csv(files[1]); b <- read.csv(files[2])
  coords <- c("Latitude", "Longitude", "variable")
  if (!identical(a[coords], b[coords])) stop("Reference output grid/order does not match.")
  ca <- as.character(read.csv(files[3])$covariates)
  cb <- as.character(read.csv(files[4])$covariates)
  ans <- do.call(rbind, lapply(unique(a$variable), function(v) {
    d <- a$value[a$variable == v] - b$value[b$variable == v]
    d <- d[is.finite(d)]
    data.frame(variable = v, comparable_points = length(d),
               mean_absolute_difference = if (length(d)) mean(abs(d)) else NA_real_,
               max_absolute_difference = if (length(d)) max(abs(d)) else NA_real_,
               root_mean_square_difference = if (length(d)) sqrt(mean(d^2)) else NA_real_,
               same_selected_covariates = setequal(ca, cb))
  }))
  atomic_write_csv(ans, destination)
  invisible(TRUE)
}

# Stack effect rows need not equal observation rows after INLA compression.
# Project effects through the stack A matrix to recover covariates in DATA order.
# The finite-weight check distinguishes structural NA padding from raster NA.
recover_covariates <- function(dat, AA, index, candidates) {
  ans <- lapply(candidates, function(nm) {
    v <- dat[[nm]]
    if (is.null(v) || length(v) != ncol(AA)) stop("Cannot recover cached covariate: ", nm)
    ok <- is.finite(v)
    weights <- as.numeric(AA %*% as.numeric(ok))[index]
    z <- as.numeric(AA %*% replace(v, !ok, 0))[index]
    z[abs(weights - 1) > 1e-8] <- NA_real_
    z
  })
  names(ans) <- candidates
  as.data.frame(ans, check.names = FALSE)
}

audit_covariates <- function(observed, predicted, tolerance) {
  fraction_names <- c("tree", "grassland", "shrub", "cropland", "built",
                      "bare", "water", "wetland", "mangrove", "tsetse")
  do.call(rbind, lapply(names(observed), function(nm) {
    x <- observed[[nm]]; p <- predicted[[nm]]
    v <- x[is.finite(x)]; q <- p[is.finite(p)]
    lo <- if (length(v)) min(v) else NA_real_
    hi <- if (length(v)) max(v) else NA_real_
    varying <- length(v) > 1L && hi - lo > tolerance * max(1, abs(lo), abs(hi))
    reason <- if (!length(v)) "all_missing_at_observations" else if (!varying)
      "constant_at_observations" else "eligible"
    unit_scale <- nm %in% fraction_names
    invalid_unit <- unit_scale && any(c(v, q) < -tolerance | c(v, q) > 1 + tolerance)
    s <- if (unit_scale || !varying) 1 else sd(v)
    data.frame(covariate = nm, eligible = varying, reason = reason,
      observed_n = length(x), observed_missing = sum(!is.finite(x)),
      observed_distinct = length(unique(v)), observed_nonzero = sum(abs(v) > tolerance),
      observed_min = lo, observed_max = hi,
      rare_support = length(v) > 0 && max(tabulate(match(v, unique(v)))) >= length(v) - 5L,
      prediction_missing = sum(!is.finite(p)),
      prediction_outside_observed = if (length(v)) sum(q < lo - tolerance | q > hi + tolerance) else NA_integer_,
      centre = if (length(v)) mean(v) else 0, scale = s,
      scale_description = if (unit_scale) "one unit (0 to 1)" else "one observed-location SD",
      invalid_fraction_units = invalid_unit)
  }))
}

transform_covariates <- function(x, audit) {
  keep <- audit$covariate[audit$eligible]
  out <- x[, keep, drop = FALSE]
  for (nm in keep) {
    a <- audit[match(nm, audit$covariate), ]
    out[[nm]] <- (out[[nm]] - a$centre) / a$scale
  }
  out
}

nearest_distances_km <- function(predicted_m, observed_m) {
  if (!nrow(observed_m) || any(!is.finite(observed_m)) || any(!is.finite(predicted_m)))
    stop("Invalid locations for distance diagnostics.")
  d2 <- rep(Inf, nrow(predicted_m))
  for (j in seq_len(nrow(observed_m))) {
    d2 <- pmin(d2, (predicted_m[, 1] - observed_m[j, 1])^2 +
                      (predicted_m[, 2] - observed_m[j, 2])^2)
  }
  sqrt(d2) / 1000
}

support_at_predictions <- function(prepared, covariates) {
  n <- nrow(prepared$prediction_lonlat)
  missing <- outside <- rep(0L, n)
  for (nm in covariates) {
    p <- prepared$raw_prediction_covariates[[nm]]
    a <- prepared$covariate_audit[match(nm, prepared$covariate_audit$covariate), ]
    missing <- missing + !is.finite(p)
    outside <- outside + (is.finite(p) & (p < a$observed_min - prepared$support_tolerance |
                                          p > a$observed_max + prepared$support_tolerance))
  }
  data.frame(missing_selected_covariates = missing, outside_selected_ranges = outside,
             valid_prediction = missing == 0L)
}

node_mesh_rows <- function(ids, n) {
  if (length(ids) != n || anyNA(ids) || anyDuplicated(ids)) stop("Unexpected spatial node IDs.")
  if (setequal(ids, 0:(n - 1L))) return(as.integer(ids) + 1L)
  if (setequal(ids, seq_len(n))) return(as.integer(ids))
  stop("Spatial IDs match neither zero-based nor one-based mesh indexing.")
}

fixed_prior_control <- function(slope_sd, intercept_sd) {
  # b0 is an explicitly named covariate, not R's implicit intercept.
  list(mean = 0, prec = list(b0 = 1 / intercept_sd^2, default = 1 / slope_sd^2))
}

# Never reuse an observed IID effect at a prediction cell, even if coordinates
# happen to coincide. The target is a new exchangeable location-level unit.
# No sum-to-zero constraint: adding prediction units must not change the prior
# on observed effects. NA in an INLA random-effect index means no contribution.
iid_indices <- function(n_observed, n_predicted) {
  stopifnot(n_observed > 0, n_predicted > 0)
  list(observed = seq_len(n_observed),
       new_location = n_observed + seq_len(n_predicted),
       surface = rep(NA_integer_, n_predicted))
}

model_formula <- function(covariates, spde, use_location_iid, iid_prior) {
  iid_hyper <- list(prec = list(prior = "pc.prec", param = iid_prior))
  terms <- c("0", "b0", covariates, "f(s, model = spde)")
  if (use_location_iid) terms <- c(terms,
    "f(site_iid, model = 'iid', constr = FALSE, hyper = iid_hyper)")
  as.formula(paste("y ~", paste(terms, collapse = " + ")), env = environment())
}

projection_table <- function(prevalence, lonlat) {
  ans <- rbind(
    data.frame(Latitude = lonlat[, 1], Longitude = lonlat[, 2],
               value = prevalence$mean, variable = "Mean"),
    data.frame(Latitude = lonlat[, 1], Longitude = lonlat[, 2],
               value = prevalence[["0.025quant"]], variable = "2.5th percentile"),
    data.frame(Latitude = lonlat[, 1], Longitude = lonlat[, 2],
               value = prevalence[["0.975quant"]], variable = "97.5th percentile"))
  ans$variable <- as.factor(ans$variable)
  ans
}

iid_sd_quantiles <- function(hyperparameters) {
  j <- which(rownames(hyperparameters) == "Precision for site_iid")
  cols <- c("0.025quant", "0.5quant", "0.975quant")
  if (length(j) != 1L || !all(cols %in% names(hyperparameters)))
    stop("Cannot identify the site_iid precision posterior. Inspect the saved INLA fit.")
  q <- as.numeric(hyperparameters[j, cols])
  if (any(!is.finite(q) | q <= 0) || any(diff(q) < 0)) stop("Invalid IID precision quantiles.")
  # The transformation is decreasing, so the two interval endpoints reverse.
  # Do not misreport 1/sqrt(mean precision) as posterior mean SD.
  data.frame(parameter = "IID logit-scale SD", lower = 1 / sqrt(q[3]),
             median = 1 / sqrt(q[2]), upper = 1 / sqrt(q[1]))
}

surface_result <- function(saved) {
  ans <- saved
  ans$prevalence <- saved$surface_prevalence
  ans$linear_predictor <- saved$surface_linear_predictor
  ans
}

uncertainty_rows <- function(saved, prepared) {
  pp <- saved$prevalence
  data.frame(Latitude = prepared$prediction_lonlat[, 1], Longitude = prepared$prediction_lonlat[, 2],
    mean = pp$mean, lower = pp[["0.025quant"]], upper = pp[["0.975quant"]],
    interval_width = pp[["0.975quant"]] - pp[["0.025quant"]],
    logit_sd = saved$linear_predictor$sd,
    nearest_observation_km = prepared$nearest_observation_km,
    support_at_predictions(prepared, saved$covariates))
}

uncertainty_summary <- function(d, draw, scenario = "primary") {
  groups <- list(all = rep(TRUE, nrow(d)), within_selected_ranges = d$outside_selected_ranges == 0,
                 outside_selected_ranges = d$outside_selected_ranges > 0)
  band <- cut(d$nearest_observation_km, c(-Inf, 10, 25, 50, 100, Inf),
              labels = c("0_10km", "10_25km", "25_50km", "50_100km", "over100km"))
  for (v in levels(band)) groups[[v]] <- band == v
  do.call(rbind, lapply(names(groups), function(nm) {
    all <- groups[[nm]]; ok <- all & is.finite(d$mean)
    data.frame(draw = draw, scenario = scenario, group = nm, grid_points = sum(all), valid_points = sum(ok),
      mean_prevalence_unweighted = if (any(ok)) mean(d$mean[ok]) else NA_real_,
      median_interval_width = if (any(ok)) median(d$interval_width[ok]) else NA_real_,
      fraction_width_above_0.8 = if (any(ok)) mean(d$interval_width[ok] > .8) else NA_real_,
      median_logit_sd = if (any(ok)) median(d$logit_sd[ok]) else NA_real_)
  }))
}


# Each parent covariate is one selection unit. Its squared term, when
# estimable, is always added/removed with its linear term (strong hierarchy).
# z uses the existing observation-only centring/scaling; the basis is z, z^2.
# Curvature coefficients have the same explicit N(0, slope_sd^2) prior.
audit_quadratic_basis <- function(observed, audit, degree = 2L) {
  z <- transform_covariates(observed, audit)
  if (!ncol(z)) return(data.frame(covariate = character(), linear_term = character(),
    quadratic_term = character(), degree = integer(), reason = character(),
    observed_distinct = integer(), quadratic_min = numeric(), quadratic_max = numeric()))
  do.call(rbind, lapply(names(z), function(nm) {
    v <- z[[nm]]; v <- v[is.finite(v)]
    distinct <- length(unique(v)); q <- v^2
    if (any(!is.finite(q))) stop("Non-finite squared observed covariate: ", nm)
    independent <- degree == 2L && distinct >= 3L &&
      qr(cbind(1, v, q), tol = 1e-10)$rank == 3L
    data.frame(covariate = nm, linear_term = nm,
      quadratic_term = if (independent) paste0(nm, "_sq") else NA_character_,
      degree = if (independent) 2L else 1L,
      reason = if (degree == 1L) "linear_requested" else if (distinct < 3L)
        "fewer_than_three_observed_values" else if (!independent)
        "square_collinear_with_intercept_and_linear" else "quadratic",
      observed_distinct = distinct, quadratic_min = min(q), quadratic_max = max(q))
  }))
}

effect_terms <- function(covariates, prepared) {
  if (!length(covariates)) return(character())
  basis <- prepared$quadratic_basis
  if (is.null(basis)) return(as.character(covariates))
  at <- match(covariates, basis$covariate)
  if (anyNA(at)) stop("Unknown parent covariate in quadratic term mapping.")
  unlist(lapply(at, function(j) c(basis$linear_term[j],
    if (basis$degree[j] == 2L) basis$quadratic_term[j])), use.names = FALSE)
}

transform_model_covariates <- function(x, prepared) {
  out <- transform_covariates(x, prepared$covariate_audit)
  basis <- prepared$quadratic_basis
  if (!is.null(basis)) for (j in which(basis$degree == 2L)) {
    nm <- basis$linear_term[j]; qnm <- basis$quadratic_term[j]
    if (qnm %in% names(out)) stop("Quadratic term name collides with a covariate: ", qnm)
    out[[qnm]] <- out[[nm]]^2
    if (any(is.finite(out[[nm]]) & !is.finite(out[[qnm]])))
      stop("Squaring overflow for covariate: ", nm)
  }
  out
}

# Candidate fits are identical across the old and new final-choice rules.
# Every other manifest field remains mandatory for checkpoint compatibility.
candidate_fit_signature <- function(signature) {
  signature$selection_rule <- NULL
  signature
}

resolve_candidate_source <- function(cfg, signature) {
  requested <- cfg$candidate_source_dir
  if (!nzchar(requested)) return("")
  automatic <- identical(requested, "auto")
  roots <- if (automatic) list.dirs(path.expand(cfg$output_parent), recursive = FALSE,
                                  full.names = TRUE) else path.expand(requested)
  if (automatic) roots <- roots[grepl("^Spatial_joint_ETH_degree", basename(roots))]
  draws <- cfg$iteration_columns - 7L
  compatible <- vapply(roots, function(root) {
    tryCatch({
      old <- readRDS(file.path(root, "run_manifest.rds"))
      # Automatic discovery uses original search runs, not derived exports.
      if (automatic && !is.null(old$selection_rule)) return(FALSE)
      if (!identical(candidate_fit_signature(old), candidate_fit_signature(signature))) return(FALSE)
      if (dir.exists(file.path(root, ".run_lock"))) return(FALSE)
      required <- c(file.path(root, "Checkpoints_ETH", "spatial_preparation.rds"),
        unlist(lapply(draws, function(draw) file.path(root, "Checkpoints_ETH",
          sprintf("draw_%04d", draw), c("candidate_scores.rds", "completed.rds")))))
      all(file.exists(required))
    }, error = function(e) FALSE)
  }, logical(1))
  roots <- roots[compatible]
  if (!length(roots)) stop("No completed, compatible candidate-source run found. Set candidate_source_dir to the previous quadratic run folder. No search has been started. For a genuinely NEW search, set candidate_source_dir = ''.")
  if (length(roots) > 1L) {
    stamp <- file.info(file.path(roots, "run_manifest.rds"))$mtime
    roots <- roots[order(stamp, roots, decreasing = TRUE)]
  }
  normalizePath(roots[1], mustWork = TRUE)
}

import_candidate_scores <- function(source, destination, context) {
  old <- readRDS(source)
  expected <- context
  expected$model_signature <- candidate_fit_signature(expected$model_signature)
  old_context <- old$context
  old_context$model_signature <- candidate_fit_signature(old_context$model_signature)
  if (!identical(old_context, expected)) stop("Imported candidate context mismatch: input, draw, response, model or candidate order differs.")
  if (!length(old$records)) stop("Candidate source has no completed fits.")
  # Validate contents before attaching the new final-selection rule.
  select_best_visited(old$records, context$candidates)
  atomic_save_rds(list(context = context, records = old$records,
    imported_from = normalizePath(source, mustWork = TRUE)), destination)
  invisible(length(old$records))
}

select_best_visited <- function(records, candidates) {
  if (!length(records) || is.null(names(records)) || anyDuplicated(names(records)))
    stop("Missing or duplicate candidate records.")
  key_for <- function(covs) paste0("m_", paste(match(covs, candidates), collapse = "_"))
  for (key in names(records)) {
    r <- records[[key]]
    if (length(r$waic) != 1L || !is.finite(r$waic) || anyDuplicated(r$covariates) ||
        !all(r$covariates %in% candidates) || !identical(key_for(r$covariates), key))
      stop("Invalid candidate record: ", key)
  }
  scores <- vapply(records, function(r) r$waic, numeric(1))
  size <- vapply(records, function(r) length(r$covariates), integer(1))
  # Lowest numerical WAIC wins. Exact ties: fewer groups, then stable key order.
  key <- names(records)[order(scores, size, names(records))[1]]
  best <- records[[key]]
  delta <- vapply(seq_along(best$covariates), function(j) {
    deleted <- records[[key_for(best$covariates[-j])]]
    if (is.null(deleted)) NA_real_ else deleted$waic - best$waic
  }, numeric(1))
  list(covariates = best$covariates, waic = best$waic, delta = delta,
       key = key, rule = "minimum_WAIC_among_visited", visited = length(records))
}

# A parent process enforces deadlines while INLA runs in killable child processes.
# setTimeLimit() alone does not reliably interrupt an external INLA executable.
elapsed_now <- function() unname(proc.time()[["elapsed"]])
pause_run <- function(message) {
  stop(structure(list(message = message, call = NULL), class = c("aat_pause", "error", "condition")))
}

# Copy INLA's own options, not unrelated R options or the user's environment.
# A clean child session otherwise loses e.g. a custom inla.call/solver setting.
capture_inla_runtime <- function() {
  # Older releases may print defaults or return an environment/NULL for a
  # bulk request. Use named getters only, and transfer supported non-NULL
  # runtime values. Missing/deprecated optional settings retain the defaults
  # of the SAME INLA version in the worker. Never copy shared fit directories.
  option_names <- c("inla.call", "inla.arg", "inla.mode", "smtp", "mkl",
    "blas.num.threads", "pardiso.license", "fmesher.call", "fmesher.arg",
    "internal.binary.mode", "internal.experimental.mode", "scale.model.default")
  opts <- list()
  for (nm in option_names) {
    warned <- FALSE
    value <- tryCatch(withCallingHandlers(INLA::inla.getOption(nm),
      warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
      error = function(e) NULL)
    if (!warned && !is.null(value)) opts[nm] <- list(value)
  }
  env_names <- c("PATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH")
  list(options = opts, environment = Sys.getenv(env_names, unset = NA_character_),
       architecture = R.version$arch, r_version = as.character(getRversion()),
       inla_version = as.character(packageVersion("INLA")))
}

restore_inla_runtime <- function(runtime) {
  if (!identical(R.version$arch, runtime$architecture) ||
      !identical(as.character(getRversion()), runtime$r_version))
    stop("Worker R architecture/version differs from the parent R session.")
  if (!identical(as.character(packageVersion("INLA")), runtime$inla_version))
    stop("Worker INLA version differs from the parent R session.")
  vars <- runtime$environment
  if (anyNA(vars)) Sys.unsetenv(names(vars)[is.na(vars)])
  if (any(!is.na(vars))) do.call(Sys.setenv, as.list(vars[!is.na(vars)]))
  do.call(INLA::inla.setOption, runtime$options)
  invisible(TRUE)
}

write_inla_failure <- function(e, cfg, covariates, work_dir = NULL) {
  destination <- file.path(cfg$diagnostics_dir, paste0("INLA_failure_worker_", Sys.getpid(), ".txt"))
  lines <- c("INLA FAILURE REPORT", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"),
    paste("R:", R.version.string), paste("Architecture:", R.version$arch),
    paste("INLA:", as.character(packageVersion("INLA"))),
    paste("INLA mode:", cfg$inla_mode), "INLA threads: 1:1",
    paste("Covariates:", paste(covariates, collapse = ", ")),
    paste("Worker library paths:", paste(.libPaths(), collapse = "; ")),
    "", "R ERROR:", conditionMessage(e), "", paste("Local fit directory:", work_dir))
  # Include complete available backend logs; no blanket environment dump.
  if (!is.null(work_dir) && dir.exists(work_dir)) {
    logs <- list.files(work_dir, pattern = "(Logfile|R_console|stderr|stdout).*\\.(txt|log)$",
                       recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    for (f in logs) {
      content <- tryCatch(readLines(f, warn = FALSE), error = function(x) conditionMessage(x))
      lines <- c(lines, "", paste("LOG:", f), content)
    }
  }
  writeLines(lines, destination)
  destination
}

make_worker_pool <- function(n, functions_file, preparation_file, cfg, deadline, log_message) {
  workers <- vector("list", n)
  closed <- FALSE
  close_pool <- function() {
    if (!closed) for (p in workers) if (!is.null(p)) {
      try(p$kill_tree(), silent = TRUE)
      try(p$close(grace = 0), silent = TRUE)
    }
    closed <<- TRUE
    invisible(NULL)
  }
  initialized <- FALSE
  on.exit(if (!initialized) close_pool(), add = TRUE)
  # No shell interpolation; paths and functions travel through callr's R API.
  for (i in seq_len(n)) {
    if (elapsed_now() >= deadline) pause_run("Session modelling budget reached during worker startup.")
    env <- c(callr::rcmd_safe_env(), cfg$runtime$environment, OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
    env <- env[!duplicated(names(env), fromLast = TRUE)]
    env <- env[!is.na(env)] # missing parent variables are unset during initialization
    opts <- callr::r_session_options(libpath = .libPaths(), env = env,
                                     user_profile = FALSE, system_profile = FALSE, supervise = TRUE)
    workers[[i]] <- callr::r_session$new(options = opts, wait = TRUE, wait_timeout = 10000)
  }
  run_jobs <- function(jobs, on_result = function(job, value, seconds) NULL) {
    if (!length(jobs)) return(invisible(list()))
    if (closed) stop("Worker pool is closed.")
    active <- vector("list", n)
    next_job <- 1L
    finished <- 0L
    values <- vector("list", length(jobs))
    last_status <- elapsed_now()
    ok <- FALSE
    on.exit(if (!ok) close_pool(), add = TRUE)
    repeat {
      # Reap completed jobs BEFORE applying deadlines, so completed scores survive.
      for (i in seq_len(n)) if (!is.null(active[[i]])) {
        p <- workers[[i]]
        ready <- p$poll_process(0)
        if (identical(ready, "ready")) {
          event <- p$read()
          if (!is.null(event)) {
            worker_log <- file.path(dirname(dirname(preparation_file)), "Diagnostics_ETH", "worker_output.log")
            lines <- c(event$stdout, event$stderr)
            lines <- lines[nzchar(lines)]
            if (length(lines) && dir.exists(dirname(worker_log)))
              cat(paste0("\n", active[[i]]$job$label, "\n"), lines, "\n", file = worker_log, append = TRUE)
          }
          if (!is.null(event) && !is.null(event$error)) stop(event$error)
          if (!is.null(event) && event$code >= 500) stop("Worker exited: ", active[[i]]$job$label)
          if (!is.null(event) && event$code == 200) {
            a <- active[[i]]
            if (is.list(event$result) && isTRUE(event$result$aat_worker_failure)) {
              report_path <- file.path(dirname(dirname(preparation_file)), "Diagnostics_ETH", "Failure_report.txt")
              source_path <- event$result$failure_report
              if (length(source_path) == 1L && file.exists(source_path))
                file.copy(source_path, report_path, overwrite = TRUE)
              if (!file.exists(report_path)) writeLines(event$result$message, report_path)
              if (length(lines)) cat("\nCALLR STDOUT/STDERR:\n", lines, "\n",
                                     file = report_path, append = TRUE)
              log_message("FAILED: ", a$job$label, ". Detailed report: ", report_path)
              stop(event$result$message, "\nRead Diagnostics_ETH/Failure_report.txt.", call. = FALSE)
            }
            seconds <- elapsed_now() - a$start
            on_result(a$job, event$result, seconds)
            values[[a$index]] <- event$result
            active[i] <- list(NULL)
            finished <- finished + 1L
          }
        } else if (!p$is_alive()) stop("Worker crashed: ", active[[i]]$job$label)
      }
      if (finished == length(jobs)) break
      now <- elapsed_now()
      if (now >= deadline) pause_run("Session modelling budget reached. Completed work is checkpointed.")
      for (a in active) if (!is.null(a) && now - a$start >= a$job$timeout) {
        pause_run(paste0(a$job$label, " exceeded its ", a$job$timeout,
                         " second limit. No failed model was scored or skipped."))
      }
      for (i in seq_len(n)) if (is.null(active[[i]]) && next_job <= length(jobs)) {
        j <- next_job
        workers[[i]]$call(jobs[[j]]$fun, args = jobs[[j]]$args)
        active[[i]] <- list(index = j, job = jobs[[j]], start = elapsed_now())
        next_job <- next_job + 1L
      }
      if (now - last_status >= 30) {
        labels <- vapply(Filter(Negate(is.null), active), function(a) a$job$label, character(1))
        log_message("RUNNING: ", paste(labels, collapse = "; "),
          " | session budget remaining ", round(max(0, deadline - now) / 60, 1), " min")
        last_status <- now
      }
      Sys.sleep(0.1)
    }
    ok <- TRUE
    values
  }
  initialize <- function(functions_file, preparation_file, cfg) {
    suppressPackageStartupMessages({library(INLA); library(Matrix); library(matrixStats)})
    funs <- readRDS(functions_file)
    for (nm in names(funs)) {
      environment(funs[[nm]]) <- .GlobalEnv
      assign(nm, funs[[nm]], envir = .GlobalEnv)
    }
    restore_inla_runtime(cfg$runtime)
    assign(".aat_prepared", readRDS(preparation_file), envir = .GlobalEnv)
    assign(".aat_cfg", cfg, envir = .GlobalEnv)
    TRUE
  }
  init_jobs <- lapply(seq_len(n), function(i) list(label = paste0("worker ", i, " initialization"),
    timeout = 60, fun = initialize,
    args = list(functions_file = functions_file, preparation_file = preparation_file, cfg = cfg)))
  run_jobs(init_jobs)
  initialized <- TRUE
  list(run = run_jobs, close = close_pool)
}

# One numerical configuration and one observed-data layout for ALL fits.
fit_observations <- function(prepared, y, covariates, cfg) {
  dat <- prepared$observation_data
  dat$y[prepared$observation_index] <- y
  trials <- dat[["numtrials"]]
  observed <- prepared$observation_index
  if (!is.numeric(trials) || length(trials) != length(dat$y))
    stop("Missing or misaligned binomial trial counts in observation_data$numtrials.")
  yy <- dat$y[observed]; nn <- trials[observed]
  if (length(y) != length(observed) || any(!is.finite(yy)) || any(!is.finite(nn)) ||
      any(nn <= 0 | yy < 0 | yy > nn) || any(abs(nn - round(nn)) > 1e-8) ||
      any(abs(yy - round(yy)) > 1e-8))
    stop("Invalid binomial data: require integer 0 <= positives <= sample size at every fitted location.")
  if (!is.null(prepared$numtrials) && !isTRUE(all.equal(as.numeric(nn),
      as.numeric(prepared$numtrials), tolerance = 0)))
    stop("Stack trial counts differ from the prepared sampled-location totals.")
  # INLA evaluates some arguments using non-standard evaluation. In this
  # workflow the formula has its own environment, where the caller's local
  # 'dat' does not exist. An inline Ntrials=dat$numtrials can consequently
  # become NULL and silently fall back to ONE trial (observed failure report).
  # Match the successful earlier script: resolve arguments first, then do.call.
  formula <- model_formula(effect_terms(covariates, prepared), prepared$spde,
                           prepared$use_location_iid, prepared$iid_prior)
  job_dir <- tempfile(pattern = paste0("fit_", Sys.getpid(), "_"), tmpdir = cfg$local_work_root)
  if (!dir.create(job_dir)) stop("Cannot create local INLA working directory.")
  console_file <- file.path(job_dir, "R_console.log")
  con <- file(console_file, open = "wt")
  original_sink <- sink.number()
  sink(con)
  successful <- FALSE
  on.exit({
    while (sink.number() > original_sink) sink()
    close(con)
    if (successful) unlink(job_dir, recursive = TRUE)
  }, add = TRUE)
  args <- list(formula = formula,
    data = dat, family = "binomial", Ntrials = as.numeric(trials),
    inla.mode = cfg$inla_mode, num.threads = "1:1",
    control.inla = cfg$control_inla,
    control.fixed = fixed_prior_control(cfg$coefficient_prior_sd, cfg$intercept_prior_sd),
    control.compute = list(waic = TRUE, config = TRUE, dic = FALSE,
                           return.marginals.predictor = FALSE),
    control.predictor = list(link = 1, compute = TRUE, A = prepared$observation_A),
    # Keep backend files until success; the working directory must not yet exist.
    working.directory = file.path(job_dir, "inla"), keep = TRUE,
    verbose = TRUE, silent = FALSE)
  res <- tryCatch(do.call(INLA::inla, args), error = function(e) {
      cat("\nINLA ERROR: ", conditionMessage(e), "\n", sep = "")
      flush(con)
      write_inla_failure(e, cfg, covariates, job_dir)
      stop(e)
    })
  if (length(res$waic$waic) != 1L || !is.finite(res$waic$waic)) stop("INLA returned non-finite WAIC.")
  successful <- TRUE
  res
}

worker_task_inner <- function(kind, covariates, y, path = NULL, expected_waic = NULL, draw = NULL) {
  p <- .aat_prepared
  cfg <- .aat_cfg
  if (kind %in% c("score", "selected_fit", "assigned_fit")) {
    res <- fit_observations(p, y, covariates, cfg)
    value <- res$waic$waic
    if (kind %in% c("selected_fit", "assigned_fit")) {
      if (kind == "selected_fit" && abs(value - expected_waic) > cfg$waic_recheck_tolerance)
        stop("Selected observed-data refit WAIC changed by ", signif(value - expected_waic, 6),
             ". Inspect numerical stability before exporting predictions.")
      atomic_save_rds(list(covariates = covariates, response = y, waic = value, fit = res), path)
      # Retain the checkpoint fit when requested; avoid a second large archive copy.
    }
    rm(res); gc(verbose = FALSE)
    return(list(waic = value))
  }
  if (kind == "sample") {
    obj <- readRDS(file.path(dirname(path), "selected_observation_fit.rds"))
    if (!identical(obj$covariates, covariates) || !identical(obj$response, y)) stop("Posterior fit mismatch.")
    joint <- sample_joint_posterior(obj$fit, p, covariates, cfg, draw)
    compact <- list(covariates = covariates, response = y, joint = joint,
      waic = obj$waic, fixed = obj$fit$summary.fixed, hyperparameters = obj$fit$summary.hyperpar,
      spatial_nodes = obj$fit$summary.random$s,
      observed_iid = if (p$use_location_iid) obj$fit$summary.random$site_iid else NULL,
      observed_prevalence = obj$fit$summary.fitted.values[p$observation_index, , drop = FALSE])
    atomic_save_rds(compact, path)
    return(list(samples = ncol(joint$beta)))
  }
  if (kind == "predict") {
    obj <- readRDS(file.path(dirname(path), "joint_posterior.rds"))
    if (!identical(obj$covariates, covariates) || !identical(obj$response, y)) stop("Joint posterior mismatch.")
    t0 <- elapsed_now()
    saved <- predict_joint(obj, p, cfg, draw)
    saved$elapsed_seconds <- elapsed_now() - t0
    atomic_save_rds(saved, path)
    return(list(seconds = saved$elapsed_seconds, mc_check = saved$mc_check))
  }
  stop("Unknown worker task: ", kind)
}

# Return a normal callr result on R-level failure, so the parent can retain the
# report before killing the pool. A failed model is never returned as a score.
worker_task <- function(kind, covariates, y, path = NULL, expected_waic = NULL, draw = NULL) {
  tryCatch(worker_task_inner(kind, covariates, y, path, expected_waic, draw), error = function(e) {
    report <- file.path(.aat_cfg$diagnostics_dir, paste0("INLA_failure_worker_", Sys.getpid(), ".txt"))
    if (!file.exists(report)) report <- write_inla_failure(e, .aat_cfg, covariates)
    list(aat_worker_failure = TRUE, message = conditionMessage(e), failure_report = report)
  })
}

make_parallel_scorer <- function(candidates, checkpoint_file, context, pool, y, cfg, report) {
  records <- list()
  if (file.exists(checkpoint_file)) {
    old <- readRDS(checkpoint_file)
    if (!identical(old$context, context)) stop("Candidate checkpoint context mismatch.")
    records <- old$records
  }
  hits <- new_fits <- 0L
  key_for <- function(covs) {
    if (anyDuplicated(covs) || !all(covs %in% candidates)) stop("Invalid candidate set.")
    paste0("m_", paste(match(covs, candidates), collapse = "_"))
  }
  score_many <- function(sets) {
    keys <- vapply(sets, key_for, character(1))
    cached <- vapply(keys, function(k) !is.null(records[[k]]), logical(1))
    hits <<- hits + sum(cached)
    pending <- which(!cached & !duplicated(keys))
    jobs <- lapply(pending, function(j) {
      report("start", keys[j], sets[[j]], NULL)
      list(label = keys[j], key = keys[j], covariates = sets[[j]],
        timeout = cfg$fit_timeout_seconds, fun = worker_task,
        args = list(kind = "score", covariates = sets[[j]], y = y))
    })
    pool$run(jobs, on_result = function(job, value, seconds) {
      if (length(value$waic) != 1L || !is.finite(value$waic)) stop("Non-finite candidate score.")
      record <- list(covariates = job$covariates, waic = value$waic, elapsed_seconds = seconds,
                      completed = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"))
      records[[job$key]] <<- record
      new_fits <<- new_fits + 1L
      atomic_save_rds(list(context = context, records = records), checkpoint_file)
      report("end", job$key, job$covariates, record)
    })
    # Return in the original formula/loop order, irrespective of finish order.
    vapply(keys, function(k) records[[k]]$waic, numeric(1), USE.NAMES = FALSE)
  }
  list(score = function(covs) score_many(list(covs))[[1]], score_many = score_many,
       records = function() records,
       stats = function() list(cache_hits = hits, new_fits = new_fits))
}

build_model_stacks <- function(prepared) {
  n <- nrow(prepared$observation_lonlat)
  xo <- transform_model_covariates(prepared$raw_observation_covariates, prepared)
  if (any(!is.finite(as.matrix(xo)))) stop("Missing observed covariates.")
  effects <- data.frame(b0 = rep(1, n), xo)
  if (prepared$use_location_iid) effects$site_iid <- seq_len(n)
  idx <- INLA::inla.spde.make.index("s", prepared$spde$n.spde)
  est <- INLA::inla.stack(tag = "est", data = list(y = rep(0, n), numtrials = prepared$numtrials),
    A = list(1, prepared$A), effects = list(effects, s = idx))
  prepared$observation_data <- INLA::inla.stack.data(est)
  prepared$observation_A <- INLA::inla.stack.A(est)
  prepared$observation_index <- INLA::inla.stack.index(est, "est")$data
  # Prediction rows do not enter an INLA fit. Ap is used on joint SPDE draws.
  for (nm in c("full_data", "full_A", "full_observation_index", "prediction_index",
                "surface_prediction_index")) prepared[[nm]] <- NULL
  prepared
}

sample_component_rows <- function(labels, component, expected_n) {
  at <- which(startsWith(labels, paste0(component, ":")))
  if (expected_n == 1L && !length(at)) at <- which(labels == component)
  if (length(at) != expected_n) stop("Posterior component mapping failed: ", component)
  if (expected_n > 1L) {
    ids <- suppressWarnings(as.integer(substring(labels[at], nchar(component) + 2L)))
    at <- at[order(node_mesh_rows(ids, expected_n))]
  }
  at
}

sample_joint_posterior <- function(res, prepared, covariates, cfg, draw) {
  if (is.null(res$misc$configs)) stop("INLA fit has no posterior configuration; config=TRUE is required.")
  fixed_names <- c("b0", effect_terms(covariates, prepared))
  components <- c("s", fixed_names, if (prepared$use_location_iid) "site_iid")
  selection <- setNames(rep(list(0L), length(components)), components)
  seed <- as.integer((as.double(cfg$seed) + 1009 * draw) %% 2147483646 + 1)
  set.seed(seed)
  samples <- INLA::inla.posterior.sample(cfg$posterior_samples, res, selection = selection,
    intern = FALSE, use.improved.mean = TRUE, add.names = FALSE, seed = seed, num.threads = 1L)
  labels <- rownames(samples[[1]]$latent)
  if (is.null(labels)) labels <- names(samples[[1]]$latent)
  if (is.null(labels) || anyDuplicated(labels)) stop("Missing/duplicate posterior latent labels.")
  b <- vapply(fixed_names, function(nm) sample_component_rows(labels, nm, 1L), integer(1))
  s <- sample_component_rows(labels, "s", prepared$spde$n.spde)
  nobs <- nrow(prepared$observation_lonlat)
  u <- if (prepared$use_location_iid) sample_component_rows(labels, "site_iid", nobs) else integer(0)
  extract <- function(rows) vapply(samples, function(z) as.numeric(z$latent)[rows], numeric(length(rows)))
  beta <- matrix(extract(b), nrow = length(b), dimnames = list(fixed_names, NULL))
  spatial <- matrix(extract(s), nrow = length(s))
  iid <- if (length(u)) matrix(extract(u), nrow = nobs) else matrix(0, nobs, length(samples))
  sigma <- rep(0, length(samples))
  if (prepared$use_location_iid) {
    hn <- rownames(samples[[1]]$hyperpar)
    if (is.null(hn)) hn <- names(samples[[1]]$hyperpar)
    j <- which(hn == "Precision for site_iid")
    if (length(j) != 1L) stop("Cannot find sampled IID precision on user scale.")
    precision <- vapply(samples, function(z) as.numeric(z$hyperpar)[j], numeric(1))
    if (any(!is.finite(precision) | precision <= 0)) stop("Invalid sampled IID precision.")
    sigma <- 1 / sqrt(precision)
  }
  if (any(!is.finite(beta)) || any(!is.finite(spatial)) || any(!is.finite(iid)))
    stop("Non-finite joint posterior samples.")
  list(beta = beta, spatial = spatial, observed_iid = iid, sigma_iid = sigma, seed = seed)
}

# Gaussian quadrature for a STANDARD NORMAL expectation (not exp(-x^2)).
normal_quadrature <- function(n) {
  J <- matrix(0, n, n)
  j <- seq_len(n - 1L)
  J[cbind(j, j + 1L)] <- J[cbind(j + 1L, j)] <- sqrt(j)
  e <- eigen(J, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}
integrate_iid <- function(eta, sigma, rule) {
  ans <- numeric(length(eta))
  for (k in seq_along(rule$z)) ans <- ans + rule$w[k] * plogis(eta + sigma * rule$z[k])
  ans
}

# Precompute logistic-normal expectations once. No probability/interval clipping.
# Small, fixed interpolation error is checked against a higher quadrature order.
make_iid_integrator <- function(sigma, tolerance = 1e-4) {
  max_sd <- max(sigma)
  if (!is.finite(max_sd) || max_sd < 0) stop("Invalid IID SD.")
  if (max_sd == 0) return(list(evaluate = function(eta, sigma) plogis(eta),
    error = 0, order = 0L, interpolation_step = 0))
  step <- .05
  sigmas <- seq(0, ceiling(max_sd / step) * step, by = step)
  etas <- seq(-25, 25, by = step)
  # Avoid an enormous table for a pathological SD posterior; fail explicitly.
  if (length(sigmas) * length(etas) > 2e6)
    stop("IID posterior is too diffuse for the integration table. Inspect IID_SD and the fitted model.")
  ee <- rep(etas, length(sigmas)); ss <- rep(sigmas, each = length(etas))
  evaluate_table <- function(eta, sigma, table, rule) {
    if (length(sigma) == 1L) sigma <- rep(sigma, length(eta))
    stopifnot(length(sigma) == length(eta), all(is.finite(eta)), all(is.finite(sigma)))
    ans <- numeric(length(eta))
    ok <- abs(eta) <= 25 & sigma >= 0 & sigma <= tail(sigmas, 1)
    if (any(ok)) {
      a <- (eta[ok] + 25) / step
      b <- sigma[ok] / step
      i <- pmin(floor(a), length(etas) - 2L); j <- pmin(floor(b), length(sigmas) - 2L)
      da <- a - i; db <- b - j
      k <- i + 1L + j * length(etas)
      ans[ok] <- (1 - da) * (1 - db) * table[k] + da * (1 - db) * table[k + 1L] +
        (1 - da) * db * table[k + length(etas)] + da * db * table[k + length(etas) + 1L]
    }
    if (any(!ok)) ans[!ok] <- integrate_iid(eta[!ok], sigma[!ok], rule)
    ans
  }
  # Test between table nodes as well as on nodes; includes the largest sampled SD.
  test_eta <- rep(seq(-24.975, 24.975, length.out = 121), 13)
  test_sd <- rep(seq(0, max_sd, length.out = 13), each = 121)
  for (order in c(31L, 61L, 121L)) {
    rule <- normal_quadrature(order)
    reference <- normal_quadrature(2L * order + 1L)
    tab <- integrate_iid(ee, ss, rule)
    err <- max(abs(evaluate_table(test_eta, test_sd, tab, rule) -
                     integrate_iid(test_eta, test_sd, reference)))
    if (err <= tolerance) break
  }
  if (err > tolerance) stop("IID integration accuracy check failed; maximum probability error ", err)
  rm(ee, ss)
  list(evaluate = function(eta, sigma) evaluate_table(eta, sigma, tab, rule),
       error = err, order = order, interpolation_step = step)
}

summarise_rows <- function(x) {
  # A one-cell diagnostic subset must still return a 1 x 3 quantile matrix.
  q <- matrixStats::rowQuantiles(x, probs = c(.025, .5, .975), type = 8,
                                useNames = FALSE, drop = FALSE)
  data.frame(mean = rowMeans(x), sd = matrixStats::rowSds(x),
    `0.025quant` = q[, 1], `0.5quant` = q[, 2], `0.975quant` = q[, 3],
    mode = NA_real_, check.names = FALSE)
}
empty_summary <- function(n) {
  as.data.frame(matrix(NA_real_, n, 6L, dimnames = list(NULL,
    c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode"))), check.names = FALSE)
}

predict_joint <- function(obj, prepared, cfg, draw) {
  joint <- obj$joint; ns <- ncol(joint$beta)
  covariates <- obj$covariates
  valid <- support_at_predictions(prepared, covariates)$valid_prediction
  m <- length(valid)
  xp <- transform_model_covariates(prepared$raw_prediction_covariates, prepared)
  fixed_terms <- effect_terms(covariates, prepared)
  X <- cbind(b0 = 1, as.matrix(xp[, fixed_terms, drop = FALSE]))
  if (!identical(colnames(X), rownames(joint$beta))) stop("Prediction fixed-effect order mismatch.")
  fun <- make_iid_integrator(joint$sigma_iid, cfg$integration_tolerance)
  main <- surface <- expected <- lin <- surf_lin <- expected_lin <- empty_summary(m)
  ok <- which(valid)
  if (!length(ok)) stop("No grid locations have all selected covariates.")
  check_at <- ok[unique(round(seq(1, length(ok), length.out = min(64L, length(ok)))))]
  half <- seq_len(ns %/% 2L)
  other_half <- setdiff(seq_len(ns), half)
  mc <- list()
  # Seed per cell makes random new-location IID draws invariant to chunk size.
  for (start in seq.int(1L, length(ok), by = cfg$prediction_chunk_size)) {
    rows <- ok[start:min(length(ok), start + cfg$prediction_chunk_size - 1L)]
    eta <- as.matrix(prepared$Ap[rows, , drop = FALSE] %*% joint$spatial) +
      X[rows, , drop = FALSE] %*% joint$beta
    if (any(!is.finite(eta))) stop("Non-finite projected latent field.")
    sigma <- rep(joint$sigma_iid, each = length(rows))
    averaged <- matrix(fun$evaluate(as.numeric(eta), sigma), nrow = length(rows))
    surface_p <- plogis(eta)
    z <- t(vapply(rows, function(i) {
      set.seed(as.integer((as.double(joint$seed) + 104729 * i) %% 2147483646 + 1))
      rnorm(ns)
    }, numeric(ns)))
    new_eta <- eta + matrix(sigma, nrow = length(rows)) * z
    new_p <- plogis(new_eta)
    main[rows, ] <- summarise_rows(new_p)
    surface[rows, ] <- summarise_rows(surface_p)
    expected[rows, ] <- summarise_rows(averaged)
    # E_theta E_u[p|theta] is the SAME mean for the first and third targets.
    main$mean[rows] <- expected$mean[rows]
    lin[rows, ] <- summarise_rows(new_eta)
    surf_lin[rows, ] <- summarise_rows(eta)
    lin$mean[rows] <- surf_lin$mean[rows]
    lin$sd[rows] <- sqrt(surf_lin$sd[rows]^2 + mean(joint$sigma_iid^2))
    # qlogis can be infinite at machine-rounded 0/1; retain NA for that logit
    # diagnostic only. Probability summaries are not truncated or altered.
    averaged_logit <- qlogis(averaged)
    finite_rows <- which(rowSums(!is.finite(averaged_logit)) == 0L)
    if (length(finite_rows)) expected_lin[rows[finite_rows], ] <-
      summarise_rows(averaged_logit[finite_rows, , drop = FALSE])
    jj <- which(rows %in% check_at)
    if (length(jj)) for (target in c("new_location", "surface_IID_zero", "expected_prevalence")) {
      xx <- switch(target, new_location = new_p, surface_IID_zero = surface_p, expected_prevalence = averaged)
      a <- summarise_rows(xx[jj, half, drop = FALSE]); b <- summarise_rows(xx[jj, other_half, drop = FALSE])
      mc[[length(mc) + 1L]] <- data.frame(cell = rows[jj], target = target, samples = ns,
        mean_half_difference_pp = 100 * abs(a$mean - b$mean),
        lower_half_difference_pp = 100 * abs(a[["0.025quant"]] - b[["0.025quant"]]),
        upper_half_difference_pp = 100 * abs(a[["0.975quant"]] - b[["0.975quant"]]))
    }
  }
  mc <- do.call(rbind, mc)
  mc$exceeds_screen <- pmax(mc$mean_half_difference_pp, mc$lower_half_difference_pp,
                             mc$upper_half_difference_pp) > cfg$mc_screen_pp
  # This is a Monte Carlo screen, not a proof of model adequacy/convergence.
  xo <- transform_model_covariates(prepared$raw_observation_covariates, prepared)
  xo <- cbind(b0 = 1, as.matrix(xo[, fixed_terms, drop = FALSE]))
  pp_obs <- summarise_rows(plogis(as.matrix(prepared$A %*% joint$spatial) +
                                  xo %*% joint$beta + joint$observed_iid))
  obs <- obj$observed_prevalence
  calibration <- data.frame(location = seq_len(nrow(obs)), inla_mean = obs$mean,
    sampled_mean = pp_obs$mean, mean_difference_pp = 100 * (pp_obs$mean - obs$mean),
    inla_lower = obs[["0.025quant"]], sampled_lower = pp_obs[["0.025quant"]],
    inla_upper = obs[["0.975quant"]], sampled_upper = pp_obs[["0.975quant"]])
  list(covariates = covariates, prevalence = main, linear_predictor = lin,
    surface_prevalence = surface, surface_linear_predictor = surf_lin,
    expected_prevalence = expected, expected_linear_predictor = expected_lin,
    use_location_iid = prepared$use_location_iid,
    iid_sd = if (prepared$use_location_iid) iid_sd_quantiles(obj$hyperparameters) else NULL,
    observed_iid = obj$observed_iid, observed_prevalence = obs, valid_prediction = valid,
    fixed = obj$fixed, hyperparameters = obj$hyperparameters, spatial_nodes = obj$spatial_nodes,
    waic_observations = obj$waic, posterior_samples = ns,
    monte_carlo = mc, mc_check = if (any(mc$exceeds_screen)) "REVIEW_MONTE_CARLO" else "SCREEN_PASSED",
    posterior_calibration = calibration,
    integration_check = data.frame(max_checked_probability_error = fun$error,
      quadrature_order = fun$order, interpolation_step = fun$interpolation_step,
      tolerance = cfg$integration_tolerance))
}

expected_result <- function(saved) {
  ans <- saved
  ans$prevalence <- saved$expected_prevalence
  ans$linear_predictor <- saved$expected_linear_predictor
  ans
}

run_aat_quick <- function(cfg) {
  if (!is.character(cfg$candidate_source_dir) || length(cfg$candidate_source_dir) != 1L ||
      is.na(cfg$candidate_source_dir)) stop("candidate_source_dir must be 'auto', a run folder, or ''.")
  if (length(cfg$covariate_degree) != 1L || is.na(cfg$covariate_degree) ||
      !cfg$covariate_degree %in% c(1L, 2L)) stop("covariate_degree must be 1 or 2.")
  required_packages <- c("INLA", "sf", "sp", "afrilearndata", "concaveman", "geodata",
                         "raster", "terra", "Matrix", "callr", "matrixStats")
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
    ". For the new runtime helpers, run install.packages(c('callr', 'matrixStats')).")
  suppressPackageStartupMessages({
    library(INLA); library(sf); library(sp); library(afrilearndata)
    library(concaveman); library(geodata); library(raster)
  })
  cfg$runtime <- capture_inla_runtime()
  if (is.null(cfg$inla_mode)) cfg$inla_mode <- INLA::inla.getOption("inla.mode")
  if (!is.character(cfg$inla_mode) || length(cfg$inla_mode) != 1L || !nzchar(cfg$inla_mode))
    stop("Cannot identify the current INLA mode. Set inla_mode to your working session's mode.")
  cfg$local_work_root <- tempfile("aat_inla_workers_")
  if (!dir.create(cfg$local_work_root)) stop("Cannot create local INLA scratch directory.")
  profiles <- list(pilot20km = c(20000, 50000), check10km = c(10000, 25000),
                   reference6km = c(6000, 15000))
  if (!cfg$mesh_profile %in% names(profiles)) stop("Unknown mesh_profile.")
  for (nm in c("coefficient_prior_sd", "intercept_prior_sd", "support_tolerance",
      "max_session_minutes", "fit_timeout_seconds", "sampling_timeout_seconds",
      "prediction_timeout_seconds", "waic_recheck_tolerance", "integration_tolerance", "mc_screen_pp",
      "waic_min_improvement", "waic_removal_tolerance")) {
    z <- cfg[[nm]]
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z <= 0) stop("Invalid ", nm)
  }
  for (nm in c("workers", "posterior_samples", "prediction_chunk_size", "seed")) {
    z <- cfg[[nm]]
    if (length(z) != 1L || !is.finite(z) || z < 1 || z != as.integer(z)) stop("Invalid ", nm)
  }
  if (cfg$posterior_samples < 1024L || cfg$posterior_samples %% 2L)
    stop("posterior_samples must be an even integer of at least 1024.")
  if (!is.logical(cfg$use_location_iid) || length(cfg$use_location_iid) != 1L || is.na(cfg$use_location_iid))
    stop("use_location_iid must be TRUE or FALSE.")
  z <- cfg$iid_prior
  if (!is.numeric(z) || length(z) != 2L || any(!is.finite(z)) || z[1] <= 0 || z[2] <= 0 || z[2] >= 1)
    stop("iid_prior must be c(positive SD, probability between 0 and 1).")
  for (nm in c("diagnostic_draws", "save_fit_draws", "save_joint_draws")) {
    z <- cfg[[nm]]
    if (!is.numeric(z) || anyNA(z) || any(z < 1 | z != as.integer(z))) stop("Invalid ", nm)
  }
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores)) cores <- 2L
  cfg$workers <- as.integer(min(cfg$workers, max(1L, cores - 1L)))
  countries_to_infer <- "Ethiopia"
  # The 14 actual selection candidates in the original Ethiopia analysis.
  candidate_covariates <- c("elevation", "precipitation", "human_fp", "pop_den", "tree", "grassland",
                            "shrub", "cropland", "built", "bare", "water", "wetland", "cattle", "tsetse")
  spatial_crs <- paste("+proj=aeqd +lat_0=9 +lon_0=40 +datum=WGS84", "+units=m +no_defs +type=crs")
  mesh_cutoff_m <- 3000
  mesh_max_edge_m <- profiles[[cfg$mesh_profile]]
  mesh_offset_m <- c(50000, 100000)
  range_prior_m <- c(50000, 0.1)
  sigma_prior <- c(1, 0.1)
  if (sf::st_is_longlat(sf::st_crs(spatial_crs)) ||
      !sf::st_crs(spatial_crs)$units_gdal %in% c("metre", "meter", "metres", "meters", "m"))
    stop("spatial_crs must use projected metres.")

  data_path <- path.expand(cfg$data_path)
  if (!file.exists(data_path)) stop("Missing input: ", data_path, ". Check the project working directory.")
  data_path <- normalizePath(data_path, mustWork = TRUE)
  message("Reading ", data_path)
  data_original <- tryCatch(read.csv(data_path), error = function(e) {
    stop("Cannot read CSV: ", data_path, "\n", conditionMessage(e),
         "\nCheck that the file contains readable data and is fully downloaded from OneDrive.", call. = FALSE)
  })
  if (!nrow(data_original) || !all(c("longitude", "latitude", "sample_size") %in% names(data_original))) {
    stop("Input must contain data rows and longitude, latitude, sample_size columns.")
  }
  data_original <- filter_ethiopia_inputs(data_original)
  if (!is.numeric(data_original$sample_size) || any(!is.finite(data_original$sample_size)) ||
      any(data_original$sample_size <= 0) ||
      any(abs(data_original$sample_size - round(data_original$sample_size)) > 1e-8)) {
    stop("sample_size must contain positive integer counts.")
  }
  iteration_columns <- cfg$iteration_columns
  if (!is.numeric(iteration_columns) || !length(iteration_columns) || anyNA(iteration_columns) ||
      any(iteration_columns < 8 | iteration_columns > ncol(data_original)) ||
      any(iteration_columns != as.integer(iteration_columns)) || anyDuplicated(iteration_columns)) {
    stop("iteration_columns must identify distinct available adjusted-count columns, starting at column 8.")
  }
  iteration_columns <- as.integer(iteration_columns)
  input_md5 <- unname(tools::md5sum(data_path))
  if (is.na(input_md5)) stop("Could not fingerprint input file.")
  package_names <- c("INLA", "sf", "sp", "afrilearndata", "concaveman", "geodata", "terra", "raster", "Matrix")
  package_versions <- vapply(package_names, function(p) as.character(packageVersion(p)), character(1))

  signature <- list(workflow_version = "aat_joint_grouped_quadratic_ETH_20260930_v1", input_md5 = input_md5,
    country = "ETH", input_filter = "ETH_bbox_original_orientation_v1",
    elevation_policy = "unmasked_30s_crs_extraction_prediction_aggregate5_na_rm_v1",
    candidate_covariates = candidate_covariates,
    selection_rule = "empirical_random_complete_set_v1",
    random_catalogue_hash = cfg$random_catalogue_hash,
    covariate_degree = cfg$covariate_degree, basis_version = "observed_scaled_z_plus_z_squared_v1",
    spatial_crs = spatial_crs, cutoff_m = mesh_cutoff_m, max_edge_m = mesh_max_edge_m,
    offset_m = mesh_offset_m, range_prior_m = range_prior_m, sigma_prior = sigma_prior,
    candidate_mode = "observations", package_versions = package_versions,
    R_version = as.character(getRversion()), num_threads = "1:1", inla_mode = cfg$inla_mode,
    inla_runtime_options = cfg$runtime$options, r_architecture = cfg$runtime$architecture,
    inla_runtime_environment = cfg$runtime$environment,
    control_inla = cfg$control_inla, coefficient_prior_sd = cfg$coefficient_prior_sd,
    intercept_prior_sd = cfg$intercept_prior_sd, use_location_iid = cfg$use_location_iid,
    iid_prior = cfg$iid_prior, support_tolerance = cfg$support_tolerance,
    thresholds = c(backward = cfg$waic_removal_tolerance, forward = -cfg$waic_min_improvement,
                   final_backward = cfg$waic_removal_tolerance),
    missing_prediction_policy = "NA_if_any_selected_covariate_missing")
  prediction_signature <- list(version = "joint_three_targets_v1", posterior_samples = cfg$posterior_samples,
    integration_tolerance = cfg$integration_tolerance, seed = cfg$seed,
    waic_recheck_tolerance = cfg$waic_recheck_tolerance, mc_screen_pp = cfg$mc_screen_pp,
    matrixStats_version = as.character(packageVersion("matrixStats")))
  # Dedicated, stable subfolder: never replace an original fitted run.
  cfg$reuse_preparation_dir <- random_preparation_source(cfg, signature)
  output_root <- cfg$random_output_root
  if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)
  output_root <- normalizePath(output_root, mustWork = TRUE)
  lock_dir <- file.path(output_root, ".run_lock")
  if (!dir.create(lock_dir, showWarnings = FALSE))
    stop("Run folder is locked. If an earlier R process was killed, remove .run_lock only after confirming it has stopped.")
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  manifest_file <- file.path(output_root, "run_manifest.rds")
  if (file.exists(manifest_file)) {
    if (!identical(readRDS(manifest_file), signature)) stop("Input/model/version mismatch. Use a new output folder.")
  } else {
    if (nzchar(cfg$resume_dir)) stop("resume_dir has no matching run manifest.")
    atomic_save_rds(signature, manifest_file)
  }
  pred_manifest_file <- file.path(output_root, "prediction_manifest.rds")
  # Each completed draw records its own prediction signature. Changed sampling
  # settings regenerate predictions, retaining model selection checkpoints.
  atomic_save_rds(prediction_signature, pred_manifest_file)
  output_dirs <- setNames(file.path(output_root, c("WAIC_increase_ETH", "Covariates_ETH",
    "Projections_ETH", "Diagnostics_ETH", "Checkpoints_ETH", "Surface_projections_ETH", "Expected_prevalence_ETH")),
    c("waic", "covariates", "projections", "diagnostics", "checkpoints", "surface", "expected"))
  for (d in output_dirs) if (!dir.exists(d) && !dir.create(d)) stop("Cannot create ", d)
  cfg$diagnostics_dir <- output_dirs[["diagnostics"]]
  log_file <- file.path(output_dirs[["diagnostics"]], "progress.log")
  log_message <- function(...) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
    cat(line, "\n"); flush.console()
    cat(line, "\n", file = log_file, append = TRUE)
  }
  log_message("RESULTS: ", output_root)
  log_message("Covariates: whole sets sampled with replacement from the frozen original ETH catalogue.")
  log_message("One observation fit per new count draw; no covariate search is run.")
  log_message("INLA mode: ", cfg$inla_mode, "; runtime options copied from this session.")
  log_message("Model time budget: ", cfg$max_session_minutes, " minutes; first-time preparation is separate.")
  if (nzchar(cfg$reuse_preparation_dir)) log_message("Reusing spatial preparation from: ", cfg$reuse_preparation_dir)
  atomic_save_rds(cfg, file.path(output_root, "settings_used.rds"))
  writeLines(c(
    "ETH RANDOM COVARIATE-SET CONTINUATION",
    "Each new adjusted-count draw receives one complete covariate set, sampled with replacement.",
    "Each original completed donor draw has equal probability; duplicate copies of the same draw are removed.",
    "Repeated sets in DIFFERENT original draws retain their empirical frequency.",
    "The catalogue is frozen; newly generated random selections never enter the donor pool.",
    "This is an approximation to repeating covariate selection, not Bayesian model averaging.",
    "The association between a count draw and its selected covariates is not retained under reassignment.",
    "No WAIC comparison is made to the donor's different response; there is no selection score for the new draw.",
    "Quadratic groups, priors, spatial model and location IID are unchanged from the supplied ETH engine.",
    "Every new draw is refitted to its own adjusted counts. No coefficients or hyperparameters are borrowed.",
    "Expected_prevalence_ETH integrates IID inside each joint posterior sample.",
    "Projections_ETH within this run contains NEW-LOCATION predictions including new IID realizations.",
    "Surface_projections_ETH contains predictions conditional on IID=0.",
    "Only Expected_prevalence_ETH is published into the parent Results_ETH folder.",
    "Joint posterior sampling and the original numerical diagnostics are retained.",
    "Saved fit checkpoints can be resampled later for correct joint area aggregation.",
    "Per-draw quantiles cannot be averaged into pooled or administrative-area credible intervals.",
    "The empirical-model approximation should be checked against full selection on a validation subset.",
    "Missing covariates remain NA and extrapolation is flagged, as in the original script.",
    "Original Latitude holds longitude; original Longitude holds latitude. This intentional CSV-label swap is unchanged.",
    "Flat merged donor folders lack full provenance; use only original ETH runs with the same model/data settings."),
    file.path(output_dirs[["diagnostics"]], "Review_notes.txt"))

  # Convert ordered longitude/latitude rows to the SAME projected CRS, in metres.
  project_to_m <- function(lonlat) {
    lonlat <- as.matrix(lonlat)[, 1:2, drop = FALSE]
    if (!nrow(lonlat) || any(!is.finite(lonlat)) ||
        any(abs(lonlat[, 1]) > 180) || any(abs(lonlat[, 2]) > 90)) {
      stop("Invalid longitude/latitude coordinates supplied for projection.")
    }
    points_ll <- sf::st_as_sf(
      data.frame(longitude = lonlat[, 1], latitude = lonlat[, 2]),
      coords = c("longitude", "latitude"), crs = "OGC:CRS84"
    )
    points_m <- sf::st_transform(points_ll, crs = spatial_crs)
    xy <- sf::st_coordinates(points_m)[, 1:2, drop = FALSE]
    if (nrow(xy) != nrow(lonlat) || any(!is.finite(xy))) {
      stop("Coordinate projection failed or changed the number of locations.")
    }
    roundtrip <- sf::st_coordinates(sf::st_transform(points_m, "OGC:CRS84"))[, 1:2, drop = FALSE]
    if (max(abs(roundtrip - lonlat)) > 1e-7) stop("Coordinate round-trip check failed.")
    xy
  }

  check_A <- function(A, n_locations, label) {
    if (nrow(A) != n_locations || any(!is.finite(A@x)) ||
        any(abs(as.numeric(Matrix::rowSums(A)) - 1) > 1e-6)) {
      stop(label, ": some locations are outside the mesh or have invalid projection weights.")
    }
  }

  # Use raster::extract for RasterLayer inputs and terra::extract for SpatRaster
  # inputs. Coordinates are transformed to the actual raster CRS, row order intact.
  extract_legacy_raster <- function(r, lonlat) {
    r_crs <- sf::st_crs(raster::crs(r))
    if (is.na(r_crs)) stop("A raster has no CRS; define its true CRS before extraction.")
    pts <- sf::st_as_sf(data.frame(lon = lonlat[, 1], lat = lonlat[, 2]),
                         coords = c("lon", "lat"), crs = "OGC:CRS84")
    xy <- sf::st_coordinates(sf::st_transform(pts, r_crs))[, 1:2, drop = FALSE]
    as.numeric(raster::extract(r, xy))
  }


  prepare_spatial_inputs <- function() {
    coord_sample_mapping <- aggregate(sample_size ~ longitude + latitude, data = data_original, FUN = sum)
    log_message("Preparing ", nrow(coord_sample_mapping), " unique locations from ", nrow(data_original), " records.")
    positive <- rep(0, nrow(coord_sample_mapping))
    kept_location_rows <- seq_len(nrow(coord_sample_mapping))
    sample_size <- coord_sample_mapping$sample_size
    # Create data object with aggregated coordinates
    data <- data.frame(
      longitude = as.numeric(coord_sample_mapping$longitude),
      latitude = as.numeric(coord_sample_mapping$latitude),
      stringsAsFactors = FALSE
    )
    #add elevation data
    #ELEVATION
    # Keep actual elevation near country borders rather than a second country's
    # raster mask. The original analysis polygon still determines the domain.
    r_elv <- geodata::elevation_30s(country = "ETH", path = "Data/Covariates", mask = FALSE)
    data$elevation <- extract_elevation_values(r_elv, data[, c("longitude", "latitude")])
    #PRECIPITATION
    r_prec<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "prec")
    data$precipitation<-(terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_1+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_2+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_3+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_4+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_5+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_6+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_7+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_8+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_9+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_10+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_11+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_prec_12)/12
    #AVERAGE TEMPERATURE
    r_tavg<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "tavg")
    data$tavg<-(terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_1+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_2+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_3+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_4+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_5+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_6+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_7+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_8+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_9+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_10+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_11+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tavg_12)/12
    #MINIMUM TEMPERATURE
    r_tmin<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "tmin")
    data$tmin<-(terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_1+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_2+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_3+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_4+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_5+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_6+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_7+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_8+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_9+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_10+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_11+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmin_12)/12
    #MAXIMUM TEMPERATURE
    r_tmax<- worldclim_country(country = "ETH",path = "Data/Covariates", var = "tmax")
    data$tmax<-(terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_1+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_2+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_3+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_4+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_5+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_6+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_7+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_8+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_9+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_10+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_11+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$ETH_wc2.1_30s_tmax_12)/12
    #POPULATION DENSITY
    r_pop_den <- population(year=2010, res=0.5, path="Data/Covariates")
    data$pop_den <-terra::extract(r_pop_den, data[, c("longitude", "latitude")])$population_density
    #LANDCOVER
    r_tree <- landcover("trees",path="Data/Covariates")
    data$tree <- terra::extract(r_tree, data[, c("longitude", "latitude")])$trees
    r_grassland <- landcover("grassland",path="Data/Covariates")
    data$grassland <- terra::extract(r_grassland, data[, c("longitude", "latitude")])$grassland
    r_shrub <- landcover("shrubs",path="Data/Covariates")
    data$shrub <- terra::extract(r_shrub, data[, c("longitude", "latitude")])$shrub
    r_cropland <- landcover("cropland",path="Data/Covariates")
    data$cropland <- terra::extract(r_cropland, data[, c("longitude", "latitude")])$cropland
    r_built <- landcover("built",path="Data/Covariates")
    data$built <- terra::extract(r_built, data[, c("longitude", "latitude")])$built
    r_bare <- landcover("bare",path="Data/Covariates")
    data$bare <- terra::extract(r_bare, data[, c("longitude", "latitude")])$bare
    r_water <- landcover("water",path="Data/Covariates")
    data$water <- terra::extract(r_water, data[, c("longitude", "latitude")])$water
    r_wetland <- landcover("wetland",path="Data/Covariates")
    data$wetland <- terra::extract(r_wetland, data[, c("longitude", "latitude")])$wetland
    r_mangrove <- landcover("mangroves",path="Data/Covariates")
    data$mangrove <- terra::extract(r_mangrove, data[, c("longitude", "latitude")])$mangroves
    r_moss <- landcover("moss",path="Data/Covariates")
    data$moss <- terra::extract(r_moss, data[, c("longitude", "latitude")])$moss
    data$tsetse_habitat <- data$tree + data$wetland
    #HUMAN FOOTPRINT
    r_human_fp<-footprint(year=2009, res=30, path="Data/Covariates")
    data$human_fp<-terra::extract(r_human_fp, data[, c("longitude", "latitude")])[,2]
    #CATTLE DATA
    r_cattle<-raster("Data/Covariates/livestock/cattle_2015/5_Ct_2015_Da.tif")
    data$cattle<-extract_legacy_raster(r_cattle, data[, c("longitude", "latitude")])
    data$cattle[is.na(data$cattle)] <- 0
    r_buffalo<-raster("Data/Covariates/livestock/buffalo_2015/5_Bf_2015_Da.tif")
    data$buffalo<-extract_legacy_raster(r_buffalo, data[, c("longitude", "latitude")])
    data$buffalo[is.na(data$buffalo)] <- 0
    r_goat<-raster("Data/Covariates/livestock/goats_2015/5_Gt_2015_Da.tif")
    data$goat<-extract_legacy_raster(r_goat, data[, c("longitude", "latitude")])
    data$goat[is.na(data$goat)] <- 0
    r_horse<-raster("Data/Covariates/livestock/horses_2015/5_Ho_2015_Da.tif")
    data$horse<-extract_legacy_raster(r_horse, data[, c("longitude", "latitude")])
    data$horse[is.na(data$horse)] <- 0
    r_pig<-raster("Data/Covariates/livestock/pigs_2015/5_Pg_2015_Da.tif")
    data$pig<-extract_legacy_raster(r_pig, data[, c("longitude", "latitude")])
    data$pig[is.na(data$pig)] <- 0
    r_sheep<-raster("Data/Covariates/livestock/sheep_2015/5_Sh_2015_Da.tif")
    data$sheep<-extract_legacy_raster(r_sheep, data[, c("longitude", "latitude")])
    data$sheep[is.na(data$sheep)] <- 0
    tsetse<-raster("Data/Covariates/tsenumbspec")
    tsetse[tsetse > 1] <- 1
    data$tsetse <-extract_legacy_raster(tsetse, data[, c("longitude", "latitude")])
    data$total_animal<-data$cattle+data$buffalo+data$goat+data$horse+data$pig+data$sheep
    latitudes<-data$latitude
    longitudes<-data$longitude
    elevations<-data$elevation
    precipitations<-data$precipitation
    tavgs<-data$tavg
    tmins<-data$tmin
    tmaxs<-data$tmax
    human_fps<-data$human_fp
    pop_dens<-data$pop_den
    trees<-data$tree
    grasslands<-data$grassland
    shrubs<-data$shrub
    croplands<-data$cropland
    builts<-data$built
    bares<-data$bare
    waters<-data$water
    wetlands<-data$wetland
    mangroves<-data$mangrove
    mosses<-data$moss
    tsetse_habitats<-data$tsetse_habitat
    cattles<-data$cattle
    buffalos<-data$buffalo
    goats<-data$goat
    horses<-data$horse
    pigs<-data$pig
    sheeps<-data$sheep
    total_animals<-data$total_animal
    tsetses<-data$tsetse
    countries_to_infer <- sort(countries_to_infer)
    #if Africa is in this list then perform the calculations for the whole continent
    if("Africa" %in% countries_to_infer){
      border <- st_coordinates(st_geometry(africontinent))
    #remove Madagascar
      border <- border[!(border[,1] > 40 & border[,1] < 60 & border[,2] > -30 & border[,2] < -10),]
      border <- border[,1:2]
      #remove Madagascar from africountries
      africountries <- africountries[!(africountries$name=="Madagascar"),]
      country <- africountries
    } else{
    v=africountries$name %in% countries_to_infer
    country <- africountries[which(v==TRUE),]
    }
    #remove all countries except the ones in countries_to_infer
    africountries_toplot <- africountries[which(africountries$name %in% countries_to_infer),]
    if (is.na(st_crs(country))) stop("Country boundary has no CRS.")
    country <- st_transform(country, "OGC:CRS84")
    country <- st_coordinates(st_geometry(country))
    poly.df <- data.frame("long"=country[,1],"lat"=country[,2])
    poly.sf <- st_as_sf(poly.df, coords = c("long","lat"))
    poly <- concaveman(poly.sf, length_threshold = 0,concavity=1.1)
    point.df <- data.frame("long"=longitudes,"lat"=latitudes)
    point.sf <- st_as_sf(point.df, coords = c("long","lat"))
    long_lat_concave <- st_coordinates(st_geometry(poly))
    border <- long_lat_concave
    inside <- lengths(st_intersects(point.sf, poly)) > 0
    if (!any(inside)) stop("No observations lie within the original country mask.")
    kept_location_rows <- kept_location_rows[inside]
    latitudes <- latitudes[inside]
    longitudes <- longitudes[inside]
    sample_size <- sample_size[inside]
    positive <- positive[inside]
    elevations <- elevations[inside]
    precipitations <- precipitations[inside]
    tavgs <- tavgs[inside]
    tmins <- tmins[inside]
    tmaxs <- tmaxs[inside]
    human_fps <- human_fps[inside]
    pop_dens <- pop_dens[inside]
    trees <- trees[inside]
    grasslands <- grasslands[inside]
    shrubs <- shrubs[inside]
    croplands <- croplands[inside]
    builts <- builts[inside]
    bares <- bares[inside]
    waters <- waters[inside]
    wetlands <- wetlands[inside]
    mangroves <- mangroves[inside]
    mosses <- mosses[inside]
    tsetse_habitats <- tsetse_habitats[inside]
    cattles <- cattles[inside]
    buffalos <- buffalos[inside]
    goats <- goats[inside]
    horses <- horses[inside]
    pigs <- pigs[inside]
    sheeps <- sheeps[inside]
    total_animals <- total_animals[inside]
    tsetses <- tsetses[inside]
    #remove NA values
    ind <- !is.na(latitudes)
    kept_location_rows <- kept_location_rows[ind]
    latitudes <- latitudes[ind]
    longitudes <- longitudes[ind]
    sample_size <- sample_size[ind]
    positive <- positive[ind]
    elevations <- elevations[ind]
    precipitations <- precipitations[ind]
    tavgs <- tavgs[ind]
    tmins <- tmins[ind]
    tmaxs <- tmaxs[ind]
    human_fps <- human_fps[ind]
    pop_dens <- pop_dens[ind]
    trees <- trees[ind]
    grasslands <- grasslands[ind]
    shrubs <- shrubs[ind]
    croplands <- croplands[ind]
    builts <- builts[ind]
    bares <- bares[ind]
    waters <- waters[ind]
    wetlands <- wetlands[ind]
    mangroves <- mangroves[ind]
    mosses <- mosses[ind]
    tsetse_habitats <- tsetse_habitats[ind]
    cattles <- cattles[ind]
    buffalos <- buffalos[ind]
    goats <- goats[ind]
    horses <- horses[ind]
    pigs <- pigs[ind]
    sheeps <- sheeps[ind]
    total_animals <- total_animals[ind]
    tsetses <- tsetses[ind]
    prevalence <- positive/sample_size
    # Geographic coordinates retained for raster extraction and CSV compatibility.
    coo_ll <- cbind(longitudes, latitudes)
    check_elevation_sites(elevations, coo_ll, sample_size, kept_location_rows,
                          output_dirs[["diagnostics"]])
    # Keep the original concave-hull mask and 200 x 200 geographic grid exactly.
    bb <- sp::bbox(border)
    x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 200)
    y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 200)
    coop <- as.matrix(expand.grid(x, y))
    ind <- sp::point.in.polygon(coop[, 1], coop[, 2], border[, 1], border[, 2])
    coop <- coop[which(ind == 1), , drop = FALSE]

    # Transform observations AND predictions before passing either to INLA.
    # Do not overwrite coop: all extraction and exports below require degrees.
    coo <- project_to_m(coo_ll)
    coop_m <- project_to_m(coop)

    # Cover the entire prediction domain, not only the observation convex hull.
    # loc.domain determines extent without inserting every grid point as a node.
    domain_m <- rbind(coo, coop_m)
    domain_m <- domain_m[chull(domain_m[, 1], domain_m[, 2]), , drop = FALSE]
    mesh <- inla.mesh.2d(
      loc = coo, loc.domain = domain_m,
      offset = mesh_offset_m, cutoff = mesh_cutoff_m,
      max.edge = mesh_max_edge_m
    )
    cat("Projected mesh:", mesh$n, "vertices;", nrow(coop), "prediction locations\n")
    spde <- inla.spde2.pcmatern(
      mesh = mesh, alpha = 2,
      prior.range = range_prior_m, prior.sigma = sigma_prior, constr = TRUE
    )
    indexs <- inla.spde.make.index("s", spde$n.spde)
    A <- inla.spde.make.A(mesh = mesh, loc = coo)
    Ap <- inla.spde.make.A(mesh = mesh, loc = coop_m)
    check_A(A, nrow(coo), "Observation projector")
    check_A(Ap, nrow(coop_m), "Prediction projector")

    # Preserve the original five-cell aggregation; ignore NoData cells within
    # partially covered blocks. Entirely missing blocks remain missing.
    ra_elv <- terra::aggregate(r_elv, fact = 5, fun = mean, na.rm = TRUE)
    dp <- as.data.frame(coop)
    #coop <- terra::crds(ra)
    dp$elevation <- extract_elevation_values(ra_elv, coop)
    dp$precipitation <- rowSums(terra::extract(r_prec, coop)[, 1:12])/12
    dp$tavg <- rowSums(terra::extract(r_tavg, coop)[, 1:12])/12
    dp$tmin <- rowSums(terra::extract(r_tmin, coop)[, 1:12])/12
    dp$tmax <- rowSums(terra::extract(r_tmax, coop)[, 1:12])/12
    dp$human_fp <- terra::extract(r_human_fp, coop)[, 1]
    dp$pop_den <- terra::extract(r_pop_den, coop)[, 1]
    dp$tree <- terra::extract(r_tree, coop)[, 1]
    dp$grassland <- terra::extract(r_grassland, coop)[, 1]
    dp$shrub <- terra::extract(r_shrub, coop)[, 1]
    dp$cropland <- terra::extract(r_cropland, coop)[, 1]
    dp$built <- terra::extract(r_built, coop)[, 1]
    dp$bare <- terra::extract(r_bare, coop)[, 1]
    dp$water <- terra::extract(r_water, coop)[, 1]
    dp$wetland <- terra::extract(r_wetland, coop)[, 1]
    dp$mangrove <- terra::extract(r_mangrove, coop)[, 1]
    dp$moss <- terra::extract(r_moss, coop)[, 1]
    dp$tsetse_habitat <- dp$tree + dp$wetland
    dp$cattle <- extract_legacy_raster(r_cattle, coop)
    dp$cattle[is.na(dp$cattle)] <- 0
    dp$buffalo <- extract_legacy_raster(r_buffalo, coop)
    dp$buffalo[is.na(dp$buffalo)] <- 0
    dp$goat <- extract_legacy_raster(r_goat, coop)
    dp$goat[is.na(dp$goat)] <- 0
    dp$horse <- extract_legacy_raster(r_horse, coop)
    dp$horse[is.na(dp$horse)] <- 0
    dp$pig <- extract_legacy_raster(r_pig, coop)
    dp$pig[is.na(dp$pig)] <- 0
    dp$sheep <- extract_legacy_raster(r_sheep, coop)
    dp$sheep[is.na(dp$sheep)] <- 0
    dp$total_animal <- dp$cattle + dp$buffalo + dp$goat + dp$horse + dp$pig + dp$sheep
    dp$tsetse <- extract_legacy_raster(tsetse, coop)
    # The SpatRaster sources above are geographic geodata products; their
    # extraction still uses geographic coo_ll/coop, never projected metres.
    # stack for estimation stk.e
    stk.e <- inla.stack(
      tag = "est",
      data = list(y = positive, numtrials = sample_size),
      A = list(1, A),
      effects = list(data.frame(b0 = rep(1, nrow(coo)), elevation=elevations,
      precipitation=precipitations, tavg=tavgs, tmin=tmins, tmax=tmaxs,
      human_fp=human_fps, pop_den=pop_dens, tree=trees, grassland=grasslands,
      shrub=shrubs, cropland=croplands, built=builts, bare=bares, water=waters, wetland=wetlands,
      mangrove=mangroves, moss=mosses, tsetse_habitat=tsetse_habitats, cattle=cattles, goat=goats,
      horse=horses, pig=pigs, buffalo=buffalos, sheep=sheeps, total_animal=total_animals, tsetse=tsetses), s = indexs)
    )
    # stack for prediction stk.p
    stk.p <- inla.stack(
      tag = "pred",
      data = list(y = NA, numtrials = NA),
      A = list(1, Ap),
      effects = list(data.frame(b0 = rep(1, nrow(coop)),elevation=dp$elevation,
      precipitation=dp$precipitation, tavg=dp$tavg, tmin=dp$tmin, tmax=dp$tmax,
      human_fp=dp$human_fp, pop_den=dp$pop_den, tree=dp$tree, grassland=dp$grassland,
      shrub=dp$shrub, cropland=dp$cropland, built=dp$built, bare=dp$bare, water=dp$water, 
      wetland=dp$wetland, mangrove=dp$mangrove, moss=dp$moss, tsetse_habitat=dp$tsetse_habitat, cattle=dp$cattle,
      goat=dp$goat, horse=dp$horse, pig=dp$pig, buffalo=dp$buffalo, sheep=dp$sheep,
       total_animal=dp$total_animal, tsetse=dp$tsetse), s = indexs)
    )
    # stk.full has stk.e and stk.p
    stk.full <- inla.stack(stk.e, stk.p)
    covariates <- candidate_covariates
    cov_temp <- character(0)
    for (i in seq_along(covariates)){
      col_number <- which(colnames(dp) == covariates[i])
      x=!is.na(dp[,col_number])
      if (any(x) && max(dp[,col_number][x]) > 0) {
        cov_temp <- c(cov_temp,covariates[i])
      }
    }
    covariates<-cov_temp
    full_covariates<-covariates


    list(
      mesh = mesh, spde = spde, A = A, Ap = Ap,
      observation_m = coo, prediction_m = coop_m,
      observation_lonlat = coo_ll, prediction_lonlat = coop,
      kept_location_rows = kept_location_rows,
      observation_data = inla.stack.data(stk.e),
      full_data = inla.stack.data(stk.full),
      observation_A = inla.stack.A(stk.e), full_A = inla.stack.A(stk.full),
      observation_index = inla.stack.index(stk.e, tag = "est")$data,
      full_observation_index = inla.stack.index(stk.full, tag = "est")$data,
      prediction_index = inla.stack.index(stk.full, tag = "pred")$data,
      candidate_covariates = as.character(covariates),
      coordinate_mapping = coord_sample_mapping
    )
  }

  remesh_prepared <- function(p, profile) {
    domain <- rbind(p$observation_m, p$prediction_m)
    domain <- domain[chull(domain[, 1], domain[, 2]), , drop = FALSE]
    p$mesh <- INLA::inla.mesh.2d(loc = p$observation_m, loc.domain = domain,
      offset = mesh_offset_m, cutoff = mesh_cutoff_m, max.edge = profiles[[profile]])
    p$spde <- INLA::inla.spde2.pcmatern(p$mesh, alpha = 2,
      prior.range = range_prior_m, prior.sigma = sigma_prior, constr = TRUE)
    p$A <- INLA::inla.spde.make.A(mesh = p$mesh, loc = p$observation_m)
    p$Ap <- INLA::inla.spde.make.A(mesh = p$mesh, loc = p$prediction_m)
    check_A(p$A, nrow(p$observation_m), "Observation projector")
    check_A(p$Ap, nrow(p$prediction_m), "Prediction projector")
    p
  }

  preparation_file <- file.path(output_dirs[["checkpoints"]], "spatial_preparation.rds")
  if (file.exists(preparation_file)) {
    log_message("Loading saved spatial preparation.")
    prepared <- readRDS(preparation_file)
  } else {
    t0 <- proc.time()[["elapsed"]]
    source_manifest <- NULL
    source_dir <- path.expand(cfg$reuse_preparation_dir)
    if (nzchar(source_dir) && dir.exists(source_dir)) {
      paths <- file.path(source_dir, c("run_manifest.rds", "Checkpoints_ETH/spatial_preparation.rds"))
      if (!all(file.exists(paths))) stop("Incomplete reuse_preparation_dir. Choose a complete run, or set it to ''.")
      source_manifest <- readRDS(paths[1])
      required <- c("country", "input_filter", "elevation_policy", "candidate_covariates", "input_md5", "spatial_crs", "cutoff_m", "offset_m", "range_prior_m", "sigma_prior", "package_versions")
      mismatch <- required[!vapply(required, function(nm) identical(source_manifest[[nm]], signature[[nm]]), logical(1))]
      if (length(mismatch)) stop("Preparation reuse mismatch: ", paste(mismatch, collapse = ", "),
                                 ". Set reuse_preparation_dir='' to prepare a fresh run.")
      log_message("Reusing ONLY spatial preparation from: ", source_dir)
      prepared <- readRDS(paths[2])
    } else {
      if (nzchar(source_dir)) log_message("Previous preparation folder unavailable; preparing from original inputs.")
      log_message("Preparing rasters, mesh and projection matrices ONCE.")
      prepared <- prepare_spatial_inputs()
    }
    pool <- candidate_covariates
    if (is.null(prepared$raw_observation_covariates)) {
      prepared$raw_observation_covariates <- recover_covariates(prepared$observation_data,
        prepared$observation_A, prepared$observation_index, pool)
      prepared$raw_prediction_covariates <- recover_covariates(prepared$full_data,
        prepared$full_A, prepared$prediction_index, pool)
      prepared$numtrials <- prepared$observation_data$numtrials[prepared$observation_index]
    }
    if (!is.null(source_manifest) && !identical(source_manifest$max_edge_m, mesh_max_edge_m)) {
      log_message("Rebuilding mesh/projectors on the cached coordinates; raster covariates are unchanged.")
      prepared <- remesh_prepared(prepared, cfg$mesh_profile)
    }
    prepared$covariate_audit <- audit_covariates(prepared$raw_observation_covariates,
                                                prepared$raw_prediction_covariates, cfg$support_tolerance)
    prepared$covariate_audit$slope_prior_sd <- cfg$coefficient_prior_sd
    prepared$support_tolerance <- cfg$support_tolerance
    atomic_write_csv(prepared$covariate_audit, file.path(output_dirs[["diagnostics"]], "Covariate_support.csv"))
    aa <- prepared$covariate_audit
    if (any(aa$invalid_fraction_units)) stop("Land-cover/tsetse values are not 0--1 fractions. Inspect Covariate_support.csv and raster units.")
    bad <- aa$covariate[aa$eligible & aa$observed_missing > 0]
    if (length(bad)) stop("Partly missing covariates at fitted sites: ", paste(bad, collapse = ", "),
                          ". See Covariate_support.csv; repair extraction rather than treating missing values as zeros.")
    prepared$candidate_covariates <- aa$covariate[aa$eligible]
    prepared$nearest_observation_km <- nearest_distances_km(prepared$prediction_m, prepared$observation_m)
    prepared$use_location_iid <- cfg$use_location_iid
    prepared$iid_prior <- cfg$iid_prior
    prepared$quadratic_basis <- audit_quadratic_basis(prepared$raw_observation_covariates,
      prepared$covariate_audit, cfg$covariate_degree)
    prepared <- build_model_stacks(prepared)
    prepared$preparation_seconds <- proc.time()[["elapsed"]] - t0
    atomic_save_rds(prepared, preparation_file)
    log_message("Preparation saved; elapsed ", round(prepared$preparation_seconds, 1), " seconds.")
  }
  atomic_write_csv(prepared$covariate_audit, file.path(output_dirs[["diagnostics"]], "Covariate_support.csv"))
  if (is.null(prepared$quadratic_basis)) stop("Missing saved quadratic basis. Use a new run folder.")
  donor_covariates <- unique(unlist(cfg$random_snapshot$catalogue$sets))
  if (!all(donor_covariates %in% prepared$candidate_covariates)) stop(
    "The donor catalogue includes covariates unavailable in this prepared ETH model: ",
    paste(setdiff(donor_covariates, prepared$candidate_covariates), collapse = ", "))

  atomic_write_csv(prepared$quadratic_basis, file.path(output_dirs[["diagnostics"]], "Quadratic_basis.csv"))
  log_message("Quadratic effects: ", sum(prepared$quadratic_basis$degree == 2L),
    "; selected as linear-plus-square groups. Remaining eligible covariates are linear.")
  excluded <- prepared$covariate_audit$covariate[!prepared$covariate_audit$eligible]
  log_message("Excluded unestimable covariates: ", if (length(excluded)) paste(excluded, collapse = ", ") else "none")
  xo <- as.matrix(transform_model_covariates(prepared$raw_observation_covariates, prepared))
  if (ncol(xo)) {
    cc <- cor(xo)
    atomic_write_csv(data.frame(covariate = rownames(cc), cc, check.names = FALSE),
                     file.path(output_dirs[["diagnostics"]], "Covariate_correlations.csv"))
    rank_x <- qr(cbind(1, xo), tol = 1e-10)$rank
    log_message("Observed design rank: ", rank_x, "/", ncol(xo) + 1L,
                ". Correlated varying covariates remain eligible; priors are explicit.")
  }
  if (any(prepared$covariate_audit$rare_support & prepared$covariate_audit$eligible))
    log_message("Rare covariates retained and flagged in Covariate_support.csv; inspect prior sensitivity.")
  log_message("Mesh vertices: ", prepared$mesh$n, "; fitted locations: ",
              length(prepared$kept_location_rows), "; prediction locations: ",
              nrow(prepared$prediction_lonlat))
  log_message("Candidates: ", paste(prepared$candidate_covariates, collapse = ", "))
  capture.output(sessionInfo(), file = file.path(output_dirs[["diagnostics"]], "sessionInfo.txt"))
  timing_path <- file.path(output_dirs[["diagnostics"]], "Timing_summary.csv")


  write_uncertainty_diagnostics <- function(saved, draw, y) {
    d <- uncertainty_rows(saved, prepared)
    sd <- uncertainty_rows(surface_result(saved), prepared)
    atomic_write_csv(uncertainty_summary(d, draw), file.path(output_dirs[["diagnostics"]],
      paste0("Uncertainty_summary_model_", draw, ".csv")))
    atomic_write_csv(uncertainty_summary(sd, draw, "surface_IID_zero"), file.path(output_dirs[["diagnostics"]],
      paste0("Uncertainty_summary_surface_model_", draw, ".csv")))
    ed <- uncertainty_rows(expected_result(saved), prepared)
    atomic_write_csv(uncertainty_summary(ed, draw, "expected_prevalence"),
      file.path(output_dirs[["diagnostics"]], paste0("Uncertainty_summary_expected_model_", draw, ".csv")))
    targets <- rbind(uncertainty_summary(ed, draw, "expected_prevalence"), uncertainty_summary(d, draw, "new_location"),
                     uncertainty_summary(sd, draw, "surface_IID_zero"))
    atomic_write_csv(targets[targets$group == "all", ], file.path(output_dirs[["diagnostics"]],
      paste0("Prediction_targets_model_", draw, ".csv")))
    if (saved$use_location_iid) atomic_write_csv(saved$iid_sd, file.path(output_dirs[["diagnostics"]],
      paste0("IID_SD_model_", draw, ".csv")))
    if (!draw %in% cfg$diagnostic_draws) return(invisible(NULL))
    atomic_write_csv(d, file.path(output_dirs[["diagnostics"]], paste0("Prediction_diagnostics_model_", draw, ".csv")))
    atomic_write_csv(sd, file.path(output_dirs[["diagnostics"]], paste0("Prediction_diagnostics_surface_model_", draw, ".csv")))
    atomic_write_csv(ed, file.path(output_dirs[["diagnostics"]], paste0("Prediction_diagnostics_expected_model_", draw, ".csv")))
    if (saved$use_location_iid) atomic_write_csv(data.frame(
      longitude_deg = prepared$observation_lonlat[, 1], latitude_deg = prepared$observation_lonlat[, 2],
      saved$observed_iid), file.path(output_dirs[["diagnostics"]], paste0("Observed_IID_effects_model_", draw, ".csv")))
    oo <- saved$observed_prevalence
    atomic_write_csv(data.frame(longitude_deg = prepared$observation_lonlat[, 1],
      latitude_deg = prepared$observation_lonlat[, 2], adjusted_positive = y, sample_size = prepared$numtrials,
      adjusted_observed_prevalence = y / prepared$numtrials, fitted_mean = oo$mean,
      fitted_lower = oo[["0.025quant"]], fitted_upper = oo[["0.975quant"]]),
      file.path(output_dirs[["diagnostics"]], paste0("Observed_fit_model_", draw, ".csv")))
    if (isTRUE(cfg$make_diagnostic_plots)) {
      tryCatch({
        grDevices::pdf(file.path(output_dirs[["diagnostics"]], paste0("Uncertainty_diagnostics_model_", draw, ".pdf")),
                        width = 10, height = 8)
        tryCatch({
          par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
          ok <- which(is.finite(d$mean))
          if (!length(ok)) stop("No valid predictions to plot.")
          ii <- ok[unique(round(seq(1, length(ok), length.out = min(3000, length(ok)))))]
          col <- ifelse(d$outside_selected_ranges[ii] > 0, "#c65c33", "#286b9c")
          plot(d$nearest_observation_km[ii], d$interval_width[ii], pch = 16, cex = .35, col = col,
            xlab = "Distance to nearest observation (km)", ylab = "95% interval width (0 to 1)", ylim = c(0, 1),
            main = "Distance and uncertainty")
          legend("bottomright", c("Within observed ranges", "Outside a selected range"),
                 col = c("#286b9c", "#c65c33"), pch = 16, cex = .7, bty = "n")
          boxplot(d$interval_width ~ cut(d$nearest_observation_km, c(-Inf, 10, 25, 50, 100, Inf),
            labels = c("0-10", "10-25", "25-50", "50-100", ">100")),
            xlab = "Distance band (km)", ylab = "95% interval width", ylim = c(0, 1), main = "All valid predictions")
          plot(d$Longitude[ii] ~ d$Latitude[ii], pch = 15, cex = .45,
            col = hcl.colors(100, "Inferno")[pmax(1, pmin(100, ceiling(100 * d$interval_width[ii])))],
            xlab = "Longitude (degrees E)", ylab = "Latitude (degrees N)", main = "Interval width: dark 0, light 1")
          points(prepared$observation_lonlat, pch = 3, cex = .4)
          plot(y / prepared$numtrials, oo$mean, pch = 16, cex = .6,
            xlim = c(0, 1), ylim = c(0, 1), xlab = "Adjusted observed prevalence", ylab = "Fitted mean prevalence",
            main = "In-sample fit (not validation)")
          abline(0, 1, col = "grey50", lty = 2)
          if (saved$use_location_iid) {
            plot(sd$nearest_observation_km[ii], sd$interval_width[ii], pch = 16, cex = .35, col = col,
              xlab = "Distance to nearest observation (km)", ylab = "95% surface interval width", ylim = c(0, 1),
              main = "Conditional surface: IID = 0")
            plot(sd$interval_width[ii], d$interval_width[ii], pch = 16, cex = .35, col = col,
              xlab = "Surface interval width (IID=0)", ylab = "New-location interval width",
              xlim = c(0, 1), ylim = c(0, 1), main = "Different prediction targets")
            abline(0, 1, col = "grey50", lty = 2)
            plot(sd$Longitude[ii] ~ sd$Latitude[ii], pch = 15, cex = .45,
              col = hcl.colors(100, "Inferno")[pmax(1, pmin(100, ceiling(100 * sd$interval_width[ii])))],
              xlab = "Longitude (degrees E)", ylab = "Latitude (degrees N)", main = "Conditional surface interval width")
            points(prepared$observation_lonlat, pch = 3, cex = .4)
            plot(sd$mean[ii], d$mean[ii], pch = 16, cex = .35, col = col,
              xlab = "Surface mean prevalence (IID=0)", ylab = "New-location mean prevalence",
              xlim = c(0, 1), ylim = c(0, 1), main = "Averaging on a nonlinear link changes means")
            abline(0, 1, col = "grey50", lty = 2)
          }
        }, finally = grDevices::dev.off())
      }, error = function(e) log_message("Diagnostic plot unavailable: ", conditionMessage(e)))
    }
  }


  pool <- NULL
  deadline <- elapsed_now() + 60 * cfg$max_session_minutes
  on.exit(if (!is.null(pool)) pool$close(), add = TRUE)
  get_pool <- function() {
    if (is.null(pool)) {
      helper_names <- c("atomic_save_rds", "elapsed_now", "fixed_prior_control", "model_formula",
        "restore_inla_runtime", "write_inla_failure", "worker_task_inner",
        "fit_observations", "sample_joint_posterior", "sample_component_rows", "node_mesh_rows",
        "transform_covariates", "transform_model_covariates", "effect_terms",
        "support_at_predictions", "normal_quadrature", "integrate_iid",
        "make_iid_integrator", "summarise_rows", "empty_summary", "predict_joint", "iid_sd_quantiles")
      functions_file <- file.path(output_dirs[["checkpoints"]], "worker_functions.rds")
      funs <- mget(helper_names, envir = environment(run_aat_quick), inherits = TRUE)
      atomic_save_rds(funs, functions_file)
      pool <<- make_worker_pool(cfg$workers, functions_file, preparation_file, cfg, deadline, log_message)
    }
    pool
  }
  stage <- function(kind, covariates, y, path, expected_waic, draw, timeout) {
    t0 <- elapsed_now()
    log_message("DRAW ", draw, " | ", toupper(kind), " START | limit ", timeout, " s")
    get_pool()$run(list(list(label = paste0("draw ", draw, " ", kind), timeout = timeout,
      fun = worker_task, args = list(kind = kind, covariates = covariates, y = y,
        path = path, expected_waic = expected_waic, draw = draw))))
    log_message("DRAW ", draw, " | ", toupper(kind), " END | ", round(elapsed_now() - t0, 1), " s")
    invisible(elapsed_now() - t0)
  }

  run_one_draw <- function(iteration) {
    draw <- iteration - 7L
    draw_start <- elapsed_now()
    draw_dir <- file.path(output_dirs[["checkpoints"]], sprintf("draw_%04d", draw))
    if (!dir.exists(draw_dir)) dir.create(draw_dir)
    done_file <- file.path(draw_dir, "completed.rds")
    result_paths <- c(
      file.path(output_dirs[["waic"]], paste0("WAIC_increase_model_", draw, ".csv")),
      file.path(output_dirs[["covariates"]], paste0("Covariates_model_", draw, ".csv")),
      file.path(output_dirs[["projections"]], paste0("Projections_model_", draw, ".csv")))
    extra_paths <- c(file.path(output_dirs[["surface"]], paste0("Projections_model_", draw, ".csv")),
      file.path(output_dirs[["expected"]], paste0("Projections_model_", draw, ".csv")),
      file.path(output_dirs[["diagnostics"]], paste0("Prediction_targets_model_", draw, ".csv")),
      file.path(output_dirs[["diagnostics"]], paste0("Monte_Carlo_model_", draw, ".csv")),
      file.path(output_dirs[["diagnostics"]], paste0("Posterior_calibration_model_", draw, ".csv")))
    if (file.exists(done_file)) {
      done <- readRDS(done_file)
      if (identical(done$prediction_signature, prediction_signature) && all(file.exists(done$files))) {
        log_message("DRAW ", draw, " already complete; reusing outputs.")
        return(done$timing)
      }
    }
    if (elapsed_now() >= deadline) pause_run("Session modelling budget reached before next draw.")
    counts <- aggregate_draw_counts(data_original, iteration, prepared$coordinate_mapping)
    y <- counts[prepared$kept_location_rows]
    stopifnot(length(y) == length(prepared$observation_index), all(y <= prepared$numtrials))
    assignment <- cfg$random_snapshot$plan[cfg$random_snapshot$plan$draw == draw, , drop = FALSE]
    if (nrow(assignment) != 1L) stop("Missing/duplicate saved random assignment for draw ", draw)
    covariates <- cfg$random_snapshot$catalogue$sets[[assignment$donor_row]]
    selected <- list(covariates = covariates, waic = NA_real_,
      delta = rep(NA_real_, length(covariates)), rule = "empirical_random_complete_set_v1",
      donor_draw = assignment$donor_draw, draw = draw)
    # An assignment is not a candidate search; no donor WAIC is copied to this response.
    atomic_save_rds(selected, file.path(draw_dir, "selected_model.rds"))
    atomic_write_csv(assignment, file.path(output_dirs[["diagnostics"]],
      paste0("Random_assignment_model_", draw, ".csv")))
    log_message("DRAW ", draw, " | RANDOM SET | donor draw ", assignment$donor_draw,
      " | ", length(covariates), " covariates: ", paste(covariates, collapse = ", "))
    context_file <- file.path(draw_dir, "prediction_context.rds")
    prediction_checkpoint <- file.path(draw_dir, "final_predictions.rds")
    joint_checkpoint <- file.path(draw_dir, "joint_posterior.rds")
    fit_checkpoint <- file.path(draw_dir, "selected_observation_fit.rds")
    if (!file.exists(context_file) || !identical(readRDS(context_file), prediction_signature)) {
      unlink(c(prediction_checkpoint, joint_checkpoint, done_file))
      atomic_save_rds(prediction_signature, context_file)
    }
    if (!file.exists(prediction_checkpoint)) {
      if (!file.exists(joint_checkpoint)) {
        if (!file.exists(fit_checkpoint)) stage("assigned_fit", covariates, y, fit_checkpoint,
          selected$waic, draw, cfg$fit_timeout_seconds)
        stage("sample", covariates, y, joint_checkpoint, selected$waic, draw, cfg$sampling_timeout_seconds)
      }
      stage("predict", covariates, y, prediction_checkpoint, selected$waic, draw, cfg$prediction_timeout_seconds)
    }
    saved <- readRDS(prediction_checkpoint)
    if (!identical(saved$covariates, covariates)) stop("Final prediction covariate mismatch.")
    if (length(saved$waic_observations) != 1L || !is.finite(saved$waic_observations))
      stop("Nonfinite fitted WAIC for the newly assigned model.")
    selected$waic <- saved$waic_observations
    atomic_save_rds(selected, file.path(draw_dir, "selected_model.rds"))
    coop <- prepared$prediction_lonlat
    atomic_write_csv(data.frame(covariates = covariates, waic_increase = selected$delta), result_paths[1])
    atomic_write_csv(data.frame(covariates = covariates), result_paths[2])
    atomic_write_csv(prepared$quadratic_basis[
      match(covariates, prepared$quadratic_basis$covariate), , drop = FALSE],
      file.path(output_dirs[["diagnostics"]], paste0("Selected_effect_terms_model_", draw, ".csv")))
    atomic_write_csv(projection_table(saved$prevalence, coop), result_paths[3])
    atomic_write_csv(projection_table(saved$surface_prevalence, coop), extra_paths[1])
    atomic_write_csv(projection_table(saved$expected_prevalence, coop), extra_paths[2])
    for (nm in c("fixed", "hyperparameters")) {
      tab <- saved[[nm]]
      stem <- if (nm == "fixed") "Fixed_effects" else "Hyperparameters"
      atomic_write_csv(data.frame(parameter = rownames(tab), tab, row.names = NULL),
        file.path(output_dirs[["diagnostics"]], paste0(stem, "_model_", draw, ".csv")))
    }
    nodes <- saved$spatial_nodes
    mesh_rows <- node_mesh_rows(nodes$ID, prepared$mesh$n)
    atomic_write_csv(data.frame(x_m = prepared$mesh$loc[mesh_rows, 1],
      y_m = prepared$mesh$loc[mesh_rows, 2], nodes),
      file.path(output_dirs[["diagnostics"]], paste0("Spatial_nodes_model_", draw, ".csv")))
    write_uncertainty_diagnostics(saved, draw, y)
    atomic_write_csv(saved$monte_carlo, extra_paths[4])
    atomic_write_csv(saved$posterior_calibration, extra_paths[5])
    atomic_write_csv(saved$integration_check, file.path(output_dirs[["diagnostics"]],
      paste0("Integration_check_model_", draw, ".csv")))
    log_message("DRAW ", draw, " | Monte Carlo: ", saved$mc_check,
      " | max observed fitted-mean reconstruction difference ",
      round(max(abs(saved$posterior_calibration$mean_difference_pp)), 3), " pp")
    if (saved$mc_check != "SCREEN_PASSED") log_message(
      "Review Monte_Carlo_model_", draw, ".csv. To refine endpoints, increase posterior_samples and resume; selection is reused.")
    med_width <- function(pp) 100 * median(pp[["0.975quant"]] - pp[["0.025quant"]], na.rm = TRUE)
    log_message("DRAW ", draw, " | median 95% widths (pp): new location ",
      round(med_width(saved$prevalence), 2), "; IID=0 surface ", round(med_width(saved$surface_prevalence), 2),
      "; expected prevalence ", round(med_width(saved$expected_prevalence), 2))
    timing <- data.frame(draw = draw, donor_draw = assignment$donor_draw,
      input_column = iteration, mesh_vertices = prepared$mesh$n,
      observation_locations = length(y), prediction_locations = nrow(coop),
      candidate_search_fits = 0L, prediction_calculation_seconds = saved$elapsed_seconds,
      session_draw_seconds = elapsed_now() - draw_start,
      selected_covariates = paste(covariates, collapse = " + "),
      waic_selection = NA_real_, waic_posterior = saved$waic_observations,
      posterior_samples = saved$posterior_samples, monte_carlo_status = saved$mc_check)
    published <- file.path(cfg$random_expected_dir, paste0("Projections_model_", draw, ".csv"))
    random_publish(extra_paths[2], published)
    # The completion marker is written last. Partial exports never count as complete.
    atomic_save_rds(list(timing = timing, prediction_signature = prediction_signature,
      files = c(result_paths, extra_paths, published)), done_file)
    if (!draw %in% cfg$save_fit_draws) unlink(fit_checkpoint)
    if (!draw %in% cfg$save_joint_draws) unlink(joint_checkpoint)
    log_message("DRAW ", draw, " COMPLETE | this session's draw work ", round(timing$session_draw_seconds / 60, 2), " min")
    log_message("Expected prevalence published: ", published)
    timing
  }

  refresh_timings <- function() {
    files <- list.files(output_dirs[["checkpoints"]], pattern = "^completed\\.rds$", recursive = TRUE, full.names = TRUE)
    if (length(files)) {
      timing <- do.call(rbind, lapply(files, function(f) readRDS(f)$timing))
      timing <- timing[order(timing$draw), , drop = FALSE]
      atomic_write_csv(timing, timing_path)
    }
  }
  for (iteration in iteration_columns) {
    outcome <- tryCatch({run_one_draw(iteration); "completed"},
      aat_pause = function(e) {
        log_message("PAUSED during draw ", iteration - 7L, ": ", conditionMessage(e))
        "paused"
      }, error = function(e) {
        log_message("STOPPED during draw ", iteration - 7L, ": ", conditionMessage(e))
        log_message("Completed predictions and assigned sets are saved. Run folder: ", output_root)
        stop(e)
      })
    refresh_timings()
    if (outcome == "paused") {
      log_message("Source this script again to continue automatically from the saved assignments.")
      log_message("No incomplete draw is marked complete. A stage timeout may need a larger stage limit.")
      return(invisible(output_root))
    }
    if (nzchar(cfg$reference_output_dir)) tryCatch(compare_run_outputs(output_root, cfg$reference_output_dir,
      iteration - 7L, file.path(output_dirs[["diagnostics"]], paste0("Reference_comparison_model_", iteration - 7L, ".csv"))),
      error = function(e) log_message("Reference comparison unavailable: ", conditionMessage(e)))
  }
  log_message("REQUESTED DRAWS COMPLETE. Review Prediction_targets, Monte_Carlo and Posterior_calibration diagnostics.")
  invisible(output_root)
}

run_random_eth(settings, random_settings)
