# ===========================================================================
# AAT ETHIOPIA — 200-KM BLOCKED CROSS-VALIDATION OF REGIONAL SITE AVERAGES
# SEPARATE analysis: retains the original 100-km SITE-LEVEL validation script.
#
# PRIMARY VALIDATION: observed mean diagnostic-adjusted prevalence across >=3
# distinct held-out sites in each ~200-km geographic block. Prediction is the
# posterior mean of the expected new-site prevalence at those SAME locations.
# A joint-posterior, new-site-IID and binomial simulation supplies the predictive
# distribution for that observed average (do NOT average site-level intervals).
# Equal weight per site is primary; cattle/sample-size weighting is sensitivity.
# These are averages of SAMPLED populations, not design-unbiased areal prevalence.
# Secondary site-level scores remain available to document fine-scale limits.
# Training folds still rebuild mesh/scaling and run training-only WAIC selection.
#
# IMPORTANT: Precomputed diagnostic-corrected counts used the full raw dataset.
# This is NOT independent end-to-end validation of the diagnostic model.
# Requires original model script and run checkpoint; run from project directory.
# ===========================================================================

cv <- list(
  # Blank paths are discovered from this script's directory and project root.
  # original_script: current ETH production engine (not Model_selection_ETH.r).
  # source_run_dir: ETH run with Checkpoints_ETH/spatial_preparation.rds.
  # Set source_run_dir explicitly to pin a run when several ETH runs exist.
  project_root = "",
  original_script = "",
  source_run_dir = "",
  data_path = "",            # Blank uses the original run's saved CSV path.
  output_dir = "CV_ETH_200km_block_mean_validation",

  # Pilot with draw 1, then extend to e.g. 1:3 without changing the folds.
  # Draw ID 1 corresponds to input CSV column 8, ID 2 to column 9, etc.
  draw_ids = 1L,
  folds = 5L,
  block_width_km = 200,       # 200-km geographic blocks in metre-based AEQD CRS
  exclusion_buffer_km = 25, # Keep the prior spatial exclusion for comparability
  min_training_sites = 30L,
  min_test_sites = 3L,
  min_evaluation_sites_per_block = 3L,  # No 1-site "regional averages"
  posterior_samples = 2048L,  # >= 1024 and even, as in original implementation
  predictive_simulations = 4096L,
  iid_quadrature_nodes = 31L,
  seed = 20261001L,

  # The primary score gives each eligible block EQUAL weight irrespective of
  # its number of sites; reports additional cattle-weighted sensitivity results.
  # Retain the site-level comparison for audit; primary CSVs now score BLOCKS.
  # Model-comparison options; selection ALWAYS uses training data only.
  compare_spatial_null = TRUE,  # Intercept + same SPDE and location IID
  compare_no_iid = FALSE,      # Optional: selected terms + SPDE, without IID
  log_to_console = TRUE
)


run_eth_cv <- function(cv) {

# The 200-km experiment must live alongside, never replace, prior analyses.
# Prefix every exported file so copied results remain identifiable outside
# their dedicated output directory.
block_output <- function(filename) {
  file.path(cv$output_dir, paste0("block200_", filename))
}
if (basename(sub("/+$", "", cv$output_dir)) %in%
    c("CV_ETH_spatial_blocks", "CV_ETH_spatial_blocks_200km")) {
  stop("200-km CV must use its own NEW output_dir, not a previous CV output folder.")
}

# Paths can be left blank for discovery, or set explicitly in cv above.
# A completed ETH production run is required, just as a completed NGA run
# was required by the supplied scripts. Model_selection_ETH.r is a data-loading
# reference, not a replacement for the production helper API/checkpoint.
eth_script_directory <- function() {
  paths <- lapply(sys.frames(), function(x) x$ofile)
  paths <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x), paths)
  if (length(paths)) return(dirname(normalizePath(tail(paths, 1L)[[1L]], mustWork = TRUE)))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1L]), mustWork = TRUE)))
  getwd()
}
eth_ancestors <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  out <- path
  while (!identical(dirname(path), path)) {
    path <- dirname(path)
    out <- c(out, path)
  }
  out
}
eth_resolve_existing <- function(path, roots, directory = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("Expected one non-empty path.")
  path <- path.expand(path)
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", path)
  candidates <- if (absolute) path else unique(c(file.path(roots, path), path))
  ok <- if (directory) dir.exists(candidates) else file.exists(candidates) & !dir.exists(candidates)
  if (!any(ok)) stop("Cannot find ", path, "\nSearched: ", paste(candidates, collapse = "\n  "))
  normalizePath(candidates[which(ok)[1L]], mustWork = TRUE)
}
eth_project_root <- function(explicit, roots) {
  if (nzchar(explicit)) return(eth_resolve_existing(explicit, roots, directory = TRUE))
  candidates <- unique(unlist(lapply(roots, eth_ancestors)))
  ok <- vapply(candidates, function(p)
    dir.exists(file.path(p, "Code", "TestSensSpec")) &&
      dir.exists(file.path(p, "Code", "Prevalence")), logical(1))
  if (!any(ok)) stop("Cannot locate the DIMCAT project root. Set cv$project_root to the folder containing Code/.")
  candidates[which(ok)[1L]]
}
eth_find_model <- function(explicit, roots) {
  if (nzchar(explicit)) return(eth_resolve_existing(explicit, roots))
  canonical <- unique(file.path(roots, "AAT_spatial_quick_joint_posterior_ETH.R"))
  canonical <- canonical[file.exists(canonical)]
  if (length(canonical)) return(normalizePath(canonical[1L], mustWork = TRUE))
  candidates <- unique(unlist(lapply(roots, function(p)
    list.files(p, pattern = "^AAT_spatial_quick_joint_posterior_ETH.*[.]r$",
               full.names = TRUE, ignore.case = TRUE))))
  if (length(candidates) != 1L)
    stop("Set cv$original_script to the current ETH production model script. Found ",
         length(candidates), " non-canonical candidates. Do not use Model_selection_ETH.r.")
  normalizePath(candidates, mustWork = TRUE)
}
eth_find_run <- function(explicit, roots) {
  complete <- function(p) all(file.exists(file.path(p,
    c("Checkpoints_ETH/spatial_preparation.rds", "run_manifest.rds", "settings_used.rds"))))
  if (nzchar(explicit)) {
    p <- eth_resolve_existing(explicit, roots, directory = TRUE)
    if (!complete(p)) stop("ETH source_run_dir must contain Checkpoints_ETH/spatial_preparation.rds, run_manifest.rds and settings_used.rds.")
    return(p)
  }
  candidates <- unique(unlist(lapply(roots, function(p)
    list.dirs(p, full.names = TRUE, recursive = FALSE))))
  candidates <- candidates[grepl("^Spatial_joint_ETH", basename(candidates))]
  candidates <- candidates[vapply(candidates, complete, logical(1))]
  if (!length(candidates)) stop("No compatible ETH preparation was found. Set cv$source_run_dir to your completed ETH run (with its preparation checkpoint).")
  # Use the most recently saved preparation, not a folder's mutable mtime.
  stamp <- file.info(file.path(candidates, "Checkpoints_ETH", "spatial_preparation.rds"))$mtime
  candidates <- candidates[order(stamp, candidates, decreasing = TRUE)]
  if (length(candidates) > 1L)
    message("Multiple ETH runs found; using newest preparation. Set cv$source_run_dir explicitly to choose another run.")
  normalizePath(candidates[1L], mustWork = TRUE)
}
eth_coordinate_key <- function(lon, lat) paste(sprintf("%.17g", lon), sprintf("%.17g", lat), sep = "|")
eth_prepare_input <- function(raw, model, mapping) {
  required_columns <- c("longitude", "latitude", "sample_size")
  if (!all(required_columns %in% names(raw))) stop("Input lacks longitude, latitude or sample_size.")
  for (nm in required_columns) raw[[nm]] <- suppressWarnings(as.numeric(as.character(raw[[nm]])))
  raw <- raw[is.finite(raw$longitude) & is.finite(raw$latitude), , drop = FALSE]
  if (!nrow(raw)) stop("No finite input coordinates.")
  if (is.function(model$filter_ethiopia_inputs)) {
    # Match the current production engine's existing orientation/filter logic.
    raw <- model$filter_ethiopia_inputs(raw)
  } else {
    # Same orientation test and ETH bbox as the supplied Model_selection_ETH.r.
    # Coordinate orientation is not otherwise changed by this CV adaptation.
    lat_ok <- sum(raw$latitude >= 3 & raw$latitude <= 15)
    lon_ok <- sum(raw$longitude >= 33 & raw$longitude <= 48)
    lat_as_lon <- sum(raw$latitude >= 33 & raw$latitude <= 48)
    lon_as_lat <- sum(raw$longitude >= 3 & raw$longitude <= 15)
    if (lat_as_lon > lon_ok && lon_as_lat > lat_ok) {
      tmp <- raw$longitude; raw$longitude <- raw$latitude; raw$latitude <- tmp
      message("Applied the ETH reference script's existing coordinate orientation rule.")
    }
  }
  keep <- is.finite(raw$longitude) & is.finite(raw$latitude) &
    raw$longitude >= 33 & raw$longitude <= 48 & raw$latitude >= 3 & raw$latitude <= 15
  raw <- raw[keep, , drop = FALSE]
  if (!nrow(raw)) stop("No records remain in the ETH input bounding box.")
  if (any(!is.finite(raw$sample_size) | raw$sample_size <= 0 |
          abs(raw$sample_size - round(raw$sample_size)) > 1e-8))
    stop("ETH input sample_size must contain positive integers.")
  # Bbox is a prefilter only. The saved production mapping/domain decides the
  # final sites, so nearby-country records cannot enter the fitted ETH dataset.
  if (!all(required_columns %in% names(mapping))) stop("Saved coordinate_mapping lacks required columns.")
  key <- eth_coordinate_key(raw$longitude, raw$latitude)
  target <- eth_coordinate_key(mapping$longitude, mapping$latitude)
  if (anyDuplicated(target)) stop("Duplicate coordinates in the saved production mapping.")
  raw <- raw[key %in% target, , drop = FALSE]
  key <- eth_coordinate_key(raw$longitude, raw$latitude)
  totals <- tapply(raw$sample_size, key, sum)
  aligned <- as.numeric(totals[target])
  if (length(aligned) != nrow(mapping) || anyNA(aligned) ||
      !isTRUE(all.equal(aligned, as.numeric(mapping$sample_size), tolerance = 1e-8)))
    stop("ETH input coordinates/sample totals do not match the saved preparation. Use its original CSV and coordinate-filtering code.")
  raw
}

# ------------------------ Preparatory checks -------------------------------
launch_roots <- unique(c(eth_script_directory(), getwd()))
launch_wd <- getwd()
on.exit(setwd(launch_wd), add = TRUE)
cv$project_root <- eth_project_root(cv$project_root, launch_roots)
analysis_dir <- file.path(cv$project_root, "Code", "Prevalence", "Bovine BCT and PCR")
search_roots <- unique(c(analysis_dir, launch_roots, cv$project_root))
cv$original_script <- eth_find_model(cv$original_script, search_roots)
cv$source_run_dir <- eth_find_run(cv$source_run_dir, unique(c(dirname(cv$original_script), search_roots)))
if (!grepl("^(/|[A-Za-z]:[/\\\\])", path.expand(cv$output_dir)))
  cv$output_dir <- file.path(analysis_dir, cv$output_dir)
cv$output_dir <- normalizePath(path.expand(cv$output_dir), mustWork = FALSE)
if (grepl("NGA|Nigeria", basename(cv$output_dir), ignore.case = TRUE))
  stop("ETH outputs must not be written to a Nigeria output folder.")
cat("\nETH model helpers: ", cv$original_script, "\nETH preparation: ", cv$source_run_dir,
    "\nETH CV outputs: ", cv$output_dir, "\n", sep = "")
setwd(cv$project_root)
required <- c("INLA", "sf", "Matrix", "matrixStats")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(INLA); library(Matrix)})
if (cv$folds < 2L || cv$block_width_km <= 0 || cv$exclusion_buffer_km < 0 ||
    cv$posterior_samples < 1024L || (cv$posterior_samples %% 2L) != 0L ||
    cv$predictive_simulations < 500L || cv$iid_quadrature_nodes < 3L)
  stop("Invalid CV configuration.")
if (anyNA(cv$draw_ids) || !length(cv$draw_ids) || any(!is.finite(cv$draw_ids)) ||
    any(cv$draw_ids < 1L | cv$draw_ids != floor(cv$draw_ids)) || anyDuplicated(cv$draw_ids))
  stop("draw_ids must be distinct positive integers.")
if (!is.null(cv$min_evaluation_sites_per_block) && cv$min_evaluation_sites_per_block < 2L)
  stop("min_evaluation_sites_per_block must be at least 2.")

# Load the production helper definitions without starting the full model run.
original_expr <- parse(file = cv$original_script, keep.source = FALSE)
last_text <- paste(deparse(tail(original_expr, 1L)[[1L]]), collapse = "")
if (!identical(gsub("[[:space:]]", "", last_text), "run_aat_quick(settings)"))
  stop("ETH production file must end in run_aat_quick(settings). Review the helper loader if its API changed.")
model <- new.env(parent = globalenv())
for (expr in head(original_expr, -1L)) eval(expr, envir = model)
needed_fns <- c("coordinate_key", "aggregate_draw_counts", "audit_covariates",
  "audit_quadratic_basis", "select_per_draw", "select_best_visited", "fit_observations",
  "sample_joint_posterior", "build_model_stacks", "transform_model_covariates",
  "effect_terms", "normal_quadrature")
missing_fns <- needed_fns[!vapply(needed_fns, function(nm) is.function(model[[nm]]), logical(1))]
if (length(missing_fns)) stop("ETH production helper API is incompatible; missing: ", paste(missing_fns, collapse = ", "))

prep_path <- file.path(cv$source_run_dir, "Checkpoints_ETH", "spatial_preparation.rds")
manifest_path <- file.path(cv$source_run_dir, "run_manifest.rds")
source_settings_path <- file.path(cv$source_run_dir, "settings_used.rds")
base <- readRDS(prep_path)
manifest <- readRDS(manifest_path)
original_settings <- readRDS(source_settings_path)
needed_prep <- c("coordinate_mapping", "kept_location_rows", "observation_lonlat",
  "observation_m", "prediction_m", "raw_observation_covariates", "numtrials")
if (!all(needed_prep %in% names(base))) stop("Cached ETH preparation lacks: ", paste(setdiff(needed_prep, names(base)), collapse = ", "))
if (!is.null(manifest$country) &&
    !all(toupper(as.character(manifest$country)) %in% c("ETH", "ETHIOPIA")))
  stop("The source manifest is not for Ethiopia.")
if (!isTRUE(manifest$use_location_iid) || !isTRUE(original_settings$use_location_iid))
  stop("Select an ETH production run with the location IID effect enabled.")
degree <- manifest$covariate_degree
if (length(degree) != 1L || is.na(degree) || !degree %in% c(1L, 2L))
  stop("Saved covariate_degree must be 1 (linear) or 2 (quadratic).")
if (!isTRUE(original_settings$covariate_degree == degree))
  stop("Saved settings and manifest disagree on covariate_degree.")
# The helper's historical name audit_quadratic_basis also supports degree=1.
# Never introduce quadratic effects into a saved linear production model.
crs <- sf::st_crs(manifest$spatial_crs)
if (is.na(crs) || isTRUE(sf::st_is_longlat(crs)) ||
    !isTRUE(tolower(crs$units_gdal) %in% c("metre", "meter", "metres", "meters", "m")))
  stop("The ETH preparation must use a projected CRS in metres.")
needed_manifest <- c("cutoff_m", "offset_m", "max_edge_m", "range_prior_m", "sigma_prior", "input_md5")
if (!all(needed_manifest %in% names(manifest))) stop("Missing saved mesh/prior/input settings.")
input_path <- if (nzchar(cv$data_path)) cv$data_path else original_settings$data_path
input_path <- eth_resolve_existing(input_path, unique(c(cv$project_root, search_roots)))
input_md5 <- unname(tools::md5sum(input_path))
if (!identical(input_md5, manifest$input_md5)) stop("Adjusted-count CSV differs from the ETH preparation's original input (MD5 mismatch).")
original_settings$data_path <- input_path
raw <- eth_prepare_input(read.csv(input_path, check.names = FALSE), model, base$coordinate_mapping)
if (any(cv$draw_ids + 7L > ncol(raw))) stop("Requested draw ID is not in the adjusted-count CSV.")
n_sites <- nrow(base$observation_m)
if (!n_sites || nrow(base$observation_lonlat) != n_sites ||
    nrow(base$raw_observation_covariates) != n_sites || length(base$numtrials) != n_sites)
  stop("Misaligned source preparation.")
ll <- as.matrix(base$observation_lonlat)
if (ncol(ll) != 2L || any(!is.finite(ll)) ||
    any(ll[, 1L] < 33 | ll[, 1L] > 48 | ll[, 2L] < 3 | ll[, 2L] > 15))
  stop("Saved observation coordinates are not in the ETH domain/orientation.")
expected_xy <- sf::st_coordinates(sf::st_transform(
  sf::st_as_sf(data.frame(lon = ll[, 1L], lat = ll[, 2L]), coords = c("lon", "lat"), crs = 4326), crs))
if (!identical(dim(expected_xy), dim(base$observation_m)) ||
    any(!is.finite(base$observation_m)) || max(abs(expected_xy - base$observation_m)) > 1)
  stop("Saved projected coordinates do not match the ETH geographic coordinates and CRS.")
kept <- base$kept_location_rows
if (length(kept) != n_sites || anyNA(kept) || anyDuplicated(kept) ||
    any(kept < 1L | kept > nrow(base$coordinate_mapping))) stop("Invalid kept_location_rows.")
mapping_ll <- base$coordinate_mapping[kept, c("longitude", "latitude"), drop = FALSE]
if (!identical(eth_coordinate_key(ll[, 1L], ll[, 2L]),
               eth_coordinate_key(mapping_ll$longitude, mapping_ll$latitude)))
  stop("Prepared ETH observation order differs from the coordinate mapping.")
if (!identical(as.numeric(base$numtrials), as.numeric(base$coordinate_mapping$sample_size[kept])))
  stop("Prepared ETH sample sizes differ from the coordinate mapping.")
cat("Verified ETH input: ", nrow(raw), " records; ", n_sites,
    " fitted sites; covariate degree ", degree, "; projected metre CRS; location IID enabled.\n", sep = "")


# ----------------- Immutable spatial block assignment -----------------------
# The original preparation has already collapsed identical coordinates. No
# exact coordinate can be assigned to different folds. Fold assignments do not
# use case counts or environmental covariate values.
xy <- base$observation_m
block_x <- floor((xy[, 1L] - min(xy[, 1L])) / (1000 * cv$block_width_km))
block_y <- floor((xy[, 2L] - min(xy[, 2L])) / (1000 * cv$block_width_km))
block_id <- paste(block_x, block_y, sep = "_")
blocks <- split(seq_len(n_sites), block_id)
if (length(blocks) < cv$folds)
  stop("Too few occupied spatial blocks. Reduce block_width_km or folds.")
set.seed(cv$seed)
block_names <- sample(names(blocks))
block_names <- block_names[order(-lengths(blocks[block_names]))]
fold_load <- integer(cv$folds)
block_fold <- setNames(integer(length(blocks)), names(blocks))
for (bn in block_names) {
  candidates <- which(fold_load == min(fold_load))
  f <- candidates[sample.int(length(candidates), 1L)]
  block_fold[bn] <- f
  fold_load[f] <- fold_load[f] + length(blocks[[bn]])
}
fold_id <- unname(block_fold[block_id])
if (anyNA(fold_id) || any(tabulate(fold_id, cv$folds) < cv$min_test_sites))
  stop("Empty/undersized fold. Change folds or block_width_km.")

# Cache identity: adding further diagnostic draws is allowed, but changing
# blocks, source model, or scoring specification requires a NEW output folder.
cv_identity <- list(version = "eth_spatial_block_mean_cv_20261023_v1",
  country = "ETH", covariate_degree = degree,
  source_run_dir = cv$source_run_dir,
  manifest_md5 = unname(tools::md5sum(manifest_path)),
  settings_md5 = unname(tools::md5sum(source_settings_path)),
  original_script_md5 = unname(tools::md5sum(cv$original_script)),
  preparation_md5 = unname(tools::md5sum(prep_path)), input_md5 = input_md5,
  block_id = block_id, fold_id = fold_id, folds = cv$folds,
  buffer_km = cv$exclusion_buffer_km, block_width_km = cv$block_width_km,
  min_training_sites = cv$min_training_sites, min_test_sites = cv$min_test_sites,
  min_evaluation_sites_per_block = cv$min_evaluation_sites_per_block,
  block_averaging = c("equal_sites", "cattle_weighted"),
  posterior_samples = cv$posterior_samples,
  predictive_simulations = cv$predictive_simulations,
  quadrature_nodes = cv$iid_quadrature_nodes, seed = cv$seed,
  compare_spatial_null = cv$compare_spatial_null,
  compare_no_iid = cv$compare_no_iid)
# An existing unrelated directory without our identity must not be overwritten.
# Completed checkpoints from exactly this analysis may still be resumed.
if (dir.exists(cv$output_dir) &&
    length(list.files(cv$output_dir, all.files = TRUE, no.. = TRUE)) > 0L &&
    !file.exists(block_output("cv_identity.rds"))) {
  stop("Output directory is not empty and has no matching 200-km identity. Choose a fresh output_dir.")
}
dir.create(cv$output_dir, showWarnings = FALSE, recursive = TRUE)
identity_file <- block_output("cv_identity.rds")
if (file.exists(identity_file)) {
  if (!identical(readRDS(identity_file), cv_identity))
    stop("This output folder belongs to a different CV specification. Use a NEW output_dir.")
} else saveRDS(cv_identity, identity_file)

block_assignment <- data.frame(site_id = seq_len(n_sites),
  longitude = base$observation_lonlat[, 1L], latitude = base$observation_lonlat[, 2L],
  x_m = xy[, 1L], y_m = xy[, 2L], block_id = block_id, fold = fold_id)
write.csv(block_assignment, block_output("spatial_fold_assignments.csv"), row.names = FALSE)

# Build one genuine training fold: entire held-out blocks excluded from ALL
# parameter estimation, mesh vertex insertion, scaling, quadratic eligibility
# AND model selection. Holdout covariates are used solely for prediction.
make_training_fold <- function(fold) {
  test <- which(fold_id == fold)
  train <- which(fold_id != fold)
  min_d2 <- rep(Inf, length(train))
  for (j in test) {
    min_d2 <- pmin(min_d2,
      (xy[train, 1L] - xy[j, 1L])^2 + (xy[train, 2L] - xy[j, 2L])^2)
  }
  excluded <- train[sqrt(min_d2) < 1000 * cv$exclusion_buffer_km]
  train <- setdiff(train, excluded)
  if (length(train) < cv$min_training_sites || length(test) < cv$min_test_sites)
    stop("Fold ", fold, ": insufficient train/test sites after buffer. Training=",
      length(train), ", test=", length(test), ". Reduce buffer/folds if justified.")

  out <- base
  out$observation_m <- base$observation_m[train, , drop = FALSE]
  out$observation_lonlat <- base$observation_lonlat[train, , drop = FALSE]
  out$raw_observation_covariates <- base$raw_observation_covariates[train, , drop = FALSE]
  out$numtrials <- base$numtrials[train]
  out$use_location_iid <- original_settings$use_location_iid
  out$iid_prior <- original_settings$iid_prior
  out$support_tolerance <- original_settings$support_tolerance

  # Mesh built from training observations and the already-defined study region.
  # The study boundary/grid contain no held-out outcomes; no held-out sites are
  # passed as observation mesh vertices.
  domain <- rbind(out$observation_m, base$prediction_m)
  domain <- domain[chull(domain[, 1L], domain[, 2L]), , drop = FALSE]
  out$mesh <- INLA::inla.mesh.2d(loc = out$observation_m, loc.domain = domain,
    cutoff = manifest$cutoff_m, offset = manifest$offset_m,
    max.edge = manifest$max_edge_m)
  out$spde <- INLA::inla.spde2.pcmatern(mesh = out$mesh, alpha = 2,
    prior.range = manifest$range_prior_m,
    prior.sigma = manifest$sigma_prior, constr = TRUE)
  out$A <- INLA::inla.spde.make.A(mesh = out$mesh, loc = out$observation_m)
  out$Ap <- NULL
  out$raw_prediction_covariates <- base$raw_observation_covariates[test, , drop = FALSE]
  out$covariate_audit <- model$audit_covariates(
    out$raw_observation_covariates, out$raw_prediction_covariates,
    original_settings$support_tolerance)
  invalid <- out$covariate_audit$invalid_fraction_units
  if (any(invalid)) stop("Invalid 0–1 fractional predictors in fold ", fold)
  partly_missing <- with(out$covariate_audit, covariate[eligible & observed_missing > 0])
  if (length(partly_missing)) stop("Missing training predictors: ", paste(partly_missing, collapse=", "))
  out$candidate_covariates <- out$covariate_audit$covariate[out$covariate_audit$eligible]
  out$quadratic_basis <- model$audit_quadratic_basis(out$raw_observation_covariates,
    out$covariate_audit, original_settings$covariate_degree)
  out <- model$build_model_stacks(out)
  # Held-out projection never includes an observed site-IID index.
  A_test <- INLA::inla.spde.make.A(mesh = out$mesh, loc = xy[test, , drop = FALSE])
  dtest <- sqrt(vapply(test, function(j) min((xy[train,1L]-xy[j,1L])^2 +
                                     (xy[train,2L]-xy[j,2L])^2), numeric(1))) / 1000
  list(prepared = out, train = train, test = test, excluded = excluded,
       A_test = A_test, test_covariates = base$raw_observation_covariates[test, , drop=FALSE],
       distance_to_training_km = dtest)
}

# ---------------------- Predictive scoring methods --------------------------
log_sum_exp <- function(z) {
  zmax <- max(z)
  if (!is.finite(zmax)) return(zmax)
  zmax + log(sum(exp(z - zmax)))
}

# The full predictive distribution is:
#   p(y_hold | train) = E_posterior E_new_IID Binomial(y_hold | n, logit^-1(eta+IID)).
# Gaussian quadrature integrates the new-site IID for binomial log-score and
# randomized PIT; separate draws approximate predictive count intervals.
score_new_site <- function(eta, sigma, y, n, rule, sims, random_seed) {
  stopifnot(length(eta) == length(sigma), n >= 1L, y >= 0L, y <= n)
  if (any(!is.finite(eta)) || any(!is.finite(sigma)) || any(sigma < 0))
    stop("Non-finite posterior predictor/IID SD.")
  ns <- length(eta)
  nodes <- length(rule$z)
  log_terms <- matrix(NA_real_, nrow=nodes, ncol=ns)
  lower_terms <- upper_terms <- numeric(nodes)
  mean_terms <- numeric(nodes)
  for (q in seq_len(nodes)) {
    p <- plogis(eta + sigma * rule$z[q])
    log_terms[q, ] <- dbinom(y, size=n, prob=p, log=TRUE) + log(rule$w[q])
    lower_terms[q] <- rule$w[q] * mean(pbinom(y-1L, n, p))
    upper_terms[q] <- rule$w[q] * mean(pbinom(y, n, p))
    mean_terms[q] <- rule$w[q] * mean(p)
  }
  log_score <- log_sum_exp(as.vector(log_terms)) - log(ns)
  cdf_less <- sum(lower_terms)
  cdf_leq <- sum(upper_terms)
  mean_p <- sum(mean_terms)

  set.seed(as.integer(random_seed %% .Machine$integer.max))
  sample_id <- sample.int(ns, size=sims, replace=TRUE)
  p_rep <- plogis(eta[sample_id] + sigma[sample_id] * rnorm(sims))
  y_rep <- rbinom(sims, n, p_rep)
  qi80 <- as.numeric(quantile(y_rep, c(.10, .90), names=FALSE, type=1))
  qi95 <- as.numeric(quantile(y_rep, c(.025, .975), names=FALSE, type=1))
  qpi95 <- as.numeric(quantile(p_rep, c(.025, .975), names=FALSE, type=8))
  list(mean_p=mean_p, log_score=log_score,
    count_lo80=qi80[1L], count_hi80=qi80[2L],
    count_lo95=qi95[1L], count_hi95=qi95[2L],
    p_lo95=qpi95[1L], p_hi95=qpi95[2L],
    pit=cdf_less + runif(1L) * max(0, cdf_leq-cdf_less))
}

# -------------------- Joint block-mean predictive scoring -------------------
# A region's observed target is the mean of adjusted y_i / n_i at its surveyed
# sites (not the national area-weighted mean). In every predictive replicate all
# sites share the SAME joint posterior draw for beta, SPDE and hyperparameters,
# but get INDEPENDENT new-location IID and binomial draws. This respects spatial
# posterior correlation and the variability of means of several sites.
# Posterior "expected" block prevalence integrates out each site's new IID;
# predictive block intervals additionally include IID and binomial noise.
# Separate simulations use the actual held-out-site locations and trial sizes.

crps_from_samples <- function(x, observed) {
  x <- sort(as.numeric(x))
  if (!length(x) || any(!is.finite(x)) || !is.finite(observed))
    stop("Invalid block CRPS inputs.")
  n <- length(x)
  mean(abs(x - observed)) -
    sum((2 * seq_len(n) - n - 1) * x) / n^2
}

block_predictions_from_joint <- function(eta_by_site, sigma, rule,
    site_indices, site_ids, y_all, draw, fold, model_name, seed_offset) {
  ns <- length(sigma)
  sims <- as.integer(cv$predictive_simulations)
  n_site <- length(site_indices)
  ids <- site_ids[site_indices]
  trials <- as.numeric(base$numtrials[ids])
  observed <- y_all[ids] / trials
  stopifnot(n_site >= cv$min_evaluation_sites_per_block,
            all(is.finite(observed)), all(trials > 0))

  # Each held-out block receives one posterior-index draw for ALL of its sites.
  # Only location-level IID and future sample counts vary independently by site.
  block_seed <- (as.double(cv$seed) + 1000003 * draw + 100003 * fold +
    sum(as.double(ids) * 271) + seed_offset) %% (.Machine$integer.max - 1)
  set.seed(as.integer(block_seed + 1))
  sampled_post <- sample.int(ns, size=sims, replace=TRUE)
  predicted_site_fraction <- matrix(NA_real_, nrow=n_site, ncol=sims)
  expected_site <- matrix(NA_real_, nrow=n_site, ncol=ns)
  for (a in seq_len(n_site)) {
    e <- eta_by_site[[site_indices[a]]]
    if (length(e) != ns || any(!is.finite(e)))
      stop("Joint prediction samples are missing for a held-out block site.")
    # Exact same quadrature method across sites; joint posterior index preserved.
    m <- numeric(ns)
    for (q in seq_along(rule$z))
      m <- m + rule$w[q] * plogis(e + sigma * rule$z[q])
    expected_site[a, ] <- m
    # New exchangeable site, NOT the training site's estimated IID effect.
    p_new <- plogis(e[sampled_post] + sigma[sampled_post] * rnorm(sims))
    predicted_site_fraction[a, ] <- rbinom(sims, size=trials[a], prob=p_new) / trials[a]
  }
  b <- block_id[ids[1L]]
  if (any(block_id[ids] != b)) stop("Cross-block aggregation is forbidden.")
  v <- list(equal_sites = rep(1 / n_site, n_site),
            cattle_weighted = trials / sum(trials))
  do.call(rbind, lapply(names(v), function(kind) {
    weight <- v[[kind]]
    pred <- as.numeric(crossprod(weight, predicted_site_fraction))
    expected <- as.numeric(crossprod(weight, expected_site))
    obs <- sum(weight * observed)
    qpred <- as.numeric(quantile(pred, c(.10, .90, .025, .975),
                                 type=8, names=FALSE))
    qexpected <- as.numeric(quantile(expected, c(.025, .975),
                                     type=8, names=FALSE))
    data.frame(draw=draw, fold=fold, model=model_name, block_id=b,
      averaging=kind, n_sites=n_site, n_cattle=sum(trials),
      site_id_list=paste(ids, collapse="+"),
      observed_block_mean=obs, posterior_mean_block_expected=mean(expected),
      block_crps=crps_from_samples(pred, obs),
      predictive_lo80=qpred[1L], predictive_hi80=qpred[2L],
      predictive_lo95=qpred[3L], predictive_hi95=qpred[4L],
      expected_prevalence_lo95=qexpected[1L],
      expected_prevalence_hi95=qexpected[2L],
      predictive_width95=qpred[4L]-qpred[3L],
      stringsAsFactors=FALSE)
  }))
}

score_blocks_for_fitted_model <- function(eta_by_site, sigma, rule,
                                         fold_data, y_all, draw, fold,
                                         model_name, seed_offset) {
  groups <- split(seq_along(fold_data$test), block_id[fold_data$test])
  eligible <- groups[lengths(groups) >= cv$min_evaluation_sites_per_block]
  if (!length(eligible))
    stop("Fold ",fold," has no held-out block with ",
      cv$min_evaluation_sites_per_block," independent sites.")
  do.call(rbind, lapply(eligible, function(local_index)
    block_predictions_from_joint(eta_by_site, sigma, rule,
      local_index, fold_data$test, y_all, draw, fold,
      model_name, seed_offset)))
}

score_blocks_for_constant <- function(fold_data, y_all, draw, fold) {
  tr <- fold_data$train
  p <- (sum(y_all[tr]) + 0.5) / (sum(base$numtrials[tr]) + 1)
  groups <- split(fold_data$test, block_id[fold_data$test])
  eligible <- groups[lengths(groups) >= cv$min_evaluation_sites_per_block]
  if (!length(eligible)) stop("No eligible blocks for constant baseline.")
  do.call(rbind, lapply(eligible, function(ids) {
    n <- as.numeric(base$numtrials[ids]); obs_i <- y_all[ids]/n
    set.seed(as.integer((as.double(cv$seed) + draw*1000003 +
       fold*100003 + sum(as.double(ids)*271) + 400) %%
       (.Machine$integer.max - 1)) + 1L)
    s <- cv$predictive_simulations
    pred_i <- vapply(seq_along(ids), function(j)
      rbinom(s, size=n[j], prob=p) / n[j], numeric(s))
    pred_i <- matrix(pred_i, nrow=s, ncol=length(ids))
    w <- list(equal_sites=rep(1/length(ids),length(ids)),
              cattle_weighted=n/sum(n))
    do.call(rbind,lapply(names(w),function(kind) {
      weight <- w[[kind]]
      pred <- as.numeric(pred_i %*% weight)
      obs <- sum(weight*obs_i)
      qpred <- as.numeric(quantile(pred,c(.10,.90,.025,.975),
                                   type=8,names=FALSE))
      data.frame(draw=draw, fold=fold, model="constant_binomial",
        block_id=block_id[ids[1]], averaging=kind,
        n_sites=length(ids), n_cattle=sum(n), site_id_list=paste(ids,collapse="+"),
        observed_block_mean=obs, posterior_mean_block_expected=p,
        block_crps=crps_from_samples(pred,obs),
        predictive_lo80=qpred[1],predictive_hi80=qpred[2],
        predictive_lo95=qpred[3],predictive_hi95=qpred[4],
        expected_prevalence_lo95=p,expected_prevalence_hi95=p,
        predictive_width95=qpred[4]-qpred[3], stringsAsFactors=FALSE)
    }))
  }))
}

# Validate one fit; chunking avoids a huge test-by-posterior prediction matrix.
score_fitted_model <- function(fit, prepared, covariates, fold_data,
                               y_all, draw, fold, model_name, unique_model_seed) {
  cfg <- original_settings
  cfg$posterior_samples <- as.integer(cv$posterior_samples)
  cfg$seed <- as.integer(cv$seed + 10111L * fold + 1009L * draw + unique_model_seed)
  post <- model$sample_joint_posterior(fit, prepared, covariates, cfg,
                                      draw = as.integer(draw))
  raw_test <- fold_data$test_covariates
  basis_test <- model$transform_model_covariates(raw_test, prepared)
  selected <- model$effect_terms(covariates, prepared)
  X <- cbind(b0=1, as.matrix(basis_test[, selected, drop=FALSE]))
  colnames(X) <- c("b0", selected)
  if (!identical(colnames(X), rownames(post$beta)))
    stop("Selected test design and posterior coefficients are misaligned.")
  valid <- apply(X, 1L, function(z) all(is.finite(z)))
  if (!all(valid))
    stop("Fold ",fold,": missing held-out predictors at site IDs ",
         paste(fold_data$test[!valid],collapse=","),". Compare common sites only after repair.")
  rule <- model$normal_quadrature(as.integer(cv$iid_quadrature_nodes))
  res <- vector("list", length(fold_data$test))
  eta_by_site <- vector("list", length(fold_data$test))
  for (j in which(valid)) {
    # Both fixed-effect and spatial effects are projected using the SAME
    # joint INLA posterior sample; new-site IID is never a trained-site IID.
    eta <- as.numeric(X[j,,drop=FALSE] %*% post$beta +
      fold_data$A_test[j,,drop=FALSE] %*% post$spatial)
    eta_by_site[[j]] <- eta
    id <- fold_data$test[j]
    y <- y_all[id]; n <- base$numtrials[id]
    v <- score_new_site(eta, post$sigma_iid, y, n, rule,
       as.integer(cv$predictive_simulations),
       as.double(cv$seed) + draw*1000000 + fold*100000 + id*10 + unique_model_seed)
    res[[j]] <- data.frame(draw=draw, fold=fold, model=model_name, site_id=id,
      block_id=block_id[id], longitude=base$observation_lonlat[id,1L],
      latitude=base$observation_lonlat[id,2L], n=n, y=y,
      adjusted_fraction=y/n, posterior_mean_p=v$mean_p,
      log_score=v$log_score, count_lo80=v$count_lo80, count_hi80=v$count_hi80,
      count_lo95=v$count_lo95, count_hi95=v$count_hi95,
      p_lo95=v$p_lo95, p_hi95=v$p_hi95, randomized_pit=v$pit,
      distance_to_training_km=fold_data$distance_to_training_km[j],
      selected_covariates=paste(covariates,collapse="+"))
  }
  if (!all(valid))
    stop("Fold ", fold, ": held-out predictors are missing for site IDs: ",
         paste(fold_data$test[!valid], collapse = ", "),
         ". Stop rather than compare models on different held-out site sets.")
  res <- Filter(Negate(is.null), res)
  if (!length(res)) stop("All test sites lack selected predictors in fold ", fold)
  list(site = do.call(rbind, res),
       block = score_blocks_for_fitted_model(eta_by_site, post$sigma_iid,
         rule, fold_data, y_all, draw, fold, model_name, unique_model_seed))
}

score_constant_baseline <- function(fold_data, y_all, draw, fold) {
  # Simple training-only constant binomial-prevalence reference; NOT an
  # uncertainty-corrected or spatial baseline. Jeffreys smoothing avoids 0/1.
  tr <- fold_data$train
  p <- (sum(y_all[tr])+.5)/(sum(base$numtrials[tr])+1)
  te <- fold_data$test
  set.seed(as.integer((as.double(cv$seed) + draw*1000000 + fold*100000 + 7)
                     %% .Machine$integer.max))
  data.frame(draw=draw, fold=fold, model="constant_binomial", site_id=te,
    block_id=block_id[te], longitude=base$observation_lonlat[te,1L],
    latitude=base$observation_lonlat[te,2L], n=base$numtrials[te], y=y_all[te],
    adjusted_fraction=y_all[te]/base$numtrials[te], posterior_mean_p=p,
    log_score=dbinom(y_all[te],base$numtrials[te],p,log=TRUE),
    count_lo80=qbinom(.10,base$numtrials[te],p),
    count_hi80=qbinom(.90,base$numtrials[te],p),
    count_lo95=qbinom(.025,base$numtrials[te],p),
    count_hi95=qbinom(.975,base$numtrials[te],p),
    p_lo95=p, p_hi95=p,
    randomized_pit=pbinom(y_all[te]-1,base$numtrials[te],p)+
      runif(length(te))*dbinom(y_all[te],base$numtrials[te],p),
    distance_to_training_km=fold_data$distance_to_training_km,
    selected_covariates="(training pooled prevalence)")
}

# --------------------- Training-only model selection ------------------------
fit_cfg <- original_settings
fit_cfg$local_work_root <- tempfile("aat_eth_cv_inla_", tmpdir=tempdir())
dir.create(fit_cfg$local_work_root)
on.exit(unlink(fit_cfg$local_work_root, recursive = TRUE), add = TRUE)
fit_cfg$diagnostics_dir <- cv$output_dir
fit_cfg$posterior_samples <- as.integer(cv$posterior_samples)
# Reuse existing prior/controls/INLA mode. No candidate or fit imported from
# the full-data run, because that WOULD contaminate cross-validation.

fit_cv <- function(fold_data, y_all, draw, fold, work) {
  p <- fold_data$prepared
  yt <- y_all[fold_data$train]
  candidates <- p$candidate_covariates
  cache_path <- file.path(work, "block200_training_candidate_waic.rds")
  records <- if (file.exists(cache_path)) readRDS(cache_path) else list()
  key_for <- function(covs) paste0("m_",paste(match(covs,candidates), collapse="_"))
  score <- function(covs) {
    key <- key_for(covs)
    if (is.null(records[[key]])) {
      cat("    training WAIC fit: ",key," [",paste(covs,collapse=","),"]\n",sep="")
      f <- model$fit_observations(p,yt,covs,fit_cfg)
      records[[key]] <<- list(covariates=covs,waic=f$waic$waic)
      saveRDS(records,cache_path)
      rm(f); gc(verbose=FALSE)
    }
    records[[key]]$waic
  }
  model$select_per_draw(candidates=candidates,score=score,
    backward_tolerance=fit_cfg$waic_removal_tolerance,
    addition_min_improvement=fit_cfg$waic_min_improvement,
    final_backward_tolerance=fit_cfg$waic_removal_tolerance)
  chosen <- model$select_best_visited(records,candidates)
  write.csv(data.frame(key=names(records),
    terms=vapply(records,function(z) paste(z$covariates,collapse="+"),character(1)),
    waic=vapply(records,function(z) z$waic,numeric(1))),
    file.path(work,"block200_training_candidate_waic.csv"),row.names=FALSE)
  cat("    Selected: ",paste(chosen$covariates,collapse="+"),
      "; WAIC=",format(chosen$waic,digits=8),"\n",sep="")
  selected_fit <- model$fit_observations(p,yt,chosen$covariates,fit_cfg)
  if (abs(selected_fit$waic$waic - chosen$waic) > fit_cfg$waic_recheck_tolerance)
    stop("Selected fit WAIC changed unexpectedly on refitting.")
  out <- list(score_fitted_model(selected_fit,p,chosen$covariates,fold_data,y_all,
    draw,fold,"selected_spatial_iid",100L))
  rm(selected_fit);gc(verbose=FALSE)

  if (isTRUE(cv$compare_spatial_null)) {
    f0 <- model$fit_observations(p,yt,character(),fit_cfg)
    out[[length(out)+1L]] <- score_fitted_model(f0,p,character(),fold_data,y_all,
      draw,fold,"intercept_spatial_iid",200L)
    rm(f0);gc(verbose=FALSE)
  }
  if (isTRUE(cv$compare_no_iid)) {
    p_no <- p
    p_no$use_location_iid <- FALSE
    p_no <- model$build_model_stacks(p_no)
    f_no <- model$fit_observations(p_no,yt,chosen$covariates,fit_cfg)
    out[[length(out)+1L]] <- score_fitted_model(f_no,p_no,chosen$covariates,
      fold_data,y_all,draw,fold,"selected_spatial_no_iid",300L)
    rm(f_no);gc(verbose=FALSE)
  }
  site_out <- do.call(rbind,lapply(out, `[[`, "site"))
  block_out <- do.call(rbind,lapply(out, `[[`, "block"))
  site_out <- rbind(site_out, score_constant_baseline(fold_data,y_all,draw,fold))
  block_out <- rbind(block_out,score_blocks_for_constant(fold_data,y_all,draw,fold))
  list(site=site_out,block=block_out)
}

# ------------------------- Fold / draw execution ----------------------------
cat("\nETH spatial CV:",n_sites,"unique sites;",length(blocks),"occupied blocks;",
    cv$folds,"folds; block width",cv$block_width_km,"km; buffer",
    cv$exclusion_buffer_km,"km.\n")
fold_preparations <- lapply(seq_len(cv$folds),make_training_fold)
fold_manifest <- do.call(rbind,lapply(seq_along(fold_preparations), function(f) {
  z <- fold_preparations[[f]]
  n_by_block <- table(block_id[z$test])
  data.frame(fold=f,train_sites=length(z$train),heldout_sites=length(z$test),
    occupied_heldout_blocks=length(n_by_block),
    eligible_blocks=sum(n_by_block>=cv$min_evaluation_sites_per_block),
    eligible_test_sites=sum(n_by_block[n_by_block>=cv$min_evaluation_sites_per_block]),
    buffered_training_sites_excluded=length(z$excluded),
    minimum_test_to_training_km=min(z$distance_to_training_km),
    median_test_to_training_km=median(z$distance_to_training_km))
}))
write.csv(fold_manifest,block_output("fold_geometry_summary.csv"),row.names=FALSE)
print(fold_manifest)
if (any(fold_manifest$eligible_blocks == 0L))
  stop("Some folds contain no test blocks with >=",
       cv$min_evaluation_sites_per_block,
       " sites. Reduce block width or number of folds.")
block_eligibility <- do.call(rbind, lapply(names(blocks), function(b)
  data.frame(block_id=b, fold=unname(block_fold[b]),
             sites=length(blocks[[b]]),
             eligible=length(blocks[[b]])>=cv$min_evaluation_sites_per_block)))
write.csv(block_eligibility,block_output("block_eligibility.csv"),row.names=FALSE)

for (draw in cv$draw_ids) {
  y_all <- model$aggregate_draw_counts(raw,as.integer(draw+7L),base$coordinate_mapping)
  y_all <- y_all[base$kept_location_rows]
  if (length(y_all)!=n_sites || any(!is.finite(y_all)) ||
      any(y_all<0 | y_all>base$numtrials | abs(y_all-round(y_all))>1e-8)) stop("Invalid adjusted counts for draw ",draw)
  for (fold in seq_len(cv$folds)) {
    work <- file.path(cv$output_dir,sprintf("draw_%04d",draw),sprintf("fold_%02d",fold))
    dir.create(work,showWarnings=FALSE,recursive=TRUE)
    done_site <- file.path(work,"block200_heldout_site_results.rds")
    done_block <- file.path(work,"block200_heldout_block_results.rds")
    if (file.exists(done_site) && file.exists(done_block)) {
      cat("Already completed draw ",draw," fold ",fold,"\n",sep="")
      next
    }
    if (xor(file.exists(done_site),file.exists(done_block)))
      stop("Partial checkpoint in ",work,". Inspect before restarting.")
    cat("\nDRAW ",draw," / FOLD ",fold," START\n",sep="")
    v <- fit_cv(fold_preparations[[fold]],y_all,draw,fold,work)
    # Write both complete fold results only after all requested fits succeed.
    saveRDS(v$site, done_site)
    saveRDS(v$block, done_block)
    write.csv(v$site,file.path(work,"block200_heldout_site_results.csv"),row.names=FALSE)
    write.csv(v$block,file.path(work,"block200_heldout_block_results.csv"),row.names=FALSE)
  }
}

# ------------------- Export results: block means PRIMARY --------------------
# The primary metrics now concern block-averaged, held-out sampled populations.
# Individual survey-location diagnostics are retained separately. ALL 200-km
# exported files have a block200_ prefix, within a dedicated output directory.
# Never append to an earlier 100-km or earlier 200-km output folder.

site_files <- list.files(cv$output_dir, "^block200_heldout_site_results[.]rds$",
                         recursive=TRUE,full.names=TRUE)
block_files <- list.files(cv$output_dir, "^block200_heldout_block_results[.]rds$",
                          recursive=TRUE,full.names=TRUE)
if (!length(site_files) || length(site_files) != length(block_files))
  stop("Missing/incomplete fold outputs; inspect per-draw checkpoints.")
all_sites <- do.call(rbind,lapply(site_files,readRDS))
all_blocks <- do.call(rbind,lapply(block_files,readRDS))
write.csv(all_sites,block_output("all_heldout_site_scores.csv"),row.names=FALSE)
write.csv(all_blocks,block_output("all_heldout_block_scores.csv"),row.names=FALSE)

# Check every included block has the same locations across model comparisons,
# and each block belongs exclusively to the fold in which it is held out.
ref_blocks <- all_blocks[all_blocks$model=="selected_spatial_iid" &
                         all_blocks$averaging=="equal_sites",]
if (!nrow(ref_blocks)) stop("Missing primary selected-model block scores.")
if (any(ref_blocks$n_sites<cv$min_evaluation_sites_per_block) ||
    any(ref_blocks$fold != unname(block_fold[ref_blocks$block_id])))
  stop("Invalid block membership or fold assignment in exported results.")
if (length(unique(ref_blocks$draw)) == 1L && anyDuplicated(ref_blocks$block_id))
  stop("Held-out block occurred more than once in one-draw CV.")

# Proper continuous predictive score for averages: CRPS (lower is better).
# Withholding larger blocks and changing target CAN change the scores; do not
# interpret improvements as evidence that local predictions improved.
primary <- all_blocks[all_blocks$averaging=="equal_sites",,drop=FALSE]
key <- interaction(primary$draw,primary$fold,primary$model,drop=TRUE)
fold_scores <- do.call(rbind,lapply(split(primary,key),function(z) data.frame(
  draw=z$draw[1],fold=z$fold[1],model=z$model[1],
  n_blocks=nrow(z),n_sites=sum(z$n_sites),
  mean_block_crps=mean(z$block_crps),
  rmse_adjusted_fraction=sqrt(mean((z$posterior_mean_block_expected-
                                      z$observed_block_mean)^2)),
  mean_absolute_error=mean(abs(z$posterior_mean_block_expected-
                                z$observed_block_mean)),
  predictive_80pct_coverage=mean(z$observed_block_mean>=z$predictive_lo80 &
                                  z$observed_block_mean<=z$predictive_hi80),
  predictive_95pct_coverage=mean(z$observed_block_mean>=z$predictive_lo95 &
                                  z$observed_block_mean<=z$predictive_hi95),
  mean_95pct_predictive_fraction_width=mean(z$predictive_width95)
)))
row.names(fold_scores) <- NULL
write.csv(fold_scores,block_output("fold_draw_scores.csv"),row.names=FALSE)
model_scores <- do.call(rbind,lapply(split(fold_scores,fold_scores$model),
  function(z) {
    # Pool BLOCKS equally across folds for the primary model-level descriptive
    # metric; separately report unweighted averages of the five fold scores.
    match_blocks <- primary[primary$model == z$model[1],,drop=FALSE]
    data.frame(model=z$model[1],completed_fold_draws=nrow(z),
      total_validated_blocks=nrow(match_blocks),
      distinct_geographic_blocks=length(unique(match_blocks$block_id)),
      pooled_equal_block_crps=mean(match_blocks$block_crps),
      pooled_equal_block_rmse=sqrt(mean((match_blocks$posterior_mean_block_expected-
                                            match_blocks$observed_block_mean)^2)),
      pooled_equal_block_mae=mean(abs(match_blocks$posterior_mean_block_expected-
                                         match_blocks$observed_block_mean)),
      pooled_equal_block_80pct_coverage=mean(
        match_blocks$observed_block_mean>=match_blocks$predictive_lo80 &
        match_blocks$observed_block_mean<=match_blocks$predictive_hi80),
      pooled_equal_block_95pct_coverage=mean(
        match_blocks$observed_block_mean>=match_blocks$predictive_lo95 &
        match_blocks$observed_block_mean<=match_blocks$predictive_hi95),
      pooled_equal_block_mean_95pct_width=mean(match_blocks$predictive_width95),
      mean_fold_block_crps=mean(z$mean_block_crps),
      mean_fold_rmse=mean(z$rmse_adjusted_fraction),
      mean_fold_mae=mean(z$mean_absolute_error),
      mean_fold_predictive_80pct_coverage=mean(z$predictive_80pct_coverage),
      mean_fold_predictive_95pct_coverage=mean(z$predictive_95pct_coverage),
      mean_fold_95pct_interval_fraction_width=mean(
        z$mean_95pct_predictive_fraction_width))
  }))
row.names(model_scores) <- NULL
write.csv(model_scores,block_output("model_summary_descriptive.csv"),row.names=FALSE)

# Secondary cattle-weighted block sensitivity. Does not generally estimate
# area-wide prevalence either; large survey sites receive greater weight.
weighted <- all_blocks[all_blocks$averaging=="cattle_weighted",,drop=FALSE]
write.csv(weighted,block_output("cattle_weighted_block_sensitivity.csv"),row.names=FALSE)

# Preserve original SITE-LEVEL fold scores under NEW descriptive filenames.
site_groups <- split(all_sites,interaction(all_sites$draw,all_sites$fold,
                                            all_sites$model,drop=TRUE))
site_scores <- do.call(rbind,lapply(site_groups,function(z) data.frame(
  draw=z$draw[1],fold=z$fold[1],model=z$model[1],n_sites=nrow(z),
  total_log_predictive_score=sum(z$log_score),
  mean_log_predictive_score=mean(z$log_score),
  rmse_adjusted_fraction=sqrt(mean((z$posterior_mean_p-z$adjusted_fraction)^2)),
  mean_absolute_error=mean(abs(z$posterior_mean_p-z$adjusted_fraction)),
  predictive_80pct_coverage=mean(z$y>=z$count_lo80 & z$y<=z$count_hi80),
  predictive_95pct_coverage=mean(z$y>=z$count_lo95 & z$y<=z$count_hi95),
  mean_95pct_predictive_fraction_width=mean((z$count_hi95-z$count_lo95)/z$n),
  mean_randomized_pit=mean(z$randomized_pit)
)))
row.names(site_scores) <- NULL
write.csv(site_scores,block_output("fold_draw_site_scores.csv"),row.names=FALSE)
site_model_scores <- do.call(rbind,lapply(split(site_scores,site_scores$model),
  function(z) data.frame(model=z$model[1],completed_fold_draws=nrow(z),
    mean_fold_log_score_per_site=mean(z$mean_log_predictive_score),
    mean_fold_rmse=mean(z$rmse_adjusted_fraction),
    mean_fold_predictive_80pct_coverage=mean(z$predictive_80pct_coverage),
    mean_fold_predictive_95pct_coverage=mean(z$predictive_95pct_coverage),
    mean_fold_95pct_interval_fraction_width=mean(
      z$mean_95pct_predictive_fraction_width)
)))
write.csv(site_model_scores,block_output("site_model_summary_descriptive.csv"),row.names=FALSE)

# Primary block diagnostics: prediction-versus-observation and accuracy by
# number of independent surveyed locations in the held-out block.
pdf(block_output("cv_diagnostics.pdf"),width=10,height=5)
old_par <- par(mfrow=c(1,2))
sel <- primary[primary$model=="selected_spatial_iid",]
if (nrow(sel)) {
  plot(sel$observed_block_mean,sel$posterior_mean_block_expected,
       xlim=c(0,1),ylim=c(0,1),pch=16,cex=.8,
       xlab="Observed mean of held-out site prevalences",
       ylab="Predicted mean at those held-out sites",
       main=paste0("ETH ",cv$block_width_km,"-km block-mean validation"))
  abline(0,1,col="red",lty=2)
  plot(sel$n_sites,abs(sel$posterior_mean_block_expected-sel$observed_block_mean),
       pch=16,cex=.8,
       xlab="Distinct held-out sites in the block",
       ylab="Absolute error in block mean prevalence",
       main="Block-mean error versus replication")
}
par(old_par);dev.off()

# Original fine-scale diagnostic plots remain available separately.
pdf(block_output("cv_site_diagnostics.pdf"),width=10,height=5)
old_par <- par(mfrow=c(1,2))
site_sel <- all_sites[all_sites$model=="selected_spatial_iid",]
if (nrow(site_sel)) {
  hist(site_sel$randomized_pit,breaks=seq(0,1,by=.1),
       xlab="Randomized PIT",main="ETH held-out single-site PIT")
  abline(h=nrow(site_sel)/10,col="red",lty=2)
  plot(site_sel$distance_to_training_km,
       abs(site_sel$posterior_mean_p-site_sel$adjusted_fraction),
       xlab="Nearest retained training site (km)",
       ylab="Absolute single-site error",main="Single-site error by distance",
       pch=16,cex=.6)
}
par(old_par);dev.off()

cat("\nETH 200-km BLOCK-MEAN CV complete for requested diagnostic draws.\n",
    "Outputs:",normalizePath(cv$output_dir),"\n",
    "Primary comparison: block200_model_summary_descriptive.csv and block200_fold_draw_scores.csv\n",
    "Fine-scale comparison: block200_site_model_summary_descriptive.csv\n")
print(model_scores)
cat("LIMITATION: block targets average surveyed locations, not unbiased regional means.\n",
    "Upstream diagnostic correction used the original data, including test sites.\n",
    "No raw-diagnostic end-to-end validation has been performed.\n",sep="")

}

run_eth_cv(cv)
