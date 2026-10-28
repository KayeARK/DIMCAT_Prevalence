# AAT Nigeria diagnostics: current metric-mesh / quadratic / IID workflow
# Updated to the current grouped-quadratic + location-IID production workflow.
#
# This diagnostic deliberately reuses the production spatial_preparation.rds.
# Therefore it uses the exact projected-metre mesh, SPDE, country mask,
# observation-only covariate scaling, quadratic basis and prediction grid from
# the current model, rather than recreating the old degree-based mesh.
#
# The current model engine is loaded for helper functions WITHOUT running its
# final run_aat_quick(settings) call.
#
# Main outputs:
#   Diagnostics.csv          pooled current-model diagnostic surfaces
#   Model_comparison.csv     spatial+IID vs no-spatial+IID WAIC by draw
#   Per_draw_summary.csv     simple draw-level summaries
# Checkpoints_latest_workflow/ allows safe resume after interruption.

diagnostic_settings <- list(
  data_path = "Code/TestSensSpec/adjusted_cases_all_sites.csv",

  # Blank = auto-find the model script. Otherwise paste the exact current file.
  model_script = "",

  # Blank = index ALL compatible current-workflow production folders beside
  # the model script.  Otherwise restrict diagnostics to one exact run folder.
  production_run_dir = "",

  # NULL = automatically use the first 20 available compatible current-workflow draws.
  # Or supply explicit model draw numbers, e.g. 83:102.
  draw_ids = NULL,
  n_draws = 20L,

  # NULL = use production prediction_manifest.rds (normally 4096).
  posterior_samples = NULL,

  output_dir = "Code/Prevalence/Bovine BCT and PCR/Diagnostics NGA",
  seed = 20260928L,
  integration_tolerance = 1e-4,
  prediction_chunk_size = 128L,
  resume = TRUE
)

country_code <- "NGA"

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(matrixStats)
  library(sf)
})

atomic_save_rds_diag <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Could not commit checkpoint: ", path)
  invisible(path)
}

atomic_write_csv_diag <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  write.csv(object, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) stop("Could not commit CSV: ", path)
  invisible(path)
}

diagnostic_script_directory <- function() {
  ofiles <- lapply(sys.frames(), function(x) x$ofile)
  ofiles <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x), ofiles)
  if (length(ofiles))
    return(dirname(normalizePath(tail(ofiles, 1)[[1]], mustWork = TRUE)))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg))
    return(dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)))
  getwd()
}

find_model_script <- function(explicit, cc) {
  if (nzchar(explicit))
    return(normalizePath(path.expand(explicit), mustWork = TRUE))

  # VS Code/R terminals commonly keep the project root as getwd() even when
  # this diagnostic script and the production model script live elsewhere.
  # Search both locations.  The regex also accepts downloaded names such as
  # AAT_spatial_quick_joint_posterior_NGA(3).R.
  roots <- unique(normalizePath(
    c(getwd(), diagnostic_script_directory()),
    mustWork = TRUE
  ))
  pattern <- paste0("^AAT_spatial_quick_joint_posterior_", cc, ".*[.]R$")
  fs <- unique(unlist(lapply(roots, function(root)
    list.files(root, pattern = pattern, full.names = TRUE, ignore.case = TRUE))))

  if (!length(fs))
    stop(
      "Cannot find the current ", cc, " model script in either:\n  ",
      paste(roots, collapse = "\n  "),
      "\nSet diagnostic_settings$model_script to its full path if it is stored elsewhere."
    )

  # Prefer the canonical unsuffixed name when present; otherwise use the newest
  # country-specific matching file.
  canonical <- fs[tolower(basename(fs)) ==
                    tolower(paste0("AAT_spatial_quick_joint_posterior_", cc, ".R"))]
  if (length(canonical))
    return(normalizePath(canonical[1], mustWork = TRUE))

  fs <- fs[order(file.info(fs)$mtime, decreasing = TRUE)]
  normalizePath(fs[1], mustWork = TRUE)
}

load_model_helpers_without_running <- function(path) {
  x <- readLines(path, warn = FALSE)
  keep <- !grepl("^\\s*run_aat_quick\\s*\\(\\s*settings\\s*\\)\\s*$", x)
  if (sum(!keep) != 1L)
    stop("Expected exactly one final run_aat_quick(settings) line in ", path,
         ". Refusing to source an ambiguous model script.")
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(x[keep], tmp)
  sys.source(tmp, envir = .GlobalEnv)
  invisible(TRUE)
}

discover_production_runs <- function(base, explicit, cc) {
  roots <- if (nzchar(explicit)) {
    normalizePath(path.expand(explicit), mustWork = TRUE)
  } else {
    x <- list.dirs(path.expand(base), recursive = FALSE, full.names = TRUE)
    x[grepl(paste0("^Spatial_joint_", cc), basename(x))]
  }

  if (!length(roots))
    stop("No Spatial_joint_", cc, " production folders were found under: ", base)

  rows <- lapply(roots, function(root) {
    manifest_file <- file.path(root, "run_manifest.rds")
    prep_file <- file.path(root, paste0("Checkpoints_", cc), "spatial_preparation.rds")
    cov_dir <- file.path(root, paste0("Covariates_", cc))
    if (!file.exists(manifest_file) || !file.exists(prep_file) || !dir.exists(cov_dir))
      return(NULL)

    man <- tryCatch(readRDS(manifest_file), error = function(e) NULL)
    if (is.null(man) ||
        !identical(as.integer(man$covariate_degree), 2L) ||
        !isTRUE(man$use_location_iid) ||
        is.null(man$workflow_version) ||
        !grepl("grouped_quadratic", man$workflow_version))
      return(NULL)

    fs <- list.files(cov_dir,
      pattern = "^Covariates_model_[0-9]+[.]csv$", full.names = TRUE)
    if (!length(fs)) return(NULL)
    ids <- suppressWarnings(as.integer(sub(
      "^Covariates_model_([0-9]+)[.]csv$", "\\1", basename(fs))))
    ok <- is.finite(ids)
    if (!any(ok)) return(NULL)

    data.frame(
      draw = ids[ok],
      cov_file = normalizePath(fs[ok], mustWork = TRUE),
      root = normalizePath(root, mustWork = TRUE),
      manifest_file = normalizePath(manifest_file, mustWork = TRUE),
      prep_file = normalizePath(prep_file, mustWork = TRUE),
      mtime = as.numeric(file.info(manifest_file)$mtime),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows))
    stop("Spatial_joint_", cc,
         " folders exist, but none is a compatible grouped-quadratic + IID production run with covariate files.")

  idx <- do.call(rbind, rows)

  # A draw can occur in more than one continuation folder.  Accept duplicates
  # only when their selected-covariate files are byte-identical; then use the
  # newest compatible copy.  A conflict is scientifically ambiguous and stops.
  out <- lapply(split(idx, idx$draw), function(z) {
    if (nrow(z) > 1L) {
      h <- unname(tools::md5sum(z$cov_file))
      if (length(unique(h)) != 1L)
        stop("Conflicting selected-covariate files found for model draw ", z$draw[1],
             " across: ", paste(z$root, collapse = "; "))
      z <- z[order(z$mtime, decreasing = TRUE), , drop = FALSE]
    }
    z[1, , drop = FALSE]
  })
  idx <- do.call(rbind, out)
  idx <- idx[order(idx$draw), , drop = FALSE]
  rownames(idx) <- NULL
  idx
}

choose_diagnostic_draws <- function(cfg, index) {
  if (is.null(cfg$draw_ids)) {
    n <- as.integer(cfg$n_draws)
    if (length(n) != 1L || is.na(n) || n < 1L || n != cfg$n_draws)
      stop("n_draws must be one positive integer when draw_ids = NULL.")
    if (nrow(index) < n)
      stop("Only ", nrow(index), " compatible current-workflow draws are available; requested ", n, ".")
    return(index$draw[seq_len(n)])
  }
  x <- as.integer(cfg$draw_ids)
  if (!length(x) || anyNA(x) || any(x < 1L) || anyDuplicated(x))
    stop("draw_ids must be distinct positive model-draw numbers, or NULL.")
  missing <- setdiff(x, index$draw)
  if (length(missing))
    stop("Requested current-workflow draw(s) not found in compatible production folders: ",
         paste(missing, collapse = ", "),
         ". Available range is ", min(index$draw), "-", max(index$draw), ".")
  x
}

assert_production_preparations_compatible <- function(reference_root, roots, cc) {
  roots <- unique(normalizePath(roots, mustWork = TRUE))
  ref_m <- readRDS(file.path(reference_root, "run_manifest.rds"))
  ref_p <- readRDS(file.path(reference_root, paste0("Checkpoints_", cc), "spatial_preparation.rds"))

  manifest_fields <- c(
    "input_md5", "country", "input_filter", "elevation_policy", "spatial_crs",
    "cutoff_m", "max_edge_m", "offset_m", "range_prior_m", "sigma_prior",
    "coefficient_prior_sd", "intercept_prior_sd", "covariate_degree",
    "use_location_iid", "iid_prior", "control_inla"
  )
  prep_fields <- c(
    "candidate_covariates", "use_location_iid", "numtrials", "kept_location_rows",
    "observation_lonlat", "prediction_lonlat", "observation_m", "prediction_m",
    "covariate_audit", "quadratic_basis", "coordinate_mapping"
  )

  for (root in setdiff(roots, reference_root)) {
    m <- readRDS(file.path(root, "run_manifest.rds"))
    p <- readRDS(file.path(root, paste0("Checkpoints_", cc), "spatial_preparation.rds"))

    bad_m <- manifest_fields[!vapply(manifest_fields, function(nm)
      identical(ref_m[[nm]], m[[nm]]), logical(1))]
    bad_p <- prep_fields[!vapply(prep_fields, function(nm)
      identical(ref_p[[nm]], p[[nm]]), logical(1))]

    # The sparse projectors and SPDE object must also represent the same setup.
    same_A <- identical(ref_p$A, p$A)
    same_Ap <- identical(ref_p$Ap, p$Ap)
    if (length(bad_m) || length(bad_p) || !same_A || !same_Ap) {
      stop(
        "Current-workflow production folders are not preparation-compatible.\nReference: ", reference_root,
        "\nOther: ", root,
        if (length(bad_m)) paste0("\nManifest differences: ", paste(bad_m, collapse = ", ")) else "",
        if (length(bad_p)) paste0("\nPreparation differences: ", paste(bad_p, collapse = ", ")) else "",
        if (!same_A) "\nA projector differs." else "",
        if (!same_Ap) "\nAp projector differs." else ""
      )
    }
  }
  invisible(TRUE)
}

find_project_root_for_data <- function(model_base, relative_data_path) {
  if (grepl("^/", path.expand(relative_data_path)))
    return(dirname(normalizePath(path.expand(relative_data_path), mustWork = TRUE)))
  p <- normalizePath(model_base, mustWork = TRUE)
  repeat {
    candidate <- file.path(p, relative_data_path)
    if (file.exists(candidate)) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  stop("Cannot locate ", relative_data_path,
       " in any ancestor of the current model-script directory: ", model_base)
}

prepare_input_like_production <- function(d, cc) {
  d$longitude <- suppressWarnings(as.numeric(as.character(d$longitude)))
  d$latitude <- suppressWarnings(as.numeric(as.character(d$latitude)))

  if (cc == "ETH") {
    # Exact orientation/bounding-box logic used by the current ETH engine.
    d <- filter_ethiopia_inputs(d)
  } else {
    keep <- is.finite(d$longitude) & is.finite(d$latitude)
    if (any(!keep))
      message("Removing ", sum(!keep), " rows with missing/non-finite coordinates.")
    d <- d[keep, , drop = FALSE]
    if (!nrow(d) ||
        any(abs(d$longitude) > 180 | abs(d$latitude) > 90))
      stop("No valid geographic coordinates remain.")
  }
  d
}

validate_production <- function(run_dir, cc, cfg) {
  manifest_file <- file.path(run_dir, "run_manifest.rds")
  prep_file <- file.path(run_dir, paste0("Checkpoints_", cc), "spatial_preparation.rds")
  manifest <- readRDS(manifest_file)
  prepared <- readRDS(prep_file)

  if (!isTRUE(manifest$covariate_degree == 2L))
    stop("Chosen production run is not quadratic.")
  if (!isTRUE(manifest$use_location_iid) ||
      !isTRUE(prepared$use_location_iid))
    stop("Chosen production run does not use the location IID effect.")
  if (is.null(prepared$quadratic_basis))
    stop("Production preparation lacks quadratic_basis.")
  if (is.null(prepared$raw_observation_covariates) ||
      is.null(prepared$raw_prediction_covariates) ||
      is.null(prepared$covariate_audit))
    stop("Production preparation lacks current covariate scaling/audit objects.")
  if (is.null(prepared$observation_m) ||
      is.null(prepared$prediction_m) ||
      is.null(prepared$A) || is.null(prepared$Ap) ||
      is.null(prepared$spde))
    stop("Production preparation lacks the current projected mesh/projectors.")

  crs <- sf::st_crs(manifest$spatial_crs)
  if (is.na(crs) || sf::st_is_longlat(crs))
    stop("Chosen production run is not using a projected metric CRS; refusing old mesh.")

  if (is.null(prepared$numtrials))
    stop("Production preparation lacks explicit Ntrials.")

  pred_manifest_file <- file.path(run_dir, "prediction_manifest.rds")
  pred_manifest <- if (file.exists(pred_manifest_file))
    readRDS(pred_manifest_file) else NULL

  # predict_joint() uses the production Monte Carlo screen threshold.
  cfg$mc_screen_pp <- if (!is.null(pred_manifest$mc_screen_pp))
    pred_manifest$mc_screen_pp else 3

  if (is.null(cfg$posterior_samples)) {
    cfg$posterior_samples <- if (!is.null(pred_manifest$posterior_samples))
      as.integer(pred_manifest$posterior_samples) else 4096L
  }
  if (cfg$posterior_samples < 1024L ||
      cfg$posterior_samples != as.integer(cfg$posterior_samples))
    stop("posterior_samples must be an integer >= 1024.")

  # fit_observations() expects these current-engine controls.
  cfg$coefficient_prior_sd <- manifest$coefficient_prior_sd
  cfg$intercept_prior_sd <- manifest$intercept_prior_sd
  cfg$control_inla <- manifest$control_inla
  cfg$inla_mode <- manifest$inla_mode
  if (is.null(cfg$inla_mode) || !nzchar(cfg$inla_mode))
    cfg$inla_mode <- INLA::inla.getOption("inla.mode")

  cfg$local_work_root <- tempfile("aat_diag_inla_")
  dir.create(cfg$local_work_root, recursive = TRUE, showWarnings = FALSE)
  cfg$diagnostics_dir <- cfg$output_dir
  cfg$runtime <- capture_inla_runtime()

  list(
    manifest = manifest,
    prepared = prepared,
    cfg = cfg,
    manifest_md5 = unname(tools::md5sum(manifest_file)),
    preparation_md5 = unname(tools::md5sum(prep_file))
  )
}

# Same selected fixed effects and same location IID, but no spatial Gaussian field.
fit_nospatial_iid <- function(prepared, y, covariates, manifest, cfg) {
  xo <- transform_model_covariates(prepared$raw_observation_covariates, prepared)
  terms <- effect_terms(covariates, prepared)
  n <- nrow(xo)

  dat <- data.frame(
    y = as.numeric(y),
    numtrials = as.numeric(prepared$numtrials),
    b0 = rep(1, n),
    xo[, terms, drop = FALSE],
    check.names = FALSE
  )
  if (prepared$use_location_iid) dat$site_iid <- seq_len(n)

  iid_hyper <- list(prec = list(prior = "pc.prec", param = manifest$iid_prior))
  rhs <- c("0", "b0", terms)
  if (prepared$use_location_iid)
    rhs <- c(rhs, "f(site_iid, model = 'iid', constr = FALSE, hyper = iid_hyper)")
  formula <- as.formula(paste("y ~", paste(rhs, collapse = " + ")),
                        env = environment())

  args <- list(
    formula = formula,
    data = dat,
    family = "binomial",
    Ntrials = dat$numtrials,
    inla.mode = cfg$inla_mode,
    num.threads = "1:1",
    control.inla = cfg$control_inla,
    control.fixed = fixed_prior_control(
      manifest$coefficient_prior_sd,
      manifest$intercept_prior_sd
    ),
    control.compute = list(
      waic = TRUE, config = TRUE, dic = FALSE,
      return.marginals.predictor = FALSE
    ),
    control.predictor = list(link = 1, compute = TRUE),
    verbose = FALSE
  )
  res <- do.call(INLA::inla, args)
  if (length(res$waic$waic) != 1L || !is.finite(res$waic$waic))
    stop("No-spatial+IID model returned non-finite WAIC.")
  res
}

sample_nospatial_iid <- function(res, prepared, covariates, cfg, draw) {
  fixed_names <- c("b0", effect_terms(covariates, prepared))
  selection <- setNames(rep(list(0L), length(fixed_names)), fixed_names)
  seed <- as.integer((as.double(cfg$seed) + 2003 * draw) %% 2147483646 + 1L)

  set.seed(seed)
  samples <- INLA::inla.posterior.sample(
    cfg$posterior_samples, res, selection = selection,
    intern = FALSE, use.improved.mean = TRUE,
    add.names = FALSE, seed = seed, num.threads = 1L
  )
  labels <- rownames(samples[[1]]$latent)
  if (is.null(labels)) labels <- names(samples[[1]]$latent)
  if (is.null(labels) || anyDuplicated(labels))
    stop("Missing/duplicate no-spatial posterior labels.")

  b <- vapply(fixed_names, function(nm)
    sample_component_rows(labels, nm, 1L), integer(1))
  beta <- matrix(
    vapply(samples, function(z) as.numeric(z$latent)[b],
           numeric(length(b))),
    nrow = length(b),
    dimnames = list(fixed_names, NULL)
  )

  sigma <- rep(0, cfg$posterior_samples)
  if (prepared$use_location_iid) {
    hn <- rownames(samples[[1]]$hyperpar)
    if (is.null(hn)) hn <- names(samples[[1]]$hyperpar)
    j <- which(hn == "Precision for site_iid")
    if (length(j) != 1L) stop("Cannot find sampled no-spatial IID precision.")
    precision <- vapply(samples, function(z) as.numeric(z$hyperpar)[j], numeric(1))
    if (any(!is.finite(precision) | precision <= 0))
      stop("Invalid no-spatial IID precision samples.")
    sigma <- 1 / sqrt(precision)
  }

  rm(samples)
  gc(verbose = FALSE)
  list(beta = beta, sigma_iid = sigma, seed = seed)
}

summarise_components <- function(prepared, covariates, joint, cfg) {
  ncell <- nrow(prepared$prediction_lonlat)
  support <- support_at_predictions(prepared, covariates)
  valid <- support$valid_prediction
  xp <- transform_model_covariates(prepared$raw_prediction_covariates, prepared)
  fixed_terms <- effect_terms(covariates, prepared)
  X <- cbind(b0 = 1, as.matrix(xp[, fixed_terms, drop = FALSE]))
  if (!identical(colnames(X), rownames(joint$beta)))
    stop("Prediction fixed-effect order mismatch.")

  out <- list(
    spatial_mean = rep(NA_real_, ncell),
    spatial_var = rep(NA_real_, ncell),
    fixed_mean = rep(NA_real_, ncell),
    fixed_var = rep(NA_real_, ncell)
  )

  all_rows <- seq_len(ncell)
  for (start in seq.int(1L, ncell, by = cfg$prediction_chunk_size)) {
    rows <- all_rows[start:min(ncell, start + cfg$prediction_chunk_size - 1L)]
    w <- as.matrix(prepared$Ap[rows, , drop = FALSE] %*% joint$spatial)
    out$spatial_mean[rows] <- rowMeans(w)
    out$spatial_var[rows] <- matrixStats::rowVars(w)

    ok <- rows[valid[rows]]
    if (length(ok)) {
      f <- X[ok, , drop = FALSE] %*% joint$beta
      out$fixed_mean[ok] <- rowMeans(f)
      out$fixed_var[ok] <- matrixStats::rowVars(f)
    }
  }
  out
}

predict_nospatial_iid <- function(prepared, covariates, joint, cfg) {
  ncell <- nrow(prepared$prediction_lonlat)
  support <- support_at_predictions(prepared, covariates)
  valid <- support$valid_prediction
  xp <- transform_model_covariates(prepared$raw_prediction_covariates, prepared)
  fixed_terms <- effect_terms(covariates, prepared)
  X <- cbind(b0 = 1, as.matrix(xp[, fixed_terms, drop = FALSE]))
  if (!identical(colnames(X), rownames(joint$beta)))
    stop("No-spatial prediction fixed-effect order mismatch.")

  expected_mean <- expected_var <- surface_mean <- surface_var <-
    rep(NA_real_, ncell)

  fun <- make_iid_integrator(joint$sigma_iid, cfg$integration_tolerance)
  ok_all <- which(valid)

  for (start in seq.int(1L, length(ok_all), by = cfg$prediction_chunk_size)) {
    rows <- ok_all[start:min(length(ok_all), start + cfg$prediction_chunk_size - 1L)]
    eta <- X[rows, , drop = FALSE] %*% joint$beta
    surface <- plogis(eta)
    surface_mean[rows] <- rowMeans(surface)
    surface_var[rows] <- matrixStats::rowVars(surface)

    sigma <- rep(joint$sigma_iid, each = length(rows))
    expected <- matrix(
      fun$evaluate(as.numeric(eta), sigma),
      nrow = length(rows)
    )
    expected_mean[rows] <- rowMeans(expected)
    expected_var[rows] <- matrixStats::rowVars(expected)
  }

  list(
    expected_mean = expected_mean,
    expected_var = expected_var,
    surface_mean = surface_mean,
    surface_var = surface_var,
    iid_sd_mean = mean(joint$sigma_iid),
    iid_sd_median = median(joint$sigma_iid),
    integration_error = fun$error
  )
}

pool_component <- function(results, mean_name, var_name) {
  ncell <- length(results[[1]][[mean_name]])
  count <- integer(ncell)
  sm <- sv <- ss <- numeric(ncell)

  for (z in results) {
    m <- z[[mean_name]]
    v <- z[[var_name]]
    ok <- is.finite(m) & is.finite(v) & v >= 0
    count[ok] <- count[ok] + 1L
    sm[ok] <- sm[ok] + m[ok]
    sv[ok] <- sv[ok] + v[ok]
    ss[ok] <- ss[ok] + v[ok] + m[ok]^2
  }

  mu <- within <- total <- between <- rep(NA_real_, ncell)
  ok <- count > 0L
  mu[ok] <- sm[ok] / count[ok]
  within[ok] <- pmax(sv[ok] / count[ok], 0)
  total[ok] <- pmax(ss[ok] / count[ok] - mu[ok]^2, 0)
  between[ok] <- pmax(total[ok] - within[ok], 0)

  list(
    mean = mu,
    sd = sqrt(total),
    within_sd = sqrt(within),
    between_sd = sqrt(between),
    n = count
  )
}

run_latest_diagnostics <- function(cfg, cc) {
  model_script <- find_model_script(cfg$model_script, cc)
  cat("Loading current model helpers from:\n", model_script, "\n\n", sep = "")

  # Production runs are written alongside the current model script.  Derive the
  # search root from that absolute path rather than from getwd(), because VS Code
  # may launch R with the working directory already inside this folder.
  base <- dirname(normalizePath(model_script, mustWork = TRUE))
  project_root <- find_project_root_for_data(base, cfg$data_path)
  if (!grepl("^/", path.expand(cfg$output_dir)))
    cfg$output_dir <- file.path(project_root, cfg$output_dir)
  cfg$output_dir <- path.expand(cfg$output_dir)

  load_model_helpers_without_running(model_script)
  required_helpers <- c(
    "fit_observations", "sample_joint_posterior", "predict_joint",
    "effect_terms", "transform_model_covariates", "support_at_predictions",
    "sample_component_rows", "make_iid_integrator", "nearest_distances_km",
    "capture_inla_runtime", "fixed_prior_control", "aggregate_draw_counts"
  )
  missing_helpers <- required_helpers[
    !vapply(required_helpers, exists, logical(1), mode = "function", inherits = TRUE)
  ]
  if (length(missing_helpers))
    stop("The selected model script is not the current workflow; missing helpers: ",
         paste(missing_helpers, collapse = ", "))
  if (!exists("settings", inherits = TRUE) ||
      !isTRUE(settings$covariate_degree == 2L) ||
      !isTRUE(settings$use_location_iid))
    stop("The selected model script is not configured for grouped quadratic covariates + location IID.")

  production_index <- discover_production_runs(base, cfg$production_run_dir, cc)
  draw_ids <- choose_diagnostic_draws(cfg, production_index)
  selected_index <- production_index[match(draw_ids, production_index$draw), , drop = FALSE]
  if (anyNA(selected_index$draw)) stop("Internal production-index matching failure.")

  reference_run <- selected_index$root[1]
  used_runs <- unique(selected_index$root)
  assert_production_preparations_compatible(reference_run, used_runs, cc)

  cat("Using compatible production run set:\n  ",
      paste(basename(used_runs), collapse = "\n  "), "\n", sep = "")
  cat("Diagnostic model draws: ", paste(draw_ids, collapse = ", "), "\n\n", sep = "")

  v <- validate_production(reference_run, cc, cfg)
  manifest <- v$manifest
  prepared <- v$prepared
  cfg <- v$cfg
  on.exit(unlink(cfg$local_work_root, recursive = TRUE), add = TRUE)

  cat(
    "Verified current workflow:\n",
    "  workflow: ", manifest$workflow_version, "\n",
    "  quadratic degree: ", manifest$covariate_degree, "\n",
    "  location IID: ", manifest$use_location_iid, "\n",
    "  mesh cutoff: ", manifest$cutoff_m / 1000, " km\n",
    "  max edges: ", paste(manifest$max_edge_m / 1000, collapse = ", "), " km\n",
    "  offsets: ", paste(manifest$offset_m / 1000, collapse = ", "), " km\n",
    "  SPDE range threshold: ", manifest$range_prior_m[1] / 1000, " km\n",
    "  posterior samples/draw: ", cfg$posterior_samples, "\n\n",
    sep = ""
  )

  data_path <- normalizePath(file.path(project_root, cfg$data_path), mustWork = TRUE)

  if (!identical(unname(tools::md5sum(data_path)), manifest$input_md5))
    stop("Input CSV differs from production run (MD5 mismatch).")

  data_original <- read.csv(data_path)
  data_original <- prepare_input_like_production(data_original, cc)

  if (!is.numeric(data_original$sample_size) ||
      any(!is.finite(data_original$sample_size)) ||
      any(data_original$sample_size <= 0) ||
      any(abs(data_original$sample_size - round(data_original$sample_size)) > 1e-8))
    stop("sample_size must contain positive integer counts.")

  if (max(draw_ids) + 7L > ncol(data_original))
    stop("Requested model draw exceeds available adjusted-count columns.")

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- file.path(cfg$output_dir, "Checkpoints_latest_workflow")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  signature <- list(
    country = cc,
    production_runs = sort(unique(normalizePath(selected_index$root))),
    draw_source_map = selected_index[, c("draw", "root", "cov_file"), drop = FALSE],
    model_script_md5 = unname(tools::md5sum(model_script)),
    manifest_md5 = v$manifest_md5,
    preparation_md5 = v$preparation_md5,
    input_md5 = manifest$input_md5,
    draw_ids = draw_ids,
    posterior_samples = cfg$posterior_samples,
    integration_tolerance = cfg$integration_tolerance,
    seed = cfg$seed
  )
  signature_file <- file.path(checkpoint_dir, "diagnostic_signature.rds")
  if (file.exists(signature_file)) {
    if (!identical(readRDS(signature_file), signature))
      stop("Existing diagnostic checkpoints use different settings. ",
           "Use a new output_dir or intentionally remove Checkpoints_latest_workflow.")
  } else {
    atomic_save_rds_diag(signature, signature_file)
  }

  results <- vector("list", length(draw_ids))
  names(results) <- as.character(draw_ids)

  for (ii in seq_along(draw_ids)) {
    draw <- draw_ids[ii]
    cp <- file.path(checkpoint_dir, sprintf("diagnostic_draw_%04d.rds", draw))

    if (isTRUE(cfg$resume) && file.exists(cp)) {
      cat("Draw ", draw, ": loading diagnostic checkpoint.\n", sep = "")
      results[[ii]] <- readRDS(cp)
      next
    }

    cat("\n========================================\n")
    cat("Diagnostic draw ", draw, " (", ii, "/", length(draw_ids), ")\n", sep = "")
    cat("========================================\n")

    cov_file <- selected_index$cov_file[match(draw, selected_index$draw)]
    if (length(cov_file) != 1L || is.na(cov_file))
      stop("Internal covariate-file index failure for draw ", draw)
    if (!file.exists(cov_file))
      stop("Missing selected-covariate file: ", cov_file)
    tab <- read.csv(cov_file, stringsAsFactors = FALSE)
    if (!"covariates" %in% names(tab))
      stop("Missing covariates column: ", cov_file)
    covariates <- trimws(as.character(tab$covariates))
    covariates <- covariates[!is.na(covariates) & nzchar(covariates)]
    if (anyDuplicated(covariates) ||
        !all(covariates %in% prepared$candidate_covariates))
      stop("Invalid selected parent covariates in ", cov_file)

    counts <- aggregate_draw_counts(
      data_original, draw + 7L, prepared$coordinate_mapping
    )
    y <- counts[prepared$kept_location_rows]
    if (length(y) != nrow(prepared$observation_lonlat) ||
        any(!is.finite(y)) || any(y < 0) ||
        any(y > prepared$numtrials))
      stop("Response/preparation mismatch for draw ", draw)

    cat("Selected parent covariates: ",
        if (length(covariates)) paste(covariates, collapse = " + ")
        else "(intercept only)", "\n", sep = "")

    # Full current spatial + IID model.
    fit_sp <- fit_observations(prepared, y, covariates, cfg)
    joint_sp <- sample_joint_posterior(fit_sp, prepared, covariates, cfg, draw)
    obj <- list(
      covariates = covariates,
      response = y,
      joint = joint_sp,
      waic = fit_sp$waic$waic,
      fixed = fit_sp$summary.fixed,
      hyperparameters = fit_sp$summary.hyperpar,
      spatial_nodes = fit_sp$summary.random$s,
      observed_iid = if (prepared$use_location_iid)
        fit_sp$summary.random$site_iid else NULL,
      observed_prevalence =
        fit_sp$summary.fitted.values[prepared$observation_index, , drop = FALSE]
    )
    pred_sp <- predict_joint(obj, prepared, cfg, draw)
    comp <- summarise_components(prepared, covariates, joint_sp, cfg)

    # Comparator: same selected quadratic fixed effects + same IID, no GP.
    fit_ng <- fit_nospatial_iid(prepared, y, covariates, manifest, cfg)
    joint_ng <- sample_nospatial_iid(fit_ng, prepared, covariates, cfg, draw)
    pred_ng <- predict_nospatial_iid(prepared, covariates, joint_ng, cfg)

    result <- list(
      draw = draw,
      covariates = covariates,

      # Current publication target.
      expected_mean = pred_sp$expected_prevalence$mean,
      expected_var = pred_sp$expected_prevalence$sd^2,

      # IID=0 latent surface.
      surface_mean = pred_sp$surface_prevalence$mean,
      surface_var = pred_sp$surface_prevalence$sd^2,
      eta_mean = pred_sp$surface_linear_predictor$mean,
      eta_var = pred_sp$surface_linear_predictor$sd^2,

      # Decomposition of eta = X beta + w.
      spatial_mean = comp$spatial_mean,
      spatial_var = comp$spatial_var,
      fixed_mean = comp$fixed_mean,
      fixed_var = comp$fixed_var,

      # No-spatial comparator, still fitted with the location IID.
      nogp_expected_mean = pred_ng$expected_mean,
      nogp_expected_var = pred_ng$expected_var,
      nogp_surface_mean = pred_ng$surface_mean,
      nogp_surface_var = pred_ng$surface_var,

      waic_spatial_iid = fit_sp$waic$waic,
      waic_nospatial_iid = fit_ng$waic$waic,
      iid_sd_spatial_median = median(joint_sp$sigma_iid),
      iid_sd_nospatial_median = median(joint_ng$sigma_iid),
      integration_error_spatial =
        pred_sp$integration_check$max_checked_probability_error[1],
      integration_error_nospatial = pred_ng$integration_error
    )

    atomic_save_rds_diag(result, cp)
    results[[ii]] <- result

    cat(
      "WAIC spatial+IID = ", round(result$waic_spatial_iid, 2),
      "; no-spatial+IID = ", round(result$waic_nospatial_iid, 2),
      "; Delta = ",
      round(result$waic_nospatial_iid - result$waic_spatial_iid, 2),
      "\n", sep = ""
    )

    rm(fit_sp, joint_sp, obj, pred_sp, comp,
       fit_ng, joint_ng, pred_ng, result)
    gc(verbose = FALSE)
  }

  expected <- pool_component(results, "expected_mean", "expected_var")
  surface <- pool_component(results, "surface_mean", "surface_var")
  eta <- pool_component(results, "eta_mean", "eta_var")
  spatial <- pool_component(results, "spatial_mean", "spatial_var")
  fixed <- pool_component(results, "fixed_mean", "fixed_var")
  nogp_expected <- pool_component(
    results, "nogp_expected_mean", "nogp_expected_var"
  )
  nogp_surface <- pool_component(
    results, "nogp_surface_mean", "nogp_surface_var"
  )

  if (!is.null(prepared$nearest_observation_km)) {
    distance_km <- prepared$nearest_observation_km
  } else {
    distance_km <- nearest_distances_km(
      prepared$prediction_m, prepared$observation_m
    )
  }

  diagnostics <- data.frame(
    lon = prepared$prediction_lonlat[, 1],
    lat = prepared$prediction_lonlat[, 2],

    # Main current expected-prevalence target.
    prevalence_mean = expected$mean,
    prevalence_sd = expected$sd,
    prevalence_within_sd = expected$within_sd,
    prevalence_between_draw_sd = expected$between_sd,

    # IID=0 surface, directly relevant to the GP-regression-to-zero comment.
    surface_prevalence_mean = surface$mean,
    surface_prevalence_sd = surface$sd,
    surface_prevalence_within_sd = surface$within_sd,
    surface_prevalence_between_draw_sd = surface$between_sd,

    # Latent scale.
    eta_mean = eta$mean,
    eta_sd = eta$sd,
    eta_within_sd = eta$within_sd,
    eta_between_draw_sd = eta$between_sd,

    spatial_mean = spatial$mean,
    spatial_sd = spatial$sd,
    spatial_within_sd = spatial$within_sd,
    spatial_between_draw_sd = spatial$between_sd,

    fixed_effect_mean = fixed$mean,
    fixed_effect_sd = fixed$sd,
    fixed_effect_within_sd = fixed$within_sd,
    fixed_effect_between_draw_sd = fixed$between_sd,

    # Comparator preserving fixed effects + location IID but removing GP.
    prevalence_no_gp = nogp_expected$mean,
    prevalence_no_gp_sd = nogp_expected$sd,
    surface_prevalence_no_gp = nogp_surface$mean,
    surface_prevalence_no_gp_sd = nogp_surface$sd,

    dist_to_data_km = distance_km,
    n_draws_prevalence = expected$n,
    n_draws_surface = surface$n,
    n_draws_spatial = spatial$n,

    spatial_abs_mean_over_sd = ifelse(
      is.finite(spatial$sd) & spatial$sd > 0,
      abs(spatial$mean) / spatial$sd,
      NA_real_
    )
  )

  model_comparison <- do.call(rbind, lapply(results, function(z)
    data.frame(
      draw = z$draw,
      selected_parent_covariates = paste(z$covariates, collapse = " + "),
      waic_spatial_iid = z$waic_spatial_iid,
      waic_nospatial_iid = z$waic_nospatial_iid,
      delta_waic_nospatial_minus_spatial =
        z$waic_nospatial_iid - z$waic_spatial_iid,
      iid_sd_spatial_median = z$iid_sd_spatial_median,
      iid_sd_nospatial_median = z$iid_sd_nospatial_median,
      integration_error_spatial = z$integration_error_spatial,
      integration_error_nospatial = z$integration_error_nospatial,
      stringsAsFactors = FALSE
    )
  ))

  per_draw_summary <- do.call(rbind, lapply(results, function(z) {
    ok <- is.finite(z$expected_mean)
    data.frame(
      draw = z$draw,
      valid_grid_cells = sum(ok),
      mean_expected_prevalence = mean(z$expected_mean[ok]),
      mean_surface_prevalence = mean(z$surface_mean[ok]),
      median_expected_posterior_sd =
        median(sqrt(z$expected_var[ok])),
      median_surface_posterior_sd =
        median(sqrt(z$surface_var[ok])),
      mean_spatial_field = mean(z$spatial_mean, na.rm = TRUE),
      sd_spatial_field_across_grid = sd(z$spatial_mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  atomic_write_csv_diag(
    diagnostics, file.path(cfg$output_dir, "Diagnostics.csv")
  )
  atomic_write_csv_diag(
    model_comparison, file.path(cfg$output_dir, "Model_comparison.csv")
  )
  atomic_write_csv_diag(
    per_draw_summary, file.path(cfg$output_dir, "Per_draw_summary.csv")
  )
  atomic_save_rds_diag(
    list(
      diagnostic_settings = cfg,
      production_runs = unique(selected_index$root),
      production_manifest = manifest,
      signature = signature
    ),
    file.path(cfg$output_dir, "Diagnostics_settings_used.rds")
  )

  cat("\n========================================\n")
  cat("DIAGNOSTICS COMPLETE\n")
  cat("========================================\n")
  cat("Mean Delta WAIC (no-spatial+IID - spatial+IID): ",
      round(mean(model_comparison$delta_waic_nospatial_minus_spatial), 2),
      "\n", sep = "")
  cat("Outputs written to: ", cfg$output_dir, "\n", sep = "")

  invisible(cfg$output_dir)
}

run_latest_diagnostics(diagnostic_settings, country_code)
