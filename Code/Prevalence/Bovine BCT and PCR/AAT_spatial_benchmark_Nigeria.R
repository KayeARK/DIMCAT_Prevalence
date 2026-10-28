# AAT Nigeria: spatial + location IID workflow; covariate selection PER DRAW.
# Run from the ORIGINAL project working directory, in a free R session.
# This is a complete standalone script, not a patch to paste into the old one.
#
# Defaults: ONE draw with complete selection, then four selected-model checks:
# spatial only, weaker IID prior, weaker slope priors, and a 10 km mesh.
# These checks hold covariates fixed; they do not certify a final mesh/model.
# Start from the original project working directory in a free R session.
#
# PRESERVED: original pool of 17 covariates; backward threshold 2, forward threshold
# -0.001, final backward threshold 0.01; separate search for EVERY draw; all
# coefficients/field/hyperparameters re-estimated in each model; geographic
# output grid and the legacy swapped Latitude/Longitude labels; three original
# CSV schemas and model filenames. Range prior remains P(rho < 50 km) = 0.1.
# Spatial sigma prior remains P(sigma > 1) = 0.1. Binomial likelihood unchanged.
# NEW: one Gaussian IID intercept per aggregated sampled coordinate. This can
# separate unstructured local variation from the SPDE field. It is a modelling
# change to assess, not a guarantee of narrower intervals. It is NOT a study
# random effect: records at the same coordinates remain aggregated as before.
# IID prior: P(sigma_iid > 1) = 0.1; sensitivity P(sigma_iid > 2) = 0.1.
#
# PREDICTION TARGETS (same grid, original-format columns in both folders):
# Projections_NGA: prevalence at a NEW location, including a fresh IID effect
# at every grid point and integrating uncertainty in its variance. These are
# intervals for prevalence, not for a future binomial sample proportion.
# Surface_projections_NGA: inverse-logit(X beta + spatial field), with IID=0.
# This conditional surface may be more precise, but is a DIFFERENT estimand.
# It is NOT prevalence averaged over IID heterogeneity. No IID term is silently
# removed from the main intervals. Adding IID may widen the main intervals.
#
# RETAINED FROM THE LAST REVISION: exclude constant/all-missing predictors; centre
# predictors using observed locations; divide continuous predictors by their
# observed SD. Land-cover fractions and binary tsetse retain unit scale (0--1).
# Slopes have explicit N(0, 2.5^2) priors; centred b0 has N(0, 5^2).
# A sensitivity fit uses slope SD 5. These are stated modelling assumptions,
# not automatic settings chosen to minimise interval widths or WAIC.
# No uncertainty is clipped, divided by sqrt(number of draws), or suppressed.
# Missing SELECTED prediction covariates yield NA summaries at the same grid
# rows. This avoids treating missing raster coverage as a measured zero.
# Partly missing observed covariates stop with an audit, rather than imputation.
#
# EFFICIENCY: prepare spatial inputs once; fit candidate models on observations;
# cache compact WAIC results WITHIN EACH DRAW; predict once for the selected
# model; checkpoint every completed candidate and draw; print elapsed time.
# The changed mesh and numerical layout require validation, not an assumption
# that they reproduce the fine-mesh selection decisions.
#
# NEXT RUNS:
# 1. Start with defaults. Read Diagnostics_NGA/Review_notes.txt and the
#    Uncertainty_summary/Sensitivity_comparison files afterwards.
# 2. A new run with mesh_profile="check10km" tests finer resolution.
#    "reference6km" reproduces the earlier corrected mesh settings.
# 3. candidate_mode="full" tests the original full prediction-stack layout
#    on the SAME mesh. Use a fresh output folder when changing either setting.
# 4. To resume/extend an unchanged run, paste its folder into resume_dir and
#    change iteration_columns as needed. Completed outputs are reused, and
#    completed candidate fits within an interrupted draw are reused too.
#    Do not run two processes against the same resume_dir.
# 5. reuse_preparation_dir imports ONLY spatial preparation from an earlier
#    compatible run. Changed priors always require new scores and model fits.
#    resume_dir is for this revised workflow only, with unchanged settings.
#    The default new folder never overwrites previous run results.
#
# API references:
# https://www.r-inla.org/learnmore/docs/reference/inla.mesh.2d.html
# https://www.r-inla.org/learnmore/docs/reference/control.predictor.html
# https://www.r-inla.org/learnmore/docs/reference/control.compute.html
# https://www.r-inla.org/learnmore/docs/reference/f.html
# https://inla.r-inla-download.org/r-inla.org/doc/prior/pc.prec.pdf

settings <- list(
  data_path = "Code/TestSensSpec/adjusted_cases_all_sites.csv",
  iteration_columns = 8L,                  # draw 1; later e.g. 8:12 or 8:1007
  mesh_profile = "pilot20km",              # pilot20km | check10km | reference6km
  candidate_mode = "observations",        # observations | full (comparison)
  output_parent = "Code/Prevalence/Bovine BCT and PCR",
  resume_dir = "",                         # blank creates a NEW timestamped folder
  reuse_preparation_dir = "Code/Prevalence/Bovine BCT and PCR/Spatial_benchmark_NGA_pilot20km_20260925_141439",
  reference_output_dir = "Code/Prevalence/Bovine BCT and PCR/Spatial_benchmark_NGA_pilot20km_20260925_141439",
  coefficient_prior_sd = 2.5,             # explicit assumption, on scales above
  intercept_prior_sd = 5,
  use_location_iid = TRUE,                # FALSE runs the preceding spatial-only model
  iid_prior = c(1, 0.1),                  # P(IID logit-scale SD > 1) = 0.1
  iid_sensitivity_prior = c(2, 0.1),       # one less restrictive IID-prior check
  support_tolerance = 1e-10,
  sensitivity_draws = 1L,                 # checks only on draw 1; no extra selection
  prior_sensitivity_sd = 5,               # one weaker-prior fit per checked draw
  mesh_check_profile = "check10km",       # runs only if finer than primary; "" disables
  diagnostic_draws = 1L,                  # detailed CSV/PDF; summary for every draw
  make_diagnostic_plots = TRUE,
  save_fit_draws = 1L,                     # save full INLA object for draw 1
  num_threads = NULL                      # NULL keeps INLA's existing default
)

atomic_save_rds <- function(object, path) {
  tmp <- tempfile(pattern = ".writing_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
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
# This reproduces the corrected earlier script's thresholds and search order.
select_per_draw <- function(candidates, score, note = function(...) NULL) {
  covs <- as.character(candidates)
  backward <- function(covs, threshold, phase) {
    current <- score(covs)
    repeat {
      if (!length(covs)) return(list(covariates = covs, waic = current, delta = numeric(0)))
      values <- vapply(seq_along(covs), function(j) score(covs[-j]), numeric(1))
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
  first <- backward(covs, 2, "backward_2")
  covs <- first$covariates
  current <- first$waic
  remaining <- setdiff(candidates, covs)
  while (length(remaining)) {
    values <- vapply(seq_along(remaining), function(j) score(c(covs, remaining[j])), numeric(1))
    delta <- values - current
    j <- which.min(delta)
    if (delta[j] >= -0.001) break
    note("forward", "add", remaining[j], delta[j])
    covs <- c(covs, remaining[j])
    remaining <- remaining[-j]
    current <- values[j]
  }
  backward(covs, 0.01, "backward_0.01")
}

# Store only compact scoring results, not hundreds of large INLA objects.
# A separate directory/context is used for each draw; never share scores across draws.
make_draw_scorer <- function(candidates, checkpoint_file, context, evaluate, report) {
  records <- list()
  if (file.exists(checkpoint_file)) {
    saved <- readRDS(checkpoint_file)
    if (!identical(saved$context, context)) stop("Candidate checkpoint context mismatch.")
    records <- saved$records
  }
  hits <- 0L
  new_fits <- 0L
  score <- function(covs) {
    if (anyDuplicated(covs) || !all(covs %in% candidates)) stop("Invalid candidate covariates.")
    # Keep formula order in the key too, avoiding numerical order changes.
    key <- paste0("m_", paste(match(covs, candidates), collapse = "_"))
    if (!is.null(records[[key]])) {
      hits <<- hits + 1L
      report("cache", key, covs, records[[key]])
      return(records[[key]]$waic)
    }
    report("start", key, covs, NULL)
    t0 <- proc.time()[["elapsed"]]
    value <- evaluate(covs)
    if (length(value) != 1L || !is.finite(value)) stop("Non-finite WAIC for ", key)
    record <- list(covariates = covs, waic = as.numeric(value),
                   elapsed_seconds = proc.time()[["elapsed"]] - t0,
                   completed = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"))
    records[[key]] <<- record
    new_fits <<- new_fits + 1L
    atomic_save_rds(list(context = context, records = records), checkpoint_file)
    report("end", key, covs, record)
    record$waic
  }
  list(score = score, records = function() records,
       stats = function() list(cache_hits = hits, new_fits = new_fits))
}

# Compare original-format output files without interpreting swapped column names.
compare_run_outputs <- function(current_dir, reference_dir, draw, destination) {
  pred_name <- paste0("Projections_model_", draw, ".csv")
  cov_name <- paste0("Covariates_model_", draw, ".csv")
  files <- c(file.path(current_dir, "Projections_NGA", pred_name),
             file.path(reference_dir, "Projections_NGA", pred_name),
             file.path(current_dir, "Covariates_NGA", cov_name),
             file.path(reference_dir, "Covariates_NGA", cov_name))
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

build_model_stacks <- function(prepared) {
  n <- nrow(prepared$observation_lonlat); m <- nrow(prepared$prediction_lonlat)
  xo <- transform_covariates(prepared$raw_observation_covariates, prepared$covariate_audit)
  xp <- transform_covariates(prepared$raw_prediction_covariates, prepared$covariate_audit)
  if (any(!is.finite(as.matrix(xo)))) stop("Missing observed covariates: inspect Covariate_support.csv.")
  # Placeholders are ONLY for fitting a common prediction stack. A prediction
  # with a missing selected covariate is exported as NA, never as this placeholder.
  for (nm in names(xp)) xp[[nm]][!is.finite(xp[[nm]])] <- 0
  idx <- INLA::inla.spde.make.index("s", prepared$spde$n.spde)
  eo <- data.frame(b0 = rep(1, n), xo)
  ep <- data.frame(b0 = rep(1, m), xp)
  es <- ep
  if (isTRUE(prepared$use_location_iid)) {
    ii <- iid_indices(n, m)
    eo$site_iid <- ii$observed; ep$site_iid <- ii$new_location; es$site_iid <- ii$surface
  }
  est <- INLA::inla.stack(tag = "est", data = list(y = rep(0, n), numtrials = prepared$numtrials),
    A = list(1, prepared$A), effects = list(eo, s = idx))
  pred <- INLA::inla.stack(tag = "pred", data = list(y = rep(NA_real_, m), numtrials = rep(NA_real_, m)),
    A = list(1, prepared$Ap), effects = list(ep, s = idx))
  if (isTRUE(prepared$use_location_iid)) {
    surf <- INLA::inla.stack(tag = "surface", data = list(y = rep(NA_real_, m), numtrials = rep(NA_real_, m)),
      A = list(1, prepared$Ap), effects = list(es, s = idx))
    full <- INLA::inla.stack(est, pred, surf)
  } else full <- INLA::inla.stack(est, pred)
  prepared$observation_data <- INLA::inla.stack.data(est)
  prepared$full_data <- INLA::inla.stack.data(full)
  prepared$observation_A <- INLA::inla.stack.A(est)
  prepared$full_A <- INLA::inla.stack.A(full)
  prepared$observation_index <- INLA::inla.stack.index(est, "est")$data
  prepared$full_observation_index <- INLA::inla.stack.index(full, "est")$data
  prepared$prediction_index <- INLA::inla.stack.index(full, "pred")$data
  prepared$surface_prediction_index <- if (isTRUE(prepared$use_location_iid))
    INLA::inla.stack.index(full, "surface")$data else prepared$prediction_index
  prepared
}

fit_model <- function(prepared, y, covariates, layout, slope_sd, intercept_sd, threads, final) {
  formula <- model_formula(covariates, prepared$spde, isTRUE(prepared$use_location_iid), prepared$iid_prior)
  dat <- if (layout == "observations") prepared$observation_data else prepared$full_data
  ind <- if (layout == "observations") prepared$observation_index else prepared$full_observation_index
  dat$y[ind] <- y
  args <- list(formula = formula, data = dat, family = "binomial", Ntrials = dat$numtrials,
    control.fixed = fixed_prior_control(slope_sd, intercept_sd),
    control.compute = list(waic = TRUE, dic = final, return.marginals.predictor = final),
    control.predictor = list(link = 1, compute = TRUE,
      A = if (layout == "observations") prepared$observation_A else prepared$full_A))
  if (!is.null(threads)) args$num.threads <- threads
  ans <- do.call(INLA::inla, args)
  if (length(ans$waic$waic) != 1L || !is.finite(ans$waic$waic)) stop("INLA returned non-finite WAIC.")
  ans
}

summarise_fit <- function(res, prepared, covariates, seconds) {
  valid <- support_at_predictions(prepared, covariates)$valid_prediction
  extract <- function(i) {
    pp <- res$summary.fitted.values[i, , drop = FALSE]
    eta <- res$summary.linear.predictor[i, , drop = FALSE]
    stopifnot(nrow(pp) == length(valid), all(is.finite(pp$mean)),
      all(is.finite(pp[["0.025quant"]])), all(is.finite(pp[["0.975quant"]])),
      all(pp$mean >= 0 & pp$mean <= 1),
      all(pp[["0.025quant"]] >= 0 & pp[["0.975quant"]] <= 1),
      all(pp[["0.025quant"]] <= pp[["0.975quant"]]),
      all(is.finite(eta$sd) & eta$sd >= 0))
    pp[!valid, ] <- NA_real_; eta[!valid, ] <- NA_real_
    list(prevalence = pp, linear_predictor = eta)
  }
  main <- extract(prepared$prediction_index)
  surface <- extract(prepared$surface_prediction_index)
  site <- iid_sd <- NULL
  if (isTRUE(prepared$use_location_iid)) {
    iid_sd <- iid_sd_quantiles(res$summary.hyperpar)
    all_site <- res$summary.random$site_iid
    at <- match(seq_len(nrow(prepared$observation_lonlat)), all_site$ID)
    if (anyNA(at)) stop("Observed IID summaries cannot be matched to sampled locations.")
    site <- all_site[at, , drop = FALSE]
  }
  list(covariates = covariates, prevalence = main$prevalence, linear_predictor = main$linear_predictor,
    surface_prevalence = surface$prevalence, surface_linear_predictor = surface$linear_predictor,
    use_location_iid = isTRUE(prepared$use_location_iid), iid_sd = iid_sd, observed_iid = site,
    valid_prediction = valid, spatial_nodes = res$summary.random$s,
    fixed = res$summary.fixed, hyperparameters = res$summary.hyperpar,
    observed_prevalence = res$summary.fitted.values[prepared$full_observation_index, , drop = FALSE],
    waic_full = res$waic$waic, elapsed_seconds = seconds)
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

prediction_comparison <- function(a, b, label) {
  pa <- a$prevalence; pb <- b$prevalence
  valid <- is.finite(pa$mean) & is.finite(pb$mean)
  da <- pa[["0.975quant"]] - pa[["0.025quant"]]
  db <- pb[["0.975quant"]] - pb[["0.025quant"]]
  sa <- if (is.null(a$surface_prevalence)) pa else a$surface_prevalence
  sb <- if (is.null(b$surface_prevalence)) pb else b$surface_prevalence
  data.frame(comparison = label, comparable_points = sum(valid),
    mean_absolute_prevalence_change = if (any(valid)) mean(abs(pb$mean[valid] - pa$mean[valid])) else NA_real_,
    maximum_absolute_prevalence_change = if (any(valid)) max(abs(pb$mean[valid] - pa$mean[valid])) else NA_real_,
    median_interval_width_primary = if (any(valid)) median(da[valid]) else NA_real_,
    median_interval_width_comparison = if (any(valid)) median(db[valid]) else NA_real_,
    median_surface_width_primary = if (any(valid)) median(sa[["0.975quant"]][valid] - sa[["0.025quant"]][valid]) else NA_real_,
    median_surface_width_comparison = if (any(valid)) median(sb[["0.975quant"]][valid] - sb[["0.025quant"]][valid]) else NA_real_,
    mean_absolute_surface_change = if (any(valid)) mean(abs(sa$mean[valid] - sb$mean[valid])) else NA_real_,
    waic_primary = a$waic_full, waic_comparison = b$waic_full,
    selection_repeated = FALSE)
}

run_aat_benchmark <- function(cfg) {
  suppressPackageStartupMessages({
    library(INLA); library(sf); library(sp); library(afrilearndata)
    library(concaveman); library(geodata); library(raster)
  })
  countries_to_infer <- "Nigeria"
  profiles <- list(pilot20km = c(20000, 50000), check10km = c(10000, 25000),
                   reference6km = c(6000, 15000))
  if (!cfg$mesh_profile %in% names(profiles)) stop("Unknown mesh_profile.")
  if (!cfg$candidate_mode %in% c("observations", "full")) stop("Unknown candidate_mode.")
  if (nzchar(cfg$mesh_check_profile) && !cfg$mesh_check_profile %in% names(profiles))
    stop("Unknown mesh_check_profile.")
  mesh_check_enabled <- nzchar(cfg$mesh_check_profile) &&
    all(profiles[[cfg$mesh_check_profile]] < profiles[[cfg$mesh_profile]])
  for (nm in c("coefficient_prior_sd", "intercept_prior_sd", "prior_sensitivity_sd", "support_tolerance")) {
    if (length(cfg[[nm]]) != 1L || !is.finite(cfg[[nm]]) || cfg[[nm]] <= 0) stop("Invalid ", nm)
  }
  if (!is.logical(cfg$use_location_iid) || length(cfg$use_location_iid) != 1L || is.na(cfg$use_location_iid))
    stop("use_location_iid must be TRUE or FALSE.")
  for (nm in c("iid_prior", "iid_sensitivity_prior")) {
    z <- cfg[[nm]]
    if (!is.numeric(z) || length(z) != 2L || any(!is.finite(z)) || z[1] <= 0 || z[2] <= 0 || z[2] >= 1)
      stop(nm, " must be c(positive SD threshold, probability strictly between 0 and 1).")
  }
  for (nm in c("sensitivity_draws", "diagnostic_draws", "save_fit_draws")) {
    v <- cfg[[nm]]
    if (!is.numeric(v) || anyNA(v) || any(v < 1 | v != as.integer(v))) stop("Invalid ", nm)
  }
  spatial_crs <- paste("+proj=aeqd +lat_0=9 +lon_0=8 +datum=WGS84",
                       "+units=m +no_defs +type=crs")
  mesh_cutoff_m <- 3000
  mesh_max_edge_m <- profiles[[cfg$mesh_profile]]
  mesh_offset_m <- c(30000, 60000) # unchanged while comparing mesh resolution
  range_prior_m <- c(50000, 0.1)
  sigma_prior <- c(1, 0.1)
  if (sf::st_is_longlat(sf::st_crs(spatial_crs)) ||
      !sf::st_crs(spatial_crs)$units_gdal %in% c("metre", "meter", "metres", "meters", "m")) {
    stop("spatial_crs must use projected metres.")
  }
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
  data_original$longitude <- as.numeric(data_original$longitude)
  data_original$latitude <- as.numeric(data_original$latitude)
  keep <- is.finite(data_original$longitude) & is.finite(data_original$latitude)
  if (any(!keep)) message("Removing ", sum(!keep), " rows with missing/non-finite coordinates.")
  data_original <- data_original[keep, , drop = FALSE]
  if (!nrow(data_original) || any(abs(data_original$longitude) > 180 | abs(data_original$latitude) > 90)) {
    stop("No valid geographic coordinates.")
  }
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
  signature <- list(workflow_version = "aat_iid_20260927_v1", input_md5 = input_md5,
                    spatial_crs = spatial_crs, cutoff_m = mesh_cutoff_m,
                    max_edge_m = mesh_max_edge_m, offset_m = mesh_offset_m,
                    range_prior_m = range_prior_m, sigma_prior = sigma_prior,
                    candidate_mode = cfg$candidate_mode, package_versions = package_versions,
                    R_version = as.character(getRversion()), num_threads = cfg$num_threads,
                    coefficient_prior_sd = cfg$coefficient_prior_sd,
                    intercept_prior_sd = cfg$intercept_prior_sd,
                    use_location_iid = cfg$use_location_iid,
                    iid_prior = cfg$iid_prior, iid_sensitivity_prior = cfg$iid_sensitivity_prior,
                    primary_prediction_target = "new_location_including_IID_if_enabled",
                    surface_prediction_target = "inverse_logit_Xbeta_plus_SPDE_IID_zero",
                    support_tolerance = cfg$support_tolerance,
                    prior_sensitivity_sd = cfg$prior_sensitivity_sd,
                    mesh_check_profile = cfg$mesh_check_profile,
                    missing_prediction_policy = "NA_if_any_selected_covariate_missing")
  # Draw numbers are deliberately excluded so an unchanged run can be extended.
  if (nzchar(cfg$resume_dir)) {
    output_root <- path.expand(cfg$resume_dir)
    if (!dir.exists(output_root)) stop("resume_dir does not exist.")
  } else {
    base <- paste0(if (cfg$use_location_iid) "Spatial_IID_NGA_" else "Spatial_only_NGA_",
                   cfg$mesh_profile, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    output_root <- file.path(cfg$output_parent, base)
    suffix <- 0L
    while (file.exists(output_root)) {
      suffix <- suffix + 1L
      output_root <- file.path(cfg$output_parent, paste0(base, "_", suffix))
    }
    if (!dir.create(output_root, recursive = TRUE)) stop("Cannot create output directory.")
  }
  output_root <- normalizePath(output_root, mustWork = TRUE)
  lock_dir <- file.path(output_root, ".run_lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) {
    stop("Run folder is locked. Ensure no process is using it. If a previous R session was killed, ",
         "remove its .run_lock directory only after confirming it is no longer running.")
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  manifest_file <- file.path(output_root, "run_manifest.rds")
  if (file.exists(manifest_file)) {
    if (!identical(readRDS(manifest_file), signature)) {
      stop("Input, model, mesh, priors or package versions differ from this run. Leave resume_dir blank. ",
           "Use reuse_preparation_dir to import geometry/covariates from an older run without reusing its fits.")
    }
  } else {
    if (nzchar(cfg$resume_dir)) stop("resume_dir has no matching run manifest.")
    atomic_save_rds(signature, manifest_file)
  }
  output_dirs <- setNames(file.path(output_root, c("WAIC_increase_NGA", "Covariates_NGA",
                          "Projections_NGA", "Diagnostics_NGA", "Checkpoints_NGA", "Surface_projections_NGA")),
                          c("waic", "covariates", "projections", "diagnostics", "checkpoints", "surface"))
  for (d in output_dirs) if (!dir.exists(d) && !dir.create(d)) stop("Cannot create ", d)
  log_file <- file.path(output_dirs[["diagnostics"]], "progress.log")
  log_message <- function(...) {
    line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
    cat(line, "\n"); flush.console()
    cat(line, "\n", file = log_file, append = TRUE)
  }
  log_message("RESULTS: ", output_root)
  log_message("Mesh profile: ", cfg$mesh_profile, "; candidate layout: ", cfg$candidate_mode,
              "; draws requested: ", paste(iteration_columns - 7L, collapse = ", "))
  log_message("Location IID: ", cfg$use_location_iid,
              ". Main projections include new-location heterogeneity; Surface_projections_NGA sets IID=0.")
  if (length(cfg$sensitivity_draws) && nzchar(cfg$mesh_check_profile) && !mesh_check_enabled)
    log_message("Additional mesh fit skipped: mesh_check_profile is not finer than the primary mesh.")
  writeLines(c("This is a revision/sensitivity run; validate sensitivity before final inference.",
               "Each diagnostic draw has its own complete covariate search and its own WAIC cache.",
               "INLA distances are metres. Output coordinates remain geographic degrees.",
               "Latitude contains longitude and Longitude contains latitude, intentionally.",
               "Resume reuses the SAVED covariates and mesh; clear reuse_preparation_dir for fresh raster inputs.",
               "Compare both selection/WAIC and predictions across meshes; assess buffer sensitivity separately.",
               paste("Projection:", spatial_crs), paste("Range prior:", paste(range_prior_m, collapse = ", "))),
             file.path(output_dirs[["diagnostics"]], "README.txt"))
  atomic_save_rds(cfg, file.path(output_root, "settings_used.rds"))
  writeLines(c(
    "LOCATION IID REVISION: INTERPRET BOTH PREDICTION TARGETS",
    "Same binomial likelihood, SPDE field, geographical grid and selection thresholds as the corrected benchmark.",
    "A full independent covariate search is performed for every primary diagnostic draw.",
    paste0("Location IID enabled: ", cfg$use_location_iid, "."),
    "When enabled: logit(p_i) = b0 + X_i beta + w(s_i) + u_i, with independent u_i ~ Normal(0, sigma_iid^2).",
    "One u_i per aggregated coordinate, not per study or original record. Existing aggregation is unchanged.",
    "IID effects have no sum-to-zero constraint; unseen prediction units do not constrain observed effects.",
    paste0("IID prior: P(sigma_iid > ", cfg$iid_prior[1], ") = ", cfg$iid_prior[2], " on the logit scale."),
    paste0("IID sensitivity prior: P(sigma_iid > ", cfg$iid_sensitivity_prior[1], ") = ", cfg$iid_sensitivity_prior[2], "."),
    "Projections_NGA retains the original schema and reports NEW-LOCATION prevalence including IID heterogeneity.",
    "Every prediction cell receives a distinct unseen IID index, independent conditional on its variance.",
    "The main intervals integrate hyperparameter uncertainty through INLA; IID variance is not plugged in at its mean.",
    "These are prevalence intervals, not predictive intervals for an observed binomial sample fraction.",
    "Surface_projections_NGA reports inverse-logit(b0 + X beta + w), conditional on location IID=0.",
    "Surface intervals retain joint coefficient/field uncertainty, including covariance through the INLA predictor.",
    "The zero-IID surface is a different estimand: it is NOT mean prevalence averaged over location heterogeneity.",
    "Narrower surface intervals must not be described as narrower new-location prevalence intervals.",
    "An IID effect may reduce spatial-field uncertainty but may widen total new-location uncertainty.",
    "IID_SD files transform posterior precision quantiles, reversing endpoints; no plug-in mean SD is reported.",
    "Predictors constant or entirely missing at fitted locations are excluded independently of the response.",
    "Excluded effects are fixed to zero by omission; this is an explicit assumption, not evidence of no effect elsewhere.",
    "Rare but varying predictors are retained and flagged. Missing values are never used as measured zeros.",
    "Observed-location means are subtracted from all eligible predictors. Continuous predictors use observed SDs.",
    "Land-cover fractions and binary tsetse use unit scale, so priors refer to a 0-to-1 contrast.",
    paste0("Slope prior: independent Normal(0, ", cfg$coefficient_prior_sd, "^2); centred b0: Normal(0, ", cfg$intercept_prior_sd, "^2)."),
    "The coefficient priors change the model; standardisation alone would not reduce uncertainty.",
    "SPDE priors unchanged: P(range < 50 km)=0.1; P(spatial SD > 1)=0.1.",
    "Fixed_effects_model files report coefficients on the centred/scaled modelling units in Covariate_support.csv.",
    "Spatial-only, IID-prior, slope-prior and mesh checks use primary selected covariates without repeating selection.",
    "The spatial-only check isolates the effect of adding IID for those covariates; it is not a fresh spatial-only search.",
    "Both new-location and zero-IID surface uncertainty are reported in sensitivity comparisons.",
    "If sensitivity is material, run the full search on the alternate configuration before final inference.",
    "No automatic choice is made from interval width or WAIC; narrower intervals alone do not validate a model.",
    "Conditional WAIC with location effects is not a substitute for validation at entirely held-out locations.",
    "Prediction rows with missing selected covariates retain their coordinates but have NA summaries.",
    "Other locations outside observed covariate ranges remain predictions and are flagged as extrapolations.",
    "Per-draw intervals condition on that adjusted-count draw and selected model. Retain within-draw uncertainty when pooling draws.",
    "Distances are approximate distances in the existing Nigeria-centred metre projection.",
    "Grid averages are unweighted descriptions, not cattle-population-weighted national prevalence.",
    "Observed_fit files are IN-SAMPLE checks, not held-out predictive validation.",
    "CSV Latitude/Longitude headings retain the pre-existing swap deliberately."),
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
    r_elv <- elevation_30s(country = "NGA", path = "Data/Covariates")
    data$elevation<-terra::extract(r_elv, data[, c("longitude", "latitude")])$NGA_elv_msk
    #PRECIPITATION
    r_prec<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "prec")
    data$precipitation<-(terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_1+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_2+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_3+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_4+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_5+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_6+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_7+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_8+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_9+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_10+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_11+
    terra::extract(r_prec, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_prec_12)/12
    #AVERAGE TEMPERATURE
    r_tavg<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tavg")
    data$tavg<-(terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_1+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_2+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_3+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_4+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_5+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_6+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_7+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_8+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_9+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_10+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_11+
    terra::extract(r_tavg, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tavg_12)/12
    #MINIMUM TEMPERATURE
    r_tmin<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tmin")
    data$tmin<-(terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_1+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_2+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_3+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_4+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_5+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_6+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_7+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_8+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_9+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_10+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_11+
    terra::extract(r_tmin, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmin_12)/12
    #MAXIMUM TEMPERATURE
    r_tmax<- worldclim_country(country = "NGA",path = "Data/Covariates", var = "tmax")
    data$tmax<-(terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_1+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_2+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_3+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_4+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_5+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_6+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_7+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_8+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_9+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_10+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_11+
    terra::extract(r_tmax, data[, c("longitude", "latitude")])$NGA_wc2.1_30s_tmax_12)/12
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

    ra_elv <- aggregate(r_elv, fact = 5, fun = mean)
    dp <- as.data.frame(coop)
    #coop <- terra::crds(ra)
    dp$elevation <- terra::extract(ra_elv, coop)[, 1]
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
    covariates<-c("elevation","precipitation","tavg","tmax","human_fp","pop_den","tree","grassland",
    "shrub","cropland","built","bare","water","wetland","mangrove","cattle","tsetse")
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
      paths <- file.path(source_dir, c("run_manifest.rds", "Checkpoints_NGA/spatial_preparation.rds"))
      if (!all(file.exists(paths))) stop("Incomplete reuse_preparation_dir. Choose a complete run, or set it to ''.")
      source_manifest <- readRDS(paths[1])
      required <- c("input_md5", "spatial_crs", "cutoff_m", "offset_m", "range_prior_m", "sigma_prior", "package_versions")
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
    pool <- c("elevation", "precipitation", "tavg", "tmax", "human_fp", "pop_den", "tree", "grassland",
              "shrub", "cropland", "built", "bare", "water", "wetland", "mangrove", "cattle", "tsetse")
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
    prepared <- build_model_stacks(prepared)
    prepared$preparation_seconds <- proc.time()[["elapsed"]] - t0
    atomic_save_rds(prepared, preparation_file)
    log_message("Preparation saved; elapsed ", round(prepared$preparation_seconds, 1), " seconds.")
  }
  atomic_write_csv(prepared$covariate_audit, file.path(output_dirs[["diagnostics"]], "Covariate_support.csv"))
  excluded <- prepared$covariate_audit$covariate[!prepared$covariate_audit$eligible]
  log_message("Excluded unestimable covariates: ", if (length(excluded)) paste(excluded, collapse = ", ") else "none")
  xo <- as.matrix(transform_covariates(prepared$raw_observation_covariates, prepared$covariate_audit))
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
    targets <- rbind(uncertainty_summary(d, draw, "new_location"),
                     uncertainty_summary(sd, draw, "surface_IID_zero"))
    atomic_write_csv(targets[targets$group == "all", ], file.path(output_dirs[["diagnostics"]],
      paste0("Prediction_targets_model_", draw, ".csv")))
    if (saved$use_location_iid) atomic_write_csv(saved$iid_sd, file.path(output_dirs[["diagnostics"]],
      paste0("IID_SD_model_", draw, ".csv")))
    if (!draw %in% cfg$diagnostic_draws) return(invisible(NULL))
    atomic_write_csv(d, file.path(output_dirs[["diagnostics"]], paste0("Prediction_diagnostics_model_", draw, ".csv")))
    atomic_write_csv(sd, file.path(output_dirs[["diagnostics"]], paste0("Prediction_diagnostics_surface_model_", draw, ".csv")))
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

  finer_prepared <- NULL
  spatial_only_prepared <- NULL
  run_sensitivity_checks <- function(saved, draw, y, draw_dir) {
    if (!draw %in% cfg$sensitivity_draws) return(invisible(NULL))
    spec_for <- function(profile = cfg$mesh_profile, sd = cfg$coefficient_prior_sd,
                         use_iid = cfg$use_location_iid, iid_prior = cfg$iid_prior)
      list(profile = profile, sd = sd, use_iid = use_iid, iid_prior = iid_prior)
    scenarios <- list()
    if (cfg$use_location_iid) {
      scenarios$spatial_only <- spec_for(use_iid = FALSE)
      scenarios$weaker_iid_prior <- spec_for(iid_prior = cfg$iid_sensitivity_prior)
    }
    scenarios$weaker_prior <- spec_for(sd = cfg$prior_sensitivity_sd)
    if (mesh_check_enabled)
      scenarios$finer_mesh <- spec_for(profile = cfg$mesh_check_profile)
    comparisons <- list()
    for (label in names(scenarios)) {
      spec <- scenarios[[label]]
      path <- file.path(draw_dir, paste0("sensitivity_", label, ".rds"))
      ctx <- list(signature = signature, draw = draw, response = y,
                  covariates = saved$covariates, scenario = spec)
      if (file.exists(path)) {
        cache <- readRDS(path)
        if (!identical(cache$context, ctx)) stop("Sensitivity checkpoint mismatch.")
        other <- cache$result
        log_message("DRAW ", draw, " | REUSE sensitivity ", label)
      } else {
        pp <- prepared
        if (spec$profile != cfg$mesh_profile) {
          if (is.null(finer_prepared)) {
            mp <- file.path(output_dirs[["checkpoints"]], "sensitivity_mesh.rds")
            if (file.exists(mp)) finer_prepared <<- readRDS(mp) else {
              log_message("Preparing sensitivity mesh: ", spec$profile, " (no raster downloads).")
              finer_prepared <<- build_model_stacks(remesh_prepared(prepared, spec$profile))
              atomic_save_rds(finer_prepared, mp)
            }
          }
          pp <- finer_prepared
        }
        if (!identical(spec$use_iid, pp$use_location_iid)) {
          if (is.null(spatial_only_prepared)) {
            pp$use_location_iid <- spec$use_iid
            spatial_only_prepared <<- build_model_stacks(pp)
          }
          pp <- spatial_only_prepared
        }
        pp$iid_prior <- spec$iid_prior
        log_message("DRAW ", draw, " | SENSITIVITY ", label, " START | selected covariates fixed; slope SD ", spec$sd)
        t0 <- proc.time()[["elapsed"]]
        fit <- fit_model(pp, y, saved$covariates, "full", spec$sd, cfg$intercept_prior_sd, cfg$num_threads, TRUE)
        other <- summarise_fit(fit, pp, saved$covariates, proc.time()[["elapsed"]] - t0)
        rm(fit)
        atomic_save_rds(list(context = ctx, result = other), path)
        log_message("DRAW ", draw, " | SENSITIVITY ", label, " END | ", round(other$elapsed_seconds, 1), " s")
      }
      comparisons[[label]] <- prediction_comparison(saved, other, label)
      comparisons[[label]]$comparison_mesh <- spec$profile
      comparisons[[label]]$comparison_slope_sd <- spec$sd
      comparisons[[label]]$comparison_has_iid <- spec$use_iid
      comparisons[[label]]$comparison_iid_prior_threshold <- if (spec$use_iid) spec$iid_prior[1] else NA_real_
      comparisons[[label]]$comparison_iid_prior_probability <- if (spec$use_iid) spec$iid_prior[2] else NA_real_
      comparisons[[label]]$primary_iid_sd_median <- if (saved$use_location_iid) saved$iid_sd$median else NA_real_
      comparisons[[label]]$comparison_iid_sd_median <- if (other$use_location_iid) other$iid_sd$median else NA_real_
      comparisons[[label]]$fit_seconds <- other$elapsed_seconds
      comparisons[[label]]$comparison_mesh_vertices <- nrow(other$spatial_nodes)
      atomic_write_csv(uncertainty_summary(uncertainty_rows(other, prepared), draw, label),
        file.path(output_dirs[["diagnostics"]], paste0("Uncertainty_summary_", label, "_model_", draw, ".csv")))
      atomic_write_csv(uncertainty_summary(uncertainty_rows(surface_result(other), prepared), draw, paste0(label, "_surface_IID_zero")),
        file.path(output_dirs[["diagnostics"]], paste0("Uncertainty_summary_surface_", label, "_model_", draw, ".csv")))
      hp <- other$hyperparameters
      atomic_write_csv(data.frame(parameter = rownames(hp), hp, row.names = NULL),
        file.path(output_dirs[["diagnostics"]], paste0("Hyperparameters_", label, "_model_", draw, ".csv")))
      if (other$use_location_iid) atomic_write_csv(other$iid_sd, file.path(output_dirs[["diagnostics"]],
        paste0("IID_SD_", label, "_model_", draw, ".csv")))
      atomic_write_csv(data.frame(Latitude = prepared$prediction_lonlat[, 1], Longitude = prepared$prediction_lonlat[, 2],
        primary_mean = saved$prevalence$mean, comparison_mean = other$prevalence$mean,
        primary_width = saved$prevalence[["0.975quant"]] - saved$prevalence[["0.025quant"]],
        comparison_width = other$prevalence[["0.975quant"]] - other$prevalence[["0.025quant"]],
        primary_surface_mean = saved$surface_prevalence$mean,
        comparison_surface_mean = other$surface_prevalence$mean,
        primary_surface_width = saved$surface_prevalence[["0.975quant"]] - saved$surface_prevalence[["0.025quant"]],
        comparison_surface_width = other$surface_prevalence[["0.975quant"]] - other$surface_prevalence[["0.025quant"]]),
        file.path(output_dirs[["diagnostics"]], paste0("Prediction_comparison_", label, "_model_", draw, ".csv")))
    }
    atomic_write_csv(do.call(rbind, comparisons), file.path(output_dirs[["diagnostics"]],
                                                         paste0("Sensitivity_comparison_model_", draw, ".csv")))
    invisible(NULL)
  }

  run_one_draw <- function(iteration) {
    draw <- iteration - 7L
    draw_start <- proc.time()[["elapsed"]]
    draw_dir <- file.path(output_dirs[["checkpoints"]], sprintf("draw_%04d", draw))
    if (!dir.exists(draw_dir) && !dir.create(draw_dir)) stop("Cannot create draw checkpoint directory.")
    done_file <- file.path(draw_dir, "completed.rds")
    result_paths <- c(
      file.path(output_dirs[["waic"]], paste0("WAIC_increase_model_", draw, ".csv")),
      file.path(output_dirs[["covariates"]], paste0("Covariates_model_", draw, ".csv")),
      file.path(output_dirs[["projections"]], paste0("Projections_model_", draw, ".csv"))
    )
    extra_paths <- c(
      file.path(output_dirs[["surface"]], paste0("Projections_model_", draw, ".csv")),
      file.path(output_dirs[["diagnostics"]], paste0("Uncertainty_summary_model_", draw, ".csv")),
      file.path(output_dirs[["diagnostics"]], paste0("Uncertainty_summary_surface_model_", draw, ".csv")),
      file.path(output_dirs[["diagnostics"]], paste0("Prediction_targets_model_", draw, ".csv")))
    if (cfg$use_location_iid) extra_paths <- c(extra_paths,
      file.path(output_dirs[["diagnostics"]], paste0("IID_SD_model_", draw, ".csv")))
    if (draw %in% cfg$diagnostic_draws) extra_paths <- c(extra_paths,
      file.path(output_dirs[["diagnostics"]], paste0("Prediction_diagnostics_model_", draw, ".csv")),
      file.path(output_dirs[["diagnostics"]], paste0("Prediction_diagnostics_surface_model_", draw, ".csv")))
    if (draw %in% cfg$sensitivity_draws) extra_paths <- c(extra_paths,
      file.path(output_dirs[["diagnostics"]], paste0("Sensitivity_comparison_model_", draw, ".csv")))
    if (file.exists(done_file) && all(file.exists(c(result_paths, extra_paths)))) {
      log_message("DRAW ", draw, " already complete; reusing saved outputs.")
      return(readRDS(done_file))
    }
    log_message("DRAW ", draw, " START: independent covariate search.")
    counts <- aggregate_draw_counts(data_original, iteration, prepared$coordinate_mapping)
    y <- counts[prepared$kept_location_rows]
    stopifnot(length(y) == length(prepared$observation_index),
              length(y) == length(prepared$full_observation_index))
    candidates <- prepared$candidate_covariates
    fit_number <- 0L
    fit_raw <- function(covs, layout, final = FALSE) {
      fit_model(prepared, y, covs, layout, cfg$coefficient_prior_sd,
                cfg$intercept_prior_sd, cfg$num_threads, final)
    }
    report <- function(event, key, covs, record) {
      if (event == "start") {
        fit_number <<- fit_number + 1L
        log_message("DRAW ", draw, " | FIT ", fit_number, " START | ", length(covs),
                    " covariates | ", key)
      } else if (event == "end") {
        log_message("DRAW ", draw, " | FIT ", fit_number, " END | ",
                    round(record$elapsed_seconds, 1), " s | WAIC ", sprintf("%.5f", record$waic))
      } else {
        log_message("DRAW ", draw, " | REUSE | ", key, " | WAIC ", sprintf("%.5f", record$waic))
      }
    }
    context <- list(model_signature = signature, input_md5 = input_md5, draw = draw, response = y,
                    candidates = candidates, candidate_mode = cfg$candidate_mode)
    scorer <- make_draw_scorer(
      candidates, file.path(draw_dir, "candidate_scores.rds"), context,
      evaluate = function(covs) fit_raw(covs, cfg$candidate_mode)$waic$waic,
      report = report
    )
    note <- function(phase, action, covariate, delta) {
      log_message("DRAW ", draw, " | ", phase, " | ", action, " ", covariate,
                  " | delta WAIC ", sprintf("%.5f", delta))
    }
    selected <- select_per_draw(candidates, scorer$score, note)
    atomic_save_rds(selected, file.path(draw_dir, "selected_model.rds"))
    covariates <- selected$covariates
    waic_increase <- selected$delta
    log_message("DRAW ", draw, " | SEARCH COMPLETE | selected ", length(covariates),
                " covariates; WAIC ", sprintf("%.5f", selected$waic))

    # Save every score for examining selection sensitivity across meshes/layouts.
    records <- scorer$records()
    score_table <- do.call(rbind, lapply(names(records), function(key) {
      r <- records[[key]]
      data.frame(key = key, covariates = paste(r$covariates, collapse = " + "),
                 waic = r$waic, elapsed_seconds = r$elapsed_seconds, completed = r$completed)
    }))
    atomic_write_csv(score_table, file.path(output_dirs[["diagnostics"]],
                                         paste0("Candidate_scores_model_", draw, ".csv")))

    # Final predictions are checkpointed before the CSV exports. If interrupted
    # during writing, a resume can finish exporting without repeating this fit.
    prediction_checkpoint <- file.path(draw_dir, "final_predictions.rds")
    if (file.exists(prediction_checkpoint)) {
      saved_predictions <- readRDS(prediction_checkpoint)
      if (!identical(saved_predictions$covariates, covariates)) stop("Selected-model checkpoint mismatch.")
      log_message("DRAW ", draw, " | REUSE final prediction summaries.")
    } else {
      log_message("DRAW ", draw, " | FINAL PREDICTION FIT START.")
      t0 <- proc.time()[["elapsed"]]
      res <- fit_raw(covariates, "full", final = TRUE)
      prediction_seconds <- proc.time()[["elapsed"]] - t0
      saved_predictions <- summarise_fit(res, prepared, covariates, prediction_seconds)
      if (draw %in% cfg$save_fit_draws) {
        atomic_save_rds(res, file.path(output_dirs[["diagnostics"]], paste0("Final_fit_model_", draw, ".rds")))
      }
      atomic_save_rds(saved_predictions, prediction_checkpoint)
      rm(res)
      log_message("DRAW ", draw, " | FINAL PREDICTION FIT END | ", round(prediction_seconds, 1), " s.")
    }
    waic_layout_difference <- saved_predictions$waic_full - selected$waic
    log_message("DRAW ", draw, " | full-stack minus selection WAIC = ",
                sprintf("%.6f", waic_layout_difference), ". Check numerical layout sensitivity.")
    if (cfg$candidate_mode == "observations" && abs(waic_layout_difference) >= 0.01) {
      log_message("CHECK: WAIC layout difference reaches 0.01. Compare candidate_mode='full' on this mesh ",
                  "before interpreting near-threshold selection decisions.")
    }

    # Original output schemas, order, labels and model-number convention.
    coop <- prepared$prediction_lonlat
    dpm <- projection_table(saved_predictions$prevalence, coop)
    atomic_write_csv(data.frame(covariates = covariates, waic_increase = waic_increase), result_paths[1])
    atomic_write_csv(data.frame(covariates = covariates), result_paths[2])
    atomic_write_csv(dpm, result_paths[3])
    atomic_write_csv(projection_table(saved_predictions$surface_prevalence, coop),
      file.path(output_dirs[["surface"]], paste0("Projections_model_", draw, ".csv")))
    hp <- saved_predictions$hyperparameters
    atomic_write_csv(data.frame(parameter = rownames(hp), hp, row.names = NULL),
                     file.path(output_dirs[["diagnostics"]], paste0("Hyperparameters_model_", draw, ".csv")))
    fixed <- saved_predictions$fixed
    atomic_write_csv(data.frame(parameter = rownames(fixed), fixed, row.names = NULL),
                     file.path(output_dirs[["diagnostics"]], paste0("Fixed_effects_model_", draw, ".csv")))
    # Spatial field uncertainty is saved at mesh nodes. Do NOT interpolate its
    # marginal SD with Ap; grid-point variance also requires cross-covariances.
    nodes <- saved_predictions$spatial_nodes
    mesh_rows <- node_mesh_rows(nodes$ID, prepared$mesh$n)
    atomic_write_csv(data.frame(x_m = prepared$mesh$loc[mesh_rows, 1],
                               y_m = prepared$mesh$loc[mesh_rows, 2], nodes),
                     file.path(output_dirs[["diagnostics"]], paste0("Spatial_nodes_model_", draw, ".csv")))
    write_uncertainty_diagnostics(saved_predictions, draw, y)
    med_width <- function(pp) median(pp[["0.975quant"]] - pp[["0.025quant"]], na.rm = TRUE)
    log_message("DRAW ", draw, " | median 95% widths: new location ",
      round(100 * med_width(saved_predictions$prevalence), 2), " percentage points; IID=0 surface ",
      round(100 * med_width(saved_predictions$surface_prevalence), 2), " percentage points. Different targets.")
    log_message("DRAW ", draw, " | valid predictions: ", sum(saved_predictions$valid_prediction), "/", nrow(coop),
                "; missing selected covariates are reported as NA.")
    run_sensitivity_checks(saved_predictions, draw, y, draw_dir)
    stats <- scorer$stats()
    timing <- data.frame(draw = draw, mesh_vertices = prepared$mesh$n,
                         observation_locations = length(y), prediction_locations = nrow(coop),
                         unique_candidate_fits = length(records), new_candidate_fits_this_session = stats$new_fits,
                         cache_hits_this_session = stats$cache_hits,
                         total_candidate_seconds = sum(vapply(records, function(r) r$elapsed_seconds, numeric(1))),
                         final_prediction_seconds = saved_predictions$elapsed_seconds,
                         session_draw_seconds = proc.time()[["elapsed"]] - draw_start,
                         selected_covariates = paste(covariates, collapse = " + "),
                         waic_selection = selected$waic, waic_full = saved_predictions$waic_full,
                         waic_layout_difference = waic_layout_difference)
    atomic_save_rds(timing, done_file)
    log_message("DRAW ", draw, " COMPLETE | ", length(records), " unique candidate fits | ",
                round(timing$total_candidate_seconds + timing$final_prediction_seconds, 1),
                " total fitting seconds (including any resumed fits).")
    timing
  }

  # Rebuild the timing summary from all completed checkpoints, including earlier
  # batches in this same run folder. Each draw's completion marker is written last.
  refresh_timings <- function() {
    files <- list.files(output_dirs[["checkpoints"]], pattern = "^completed\\.rds$",
                        recursive = TRUE, full.names = TRUE)
    if (length(files)) {
      timing <- do.call(rbind, lapply(files, readRDS))
      timing <- timing[order(timing$draw), , drop = FALSE]
      atomic_write_csv(timing, timing_path)
    }
  }
  for (iteration in iteration_columns) {
    tryCatch(run_one_draw(iteration), error = function(e) {
      log_message("STOPPED during draw ", iteration - 7L, ": ", conditionMessage(e))
      log_message("Completed candidate scores are saved. Resume this folder with unchanged settings.")
      stop(e)
    })
    draw <- iteration - 7L
    if (nzchar(cfg$reference_output_dir)) {
      matched <- tryCatch(compare_run_outputs(output_root, path.expand(cfg$reference_output_dir), draw,
                          file.path(output_dirs[["diagnostics"]], paste0("Reference_comparison_model_", draw, ".csv"))),
                          error = function(e) { log_message("Reference comparison: ", conditionMessage(e)); FALSE })
      if (!isTRUE(matched)) log_message("Reference comparison unavailable; fit outputs are saved.")
    }
    refresh_timings()
  }
  log_message("REQUESTED DRAWS COMPLETE. Timing summary: ", timing_path)
  log_message("Inspect Covariate_support.csv, Uncertainty_summary and Sensitivity_comparison in Diagnostics_NGA. ",
              "The sensitivity fits hold selected covariates fixed. Material changes warrant a full-search comparison.")
  log_message("Prediction_targets_model files compare new-location and conditional-surface uncertainty. ",
              "Use the surface output only when reporting that different estimand; inspect finer-mesh sensitivity before scaling up.")
  invisible(output_root)
}

run_aat_benchmark(settings)
