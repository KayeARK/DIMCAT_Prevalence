# ===========================================================================
# AAT NIGERIA — SPATIALLY BLOCKED, TRAINING-ONLY CROSS-VALIDATION
# Based on: AAT_spatial_quick_joint_posterior(1).R (29 September 2026).
#
# WHAT THIS VALIDATES:
#  - Out-of-block prediction of the *diagnostic-adjusted binomial counts*.
#  - Model structure, training-only WAIC selection, fixed effects, SPDE and IID.
#  - All tests at the same exact spatial coordinate stay in the same fold.
#  - Predictions are for a NEW site: its IID is integrated/sampled, NEVER
#    borrowed from an observed training site.
#  - Held-out raw observations are NOT an independent diagnostic validation:
#    the upstream diagnostic draws were constructed using all original data.
#    Truly end-to-end CV would re-estimate test properties without held-out data.
#
# The script does NOT change the original script, saved INLA runs or results.
# Requires the original model script and its spatial_preparation.rds checkpoint.
# Run from the SAME project working directory as the original model script.
# ===========================================================================

cv <- list(
  # EDIT THESE THREE PATHS FIRST. source_run_dir is a COMPATIBLE completed
  # Nigeria run containing Checkpoints_NGA/spatial_preparation.rds.
  original_script = "Code/Prevalence/Bovine BCT and PCR/AAT_spatial_quick_joint_posterior_NGA.R",
  source_run_dir = "Code/Prevalence/Bovine BCT and PCR/REPLACE_WITH_YOUR_RUN_FOLDER",
  output_dir = "Code/Prevalence/Bovine BCT and PCR/CV_NGA_spatial_blocks",

  # Pilot with draw 1, then extend to e.g. 1:3 without changing the folds.
  # Draw ID 1 corresponds to input CSV column 8, ID 2 to column 9, etc.
  draw_ids = 1L,
  folds = 5L,
  block_width_km = 100,       # Grid blocks in Nigeria's METRE-based AEQD CRS
  exclusion_buffer_km = 25, # Remove training sites this close to test sites
  min_training_sites = 30L,
  min_test_sites = 3L,
  posterior_samples = 2048L,  # >= 1024 and even, as in original implementation
  predictive_simulations = 4096L,
  iid_quadrature_nodes = 31L,
  seed = 20261001L,

  # Model-comparison options; selection ALWAYS uses training data only.
  compare_spatial_null = TRUE,  # Intercept + same SPDE and location IID
  compare_no_iid = FALSE,      # Optional: selected terms + SPDE, without IID
  log_to_console = TRUE
)

# ------------------------ Preparatory safety checks -------------------------
required <- c("INLA", "sf", "Matrix", "matrixStats")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(INLA); library(Matrix)})

for (nm in c("original_script", "source_run_dir", "output_dir"))
  if (!is.character(cv[[nm]]) || length(cv[[nm]]) != 1L || !nzchar(cv[[nm]]))
    stop("Set cv$", nm)
if (!file.exists(cv$original_script)) stop("Cannot find original_script: ", cv$original_script)
if (!dir.exists(cv$source_run_dir)) stop("Cannot find source_run_dir: ", cv$source_run_dir)
if (cv$folds < 2L || cv$block_width_km <= 0 || cv$exclusion_buffer_km < 0 ||
    cv$posterior_samples < 1024L || (cv$posterior_samples %% 2L) != 0L ||
    cv$predictive_simulations < 500L || cv$iid_quadrature_nodes < 3L)
  stop("Invalid CV configuration.")
if (anyNA(cv$draw_ids) || any(cv$draw_ids < 1L) || anyDuplicated(cv$draw_ids))
  stop("draw_ids must be distinct positive integers.")

# Read FUNCTION DEFINITIONS from the original file WITHOUT running its final
# run_aat_quick(settings) call. Its internal model and sampling functions are
# then reused unchanged. The original file is expected to end with this call.
original_expr <- parse(file = cv$original_script, keep.source = FALSE)
last_text <- paste(deparse(tail(original_expr, 1L)[[1L]]), collapse = "")
if (!grepl("^run_aat_quick\\(settings\\)$", gsub("[[:space:]]", "", last_text)))
  stop("Original file no longer ends in run_aat_quick(settings). Review loader.")
model <- new.env(parent = globalenv())
for (expr in head(original_expr, -1L)) eval(expr, envir = model)
needed_fns <- c("coordinate_key", "aggregate_draw_counts", "audit_covariates",
  "audit_quadratic_basis", "select_per_draw", "select_best_visited", "fit_observations",
  "sample_joint_posterior", "build_model_stacks", "transform_model_covariates",
  "effect_terms", "normal_quadrature")
if (!all(vapply(needed_fns, function(nm) is.function(model[[nm]]), logical(1))))
  stop("Original model helper API changed; compare original_script to supplied version.")

prep_path <- file.path(cv$source_run_dir, "Checkpoints_NGA", "spatial_preparation.rds")
manifest_path <- file.path(cv$source_run_dir, "run_manifest.rds")
source_settings_path <- file.path(cv$source_run_dir, "settings_used.rds")
if (!all(file.exists(c(prep_path, manifest_path, source_settings_path))))
  stop("source_run_dir must contain spatial preparation, run manifest and saved settings.")
base <- readRDS(prep_path)
manifest <- readRDS(manifest_path)
original_settings <- readRDS(source_settings_path)

needed_prep <- c("coordinate_mapping", "kept_location_rows", "observation_lonlat",
  "observation_m", "prediction_m", "raw_observation_covariates", "numtrials")
if (!all(needed_prep %in% names(base)))
  stop("Cached spatial preparation lacks fields: ",
       paste(setdiff(needed_prep, names(base)), collapse = ", "))
if (!isTRUE(manifest$use_location_iid))
  message("Note: original run did not use location IID; CV will match that setting.")
if (!identical(manifest$covariate_degree, 2L))
  stop("Expected the original quadratic model; check the original run manifest.")

input_path <- original_settings$data_path
if (!file.exists(input_path)) stop("Cannot find adjusted-count input: ", input_path)
input_md5 <- unname(tools::md5sum(input_path))
if (!identical(input_md5, manifest$input_md5))
  stop("Adjusted-count input differs from the cached preparation's input. Aborting.")
raw <- read.csv(input_path, check.names = FALSE)
raw$longitude <- as.numeric(raw$longitude)
raw$latitude <- as.numeric(raw$latitude)
raw <- raw[is.finite(raw$longitude) & is.finite(raw$latitude), , drop = FALSE]
if (any(cv$draw_ids + 7L > ncol(raw))) stop("A requested draw ID is not in input data.")
n_sites <- nrow(base$observation_m)
if (nrow(base$observation_lonlat) != n_sites || nrow(base$raw_observation_covariates) != n_sites ||
    length(base$numtrials) != n_sites) stop("Misaligned source preparation.")
if (!identical(as.numeric(base$numtrials),
               as.numeric(base$coordinate_mapping$sample_size[base$kept_location_rows])))
  stop("Mismatch between prepared sample sizes and original coordinate mapping.")

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
  f <- sample(candidates, 1L)
  block_fold[bn] <- f
  fold_load[f] <- fold_load[f] + length(blocks[[bn]])
}
fold_id <- unname(block_fold[block_id])
if (anyNA(fold_id) || any(tabulate(fold_id, cv$folds) < cv$min_test_sites))
  stop("Empty/undersized fold. Change folds or block_width_km.")

# Cache identity: adding further diagnostic draws is allowed, but changing
# blocks, source model, or scoring specification requires a NEW output folder.
cv_identity <- list(version = "nga_spatial_block_cv_20261001_v1",
  original_script_md5 = unname(tools::md5sum(cv$original_script)),
  preparation_md5 = unname(tools::md5sum(prep_path)), input_md5 = input_md5,
  block_id = block_id, fold_id = fold_id, folds = cv$folds,
  buffer_km = cv$exclusion_buffer_km, block_width_km = cv$block_width_km,
  min_training_sites = cv$min_training_sites, min_test_sites = cv$min_test_sites,
  posterior_samples = cv$posterior_samples,
  predictive_simulations = cv$predictive_simulations,
  quadrature_nodes = cv$iid_quadrature_nodes, seed = cv$seed,
  compare_spatial_null = cv$compare_spatial_null,
  compare_no_iid = cv$compare_no_iid)
dir.create(cv$output_dir, showWarnings = FALSE, recursive = TRUE)
identity_file <- file.path(cv$output_dir, "cv_identity.rds")
if (file.exists(identity_file)) {
  if (!identical(readRDS(identity_file), cv_identity))
    stop("This output folder belongs to a different CV specification. Use a NEW output_dir.")
} else saveRDS(cv_identity, identity_file)

block_assignment <- data.frame(site_id = seq_len(n_sites),
  longitude = base$observation_lonlat[, 1L], latitude = base$observation_lonlat[, 2L],
  x_m = xy[, 1L], y_m = xy[, 2L], block_id = block_id, fold = fold_id)
write.csv(block_assignment, file.path(cv$output_dir, "spatial_fold_assignments.csv"), row.names = FALSE)

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
  for (j in which(valid)) {
    # Both fixed-effect and spatial effects are projected using the SAME
    # joint INLA posterior sample; new-site IID is never a trained-site IID.
    eta <- as.numeric(X[j,,drop=FALSE] %*% post$beta +
      fold_data$A_test[j,,drop=FALSE] %*% post$spatial)
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
  do.call(rbind, res)
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
fit_cfg$local_work_root <- tempfile("aat_cv_inla_", tmpdir=tempdir())
dir.create(fit_cfg$local_work_root)
fit_cfg$diagnostics_dir <- cv$output_dir
fit_cfg$posterior_samples <- as.integer(cv$posterior_samples)
# Reuse existing prior/controls/INLA mode. No candidate or fit imported from
# the full-data run, because that WOULD contaminate cross-validation.

fit_cv <- function(fold_data, y_all, draw, fold, work) {
  p <- fold_data$prepared
  yt <- y_all[fold_data$train]
  candidates <- p$candidate_covariates
  cache_path <- file.path(work, "training_candidate_waic.rds")
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
    file.path(work,"training_candidate_waic.csv"),row.names=FALSE)
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
  out[[length(out)+1L]] <- score_constant_baseline(fold_data,y_all,draw,fold)
  do.call(rbind,out)
}

# ------------------------- Fold / draw execution ----------------------------
cat("\nSpatial CV:",n_sites,"unique sites;",length(blocks),"occupied blocks;",
    cv$folds,"folds; block width",cv$block_width_km,"km; buffer",
    cv$exclusion_buffer_km,"km.\n")
fold_preparations <- lapply(seq_len(cv$folds),make_training_fold)
fold_manifest <- do.call(rbind,lapply(seq_along(fold_preparations), function(f) {
  z <- fold_preparations[[f]]
  data.frame(fold=f,train_sites=length(z$train),heldout_sites=length(z$test),
    buffered_training_sites_excluded=length(z$excluded),
    minimum_test_to_training_km=min(z$distance_to_training_km),
    median_test_to_training_km=median(z$distance_to_training_km))
}))
write.csv(fold_manifest,file.path(cv$output_dir,"fold_geometry_summary.csv"),row.names=FALSE)
print(fold_manifest)

for (draw in cv$draw_ids) {
  y_all <- model$aggregate_draw_counts(raw,as.integer(draw+7L),base$coordinate_mapping)
  y_all <- y_all[base$kept_location_rows]
  if (length(y_all)!=n_sites || any(!is.finite(y_all)) ||
      any(y_all<0 | y_all>base$numtrials)) stop("Invalid adjusted counts for draw ",draw)
  for (fold in seq_len(cv$folds)) {
    work <- file.path(cv$output_dir,sprintf("draw_%04d",draw),sprintf("fold_%02d",fold))
    dir.create(work,showWarnings=FALSE,recursive=TRUE)
    done <- file.path(work,"heldout_results.rds")
    if (file.exists(done)) {
      cat("Already completed draw ",draw," fold ",fold,"\n",sep="")
      next
    }
    cat("\nDRAW ",draw," / FOLD ",fold," START\n",sep="")
    v <- fit_cv(fold_preparations[[fold]],y_all,draw,fold,work)
    # Write the complete fold result only after all requested fits succeed.
    saveRDS(v,done)
    write.csv(v,file.path(work,"heldout_results.csv"),row.names=FALSE)
  }
}

# Aggregate available completed runs. The units of independent geographical
# validation are SPATIAL BLOCKS, not every draw-site pair. Do not use these
# rows as independent observations in inferential confidence intervals.
all_result_files <- list.files(cv$output_dir,"^heldout_results[.]rds$",
                               recursive=TRUE,full.names=TRUE)
all_results <- do.call(rbind,lapply(all_result_files,readRDS))
write.csv(all_results,file.path(cv$output_dir,"all_heldout_site_scores.csv"),row.names=FALSE)
keys <- interaction(all_results$draw,all_results$fold,all_results$model,drop=TRUE)
groups <- split(all_results,keys)
fold_scores <- do.call(rbind,lapply(groups,function(z) data.frame(
  draw=z$draw[1L],fold=z$fold[1L],model=z$model[1L],n_sites=nrow(z),
  total_log_predictive_score=sum(z$log_score),
  mean_log_predictive_score=mean(z$log_score),
  rmse_adjusted_fraction=sqrt(mean((z$posterior_mean_p-z$adjusted_fraction)^2)),
  mean_absolute_error=mean(abs(z$posterior_mean_p-z$adjusted_fraction)),
  predictive_80pct_coverage=mean(z$y>=z$count_lo80 & z$y<=z$count_hi80),
  predictive_95pct_coverage=mean(z$y>=z$count_lo95 & z$y<=z$count_hi95),
  mean_95pct_predictive_fraction_width=mean((z$count_hi95-z$count_lo95)/z$n),
  mean_randomized_pit=mean(z$randomized_pit)
)))
row.names(fold_scores) <- NULL
write.csv(fold_scores,file.path(cv$output_dir,"fold_draw_scores.csv"),row.names=FALSE)
model_scores <- do.call(rbind,lapply(split(fold_scores,fold_scores$model),function(z) data.frame(
  model=z$model[1L],completed_fold_draws=nrow(z),
  mean_fold_log_score_per_site=mean(z$mean_log_predictive_score),
  mean_fold_rmse=mean(z$rmse_adjusted_fraction),
  mean_fold_predictive_80pct_coverage=mean(z$predictive_80pct_coverage),
  mean_fold_predictive_95pct_coverage=mean(z$predictive_95pct_coverage),
  mean_fold_95pct_interval_fraction_width=mean(z$mean_95pct_predictive_fraction_width)
)))
row.names(model_scores) <- NULL
write.csv(model_scores,file.path(cv$output_dir,"model_summary_descriptive.csv"),row.names=FALSE)

# Simple diagnostics: PIT should be approximately uniform (allowing finite
# samples and dependence); calibration should be checked along with sharpness.
pdf(file.path(cv$output_dir,"cv_diagnostics.pdf"),width=10,height=5)
old_par <- par(mfrow=c(1,2))
selected <- all_results[all_results$model=="selected_spatial_iid",]
if (nrow(selected)) {
  hist(selected$randomized_pit,breaks=seq(0,1,by=.1),xlab="Randomized PIT",
       main="Held-out selected-model PIT")
  abline(h=nrow(selected)/10,col="red",lty=2)
  plot(selected$distance_to_training_km,
       abs(selected$posterior_mean_p-selected$adjusted_fraction),
       xlab="Nearest retained training site (km)",
       ylab="Absolute error (adjusted fraction)",
       main="Held-out errors versus distance",pch=16,cex=.6)
}
par(old_par);dev.off()

cat("\nCV complete for requested draws. Output:",normalizePath(cv$output_dir),"\n")
print(model_scores)
cat("Interpretation: validation is CONDITIONAL on precomputed diagnostic-adjusted counts.\n",
    "A true held-out raw-test likelihood requires an end-to-end diagnostic refit.\n",sep="")
