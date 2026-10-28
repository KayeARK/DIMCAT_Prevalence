# ============================================================================
# AAT: selected-covariate + spatial + IID validation (two publication figures)
# Sources: non-aggregated, 100-km means, and 200-km means (six output folders).
# NO fitting, no INLA dependency, no changes to source CSVs/checkpoints.
# Requires ggplot2; install once with install.packages("ggplot2") if necessary.
#
# These are POSTERIOR PREDICTIVE intervals for corrected survey observations.
# They are NOT credible intervals for the mapped IID-integrated expected mean.
# Coverage is calculated from the saved interval endpoints, not from PIT values,
# not from averaging fold percentages, and not from averaging site quantiles.
# Repeated diagnostic draws are repeated evaluations, not independent new sites.
# ============================================================================

plot_settings <- list(
  # Blank finds the six folders beside this script or in the project analysis dir.
  base_dir = "",
  output_dir = "CV_selected_model_diagnostics_ordered",
  # Optional: replot an existing Interval_coverage_audit.csv instead of raw folders.
  # Leave blank for the normal six-folder workflow. No fitting is performed.
  audit_csv = "",
  expected_folds = 5L,
  # Optional explicit paths (absolute, or relative to base_dir).
  input_dirs = c(
    NGA_site = "CV_NGA_spatial_blocks",
    NGA_100 = "CV_NGA_100km_block_mean_validation",
    NGA_200 = "CV_NGA_200km_block_mean_validation",
    ETH_site = "CV_ETH_spatial_blocks",
    ETH_100 = "CV_ETH_100km_block_mean_validation",
    ETH_200 = "CV_ETH_200km_block_mean_validation"
  ),
  require_all_six = TRUE,  # FALSE permits a clearly documented partial report.
  draw_ids = NULL,        # NULL uses draws common to all available analyses.
  # 80% and 95% are the levels exported by all six generation scripts.
  interval_levels = c(80L, 95L),
  representative_draw = NULL, # NULL: first included draw for detailed point plots.
  dpi = 600L,
  width_inches = 10.5,
  height_inches = 7.1,
  base_font_size = 10,
  save_png = TRUE,
  save_pdf = TRUE
)

run_combined_cv_diagnostics <- function(cfg) {
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop('Install ggplot2 first: install.packages("ggplot2")')
  # Keep a private null device available for text/grob measurement in headless R.
  # This prevents a graphics package from reopening the previous export device.
  grDevices::pdf(file=NULL)
  measurement_device <- grDevices::dev.cur()
  on.exit({
    if(measurement_device %in% grDevices::dev.list()) grDevices::dev.off(measurement_device)
  },add=TRUE)
  # Imports are local: source() does not replace the user's data objects.
  script_dir <- function() {
    p <- Filter(function(x) is.character(x) && length(x)==1L && nzchar(x),
                lapply(sys.frames(), function(x) x$ofile))
    if (length(p)) return(dirname(normalizePath(tail(p,1L)[[1L]],mustWork=TRUE)))
    a <- grep("^--file=",commandArgs(FALSE),value=TRUE)
    if (length(a)) return(dirname(normalizePath(sub("^--file=","",a[1L]),mustWork=TRUE)))
    getwd()
  }
  absolute <- function(p) grepl("^(/|[A-Za-z]:[/\\\\])",path.expand(p))
  ancestors <- function(p) {
    out <- normalizePath(p,mustWork=TRUE)
    while (dirname(tail(out,1L)) != tail(out,1L)) out <- c(out,dirname(tail(out,1L)))
    out
  }
  launch <- unique(c(script_dir(),getwd()))
  audit_mode <- !is.null(cfg$audit_csv) && nzchar(cfg$audit_csv)
  if (nzchar(cfg$base_dir)) {
    root <- normalizePath(path.expand(cfg$base_dir),mustWork=TRUE)
  } else if (audit_mode) {
    root <- script_dir()
  } else {
    parents <- unique(unlist(lapply(launch,ancestors)))
    candidates <- unique(c(launch,file.path(parents,"Code","Prevalence","Bovine BCT and PCR")))
    counts <- vapply(candidates,function(p) sum(dir.exists(file.path(p,cfg$input_dirs))),integer(1))
    if (!any(counts>0)) stop("Cannot find CV output folders. Set plot_settings$base_dir to their parent directory.")
    root <- candidates[which.max(counts)]
  }
  if (!all(cfg$interval_levels %in% c(80L,95L)) || !length(cfg$interval_levels))
    stop("Existing exports support 80% and 95% only; other levels require predictive samples.")
  cfg$interval_levels <- sort(unique(as.integer(cfg$interval_levels)))
  out <- if (absolute(cfg$output_dir)) path.expand(cfg$output_dir) else file.path(root,cfg$output_dir)
  input_paths <- vapply(cfg$input_dirs,function(p) if (absolute(p)) path.expand(p) else file.path(root,p),character(1))
  if (normalizePath(out,mustWork=FALSE) %in% normalizePath(input_paths,mustWork=FALSE))
    stop("The diagnostics output folder must differ from all input folders.")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  read_csv <- function(p) read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
  bind <- function(z) { z <- Filter(function(x) !is.null(x) && nrow(x)>0L,z)
    if (!length(z)) return(data.frame()); rownames_out <- do.call(rbind,z); rownames(rownames_out)<-NULL; rownames_out }
  require_columns <- function(z,nms,label) {
    if (!all(nms %in% names(z))) stop(label,": missing columns ",paste(setdiff(nms,names(z)),collapse=", "))
  }
  finite <- function(z,nms,label) {
    if (any(!vapply(z[nms],is.numeric,logical(1))) || any(!is.finite(as.matrix(z[nms]))))
      stop(label,": non-numeric or non-finite values in ",paste(nms,collapse=", "))
  }
  group_apply <- function(z,by,fun) {
    if (!nrow(z)) return(data.frame())
    key <- do.call(interaction,c(z[by],list(drop=TRUE,lex.order=TRUE)))
    bind(lapply(split(z,key),function(g) {
      ans <- fun(g)
      cbind(g[rep(1L,nrow(ans)),by,drop=FALSE],ans)
    }))
  }
  pct <- function(x) paste0(format(round(100*x),trim=TRUE),"%")
  model_names <- c(selected_spatial_iid="Selected covariates + spatial + IID")
  target_levels <- c("Survey sites","100 km means","200 km means")
  country_names <- c(NGA="Nigeria",ETH="Ethiopia")
  records <- list(); inventory <- list(); geometry <- list()
  if (!audit_mode) for (id in names(input_paths)) {
    path <- input_paths[[id]]
    country <- strsplit(id,"_",fixed=TRUE)[[1L]][1L]
    kind <- sub("^[^_]+_","",id)
    prefix <- if (kind=="site") "" else paste0("block",kind,"_")
    filename <- paste0(prefix,if (kind=="site") "all_heldout_site_scores.csv" else "all_heldout_block_scores.csv")
    file <- file.path(path,filename)
    # Fall back to complete per-fold CSVs if the combined export is absent.
    # Do not read summary means or synthesize missing interval endpoints.
    chunks <- character()
    if (!file.exists(file) && dir.exists(path)) {
      suffix <- if (kind=="site") "heldout_results.csv" else "heldout_block_results.csv"
      chunks <- list.files(path,pattern=paste0("^",prefix,suffix,"$"),recursive=TRUE,full.names=TRUE)
      if (kind=="site" && !length(chunks)) chunks <- list.files(path,
        pattern="^heldout_site_results[.]csv$",recursive=TRUE,full.names=TRUE)
    }
    available <- file.exists(file) || length(chunks)>0L
    inventory[[id]] <- data.frame(analysis=id,country=country,target=kind,
       directory=path,available=available,source=if (file.exists(file)) file else paste(chunks,collapse=" | "))
    if (!available) next
    z <- if (file.exists(file)) read_csv(file) else bind(lapply(chunks,read_csv))
    if (!nrow(z)) stop("Empty result file for ",id)
    require_columns(z,c("draw","fold","model","block_id"),id)
    z <- z[!is.na(z$model) & z$model=="selected_spatial_iid",,drop=FALSE]
    if (!nrow(z)) stop("No selected + spatial + IID results in ",id)
    if (kind!="site") {
      require_columns(z,"averaging",id)
      z <- z[!is.na(z$averaging) & z$averaging=="equal_sites",,drop=FALSE]
      if (!nrow(z)) stop("No equal-site aggregated results in ",id)
    }
    if (anyNA(z[c("draw","fold","model","block_id")])) stop("Missing identifiers in ",id)
    if (any(!z$model %in% names(model_names))) stop("Unknown model label in ",id,": ",paste(setdiff(z$model,names(model_names)),collapse=", "))
    if (any(!is.finite(z$draw) | z$draw<1 | z$draw!=floor(z$draw)) ||
        any(!is.finite(z$fold) | z$fold<1 | z$fold!=floor(z$fold))) stop("Invalid draw/fold IDs in ",id)
    if (kind=="site") {
      need <- c("site_id","y","n","adjusted_fraction","posterior_mean_p","count_lo80","count_hi80","count_lo95","count_hi95")
      require_columns(z,need,id); finite(z,setdiff(need,"site_id"),id)
      if (any(z$n<=0 | z$y<0 | z$y>z$n | abs(z$y-round(z$y))>1e-8 | abs(z$n-round(z$n))>1e-8) ||
          any(abs(z$adjusted_fraction-z$y/z$n)>1e-8)) stop("Invalid counts/fractions in ",id)
      w <- data.frame(draw=z$draw,fold=z$fold,model=z$model,unit_id=as.character(z$site_id),
        block_id=as.character(z$block_id),averaging="equal_sites",observed=z$y/z$n,
        predicted=z$posterior_mean_p,lo80=z$count_lo80/z$n,hi80=z$count_hi80/z$n,
        lo95=z$count_lo95/z$n,hi95=z$count_hi95/z$n,n_sites=1,n_animals=z$n,
        pit=NA_real_,distance_km=NA_real_,
        log_score=NA_real_,crps=NA_real_)
    } else {
      need <- c("averaging","n_sites","n_cattle","observed_block_mean","posterior_mean_block_expected",
                "predictive_lo80","predictive_hi80","predictive_lo95","predictive_hi95")
      require_columns(z,need,id); finite(z,setdiff(need,"averaging"),id)
      if (any(!z$averaging %in% c("equal_sites","cattle_weighted")) || any(z$n_sites<2 | z$n_cattle<=0))
        stop("Invalid block averaging/support in ",id)
      w <- data.frame(draw=z$draw,fold=z$fold,model=z$model,unit_id=as.character(z$block_id),
        block_id=as.character(z$block_id),averaging=z$averaging,observed=z$observed_block_mean,
        predicted=z$posterior_mean_block_expected,lo80=z$predictive_lo80,hi80=z$predictive_hi80,
        lo95=z$predictive_lo95,hi95=z$predictive_hi95,n_sites=z$n_sites,n_animals=z$n_cattle,
        pit=NA_real_,distance_km=NA_real_,log_score=NA_real_,crps=NA_real_)
    }
    w$analysis <- id; w$country <- country_names[[country]]
    w$target <- if (kind=="site") target_levels[1L] else paste(kind,"km means")
    if (any(as.matrix(w[c("observed","predicted","lo80","hi80","lo95","hi95")]) < -1e-8) ||
        any(as.matrix(w[c("observed","predicted","lo80","hi80","lo95","hi95")]) > 1+1e-8) ||
        any(w$lo95>w$lo80 | w$lo80>w$hi80 | w$hi80>w$hi95)) stop("Invalid/non-nested intervals in ",id)
    key <- paste(w$draw,w$model,w$unit_id,w$averaging,sep="|")
    if (anyDuplicated(key)) stop("Duplicate held-out unit/model/draw in ",id,"; do not merge duplicate runs.")
    records[[id]] <- w
    gf <- file.path(path,paste0(prefix,"fold_geometry_summary.csv"))
    if (file.exists(gf)) {
      g <- read_csv(gf); g$analysis <- id
      # Keep fields common to site and aggregated geometry tables.
      keep <- intersect(c("analysis","fold","train_sites","heldout_sites","buffered_training_sites_excluded","minimum_test_to_training_km","median_test_to_training_km"),names(g))
      geometry[[id]] <- g[keep]
    }
  }
  if (audit_mode) {
    af <- if (absolute(cfg$audit_csv)) path.expand(cfg$audit_csv) else file.path(root,cfg$audit_csv)
    if (!file.exists(af)) stop("Cannot find audit_csv: ",af)
    if (normalizePath(af,mustWork=TRUE)==normalizePath(file.path(out,"Interval_coverage_audit.csv"),mustWork=FALSE))
      stop("Use a different output directory when reading a previous audit CSV.")
    a <- read_csv(af)
    need <- c("analysis","country","target","model","averaging","nominal","draw","fold",
              "unit_id","observed","predicted","lower","upper")
    require_columns(a,need,af)
    a <- a[!is.na(a$model) & a$model=="selected_spatial_iid" &
             !is.na(a$averaging) & a$averaging=="equal_sites",,drop=FALSE]
    if (!nrow(a)) stop("No selected-model/equal-site evaluations in audit CSV.")
    if (anyNA(a[need])) stop("Missing selected-model audit values.")
    finite(a,c("nominal","draw","fold","observed","predicted","lower","upper"),af)
    if (any(!a$analysis %in% names(input_paths))) stop("Unrecognised analysis IDs in audit CSV.")
    if (anyDuplicated(a[c("analysis","draw","unit_id","nominal")])) stop("Duplicate audit evaluations.")
    for (id in names(input_paths)) {
      z <- a[a$analysis==id,,drop=FALSE]
      country <- strsplit(id,"_",fixed=TRUE)[[1L]][1L]
      kind <- sub("^[^_]+_","",id)
      target <- if(kind=="site") target_levels[1L] else paste(kind,"km means")
      inventory[[id]] <- data.frame(analysis=id,country=country,target=kind,
        directory=dirname(af),available=nrow(z)>0L,source=af)
      if (!nrow(z)) next
      if (any(z$country!=country_names[[country]] | z$target!=target)) stop("Audit country/target mismatch: ",id)
      hi <- z[abs(z$nominal-.95)<1e-10,,drop=FALSE]
      lo <- z[abs(z$nominal-.80)<1e-10,,drop=FALSE]
      ix <- match(paste(hi$draw,hi$unit_id),paste(lo$draw,lo$unit_id))
      if (!nrow(hi) || nrow(hi)!=nrow(lo) || anyNA(ix) || nrow(z)!=nrow(hi)+nrow(lo))
        stop("Audit CSV requires paired 80% and 95% intervals for each unit: ",id)
      lo <- lo[ix,,drop=FALSE]
      if (any(hi$fold!=lo$fold | abs(hi$observed-lo$observed)>1e-10 |
              abs(hi$predicted-lo$predicted)>1e-10)) stop("Inconsistent audit interval pairs: ",id)
      records[[id]] <- data.frame(draw=hi$draw,fold=hi$fold,model=hi$model,
        unit_id=as.character(hi$unit_id),block_id=NA_character_,averaging="equal_sites",
        observed=hi$observed,predicted=hi$predicted,lo80=lo$lower,hi80=lo$upper,
        lo95=hi$lower,hi95=hi$upper,n_sites=NA_real_,n_animals=NA_real_,
        pit=NA_real_,distance_km=NA_real_,log_score=NA_real_,crps=NA_real_,
        analysis=id,country=country_names[[country]],target=target)
    }
  }
  inventory <- bind(inventory)
  write.csv(inventory,file.path(out,"Input_inventory.csv"),row.names=FALSE)
  absent <- inventory$analysis[!inventory$available]
  if (length(absent) && isTRUE(cfg$require_all_six)) stop("Missing raw outputs for: ",paste(absent,collapse=", "),
    ". See Input_inventory.csv; run these analyses or set require_all_six=FALSE for a labelled partial report.")
  if (!length(records)) stop("No usable raw validation outputs found.")
  if (length(absent)) warning("PARTIAL report: missing ",paste(absent,collapse=", "))
  common_draws <- Reduce(intersect,lapply(records,function(x) unique(x$draw)))
  draw_ids <- if (is.null(cfg$draw_ids)) sort(common_draws) else sort(unique(cfg$draw_ids))
  if (!length(draw_ids) || any(!draw_ids %in% common_draws)) stop("Requested draws are not present in every available analysis. Choose common completed draw IDs.")
  inventory$included_draws <- paste(draw_ids,collapse=",")
  inventory$available_draws <- vapply(inventory$analysis,function(id)
    if (id %in% names(records)) paste(sort(unique(records[[id]]$draw)),collapse=",") else "",character(1))
  write.csv(inventory,file.path(out,"Input_inventory.csv"),row.names=FALSE)
  dat <- bind(records); dat <- dat[dat$draw %in% draw_ids,,drop=FALSE]
  dat$target <- factor(dat$target,levels=target_levels)
  dat$country <- factor(dat$country,levels=c("Nigeria","Ethiopia"))
  dat$model <- factor(dat$model,levels=names(model_names))
  dat$error <- dat$predicted-dat$observed
  # Shared checks also apply to the optional audit-CSV route.
  finite(dat,c("draw","fold","observed","predicted","lo80","hi80","lo95","hi95"),"Selected results")
  if (any(dat$draw<1 | dat$draw!=floor(dat$draw) | dat$fold<1 | dat$fold!=floor(dat$fold)))
    stop("Invalid draw/fold identifiers.")
  if (anyNA(dat$unit_id) || any(!nzchar(dat$unit_id)) ||
      anyDuplicated(dat[c("analysis","draw","unit_id")])) stop("Missing or duplicate held-out units.")
  if (any(as.matrix(dat[c("observed","predicted","lo80","hi80","lo95","hi95")]) < -1e-8) ||
      any(as.matrix(dat[c("observed","predicted","lo80","hi80","lo95","hi95")]) > 1+1e-8) ||
      any(dat$lo95>dat$lo80 | dat$lo80>dat$hi80 | dat$hi80>dat$hi95))
    stop("Invalid prevalence values or non-nested predictive intervals.")
  if (length(cfg$expected_folds)!=1L || !is.finite(cfg$expected_folds) ||
      cfg$expected_folds<2L || cfg$expected_folds!=floor(cfg$expected_folds))
    stop("expected_folds must match the generation scripts (default: 5).")
  for (a in unique(dat$analysis)) for (d in draw_ids) {
    got <- sort(unique(dat$fold[dat$analysis==a & dat$draw==d]))
    if (!identical(as.integer(got),seq_len(cfg$expected_folds)))
      stop("Incomplete/unexpected folds in ",a," draw ",d,"; do not report a partial validation run.")
  }
  # Require completed, matching folds when geometry files are available.
  for (a in names(geometry)) for (d in draw_ids) {
    expected <- sort(unique(geometry[[a]]$fold))
    got <- sort(unique(dat$fold[dat$analysis==a & dat$draw==d & dat$model=="selected_spatial_iid"]))
    if (!identical(as.numeric(got),as.numeric(expected))) stop("Incomplete folds in ",a," draw ",d)
  }
  primary <- dat[dat$averaging=="equal_sites",,drop=FALSE]
  sel <- primary[primary$model=="selected_spatial_iid",,drop=FALSE]
  rep_draw <- if (is.null(cfg$representative_draw)) draw_ids[1L] else cfg$representative_draw
  if (length(rep_draw)!=1 || !rep_draw %in% draw_ids) stop("Invalid representative_draw.")
  representative <- sel[sel$draw==rep_draw,,drop=FALSE]
  # Keep original endpoint-based binary indicators; no normal approximation.
  long <- bind(lapply(cfg$interval_levels,function(level) {
    z <- dat
    z$nominal <- level/100
    z$lower <- z[[paste0("lo",level)]]; z$upper <- z[[paste0("hi",level)]]
    # Counts were divided by the same n on both sides; inclusive endpoints.
    z$covered <- z$observed>=z$lower & z$observed<=z$upper
    z$below <- z$observed<z$lower; z$above <- z$observed>z$upper
    z$width <- z$upper-z$lower
    alpha <- 1-z$nominal
    z$interval_score <- z$width + 2/alpha*pmax(z$lower-z$observed,0) + 2/alpha*pmax(z$observed-z$upper,0)
    z
  }))
  by <- c("analysis","country","target","model","averaging","nominal")
  coverage_fn <- function(z) data.frame(n_evaluations=nrow(z),n_distinct_units=length(unique(z$unit_id)),
    n_draws=length(unique(z$draw)),covered=sum(z$covered),below=sum(z$below),above=sum(z$above),
    coverage=mean(z$covered),mean_width_pp=100*mean(z$width),mean_interval_score_pp=100*mean(z$interval_score))
  pooled <- group_apply(long,by,coverage_fn)
  per_draw <- group_apply(long,c(by,"draw"),coverage_fn)
  per_fold <- group_apply(long,c(by,"draw","fold"),coverage_fn)
  write.csv(pooled,file.path(out,"Coverage_pooled.csv"),row.names=FALSE)
  write.csv(per_draw,file.path(out,"Coverage_by_draw.csv"),row.names=FALSE)
  write.csv(per_fold,file.path(out,"Coverage_by_fold_draw.csv"),row.names=FALSE)
  metrics <- group_apply(dat,c("analysis","country","target","model","averaging"),function(z)
    data.frame(n_evaluations=nrow(z),n_distinct_units=length(unique(z$unit_id)),n_draws=length(unique(z$draw)),
      RMSE_pp=100*sqrt(mean(z$error^2)),MAE_pp=100*mean(abs(z$error)),bias_pp=100*mean(z$error),
      mean_95_width_pp=100*mean(z$hi95-z$lo95),
      median_95_width_pp=100*median(z$hi95-z$lo95)))
  write.csv(metrics,file.path(out,"Prediction_metrics.csv"),row.names=FALSE)
  if (length(geometry)) write.csv(bind(geometry),file.path(out,"Fold_geometry.csv"),row.names=FALSE)
  write.csv(long[,c(by,"draw","fold","unit_id","observed","predicted","lower","upper","covered","width","interval_score")],
            file.path(out,"Interval_coverage_audit.csv"),row.names=FALSE)

  theme <- ggplot2::theme_bw(base_size=cfg$base_font_size,base_family="sans") +
    ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_line(colour="#E8E8E8",linewidth=.25),
      strip.background=ggplot2::element_rect(fill="#F1F3F5",colour=NA),
      strip.text=ggplot2::element_text(face="bold",size=cfg$base_font_size),
      legend.position="bottom",legend.title=ggplot2::element_blank(),
      plot.title=ggplot2::element_text(face="bold",size=cfg$base_font_size+2),
      plot.subtitle=ggplot2::element_text(size=cfg$base_font_size-1,colour="#444444"),
      plot.caption=ggplot2::element_text(hjust=0,size=cfg$base_font_size-2,colour="#444444"),
      plot.margin=ggplot2::margin(9,12,9,9),panel.spacing=grid::unit(.9,"lines"))
  facet <- function() ggplot2::facet_grid(country~target,drop=FALSE)
  manifest_plots <- list()
  save_plot <- function(p,name,caption,height=cfg$height_inches) {
    p <- p+theme+ggplot2::labs(caption=caption)
    temporary_exports <- character()
    on.exit(unlink(temporary_exports),add=TRUE)
    publish <- function(tmp,destination) {
      if (!file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size <= 0)
        stop("Figure export is empty: ",destination)
      if (!file.rename(tmp,destination)) {
        if(!file.copy(tmp,destination,overwrite=TRUE)) stop("Cannot save figure: ",destination)
        unlink(tmp)
      }
    }
    if (isTRUE(cfg$save_pdf)) {
      tmp <- tempfile(paste0(name,"_"),tmpdir=out,fileext=".pdf")
      temporary_exports <- c(temporary_exports,tmp)
      grDevices::pdf(tmp,width=cfg$width_inches,height=height,useDingbats=FALSE)
      device_id <- grDevices::dev.cur()
      tryCatch(print(p),finally=grDevices::dev.off(device_id))
      publish(tmp,file.path(out,paste0(name,".pdf")))
    }
    if (isTRUE(cfg$save_png)) {
      tmp <- tempfile(paste0(name,"_"),tmpdir=out,fileext=".png")
      temporary_exports <- c(temporary_exports,tmp)
      grDevices::png(tmp,width=cfg$width_inches,height=height,units="in",res=cfg$dpi,
        bg="white",type=if(capabilities("cairo")) "cairo" else getOption("bitmapType"))
      device_id <- grDevices::dev.cur()
      tryCatch(print(p),finally=grDevices::dev.off(device_id))
      publish(tmp,file.path(out,paste0(name,".png")))
    }
    manifest_plots[[name]] <<- data.frame(figure=name,caption=caption)
    invisible(p)
  }
  draw_label <- if(length(draw_ids)<=8L) paste(draw_ids,collapse=", ") else
    paste0(length(draw_ids)," common draws (",min(draw_ids)," to ",max(draw_ids),"; see inventory)")
  common_note <- paste0("Draws: ",draw_label,". Intervals predict corrected survey observations; they do not validate mapped mean-prevalence CrIs.",
    if(length(absent)) paste0("\nPARTIAL REPORT: missing ",paste(absent,collapse=", "),".") else "")
  pp <- pooled[pooled$averaging=="equal_sites" & pooled$model=="selected_spatial_iid",,drop=FALSE]
  ff <- per_fold[per_fold$averaging=="equal_sites" & per_fold$model=="selected_spatial_iid",,drop=FALSE]
  pp$label <- sprintf("%.1f%%\n(%d/%d)",100*pp$coverage,pp$covered,pp$n_evaluations)
  p <- ggplot2::ggplot(pp,ggplot2::aes(x=nominal,y=coverage))+
    ggplot2::geom_abline(slope=1,intercept=0,linetype=2,colour="#666666",linewidth=.45)+
    ggplot2::geom_point(data=ff,shape=1,size=1.8,colour="#777777",position=ggplot2::position_jitter(width=.003,height=0,seed=19))+
    ggplot2::geom_point(size=2.7,colour="#0072B2")+
    ggplot2::geom_text(ggplot2::aes(y=.075,label=label),size=2.8,colour="#004B76")+
    ggplot2::scale_x_continuous(breaks=cfg$interval_levels/100,labels=pct,limits=c(.74,1.035))+
    ggplot2::scale_y_continuous(labels=pct,limits=c(0,1.04),breaks=seq(0,1,.2))+
    facet()+ggplot2::labs(title="Predictive interval coverage",subtitle="Selected covariates + spatial + IID; blue = overall coverage; open circles = fold/draw coverage",
      x="Nominal predictive interval probability",y="Held-out observations inside interval")
  save_plot(p,"01_Coverage_selected_model",paste0("Dashed line: nominal coverage. Counts show observations inside intervals / total observations.\n",common_note,"\nRepeated draws reuse the same sites/cells. Fold dots are descriptive, not confidence intervals."))

  representative_metrics <- group_apply(representative,c("analysis","country","target"),function(z)
    data.frame(draw=rep_draw,n=nrow(z),RMSE_pp=100*sqrt(mean(z$error^2)),
               MAE_pp=100*mean(abs(z$error)),mean_95_width_pp=100*mean(z$hi95-z$lo95)))
  # Order independently within each country/target; ties use the unit ID.
  # Ordering changes display only, never observations, intervals or metrics.
  panel_ids <- c("NGA_site","NGA_100","NGA_200","ETH_site","ETH_100","ETH_200")
  panel_ids <- panel_ids[panel_ids %in% representative$analysis]
  ranked <- bind(lapply(panel_ids,function(id) {
    z <- representative[representative$analysis==id,,drop=FALSE]
    z <- z[order(z$predicted,z$unit_id),,drop=FALSE]
    z$rank <- seq_len(nrow(z))
    z$outside95 <- z$observed<z$lo95 | z$observed>z$hi95
    z
  }))
  pm <- representative_metrics[match(panel_ids,representative_metrics$analysis),,drop=FALSE]
  panel_labels <- sprintf("%s | %s\nn = %d; RMSE = %.1f pp; MAE = %.1f pp",
    pm$country,pm$target,pm$n,pm$RMSE_pp,pm$MAE_pp)
  ranked$panel <- factor(ranked$analysis,levels=panel_ids,labels=panel_labels)
  write.csv(representative_metrics,file.path(out,"Predictive_plot_metrics.csv"),row.names=FALSE)
  write.csv(ranked[,c("analysis","country","target","draw","fold","unit_id","rank","observed",
    "predicted","lo80","hi80","lo95","hi95","outside95")],
    file.path(out,"Predictive_plot_data.csv"),row.names=FALSE)
  integer_breaks <- function(limits) {
    b <- pretty(limits,n=4)
    b <- b[b>=1 & b<=limits[2] & abs(b-round(b))<1e-8]
    sort(unique(c(1,round(b))))
  }
  p <- ggplot2::ggplot(ranked,ggplot2::aes(rank,predicted))+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lo95,ymax=hi95,fill="95% interval"))+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lo80,ymax=hi80,fill="80% interval"))+
    ggplot2::geom_line(colour="#00689D",linewidth=.7)+
    ggplot2::geom_point(data=ranked[!ranked$outside95,,drop=FALSE],
      ggplot2::aes(y=observed),colour="#222222",size=1.05,alpha=.6)+
    ggplot2::geom_point(data=ranked[ranked$outside95,,drop=FALSE],
      ggplot2::aes(y=observed),colour="#C4451C",shape=17,size=1.7,alpha=.9)+
    ggplot2::scale_fill_manual(values=c("95% interval"="#DCEBF5","80% interval"="#9EC9E2"),
      breaks=c("95% interval","80% interval"),drop=FALSE)+
    # facet_wrap gives EACH panel its own rank axis, including 11-cell ETH.
    ggplot2::facet_wrap(~panel,ncol=3,scales="free_x",drop=TRUE)+
    ggplot2::scale_x_continuous(breaks=integer_breaks,expand=ggplot2::expansion(mult=c(.015,.015)))+
    ggplot2::scale_y_continuous(labels=pct,breaks=seq(0,1,.25),limits=c(0,1),
      expand=ggplot2::expansion(mult=c(.025,.025)))+
    ggplot2::labs(title="Ordered predictive intervals and held-out observations",
      subtitle=paste("Selected covariates + spatial + IID; diagnostic-correction draw",rep_draw),
      x="Held-out sites / cells ordered by posterior mean prediction",y="Prevalence")
  save_plot(p,"02_Ordered_predictive_intervals",paste0(
    "Blue line: posterior mean. Black dots: corrected observations inside the 95% interval; orange triangles: outside.\n",
    "Bands are pointwise predictive intervals. Ordering is not geographical and does not demonstrate predictive discrimination.\n",
    "Cell means give sampled sites equal weight; pp: percentage points. ",common_note),height=8)

  write.csv(bind(manifest_plots),file.path(out,"Figure_captions.csv"),row.names=FALSE)
  cat("\nSelected-model predictive coverage (equal sites):\n")
  print(pp[,c("country","target","nominal","covered","n_evaluations","coverage","mean_width_pp")],row.names=FALSE)
  writeLines(c("AAT selected-covariate + spatial + IID validation: two figures",paste("Generated:",Sys.time()),
    paste("Draw IDs:",paste(draw_ids,collapse=", ")),paste("Detailed point-plot draw:",rep_draw),
    paste("Missing analyses:",if(length(absent)) paste(absent,collapse=", ") else "none"),
    "Coverage_pooled.csv uses the ratio of all covered evaluations to all evaluations, not the unweighted mean of fold percentages.",
    "The primary aggregation gives each eligible cell equal weight and each surveyed site equal weight within a cell.",
    "Only selected_spatial_iid and equal_sites evaluations are included; all comparator models and sample-size-weighted results are excluded.",
    "The script reads saved predictions and performs no fitting or new posterior/predictive simulation.",
    "All intervals scored against observations are posterior predictive intervals (new-site IID plus binomial sampling).",
    "expected_prevalence_lo95/hi95 and p_lo95/hi95 are deliberately NOT scored against observed counts/fractions.",
    "Coverage is conditional on the included diagnostic-correction draws. Repeated draws are not independent replicates.",
    "No binomial/Wilson confidence intervals are assigned to coverage: spatial and shared-fit dependence would invalidate a naive independent-trials interpretation.",
    "Fold/draw points describe heterogeneity; they are not inferential confidence bounds.",
    "The 100-km scripts retain 200-km withholding blocks. Original non-aggregated scripts typically use 100-km withholding blocks.",
    "Thus the three primary targets do not automatically share the same CV design; consult Fold_geometry.csv and the saved assignments.",
    "Grid-cell means validate sampled populations, not unbiased areal mean prevalence. Finer cells may exclude more sparse areas.",
    "80% and 95% are the only common saved interval levels; extra levels cannot be recovered from endpoint summaries.",
    "Coverage and errors do not by themselves establish accurate spatial discrimination or calibrate mapped expected-prevalence credible intervals.",
    "Coverage pools the included draws; the predictive plot uses representative_draw, with metrics computed only from its displayed points.",
    "Source intervals, observations and coverage flags are exported in Interval_coverage_audit.csv.",
    "Predictive_plot_data.csv records the within-panel rank and outside95 flag for each displayed unit.",
    "The ordered plot joins pointwise intervals for display; it is not a simultaneous band or a spatial/temporal trajectory.",
    "Vector PDFs are suitable for layout/editing; PNGs are 600 dpi by default.",
    "Input_inventory.csv records the exact source paths and available/included draws."),file.path(out,"README_diagnostics.txt"))
  capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
  cat("\nDiagnostics saved in: ",normalizePath(out,mustWork=TRUE),"\n",sep="")
  invisible(list(output_dir=out,coverage=pooled,metrics=metrics,inventory=inventory))
}

run_combined_cv_diagnostics(plot_settings)
