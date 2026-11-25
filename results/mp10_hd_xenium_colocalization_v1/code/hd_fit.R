source(file.path(Sys.getenv("OUT"),"code/io.R"))
library(spacexr)
s <- commandArgs(TRUE)[1]; d <- file.path(out,"hd",s)
helper <- any(c("--reuse-sigma-tail","--reuse-sigma-center","--reuse-sigma-quarter") %in% commandArgs(TRUE))
if(!helper) lock <- filelock::lock(file.path(d,"fit.lock"))
if(file.exists(file.path(d,"fit_complete.txt")))quit(save="no")
ref <- readRDS(file.path(root,"highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1/models/fine_reference.rds"))
labs <- factor(gsub("/","_",ref$metadata$label,fixed=TRUE));names(labs)<-colnames(ref$counts)
refobj <- Reference(ref$counts,labs,nUMI=colSums(ref$counts),n_max_cells=5000,min_UMI=100)
x <- sparse_input(file.path(d,"counts")); u <- fread(file.path(d,"units.tsv.gz")); stopifnot(identical(colnames(x),u$region_id))
x <- x[intersect(rownames(x),rownames(ref$counts)),]; totals <- colSums(x); valid <- totals>=100
sp <- SpatialRNA(data.frame(x=u$x_um[valid],y=u$y_um[valid],row.names=u$region_id[valid]),x[,valid],totals[valid])
rm(x,ref);gc()
fit <- create.RCTD(sp,refobj,max_cores=1,UMI_min=100,UMI_min_sigma=300,counts_MIN=10)
rm(sp,refobj);gc()
fit@config$RCTDmode <- "full"
fit@config$N_fit <- 200
fit <- spacexr::fitBulk(fit)
cfgfile <- file.path(d,"RCTD_config.json")
if(helper){
  cfg<-jsonlite::read_json(cfgfile,simplifyVector=TRUE)
  stopifnot(cfg$version==as.character(packageVersion("spacexr")),cfg$RCTDmode=="full",cfg$N_fit==200)
  key<-as.character(round(cfg$sigma*100))
  for(i in 1:5){
    collection<-readRDS(system.file("extdata",paste0("Qmat/Q_mat_",i,".rds"),package="spacexr"))
    if(key %in% names(collection)){fit@internal_vars$Q_mat<-collection[[key]];break}
  }
  stopifnot(!is.null(fit@internal_vars$Q_mat))
  fit@internal_vars$X_vals<-readRDS(system.file("extdata","Qmat/X_vals.rds",package="spacexr"))
  fit@internal_vars$sigma<-cfg$sigma
  rm(collection);gc()
}else{
  fit <- spacexr::choose_sigma_c(fit)
  jsonlite::write_json(c(fit@config,list(sigma=fit@internal_vars$sigma,version=as.character(packageVersion("spacexr")))),cfgfile,auto_unbox=TRUE,digits=15)
}
# One slide calibration, native full-mode pixel fit in bounded chunks.
puck <- fit@spatialRNA
ids <- colnames(puck@counts)
batches <- split(seq_along(ids),ceiling(seq_along(ids)/256))
dir.create(file.path(d,"weights"),showWarnings=FALSE)
fit_batch <- function(i,verify=FALSE) {
  dest <- file.path(d,"weights",sprintf("%05d.tsv.gz",i)); if(!verify && file.exists(dest))return(TRUE)
  j <- batches[[i]]; sub <- fit
  sub@spatialRNA <- puck
  sub@spatialRNA@counts <- puck@counts[,j,drop=FALSE]; sub@spatialRNA@nUMI <- puck@nUMI[j]; sub@spatialRNA@coords <- puck@coords[j,,drop=FALSE]
  sub <- spacexr::fitPixels(sub,doublet_mode="full")
  w <- as.matrix(sub@results$weights); stopifnot(identical(rownames(w),ids[j]),all(is.finite(w)),all(w>=0))
  if(verify){
    old<-fread(dest,check.names=FALSE);stopifnot(identical(old$region_id,rownames(w)),identical(names(old)[-1],colnames(w)),max(abs(as.matrix(old[,-1])-w))<1e-8);return(TRUE)
  }
  tmp <- paste0(dest,".",Sys.getpid(),".tmp.gz");fwrite(data.frame(region_id=rownames(w),w,check.names=FALSE),tmp,sep="\t");file.rename(tmp,dest)
}
work <- seq_along(batches)
if(helper){
  stopifnot(fit_batch(1,verify=TRUE));cat("256 native weight vectors verified; starting reverse batch order\n");flush.console();work<-rev(work)
  if("--reuse-sigma-center" %in% commandArgs(TRUE))work<-work[order(abs(work-length(work)/2))]
  if("--reuse-sigma-quarter" %in% commandArgs(TRUE))work<-work[order(abs(work-.30*length(work)))]
}
ok <- parallel::mclapply(work,fit_batch,mc.cores=as.integer(Sys.getenv("HD_WORKERS","6")),mc.preschedule=TRUE)
stopifnot(all(vapply(ok,isTRUE,logical(1))))
if(!helper)writeLines(as.character(length(ids)),file.path(d,"fit_complete.txt"))

