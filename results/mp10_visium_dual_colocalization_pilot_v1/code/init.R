.libPaths(unique(c(path.expand("~/R/x86_64-pc-linux-gnu-library/4.5"),path.expand("~/R/x86_64-pc-linux-gnu-library/4.4"),"/refdir/Rlib",.libPaths())))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
if(requireNamespace("RhpcBLASctl",quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
options(Ncpus=1,mc.cores=1)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(Seurat)})
setDTthreads(1)
options(future.globals.maxSize=8*1024^3)
future::plan("sequential")
seed <- 20261125L
read_matrix <- function(name) {
 p <- file.path(tmp,name)
 load_int <- function(k) readBin(file.path(p,k),integer(),n=file.size(file.path(p,k))/4,endian="little")
 genes <- readLines(file.path(p,"genes")); ids <- readLines(file.path(p,"ids"))
 # Preserve total UMI without allowing Seurat symbol rewriting to merge identifiers.
 bad <- grepl("[_|]",genes)
 genes[bad] <- paste0("PilotUnusedFeature",which(bad))
 stopifnot(!anyDuplicated(genes),!anyDuplicated(ids))
 new("dgCMatrix",i=load_int("i"),p=load_int("p"),x=as.numeric(load_int("x")),Dim=as.integer(c(length(genes),length(ids))),Dimnames=list(genes,ids))
}
