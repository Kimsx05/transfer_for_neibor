Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_THREAD_LIMIT=1,RCPP_PARALLEL_NUM_THREADS=1)
P<-normalizePath('.');root<-file.path(P,'highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1');A<-dirname(root)
.libPaths(c(file.path(A,'hd_bin8_comparison_v1/Rlib'),file.path(A,'hd_annotation_diagnosis_singler_v1/Rlib'),'/refdir/Rlib','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
suppressPackageStartupMessages({library(Matrix);library(data.table);library(jsonlite);library(Seurat);library(SingleR);library(BiocParallel)})
if(requireNamespace('RhpcBLASctl',quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
if(requireNamespace('RcppParallel',quietly=TRUE))RcppParallel::setThreadOptions(numThreads=1)
setDTthreads(1);options(mc.cores=4);set.seed(123);future::plan('sequential')
put<-function(x,n)fwrite(x,file.path(root,'tables',n),sep='\t',na='NA');rd<-function(n)fread(file.path(root,'tables',n),data.table=FALSE)
lognorm<-function(x){x<-as(x,'dgCMatrix');x@x<-log1p(x@x*rep(10000/pmax(colSums(x),1),diff(x@p)));x}
gate<-function(){while(as.integer(format(Sys.time(),'%H',tz='Asia/Shanghai'))>=7 && as.integer(format(Sys.time(),'%H',tz='Asia/Shanghai'))<23)Sys.sleep(30)}
