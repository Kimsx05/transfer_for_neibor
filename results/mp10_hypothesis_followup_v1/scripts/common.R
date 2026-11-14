Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_THREAD_LIMIT=1,RCPP_PARALLEL_NUM_THREADS=1,NUMEXPR_NUM_THREADS=1,NUMBA_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,BLIS_NUM_THREADS=1)
.libPaths(unique(c(file.path(getwd(),'mp10_immune_endothelial_mechanisms_v1/Rlib'),'/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4','/refdir/Rlib',.libPaths())))
if(requireNamespace('RhpcBLASctl',quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
if(requireNamespace('RcppParallel',quietly=TRUE))RcppParallel::setThreadOptions(numThreads=1)
options(Ncpus=1,mc.cores=1,warn=1)
suppressPackageStartupMessages({library(data.table);library(Matrix);library(qs);library(jsonlite);library(digest)})
setDTthreads(1)
P<-normalizePath('.');B<-file.path(P,'results/mp10_hypothesis_followup_v1');S<-file.path(P,'mp10_hypothesis_followup_v1/scripts');OLD<-file.path(P,'results/mp10_immune_endothelial_mechanisms_v1');REF<-file.path(P,'resources/mp10_regulatory_priors');objpath<-file.path(P,'subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs')
rd<-function(n)fread(file.path(B,'tables',n),sep='\t');old<-function(n)fread(file.path(OLD,'tables',n),sep='\t');put<-function(x,n)fwrite(x,file.path(B,'tables',n),sep='\t',na='NA',compress=if(grepl('gz$',n))'gzip'else'none')
M20<-scan(file.path(P,'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'),what='character',quiet=TRUE)
thread_audit<-function(stage){d<-data.table(stage=stage,pid=Sys.getpid(),time=as.character(Sys.time()),BLAS=if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::blas_get_num_procs()else NA_integer_,OMP=if(requireNamespace('RhpcBLASctl',quietly=TRUE))RhpcBLASctl::omp_get_max_threads()else NA_integer_,data_table=getDTthreads());fwrite(d,file.path(B,'logs/thread_audit.tsv'),sep='\t',append=file.exists(file.path(B,'logs/thread_audit.tsv')));stopifnot(d$BLAS==1,d$OMP==1,d$data_table==1)}
thread_audit(basename(commandArgs()[grepl('^--file=',commandArgs())]))
