# Site Rprofile sets 24 threads on this host; explicitly override it before any analysis.
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_THREAD_LIMIT=1,NUMEXPR_NUM_THREADS=1,NUMBA_NUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
if(requireNamespace('RhpcBLASctl',quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
.libPaths(unique(c(file.path(getwd(),"mp10_immune_endothelial_mechanisms_v1/Rlib"),'/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4','/refdir/Rlib',.libPaths())))
if(requireNamespace('RcppParallel',quietly=TRUE))RcppParallel::setThreadOptions(numThreads=1)
options(Ncpus=1,mc.cores=1)
Sys.setenv(VECLIB_MAXIMUM_THREADS=1,BLIS_NUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(Matrix);library(qs);library(jsonlite);library(digest)})
setDTthreads(1);options(warn=1,mc.cores=1);Sys.setenv(RCPP_PARALLEL_NUM_THREADS=1)
P<-normalizePath('.');B<-file.path(P,'results/mp10_immune_endothelial_mechanisms_v1');S<-file.path(P,'mp10_immune_endothelial_mechanisms_v1/scripts');REF<-file.path(P,'resources/mp10_regulatory_priors');BASE<-file.path(P,'mp10_sc_c2l_final_annotation_v1');HDEF<-file.path(P,'mp10_cellchat_high_definition_sensitivity_v1')
tb<-function(n)file.path(B,'tables',n);put<-function(x,n)fwrite(x,tb(n),sep='\t',na='NA',compress=if(grepl('gz$',n))'gzip'else'none');rd<-function(n)fread(tb(n),sep='\t');objpath<-file.path(P,'subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs')
