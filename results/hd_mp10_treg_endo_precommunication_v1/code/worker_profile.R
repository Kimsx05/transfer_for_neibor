Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_THREAD_LIMIT=1,RCPP_PARALLEL_NUM_THREADS=1)
if(requireNamespace('RhpcBLASctl',quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
options(mc.cores=1)
.libPaths(c(file.path(getwd(),'highres_VisiumHD/analysis/hd_bin8_comparison_v1/Rlib'),file.path(getwd(),'highres_VisiumHD/analysis/hd_annotation_diagnosis_singler_v1/Rlib'),'/refdir/Rlib','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
if(requireNamespace('doSNOW',quietly=TRUE))loadNamespace('doSNOW')
