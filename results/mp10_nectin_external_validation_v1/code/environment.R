source('results/mp10_nectin_external_validation_v1/code/common.R')
vers<-c('GSVA','limma','edgeR','org.Hs.eg.db','AnnotationDbi','BiocParallel','sandwich','readxl','data.table','ggplot2','ragg','matrixStats');put(data.table(package=vers,version=vapply(vers,function(x)as.character(packageVersion(x)),'')),file.path(D,'config/R_versions.tsv'));print(sessionInfo());print(Sys.getenv(c('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','OMP_THREAD_LIMIT')));system('python --version');system('uname -sr')

