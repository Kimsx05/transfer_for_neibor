# Launch with Rscript --vanilla and single-thread environment (run.sh).
.libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4','/refdir/Rlib',.libPaths())))
options(mc.cores=1,Ncpus=1)
suppressPackageStartupMessages({library(data.table);library(Matrix);library(qs);library(edgeR);library(limma)})
setDTthreads(1)
if(requireNamespace("RhpcBLASctl",quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
out<-Sys.getenv('B_OUT','results/mp10_treg_metabolic_mechanisms_v1')
old<-'results/mp10_immune_endothelial_mechanisms_v1'
v2<-'results/mp10_hypothesis_followup_v2'
readt<-function(root,name)fread(file.path(root,'tables',name),sep='\t')
saveT<-function(x,name)fwrite(x,file.path(out,'tables',name),sep='\t',na='NA',nThread=1)
top<-scan('neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt',what='',quiet=TRUE)
stopifnot(length(top)==20,'ACAT2'%in%top,!('SOAT2'%in%top))
writeLines(capture.output(sessionInfo()),file.path(out,'logs','R-session.txt'))
