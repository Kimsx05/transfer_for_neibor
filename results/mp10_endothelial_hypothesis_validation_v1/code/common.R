options(mc.cores=1,Ncpus=1)
.libPaths(unique(c("/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5","/refdir/Rlib",.libPaths())))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(edgeR);library(limma)})
setDTthreads(1)
RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)
stopifnot(getDTthreads()==1,RhpcBLASctl::blas_get_num_procs()==1)
O="results/mp10_endothelial_hypothesis_validation_v1"
C="cache/mp10_public_resources_v1"
V="results/mp10_hypothesis_followup_v2"
B="results/mp10_immune_endothelial_mechanisms_v1"
R="resources/mp10_regulatory_priors"
tab=function(root,n)fread(file.path(root,"tables",n),sep="\t")
put=function(x,n,dir="tables")fwrite(x,file.path(O,dir,n),sep="\t",na="NA",nThread=1)
gettab=function(n)tab(O,n)
M20=scan("neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt",what="",quiet=TRUE)
stopifnot(length(M20)==20,"ACAT2"%in%M20)

