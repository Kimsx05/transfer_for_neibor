source("code/direction.R");library(data.table);library(Matrix);setDTthreads(1);RhpcBLASctl::blas_set_num_threads(1)
set.seed(123);x=matrix(rpois(300*400,lambda=.8),300);rownames(x)=paste0("g",1:300);x=as(x,"dgCMatrix")
ks=rownames(x)[1:5];d=sapply(c(FALSE,TRUE),function(depth){
 a=scTenifoldKnk::knockoutDirection(x,as.list(ks),regressLibSize=depth)
 b=gram_direction(x,rownames(x),ks,depth=depth)
 max(abs(a-b))
});stopifnot(max(d)<1e-10);fwrite(data.table(depth=c(FALSE,TRUE),max_absolute_difference=d),"qc/direction_equivalence.tsv",sep="\t")

