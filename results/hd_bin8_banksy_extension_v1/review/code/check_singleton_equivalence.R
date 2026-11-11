source(file.path(commandArgs(TRUE)[1],'code/common.R'));source(file.path(root,'code/fast_singletons.R'));set.seed(123)
for(i in 1:12){n<-300;z<-rsparsematrix(n,n,density=.03);z@x<-abs(z@x);z<-as((z+t(z))/2,'dgCMatrix');diag(z)<-1;dimnames(z)<-list(paste0('bin',1:n),paste0('bin',1:n));z<-as(z,'Graph');id<-sample(c(rep(0:5,each=45),6:35));names(id)<-colnames(z);if(i%%2==0){z[1:10,]<-0;z[,1:10]<-0};a<-Seurat:::GroupSingletons(id,z,verbose=FALSE);b<-fast_group_singletons(id,z,verbose=FALSE);stopifnot(identical(a,b));cat('MATCH',i,'\n')}
# No-singleton case retains original storage type/levels.
id<-rep(0:5,each=50);names(id)<-colnames(z);stopifnot(identical(Seurat:::GroupSingletons(id,z,verbose=FALSE),fast_group_singletons(id,z,verbose=FALSE)));cat('ALL_NATIVE_EQUIVALENCE_CHECKS_PASSED\n')
