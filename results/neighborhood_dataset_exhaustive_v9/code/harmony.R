a<-commandArgs(TRUE); suppressPackageStartupMessages(library(harmony));library(jsonlite)
x<-as.matrix(read.table(a[1],sep='\t',header=FALSE));m<-read.delim(a[2],check.names=FALSE);rownames(m)<-m$spot_uid;rownames(x)<-m$spot_uid
key<-a[3];theta<-as.numeric(a[4]);lambda<-as.numeric(a[5]);out<-a[6]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
run<-function(maxit){set.seed(123);RunHarmony(data_mat=x,meta_data=m,vars_use=key,theta=theta,lambda=lambda,sigma=.1,nclust=100,max_iter=maxit,early_stop=TRUE,ncores=1,return_object=TRUE,verbose=TRUE,.options=harmony_options(max.iter.cluster=20,epsilon.cluster=1e-3,epsilon.harmony=1e-4))}
h<-run(20); obj<-h$objective_harmony; iterations<-length(obj)-1L; extended<-FALSE
if(iterations>=20){h<-run(50);obj<-h$objective_harmony;iterations<-length(obj)-1L;extended<-TRUE}
z<-t(h$getZcorr());stopifnot(nrow(z)==nrow(x),ncol(z)==ncol(x),all(is.finite(z)))
write.table(z,file.path(out,'embedding.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(data.frame(iteration=seq_along(obj)-1,objective=obj),file.path(out,'objective.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write_json(list(implementation='R harmony',version=as.character(packageVersion('harmony')),batch_key=key,theta=theta,lambda=lambda,sigma=.1,nclust=100,max_iter=if(extended)50 else 20,early_stop=TRUE,ncores=1,seed=123,max_iter_cluster=20,epsilon_cluster=.001,epsilon_harmony=.0001,iterations=iterations,extended=extended,convergence_status=if(iterations < if(extended)50 else 20)'converged' else 'iteration_limit',fit_n=nrow(x),fit_scope='all anchors joint; not independent holdout representation'),file.path(out,'harmony_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
