.libPaths(c(file.path(commandArgs(trailingOnly=TRUE)[1],"Rlib"),.libPaths()))
suppressPackageStartupMessages({library(SingleR);library(Matrix);library(rhdf5);library(BiocParallel);library(S4Vectors)})
options(mc.cores=8);set.seed(123);Sys.setenv(OMP_NUM_THREADS=8,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
argv<-commandArgs(trailingOnly=TRUE);root<-normalizePath(argv[1]);bp<-MulticoreParam(workers=8,progressbar=FALSE,RNGseed=123)
read_matrix<-function(file){s<-h5read(file,'shape');x<-new('dgCMatrix',i=as.integer(h5read(file,'indices')),p=as.integer(h5read(file,'indptr')),x=as.double(h5read(file,'data')),Dim=as.integer(rev(s)));rownames(x)<-as.character(h5read(file,'genes'));colnames(x)<-as.character(h5read(file,'cells'));x@x<-log1p(x@x*rep(1e4/pmax(Matrix::colSums(x),1),diff(x@p)));x}
ref<-read_matrix(file.path(root,'temporary/reference_shared_counts.h5'));meta<-read.delim(file.path(root,'metadata/SingleR_reference_selection.tsv'),check.names=FALSE,stringsAsFactors=FALSE);stopifnot(identical(meta$cell_id,colnames(ref)),!anyDuplicated(rownames(ref)))
writeLines(c(R.version.string,paste('SingleR',packageVersion('SingleR')),'counts summed by approved symbols, same14327genes; log1p10k; seed123; de.method=wilcox; de.n=20 per pair; aggr.ref=TRUE; fixed training per hierarchy; num.threads8; chunks5000; prune per entire section'),file.path(root,'metadata/SingleR_method.txt'))
for(head in c('major','subtype')){
 labels<-if(head=='major')meta$major_v21 else meta$reference_celltype_v2;model_file<-file.path(root,'models',paste0('SingleR_',head,'.rds'))
 if(file.exists(model_file)){trained<-rebuildIndex(readRDS(model_file),num.threads=8)}else{cat('TRAIN',head,length(unique(labels)),ncol(ref),'\n');flush.console();trained<-trainSingleR(ref,labels=labels,test.genes=rownames(ref),de.method='wilcox',de.n=20,aggr.ref=TRUE,num.threads=8,BPPARAM=bp);saveRDS(trained,model_file,compress=TRUE)}
 for(section in c('Cervilla','D1_Johnson')){
  outfile<-file.path(root,'metadata',paste0(section,'_SingleR_',head,'.tsv.gz'));if(file.exists(outfile))next
  query<-read_matrix(file.path(root,'temporary',paste0(section,'_shared_counts.h5')));stopifnot(identical(rownames(query),rownames(ref)));all<-list()
  for(start in seq(1,ncol(query),by=5000)){
   end<-min(ncol(query),start+4999);cat('CLASSIFY',head,section,start,end,'\n');flush.console();all[[length(all)+1]]<-classifySingleR(query[,start:end,drop=FALSE],trained,prune=FALSE,num.threads=8,BPPARAM=SerialParam());
  }
  z<-do.call(rbind,all);rm(all);gc();stopifnot(nrow(z)==ncol(query));bad<-pruneScores(z);delta<-getDeltaFromMedian(z);thresholds<-pruneScores(z,get.thresholds=TRUE);idx<-match(z$labels,colnames(z$scores));selected<-z$scores[cbind(seq_len(nrow(z)),idx)];ord<-t(apply(z$scores,1,order,decreasing=TRUE));second<-colnames(z$scores)[ord[,2]];best<-colnames(z$scores)[ord[,1]]
  out<-data.frame(cell_id=colnames(query),candidate=z$labels,pruned_label=ifelse(bad,'PRUNED_UNCERTAIN',z$labels),pruned=bad,delta_median=delta,delta_next=z$delta.next,selected_score=selected,top_initial=best,second_initial=second,top_initial_score=z$scores[cbind(seq_len(nrow(z)),ord[,1])],second_initial_score=z$scores[cbind(seq_len(nrow(z)),ord[,2])],run_status='ran_cell_level',stringsAsFactors=FALSE)
  write.table(out,gzfile(outfile),sep='\t',quote=FALSE,row.names=FALSE);write.table(data.frame(label=names(thresholds),threshold=as.numeric(thresholds)),file.path(root,'metadata',paste0(section,'_SingleR_',head,'_pruning.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  # Compact diagnostic candidates within each biological major from initial fine scores, not a separate fine-tuned prediction.
  if(head=='subtype'){
   map<-setNames(meta$major_v21,meta$reference_celltype_v2);res<-data.frame(cell_id=colnames(query))
   for(major in unique(meta$major_v21)){
    cols<-which(map[colnames(z$scores)]==major);v<-z$scores[,cols,drop=FALSE];k<-max.col(v,ties.method='first');tag<-gsub('[^A-Za-z0-9]','_',major);res[[paste0(tag,'_candidate')]]<-colnames(v)[k];res[[paste0(tag,'_score')]]<-v[cbind(seq_len(nrow(v)),k)];res[[paste0(tag,'_margin')]]<-if(length(cols)>1)apply(v,1,function(x){y<-sort(x,decreasing=TRUE);y[1]-y[2]})else NA_real_
   }
   write.table(res,gzfile(file.path(root,'metadata',paste0(section,'_SingleR_within_major_initial.tsv.gz'))),sep='\t',quote=FALSE,row.names=FALSE)
  }
  cat('DONE',head,section,'cells',nrow(out),'pruned',sum(bad),'\n');flush.console();rm(query,z,out);gc()
 }
}
writeLines(capture.output(sessionInfo()),file.path(root,'metadata/SingleR_session.txt'));cat('SINGLER_COMPLETE\n')
