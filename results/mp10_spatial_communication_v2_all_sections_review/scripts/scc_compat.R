# In-process compatibility only. Installed SpatialCellChat and v1 vendor stay unchanged.
# Native sparseMatrix() inferred dimensions from nonzero coordinates, truncating a
# trailing isolated spot. Explicit dimensions retain all locations, with identical
# distances and native contact/self-diagonal rules.
scc_make_computeCellDistance <- function() {
 fn<-SpatialCellChat::computeCellDistance;hits<-0L
 replace<-function(x){
  if(!is.call(x))return(x)
  if(identical(x[[1]],quote(Matrix::sparseMatrix))){
   stopifnot(!'dims'%in%names(as.list(x)));hits<<-hits+1L
   return(as.call(c(as.list(x),list(dims=quote(c(NC,NC))))))
  }
  for(i in seq_along(x)[-1])x[[i]]<-replace(x[[i]])
  x
 }
 body(fn)<-replace(body(fn));stopifnot(hits==1L);fn
}
scc_make_computeCommunProb <- function(){
 fn<-SpatialCellChat::computeCommunProb;env<-new.env(parent=environment(fn))
 env$computeCellDistance<-scc_make_computeCellDistance();environment(fn)<-env;fn
}
scc_checkpoint_signature <- function(IN,genes,spot_id) {
 files<-file.path(IN,c('counts.mtx.gz','coords.tsv','genes.tsv'))
 stopifnot(all(file.exists(files)));fi<-file.info(files)
 list(input_files=normalizePath(files),input_sizes=fi$size,input_mtime=as.numeric(fi$mtime),genes=genes,spot_id=spot_id,version=as.character(utils::packageVersion('SpatialCellChat')),feature_method='meringue',feature_p=.05,variable_both=FALSE,expression='normalizeData all genes library10000 log1p once',seed=666)
}
scc_save_feature_checkpoint <- function(chat,signature,path,origin){
 checkpoint<-list(signature=signature,var.features=chat@var.features,LR=chat@LR,signaling_genes=rownames(chat@data.signaling),origin=origin)
 tmp<-paste0(path,'.tmp.',Sys.getpid());saveRDS(checkpoint,tmp,compress=TRUE);stopifnot(file.rename(tmp,path))
}
scc_restore_features <- function(chat,IN,OUT,genes,spot_id){
 signature<-scc_checkpoint_signature(IN,genes,spot_id);path<-file.path(OUT,'feature_selection_checkpoint.rds')
 if(file.exists(path)){
  cp<-readRDS(path);stopifnot(identical(cp$signature,signature),identical(cp$signaling_genes,rownames(chat@data.signaling)))
  chat@var.features<-cp$var.features;chat@LR<-cp$LR;cat('RESTORED MERINGUE/LR checkpoint; no feature re-fit\n');return(list(chat=chat,restored=TRUE,signature=signature))
 }
 # Earlier failed runs already wrote full feature statistics and exact selected LR.
 # Rebuild small state from these tables and verify official LR selection exactly.
 fp<-file.path(OUT,'spatial_variable_gene_test.tsv');lp<-file.path(OUT,'LR_selected.tsv');par<-file.path(OUT,'parameters.tsv')
 if(file.exists(fp)&&file.exists(lp)&&file.exists(par)){
  old<-data.table::fread(par);expected<-c(package_version='0.1.0',feature_selection='MERINGUE',feature_adjP='0.05',variable_both='FALSE',expression='log1p library-normalized counts 10000')
  for(k in names(expected))stopifnot(identical(as.character(old[parameter==k,value]),unname(expected[k])))
  # Frozen input must predate completed old selection tables; no large-file hash.
  stopifnot(max(signature$input_mtime)<=min(as.numeric(file.info(c(fp,lp))$mtime)))
  z<-as.data.frame(data.table::fread(fp));stopifnot(!anyDuplicated(z$gene),identical(z$gene,rownames(chat@data.signaling)),all(c('p.adj','observed')%in%names(z)))
  rownames(z)<-z$gene;z$gene<-NULL;features<-rownames(z)[which(z$p.adj<.05 & !is.nan(z$observed))]
  chat@var.features<-list(features=features,features.info=z)
  ids<-data.table::fread(lp,select='interaction')$interaction
  chat<-SpatialCellChat::identifyOverExpressedInteractions(chat,variable.both=FALSE)
  stopifnot(!anyDuplicated(ids),setequal(rownames(chat@LR$LRsig),ids))
  # Preserve the exact saved native order (parallel LR-selection enumeration can differ).
  chat@LR$LRsig<-chat@LR$LRsig[ids,,drop=FALSE];stopifnot(identical(rownames(chat@LR$LRsig),ids))
  scc_save_feature_checkpoint(chat,signature,path,'validated frozen prior feature/LR TSV; official LR selection replay matched exactly')
  cat('RESTORED prior MERINGUE TSV:',length(features),'features;',length(ids),'exact native LR; no MERINGUE re-fit\n')
  return(list(chat=chat,restored=TRUE,signature=signature))
 }
 list(chat=chat,restored=FALSE,signature=signature)
}
