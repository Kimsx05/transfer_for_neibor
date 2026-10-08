source('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1/code/00_common.R')
stopifnot(file.exists(out('02_ucell_scores/COMPLETE')))
future::plan('sequential')
M <- qread(out('02_ucell_scores/M_RAW.qs'));Z <- qread(out('02_ucell_scores/M_Z.qs'))
instances <- data.table(instance=c('FULL_Z','FULL_RAW',as.vector(rbind(paste0('BAL_Z_seed',42:44),paste0('BAL_RAW_seed',42:44)))))
instances[,directory:=ifelse(instance=='FULL_Z','03_full_z',ifelse(instance=='FULL_RAW','04_full_raw',paste0('05_balanced/',instance)))]
write_tsv(instances,out('instance_manifest.tsv'))
for(ii in seq_len(nrow(instances))) {
 nm<-instances$instance[ii];dir<-out(instances$directory[ii]);dir.create(dir,recursive=TRUE,showWarnings=FALSE)
 ids<-if(startsWith(nm,'FULL')) rownames(M) else qread(out('05_balanced',paste0('seed',sub('.*seed','',nm),'_cell_ids.qs')))
 X<-if(grepl('_RAW',nm)) M[ids,mp,drop=FALSE] else Z[ids,mp,drop=FALSE]
 stopifnot(ncol(X)==14,identical(rownames(X),ids),all(is.finite(X)))
 id<-list(config_hash=config_hash,matrix_hash=digest(X,'sha256'),instance=nm)
 ip<-file.path(dir,'instance_identity.qs');if(file.exists(ip)) stopifnot(identical(id,qread(ip))) else qsave_atomic(id,ip)
 gp<-file.path(dir,'SNN.qs')
 if(file.exists(gp)) {snn<-qread(gp);logmsg(nm,'reuse SNN')} else {
  logmsg(nm,'FindNeighbors cells',nrow(X),'features',ncol(X));set.seed(42)
  nn<-FindNeighbors(object=X,k.param=20,nn.method='annoy',n.trees=50,annoy.metric='euclidean',l2.norm=FALSE,compute.SNN=TRUE,prune.SNN=1/15,verbose=TRUE)
  snn<-nn$snn;stopifnot(identical(rownames(snn),ids),identical(colnames(snn),ids),inherits(snn,'sparseMatrix'))
  qsave_atomic(snn,gp);rm(nn)
 }
 rr<-if(nm=='FULL_Z') c(.4,.2,.6) else .4
 for(res in rr) {
  cp<-file.path(dir,paste0('clusters_res',res,'.tsv.gz'))
  if(!file.exists(cp)) {
   logmsg(nm,'Leiden resolution',res);set.seed(42)
   cl<-FindClusters(object=snn,algorithm=4,leiden_method='leidenbase',modularity.fxn=1,resolution=res,random.seed=42,n.iter=10,group.singletons=TRUE,verbose=TRUE)
   stopifnot(setequal(rownames(cl),ids),ncol(cl)==1)
   labels<-data.table(cell_id=rownames(cl),cluster=paste0('C',cl[[1]]));labels<-labels[match(ids,cell_id)]
   write_tsv(labels,cp)
  }
 }
 up<-file.path(dir,'umap.qs')
 if(!file.exists(up)) {
  logmsg(nm,'UMAP direct 14D');set.seed(42)
  # Seurat matrix method calls uwot::umap without PCA; use 1 thread for reproducibility.
  old<-options(mc.cores=1)
  u<-RunUMAP(object=X,umap.method='uwot',n.components=2L,metric='euclidean',n.neighbors=30L,min.dist=.3,seed.use=42L,uwot.sgd=FALSE,verbose=TRUE)
  options(old)
  emb<-Embeddings(u);stopifnot(setequal(rownames(emb),ids),ncol(emb)==2)
  emb<-emb[ids,,drop=FALSE];colnames(emb)<-c('UMAP_1','UMAP_2');qsave_atomic(emb,up)
  write_tsv(as.data.table(emb,keep.rownames='cell_id'),file.path(dir,'umap_coordinates.tsv.gz'))
 }
 writeLines('complete',file.path(dir,'COMPLETE'));logmsg(nm,'complete');rm(snn,X);gc()
}
capture.output(sessionInfo(),file=out('logs/sessionInfo_graphs.txt'))
writeLines('complete',out('logs/GRAPHS_COMPLETE'))
