.libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(qs);library(Matrix);library(CellChat)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1'
source(file.path(root,'scripts/fast_trimean.R'));stopifnot(file.exists(file.path(root,'tables/triMean_acceleration_validation.tsv')));setDTthreads(1);future::plan('sequential');options(future.globals.maxSize=8*1024^3)
s<-fread(file.path(root,'tables/samples.tsv'));m<-fread(file.path(root,'tables/cells.tsv.gz'));x<-qread(file.path(root,'objects/signaling_RNA_data.qs'),nthreads=1)
m<-merge(m,s[,.(sample_uid,patient_uid_recorded,tissue_status,T_category,MP10_group)],by='sample_uid',sort=FALSE);stopifnot(all(m$barcode%in%colnames(x)));x<-x[,m$barcode]
DB<-CellChatDB.human;fwrite(as.data.table(DB$interaction),file.path(root,'tables/CellChat_DB_interactions.tsv'),sep='\t');fwrite(as.data.table(DB$complex,keep.rownames='complex'),file.path(root,'tables/CellChat_DB_complexes.tsv'),sep='\t')
m[,aggregate_id:=paste(sample_uid,communication_type,sep='|||')];groups<-unique(m$aggregate_id);gg<-match(m$aggregate_id,groups);gmat<-sparseMatrix(i=seq_len(nrow(m)),j=gg,x=1,dims=c(nrow(m),length(groups)));cnt<-tabulate(gg,nbins=length(groups))
# Cache all database gene means/detection fractions once; downstream shortlist uses the same scale.
means<-x%*%gmat%*%Diagonal(x=1/cnt);detection<-(x>0)%*%gmat%*%Diagonal(x=1/cnt);dimnames(means)<-dimnames(detection)<-list(rownames(x),groups)
gmeta<-m[,.(n_cells=.N,dataset=dataset[1],sample_uid=sample_uid[1],patient_uid=patient_uid_recorded[1],major=major[1],cell_type=communication_type[1]),by=aggregate_id]
qsave(list(mean=means,detection=detection,metadata=gmeta),file.path(root,'objects/sample_expression_support.qs'),preset='fast');fwrite(gmeta,file.path(root,'tables/sample_subtype_counts.tsv'),sep='\t');rm(means,detection,gmat);gc()
# Stratified cap protects rare populations; randomness is reproducible and never a biological replicate.
set.seed(123);sel<-m[,.I[if(.N>100)sample.int(.N,100) else seq_len(.N)],by=.(sample_uid,communication_type)]$V1
m<-m[sel];x<-x[,sel];m[,aggregate_id:=NULL]
fwrite(m[,.(n_sampled=.N),by=.(dataset,sample_uid,patient_uid_recorded,communication_type)],file.path(root,'tables/CellChat_cell_contributions.tsv'),sep='\t')
specs<-list(Normal_Tumor=list(col='tissue_status',groups=c('Normal','Tumor')),T1_T2=list(col='T_category',groups=c('T1','T2')),MP10_low_high=list(col='MP10_group',groups=c('low','high')))
tasks<-list();for(ds in c('ALL',sort(unique(m$dataset))))for(co in names(specs))for(gr in specs[[co]]$groups){ids<-which((ds=='ALL'|m$dataset==ds)&!is.na(m[[specs[[co]]$col]])&m[[specs[[co]]$col]]==gr);if(length(ids))tasks[[length(tasks)+1]]<-list(scope=ds,comparison=co,condition=gr,ids=ids)}
dir.create(file.path(root,'objects/networks'),showWarnings=FALSE)
worker<-function(task){
 key<-paste(task$scope,task$comparison,task$condition,sep='__');fp<-file.path(root,'objects/networks',paste0(key,'.qs'));logfile<-file.path(root,'logs',paste0('CellChat_',key,'.log'))
 if(file.exists(fp))return(data.table(network_id=key,status='reused'))
 sink(logfile);on.exit(sink(),add=TRUE)
 tryCatch({
  md<-copy(m[task$ids]);counts<-md[,.N,by=communication_type];keep<-counts[N>=30,communication_type];md<-md[communication_type%in%keep]
  if(length(keep)<2)stop('fewer_than_2_evaluable_cell_types')
  cat('START',key,format(Sys.time()),nrow(md),'cells',length(keep),'types\n')
  meta<-data.frame(labels=factor(md$communication_type),row.names=md$barcode)
  cc<-createCellChat(object=x[,md$barcode],meta=meta,group.by='labels');cc@DB<-DB
  cc<-subsetData(cc);cc<-identifyOverExpressedGenes(cc);cc<-identifyOverExpressedInteractions(cc)
  cc<-computeCommunProb_fast(cc,type='triMean',population.size=FALSE,raw.use=TRUE,nboot=100,seed.use=123)
  cc<-filterCommunication(cc,min.cells=30);cc<-computeCommunProbPathway(cc);cc<-aggregateNet(cc)
  # Small numerical network cache only, not full expression-bearing CellChat objects.
  net<-list(net=cc@net,netP=cc@netP,LR=cc@LR$LRsig,counts=counts,contributions=md[,.(n_cells=.N),by=.(dataset,sample_uid,patient_uid_recorded,communication_type)],scope=task$scope,comparison=task$comparison,condition=task$condition,params=list(version=as.character(packageVersion('CellChat')),type='triMean',population.size=FALSE,min.cells=30,nboot=100,seed=123,sample_type_cap=100,data_max=max(cc@data.signaling)))
  qsave(net,fp,preset='fast');cat('DONE',format(Sys.time()),'\n');data.table(network_id=key,status='completed',n_cells=nrow(md),n_types=length(keep),n_LR=nrow(cc@LR$LRsig))
 },error=function(e){cat('ERROR',conditionMessage(e),'\n');data.table(network_id=key,status='failed',error=conditionMessage(e))})
}
# Two CPU processes, each with one BLAS thread; independent condition networks.
if(Sys.getenv('MP10_ALL_ONLY')=='1')tasks<-Filter(function(z)z$scope=='ALL',tasks)
if(Sys.getenv('MP10_DATASETS_ONLY')=='1')tasks<-Filter(function(z)z$scope!='ALL',tasks)
if(Sys.getenv('MP10_SKIP_RUNNING_NORMAL')=='1')tasks<-Filter(function(z)!(z$scope=='ALL' && z$comparison=='Normal_Tumor'),tasks)
res<-parallel::mclapply(tasks,worker,mc.cores=as.integer(Sys.getenv('MP10_WORKERS','4')),mc.preschedule=FALSE)
fwrite(rbindlist(res,fill=TRUE),file.path(root,'tables/CellChat_network_status.tsv'),sep='\t');cat('NETWORKS DONE\n')
