#!/usr/local/bin/Rscript
# Full script is parsed as one function before inference to keep running jobs immutable.
main <- function() {
 REQUESTED_R_LIBS<-c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4')
 PROJECT<-'/home/data/t070721/codex_workspace/Bladder Metabolism';V1<-file.path(PROJECT,'mp10_spatial_communication_v1');BASE<-file.path(PROJECT,'mp10_spatial_communication_v2_all_sections')
 .libPaths(unique(c(file.path(V1,'Rlib'),REQUESTED_R_LIBS,.libPaths())))
 suppressPackageStartupMessages({library(data.table);library(Matrix);library(SpatialCellChat)})
 source(file.path(BASE,'scripts','scc_lib.R'),local=FALSE)
 source(file.path(BASE,'scripts','scc_compat.R'),local=FALSE);compute_prob<-scc_make_computeCommunProb()
 setDTthreads(1);options(future.globals.maxSize=12*1024^3,stringsAsFactors=FALSE);workers<-as.integer(Sys.getenv('SCC_WORKERS','6'));stopifnot(workers>=1,workers<=8)
 future::plan(future::multicore,workers=workers);on.exit(future::plan(future::sequential),add=TRUE)
 args<-commandArgs(TRUE);if(!length(args))stop('Usage: Rscript scc_run.R SECTION [--aggregate-from POSITION_NETWORK_RDS] [--output OUTPUT_DIR]')
 section<-args[1];opt<-function(k,default=NULL){j<-match(k,args);if(is.na(j))default else args[j+1]}
 OUT<-opt('--output',file.path(BASE,'spatialcellchat',section,'main'));dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
 if(file.exists(file.path(OUT,'COMPLETED.txt'))){cat('REUSE completed:',OUT,'\n');return(invisible(TRUE))}
 IN<-file.path(BASE,'input_review','sections',section);t0<-Sys.time();set.seed(666)
 meta<-fread(file.path(IN,'meta.tsv'));coords<-fread(file.path(IN,'coords.tsv'));genes<-fread(file.path(IN,'genes.tsv'))$gene
 stopifnot(identical(coords$spot_id,meta$spot_id),!anyDuplicated(meta$spot_id),!anyDuplicated(genes))
 Qlist<-scc_load_Q(IN,meta$spot_id);support<-scc_support(BASE,IN,section,Qlist)
 epgood<-support[level=='major'&node%in%c('Epi_MP10_high','Epi_MP10_low'),all(type_evaluable)]
 hlstatus<-if(epgood)'comparable' else 'network_only';if('hl_status'%in%names(meta))hlstatus<-unique(meta$hl_status)[1]
 xy<-as.matrix(coords[,.(x_um,y_um)]);rownames(xy)<-meta$spot_id;stopifnot(all(is.finite(xy)),nrow(unique(coords[,.(x_um,y_um)]))==nrow(coords))
 nearest_nn<-FNN::get.knn(xy,k=1)$nn.dist[,1];observed_nn<-median(nearest_nn);dNN<-100;scale.distance<-.0101
 connected<-which(nearest_nn<=200+1e-6);stopifnot(length(connected)>0);native_inferred_dim<-max(connected);dimension_padding<-nrow(meta)-native_inferred_dim
 cat('Distance shape audit: full spots',nrow(meta),'native inferred dimension',native_inferred_dim,'trailing isolated spots retained',dimension_padding,'\n')
 params<-data.table(parameter=c('method','package_version','commit','database','categories','expression','ALRA','feature_selection','feature_adjP','variable_both','position_filter','position_nboot','position_seed','position_quantile','min_links','min_source_positions','min_target_positions','group_P','dNN_um','observed_median_nn_um','interaction_range_um','contact_range_um','tol_um','scale_distance','Kh','Hill_n','use_AGAN','group_aggregation','background','threads','type_support_rule'),value=as.character(c('SpatialCellChat',as.character(packageVersion('SpatialCellChat')),'be88f300464cf970b36a6165ea7d2e41a47a7511','CellChatDB.human bundled','Secreted Signaling;ECM-Receptor;Cell-Cell Contact','log1p library-normalized counts 10000','not applied; measured expression retained','MERINGUE',.05,FALSE,'native filterProbability sampled observed probability quantile; NOT permutation patient P',100,666,.95,10,10,10,'not computed or used',100,observed_nn,200,105,1e-6,scale.distance,.5,1,TRUE,'official Visium Avg helper sum per-LR B; cross-position only','complete slice; Q full model proportions',workers,'at least 10 spots posterior mean abundance>=0.01 and total abundance>1; no spot-Q truncation')))
 stopifnot(as.character(packageVersion('SpatialCellChat'))=='0.1.0');params<-rbind(params,data.table(parameter=c('dimension_compatibility','native_inferred_distance_dimension','trailing_isolated_spots_retained','dimension_compatibility_triggered'),value=c('native distance sparseMatrix explicit dims=c(NC,NC); all spots and native diagonal rules retained',native_inferred_dim,dimension_padding,dimension_padding>0)))
 # Preserve prior feature-selection parameters until old TSV checkpoint validation.
 if(!file.exists(file.path(OUT,'parameters.tsv')))scc_put(params,'parameters',OUT)
 cache<-opt('--aggregate-from',file.path(OUT,'position_network_native.rds'));reused_position<-file.exists(cache)
 if(reused_position){
  cat('Resume frozen native P aggregation:',cache,'\n');net<-readRDS(cache);stopifnot(identical(net$spot_id,meta$spot_id));Plist<-net$prob_cell;LR<-net$LR
  oldpar<-net$parameters
  if(!is.null(oldpar)){oldpar<-as.data.table(oldpar);for(k in c('interaction_range_um','contact_range_um','tol_um','scale_distance','position_seed','position_quantile','Kh','Hill_n')){old<-as.numeric(oldpar[parameter==k,value]);new<-as.numeric(params[parameter==k,value]);stopifnot(length(old)==1,length(new)==1,abs(old-new)<1e-10)}}
  inference_seconds<-if(!is.null(net$inference_seconds))net$inference_seconds else NA_real_;rm(net);gc()
 }else{
  X<-as(readMM(gzfile(file.path(IN,'counts.mtx.gz'))),'CsparseMatrix');stopifnot(ncol(X)==nrow(meta),nrow(X)==length(genes),all(is.finite(X@x)),all(X@x>=0),all(colSums(X)>0));dimnames(X)<-list(genes,meta$spot_id)
  DB<-SpatialCellChat::subsetDB(SpatialCellChat::CellChatDB.human,search=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact'),non_protein=FALSE)
  subunits<-function(x){if(x%in%rownames(DB$complex)){v<-as.character(unlist(DB$complex[x,,drop=FALSE]));v[nzchar(v)]}else x}
  cov<-rbindlist(lapply(seq_len(nrow(DB$interaction)),function(i){z<-DB$interaction[i,,drop=FALSE];l<-subunits(z$ligand);r<-subunits(z$receptor);data.table(interaction=rownames(z),ligand=z$ligand,receptor=z$receptor,category=z$annotation,ligand_subunits=paste(l,collapse=';'),receptor_subunits=paste(r,collapse=';'),all_subunits_measured=all(c(l,r)%in%genes),missing_subunits=paste(setdiff(c(l,r),genes),collapse=';'))}))
  md<-data.frame(labels=factor(rep('whole_slice',nrow(meta))),samples=section,row.names=meta$spot_id)
  chat<-createSpatialCellChat(normalizeData(X),meta=md,group.by='labels',datatype='spatial',coordinates=xy,spatial.factors=list(ratio=1,tol=1e-6));chat@DB<-DB;chat<-subsetData(chat);rm(X);gc()
  restore<-scc_restore_features(chat,IN,OUT,genes,meta$spot_id);chat<-restore$chat
  if(!restore$restored){chat<-identifyOverExpressedGenes(chat,selection.method='meringue',do.grid=FALSE,thresh.p=.05);chat<-identifyOverExpressedInteractions(chat,variable.both=FALSE);scc_save_feature_checkpoint(chat,restore$signature,file.path(OUT,'feature_selection_checkpoint.rds'),'native feature and LR selection')}
  LR<-chat@LR$LRsig;rm(restore);scc_put(params,'parameters',OUT)
  scc_put(data.table(gene=rownames(chat@var.features$features.info),chat@var.features$features.info),'spatial_variable_gene_test',OUT);scc_put(data.table(interaction=rownames(LR),LR),'LR_selected',OUT)
  cov[,selected_for_inference:=interaction%in%rownames(LR)];scc_put(cov,'gene_lr_coverage',OUT)
  if(nrow(LR)>0){
   chat<-compute_prob(chat,LR.use=LR,distance.use=TRUE,interaction.range=200,contact.range=105,tol=1e-6,scale.distance=.0101,contact.dependent=TRUE,use.AGAN=TRUE,Kh=.5,n=1)
   chat<-filterProbability(chat,nboot=100,seed.use=666,thresh=.05);chat@net$tmp$prob.cell<-lapply(chat@net$tmp$prob.cell,drop0);chat@net$prob.cell<-SpatialCellChat:::my_as_sparse3Darray(chat@net$tmp$prob.cell);dimnames(chat@net$prob.cell)<-list(meta$spot_id,meta$spot_id,names(chat@net$tmp$prob.cell));chat<-filterCommunication(chat,min.cells=NULL,min.links=10,min.cells.sr=10);Plist<-chat@net$tmp$prob.cell
  }else{Plist<-list();cat('No selected LR: computed empty supported network; no artificial positive edges\n')}
  inference_seconds<-as.numeric(difftime(Sys.time(),t0,units='secs'));tmp<-file.path(OUT,paste0('position_network_native.rds.tmp.',Sys.getpid()));saveRDS(list(prob_cell=Plist,spot_id=meta$spot_id,LR=LR,parameters=params,inference_seconds=inference_seconds),tmp,compress=TRUE);stopifnot(file.rename(tmp,file.path(OUT,'position_network_native.rds')));rm(chat);gc()
 }
 agg<-scc_aggregate(Plist,meta$spot_id,Qlist,support,OUT,section)
 scc_put(data.table(section=section,range_multiplier=2,elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units='secs')),aggregation_seconds=agg$aggregation_seconds,inference_and_position_QC_seconds=inference_seconds,n_spots=nrow(meta),n_input_genes=length(genes),n_LR=length(Plist),R_version=R.version.string,status='complete',hl_status=hlstatus,reused_position_cache=reused_position,position_cache=normalizePath(cache,mustWork=FALSE)),'run_status',OUT)
 writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'));scc_atomic_lines('COMPLETE',file.path(OUT,'COMPLETED.txt'));cat('SCC COMPLETE',section,'\n')
}
main()
