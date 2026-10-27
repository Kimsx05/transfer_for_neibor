#!/usr/local/bin/Rscript
REQUESTED_R_LIBS <- c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4')
BASE <- normalizePath('mp10_spatial_communication_v1')
.libPaths(unique(c(file.path(BASE,'Rlib'),REQUESTED_R_LIBS,.libPaths())))
suppressPackageStartupMessages({library(data.table); library(Matrix); library(SpatialCellChat)})
setDTthreads(1); options(future.globals.maxSize=12*1024^3,stringsAsFactors=FALSE)
future::plan(future::multicore,workers=as.integer(Sys.getenv('SCC_WORKERS','8')))
# Nested BLAS is externally restricted to one thread in scc_batch.sh.
args<-commandArgs(TRUE); section<-args[1]; mult<-if(length(args)>1)as.numeric(args[2]) else 2
OUT<-file.path(BASE,'spatialcellchat',section,if(mult==2)'main' else 'range3dNN')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(OUT,'COMPLETED.txt'))){cat('Reuse complete frozen output:',OUT,'\n');quit(save='no',status=0)}
IN<-file.path(BASE,'input_review','sections',section)
put<-function(x,name)fwrite(x,file.path(OUT,paste0(name,'.tsv')),sep='\t',na='NA')
t0<-Sys.time(); set.seed(666)
meta<-fread(file.path(IN,'meta.tsv')); coords<-fread(file.path(IN,'coords.tsv'))
genes<-fread(file.path(IN,'genes.tsv'))$gene
X<-as(readMM(gzfile(file.path(IN,'counts.mtx.gz'))),'dgCMatrix'); stopifnot(ncol(X)==nrow(meta),nrow(X)==length(genes),all(X@x>=0),!anyDuplicated(genes),identical(coords$spot_id,meta$spot_id))
dimnames(X)<-list(genes,meta$spot_id)
Qlist<-lapply(c('major','fine'),function(level){z<-fread(file.path(IN,paste0('Q_',level,'.tsv')));stopifnot(identical(z$spot_id,meta$spot_id));q<-as.matrix(z[,-1]);rownames(q)<-z$spot_id;stopifnot(all(q>=0),max(abs(rowSums(q)-1))<1e-6);q});names(Qlist)<-c('major','fine')
xy<-as.matrix(coords[,.(x_um,y_um)]);rownames(xy)<-meta$spot_id
nn<-FNN::get.knn(xy,k=1)$nn.dist[,1]; dNN<-round(median(nn),8); scale.distance<-1.01/min(nn)
# Fixed before inspecting H/L: measured log1p(CP10k), no ALRA; MERINGUE adjusted P<.05; LR either component spatially variable; all assayed complex subunits.
params<-data.table(parameter=c('method','package_version','commit','database','categories','expression','ALRA','feature_selection','feature_adjP','variable_both','position_filter','position_nboot','position_seed','position_quantile','min_links','min_source_positions','min_target_positions','group_P','dNN_um','interaction_range_um','contact_range_um','tol_um','scale_distance','Kh','Hill_n','use_AGAN','group_aggregation','background','threads'),value=c('SpatialCellChat',as.character(packageVersion('SpatialCellChat')),'be88f300464cf970b36a6165ea7d2e41a47a7511','CellChatDB.human bundled','Secreted Signaling;ECM-Receptor;Cell-Cell Contact','log1p library-normalized counts 10000','not applied; measured expression retained','MERINGUE',.05,FALSE,'native filterProbability sampled observed probability quantile; NOT permutation patient P',100,666,.95,10,10,10,'not computed or used',dNN,mult*dNN,1.05*dNN,1e-6,scale.distance,.5,1,TRUE,'official Visium Avg helper sum per-LR B; cross-position only','complete slice; Q full model proportions',future::nbrOfWorkers()))
put(params,'parameters')
# Verify official helper against explicit F/K with zeros cleaned.
p<-sparseMatrix(i=c(1,1,2,3),j=c(1,2,3,1),x=c(.1,.2,.3,.4),dims=c(3,3));diag(p)<-0;p<-drop0(p)
q<-matrix(c(.7,.3,0,.2,.8,0,.4,.5,.1),3,3,byrow=TRUE)
f<-t(q)%*%p%*%q;k<-t(1*(q>0))%*%(p>0)%*%(1*(q>0)); b<-as.matrix(f/k);b[is.nan(b)]<-0
check<-max(abs(b-SpatialCellChat::computeAvgCommunProb_LR_Visium_Avg(p,q,1*(q>0))));stopifnot(check<1e-12)
put(data.table(check='official_helper_vs_F_over_K',max_abs_error=check,status='PASS'),'helper_check')
DB<-SpatialCellChat::subsetDB(SpatialCellChat::CellChatDB.human,search=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact'),non_protein=FALSE)
# Full coverage inventory, including missing subunits; absent genes are not treated as measured zero.
subunits<-function(x){if(x%in%rownames(DB$complex)){v<-as.character(unlist(DB$complex[x,,drop=FALSE]));v[nzchar(v)]}else x}
cov<-rbindlist(lapply(seq_len(nrow(DB$interaction)),function(i){z<-DB$interaction[i,,drop=FALSE];l<-subunits(z$ligand);r<-subunits(z$receptor);data.table(interaction=rownames(z),ligand=z$ligand,receptor=z$receptor,category=z$annotation,ligand_subunits=paste(l,collapse=';'),receptor_subunits=paste(r,collapse=';'),all_subunits_measured=all(c(l,r)%in%genes),missing_subunits=paste(setdiff(c(l,r),genes),collapse=';'))}))
# A single background idents label is bookkeeping; no dominant-type hard assignment, no unused levels.
md<-data.frame(labels=factor(rep('whole_slice',nrow(meta))),samples=section,row.names=meta$spot_id)
chat<-createSpatialCellChat(normalizeData(X),meta=md,group.by='labels',datatype='spatial',coordinates=xy,spatial.factors=list(ratio=1,tol=1e-6))
chat@DB<-DB;chat<-subsetData(chat);rm(X);gc()
chat<-identifyOverExpressedGenes(chat,selection.method='meringue',do.grid=FALSE,thresh.p=.05)
chat<-identifyOverExpressedInteractions(chat,variable.both=FALSE)
# Passing LR.use sorts categories in the native API, essential for contact-vs-diffusible index.
LR<-chat@LR$LRsig
put(data.table(gene=rownames(chat@var.features$features.info),chat@var.features$features.info),'spatial_variable_gene_test')
put(data.table(interaction=rownames(LR),LR),'LR_selected')
cov[,selected_for_inference:=interaction%in%rownames(LR)];put(cov,'gene_lr_coverage')
chat<-computeCommunProb(chat,LR.use=LR,distance.use=TRUE,interaction.range=mult*dNN,contact.range=1.05*dNN,tol=1e-6,scale.distance=scale.distance,contact.dependent=TRUE,use.AGAN=TRUE)
chat<-filterProbability(chat,nboot=100,seed.use=666,thresh=.05)
# Clean explicit zeros before native support counts; synchronize both stores.
chat@net$tmp$prob.cell<-lapply(chat@net$tmp$prob.cell,drop0)
chat@net$prob.cell<-SpatialCellChat:::my_as_sparse3Darray(chat@net$tmp$prob.cell)
dimnames(chat@net$prob.cell)<-list(meta$spot_id,meta$spot_id,names(chat@net$tmp$prob.cell))
chat<-filterCommunication(chat,min.cells=NULL,min.links=10,min.cells.sr=10)
Plist<-chat@net$tmp$prob.cell
saveRDS(list(prob_cell=Plist,spot_id=meta$spot_id,LR=LR,parameters=params),file.path(OUT,'position_network_native.rds'),compress=TRUE)
rm(chat);gc()
for(level in names(Qlist)){
  Q<-Qlist[[level]]; H<-1*(Q>0); n<-ncol(Q); N<-colnames(Q)
  total<-flow<-links<-diagflow<-matrix(0,n,n,dimnames=list(N,N)); rows<-vector('list',length(Plist)); diags<-vector('list',length(Plist)); lrstat<-vector('list',length(Plist))
  for(i in seq_along(Plist)){
    P<-drop0(Plist[[i]]);D<-Diagonal(x=diag(P));diag(P)<-0;P<-drop0(P)
    F<-as.matrix(crossprod(Q,P%*%Q));K<-as.matrix(crossprod(H,(P>0)%*%H))
    B<-SpatialCellChat::computeAvgCommunProb_LR_Visium_Avg(P,Q,H)
    expected<-F/K;expected[is.nan(expected)]<-0;stopifnot(max(abs(expected-B))<1e-10)
    DF<-as.matrix(crossprod(Q,D%*%Q))
    total<-total+B;flow<-flow+F;links<-links+K;diagflow<-diagflow+DF
    rows[[i]]<-data.table(source=rep(N,n),target=rep(N,each=n),interaction=names(Plist)[i],mean_strength=as.vector(B),flow=as.vector(F),supported_position_pair_denominator=as.vector(K))
    diags[[i]]<-data.table(source=rep(N,n),target=rep(N,each=n),interaction=names(Plist)[i],same_spot_flow=as.vector(DF))
    lrstat[[i]]<-data.table(interaction=names(Plist)[i],cross_spot_links=length(P@x),cross_spot_flow=sum(P),same_spot_links=sum(diag(D)>0),same_spot_flow=sum(D),source_positions=sum(rowSums(P)>0),target_positions=sum(colSums(P)>0))
  }
  evaluable<-colSums(Q)>0 & colSums(H)>0
  for(matname in c('total','flow','links','diagflow')){z<-get(matname);z[!evaluable,]<-NA_real_;z[,!evaluable]<-NA_real_;assign(matname,z)}
  edges<-data.table(source=rep(N,n),target=rep(N,each=n),total_mean_strength=as.vector(total),total_flow=as.vector(flow),summed_LR_supported_pair_denominators=as.vector(links),same_spot_flow_excluded=as.vector(diagflow),source_effective_abundance=rep(colSums(Q),n),target_effective_abundance=rep(colSums(Q),each=n),source_support_spots=rep(colSums(H),n),target_support_spots=rep(colSums(H),each=n),evaluable=rep(evaluable,n)&rep(evaluable,each=n),n_LR=length(Plist),group_P_test='not performed; not used')
  put(edges,paste0(level,'_edges'));fwrite(rbindlist(rows),file.path(OUT,paste0(level,'_per_LR_aggregation.tsv.gz')),sep='\t');fwrite(rbindlist(diags),file.path(OUT,paste0(level,'_same_spot_per_LR.tsv.gz')),sep='\t')
  put(data.table(node=N,total,evaluable=evaluable),paste0(level,'_total_mean_strength_matrix'));put(data.table(node=N,flow),paste0(level,'_total_flow_matrix'));put(rbindlist(lrstat),'position_LR_support')
}
Ptotal<-Reduce('+',Plist);sp<-data.table(spot_id=meta$spot_id,outgoing_including_same_spot=rowSums(Ptotal),incoming_including_same_spot=colSums(Ptotal),same_spot=diag(Ptotal));sp[,`:=`(outgoing_cross_spot=outgoing_including_same_spot-same_spot,incoming_cross_spot=incoming_including_same_spot-same_spot)];put(sp,'spot_sender_receiver')
put(data.table(section=section,range_multiplier=mult,elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units='secs')),n_spots=nrow(meta),n_input_genes=length(genes),n_LR=length(Plist),R_version=R.version.string,status='complete'),'run_status')
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
writeLines('COMPLETE',file.path(OUT,'COMPLETED.txt'))
