REQUESTED_R_LIBS <- c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4')
BASE <- normalizePath('mp10_spatial_communication_v1')
.libPaths(unique(c(file.path(BASE,'Rlib'),REQUESTED_R_LIBS,.libPaths())))
suppressPackageStartupMessages({library(data.table); library(Matrix); library(SpatialCellChat)})
setDTthreads(1);args<-commandArgs(TRUE);section<-args[1];mult<-as.numeric(args[2])
OUT<-file.path(BASE,'spatialcellchat',section,if(mult==2)'main' else 'range3dNN')
IN<-file.path(BASE,'input_review','sections',section)
put<-function(x,name)fwrite(x,file.path(OUT,paste0(name,'.tsv')),sep='\t',na='NA')
t0<-Sys.time();meta<-fread(file.path(IN,'meta.tsv'));genes<-fread(file.path(IN,'genes.tsv'))$gene
Qlist<-lapply(c('major','fine'),function(level){z<-fread(file.path(IN,paste0('Q_',level,'.tsv')));stopifnot(identical(z$spot_id,meta$spot_id));q<-as.matrix(z[,-1]);rownames(q)<-z$spot_id;q});names(Qlist)<-c('major','fine')
net<-readRDS(file.path(OUT,'position_network_native.rds'));Plist<-net$prob_cell
stopifnot(identical(net$spot_id,meta$spot_id));rm(net);gc()
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
