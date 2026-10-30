#!/usr/local/bin/Rscript
# Post-hoc sums of previously computed per-LR B/F/K. Native probability model is not refitted.
main <- function(){
 PROJECT<-'/home/data/t070721/codex_workspace/Bladder Metabolism';BASE<-file.path(PROJECT,'mp10_spatial_communication_v2_all_sections');V1<-file.path(PROJECT,'mp10_spatial_communication_v1')
 .libPaths(unique(c(file.path(V1,'Rlib'),'/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
 suppressPackageStartupMessages(library(data.table));setDTthreads(1);source(file.path(BASE,'scripts','scc_lib.R'))
 args<-commandArgs(TRUE);if(length(args)<3)stop('Usage: scc_common_lr.R SECTION KEEP_TSV OUTPUT_DIR [--source MAIN_DIR]')
 sec<-args[1];keepfile<-args[2];out<-args[3];j<-match('--source',args);src<-if(is.na(j))file.path(BASE,'spatialcellchat',sec,'main') else args[j+1]
 if(normalizePath(src,mustWork=TRUE)==normalizePath(out,mustWork=FALSE))stop('Sensitivity output must differ from original source')
 dir.create(out,recursive=TRUE,showWarnings=FALSE);keep<-fread(keepfile);stopifnot('interaction'%in%names(keep),!anyNA(keep$interaction),!anyDuplicated(keep$interaction))
 native<-fread(file.path(src,'position_LR_support.tsv'));available<-native$interaction;selected<-intersect(keep$interaction,available)
 manifest<-copy(keep);manifest[,present_in_native_inference:=interaction%in%available];scc_put(manifest,'requested_LR_manifest',out);scc_put(native[interaction%in%selected],'selected_LR_position_support',out)
 scc_put(data.table(section=sec,n_original_LR=length(available),n_requested_LR=nrow(keep),n_retained_LR=length(selected),fraction_original_retained=if(length(available))length(selected)/length(available) else NA_real_,fraction_requested_present=if(nrow(keep))length(selected)/nrow(keep) else NA_real_,status=if(length(selected))'computed_posthoc' else 'NE_empty_common_LR',original_inference=normalizePath(src),retained_manifest=normalizePath(keepfile),refitted=FALSE,aggregation='sum per-LR F/K; no pseudocount, no group P screening'),'coverage',out)
 for(lev in c('major','fine')){
  orig<-fread(file.path(src,paste0(lev,'_edges.tsv')));p<-fread(file.path(src,paste0(lev,'_per_LR_aggregation.tsv.gz')));p<-p[interaction%in%selected]
  diagfile<-file.path(src,paste0(lev,'_same_spot_per_LR.tsv.gz'));dg<-fread(diagfile);dg<-dg[interaction%in%selected]
  if(length(selected)){
   z<-p[,.(total_mean_strength=sum(mean_strength),total_flow=sum(flow),summed_LR_supported_pair_denominators=sum(supported_position_pair_denominator)),by=.(source,target)]
   d<-dg[,.(same_spot_flow_excluded=sum(same_spot_flow)),by=.(source,target)];z<-merge(z,d,by=c('source','target'),all=TRUE)
   original.order<-paste(orig$source,orig$target,sep='\r');idx<-match(original.order,paste(z$source,z$target,sep='\r'));stopifnot(!anyNA(idx),nrow(z)==nrow(orig));z<-z[idx]
   for(k in setdiff(names(z),c('source','target')))set(orig,j=k,value=z[[k]])
  }else{
   for(k in c('total_mean_strength','total_flow','summed_LR_supported_pair_denominators','same_spot_flow_excluded'))set(orig,j=k,value=rep(NA_real_,nrow(orig)))
  }
  good<-if('evaluable'%in%names(orig))as.logical(orig$evaluable) else rep(TRUE,nrow(orig));if(all(c('source_eligible','target_eligible')%in%names(orig)))good<-good & orig$source_eligible & orig$target_eligible
  for(k in c('total_mean_strength','total_flow','summed_LR_supported_pair_denominators','same_spot_flow_excluded'))orig[!good,(k):=NA_real_]
  orig[,`:=`(n_LR=length(selected),evaluable=good & length(selected)>0,sensitivity='posthoc_common_LR_original_native_network',source_native_dir=normalizePath(src))]
  scc_put(orig,paste0(lev,'_edges'),out)
  N<-unique(orig$source);stopifnot(nrow(orig)==length(N)^2)
  for(metric in c('total_mean_strength','total_flow')){mat<-dcast(orig,source~target,value.var=metric);mat<-mat[match(N,source),c('source',N),with=FALSE];setnames(mat,'source','node');scc_put(mat,paste0(lev,'_',metric,'_matrix'),out)}
 }
 scc_atomic_lines('COMPLETE_POSTHOC_NO_REFIT',file.path(out,'COMPLETED.txt'));cat('SCC COMMON LR COMPLETE',sec,length(selected),'LR\n')
}
main()
