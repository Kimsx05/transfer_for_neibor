source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R');source(file.path(S,'model_helpers.R'))
suppressPackageStartupMessages(library(nichenetr));source(file.path(S,'nichenet_native_metrics.R'));stopifnot(all(rd('nichenet_native_metrics_validation.tsv')$pass),all(rd('nichenet_sparse_target_validation.tsv')$pass));b<-qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1);pr<-b$profiles;plan<-rd('model_plan.tsv');plan<-plan[order(kind!='receiver',!primary,model_id)];lr<-rd('nichenet_lr_prior.tsv.gz');lt<-readRDS(file.path(REF,'ligand_target_matrix_nsga2r_final.rds'));m20<-scan(file.path(P,'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'),what='character',quiet=TRUE)
dir.create(tb('nichenet'),showWarnings=FALSE);lt_hash<-digest(file.path(REF,'ligand_target_matrix_nsga2r_final.rds'),file=TRUE,algo='sha256');activities<-list();states<-list();bgout<-list();links<-list();detection<-list()
median_detect<-function(pp,genes){genes<-intersect(genes,rownames(b$counts));if(!nrow(pp))return(data.table(gene=genes,median_detection=NA_real_,n_profiles=0,n_detect10=0,n_detect05=0));d<-as.matrix(b$detected_cells[genes,pp$profile_id,drop=FALSE]);d<-sweep(d,2,pp$n_cells,'/');data.table(gene=genes,median_detection=apply(d,1,median),n_profiles=nrow(pp),n_detect10=rowSums(d>=.1),n_detect05=rowSums(d>=.05))}
# All prior ligands, receptors and source coverage; never restricted to CellChat candidates.
allgenes<-union(lr$from,lr$to)
for(pt in unique(pr$profile_type)){
 pp<-pr[profile_type==pt&main_patient&n_cells>=30];dd<-median_detect(pp,allgenes);dd[,profile_type:=pt];detection[[pt]]<-dd
}
put(rbindlist(detection),'nichenet_expression_coverage_by_celltype.tsv.gz')
for(ii in seq_len(nrow(plan))){p<-plan[ii];id<-p$model_id;f<-tb(paste0('de/',id,'.tsv.gz'));if(p$only_dataset!='ALL')next
 if(!file.exists(f)){st<-CJ(target_direction=c('positive','negative','union'),detection_threshold=c(.1,.05));st[,`:=`(model_id=id,direction=if(p$kind=='Epi')'TME_to_Epi'else'Epi_to_TME',receiver=if(p$kind=='Epi')'Epi'else p$target,mask='Full',n_background=NA_integer_,n_targets=NA_integer_,n_sender_agnostic=NA_integer_,n_focused=NA_integer_,status='not_evaluable',reason='insufficient_sample_coverage_or_model_not_evaluable',n_specimens=NA_integer_)];states[[id]]<-st;next};s<-select_model(p,pr);if(is.null(s))next
 # Dataset-only Epi models support descriptive sender DE; NicheNet stability uses pooled LODO instead.
 if(p$only_dataset!='ALL')next
 de<-fread(f);ss<-unique(s$profiles$sample_uid);receiver<-if(p$kind=='Epi')'Epi'else p$target;direction<-if(p$kind=='Epi')'TME_to_Epi'else'Epi_to_TME'
 recpp<-pr[sample_uid%in%ss&profile_type==if(p$kind=='Epi')'Epi_all'else p$target]
 recdet<-median_detect(recpp,unique(lr$to));epidet<-median_detect(pr[sample_uid%in%ss&profile_type=='Epi_all'],unique(lr$from))
 tmep<-pr[sample_uid%in%ss&!profile_type%in%c('Epi_all','H50','H10','H05','L50')&n_cells>=p$min_cells]
 tmed<-rbindlist(lapply(unique(tmep$profile_type),function(pt){z<-median_detect(tmep[profile_type==pt],unique(lr$from));z[,sender_type:=pt];z}),fill=TRUE)
 for(mask in if(p$kind=='Epi')c('Full','M20-masked')else'Full'){
 bg<-intersect(de$gene,rownames(lt));if(mask=='M20-masked')bg<-setdiff(bg,m20);bgout[[paste(id,mask)]]<-data.table(model_id=id,mask=mask,gene=bg)
 for(signset in c('positive','negative','union')){
 target<-de[FDR<.05 & if(signset=='positive')logFC>0 else if(signset=='negative')logFC<0 else TRUE,gene];target<-intersect(target,bg);stopifnot(all(target%in%bg));if(mask=='M20-masked')stopifnot(!any(c(bg,target)%in%m20))
 for(threshold in c(.1,.05)){
 rec<-recdet[median_detection>=threshold,gene];pot<-intersect(unique(lr[to%in%rec,from]),colnames(lt));focused<-if(p$kind=='receiver')intersect(pot,epidet[median_detection>=threshold,gene])else intersect(pot,tmed[median_detection>=threshold,gene])
 tag<-paste(id,mask,signset,threshold,sep='_');st<-data.table(model_id=id,direction=direction,receiver=receiver,mask=mask,target_direction=signset,detection_threshold=threshold,n_background=length(bg),n_targets=length(target),n_sender_agnostic=length(pot),n_focused=length(focused),status='not_evaluable',reason='',n_specimens=length(ss))
 if(!length(target)||length(bg)<2){st[,reason:='insufficient_target_signal'];states[[tag]]<-st;next};if(!length(pot)){st[,reason:='receptor_not_detected_or_no_ligand_target_prior'];states[[tag]]<-st;next}
 of<-tb(paste0('nichenet/',tag,'.tsv.gz'));key<-digest(list(bg,target,pot,'nichenetr2.2.1.1',lt_hash),algo='sha256')
 # Native activity scores are computed once for the receptor-qualified candidate universe; focused rankings are a subset.
 allpot<-intersect(unique(lr[to%in%recdet[median_detection>=.05,gene],from]),colnames(lt))
 nativekey<-digest(list(sort(bg),sort(target),sort(allpot),lt_hash,'nichenetr2.2.1.1'),algo='sha256');nativefile<-tb(paste0('nichenet/native_',nativekey,'.tsv.gz'))
 if(file.exists(nativefile)&&file.exists(paste0(nativefile,'.sha256'))&&readLines(paste0(nativefile,'.sha256'))==nativekey){native<-fread(nativefile)}else{
 cat(tag,'targets',length(target),'union ligands',length(allpot),'\n');native<-nichenet_four_metrics(geneset=target,background_expressed_genes=bg,ligand_target_matrix=lt,potential_ligands=allpot);fwrite(native,nativefile,sep='\t',compress='gzip');writeLines(nativekey,paste0(nativefile,'.sha256'))}
 act<-copy(native[test_ligand%in%pot]);act[,`:=`(model_id=id,direction=direction,receiver=receiver,mask=mask,target_direction=signset,detection_threshold=threshold)];fwrite(act,of,sep='\t',compress='gzip');writeLines(key,paste0(of,'.sha256'))

 setorder(act,-aupr_corrected,test_ligand);act[,rank_sender_agnostic:=seq_len(.N)];act[,rank_focused:=NA_integer_];act[test_ligand%in%focused,rank_focused:=seq_len(.N)];act[,low_information:=length(target)<10];activities[[tag]]<-act
 st[,`:=`(status='completed',reason=if(length(target)<10)'low_information_fewer_than_10_targets'else'')];states[[tag]]<-st
 # Main-model target decomposition for top activity candidates and all named hypotheses; prior weights are not expression effects.
 if(p$primary && threshold==.1){ll<-union(head(act$test_ligand,10),intersect(c('NECTIN2','PVR','VEGFA','VEGFB','VEGFC','PGF'),colnames(lt)));z<-rbindlist(lapply(ll,function(l)data.table(ligand=l,target=target,prior_weight=lt[target,l])));z<-merge(z,de[,.(target=gene,target_logFC=logFC,target_p=PValue,target_q=FDR)],by='target');z[,`:=`(model_id=id,direction=direction,receiver=receiver,mask=mask,target_direction=signset,detection_threshold=threshold)];links[[tag]]<-z}
 }
 }
 }
 put(rbindlist(states,fill=TRUE),'nichenet_analysis_status.tsv');put(rbindlist(activities,fill=TRUE),'nichenet_ligand_activities.tsv.gz')
}
put(rbindlist(states,fill=TRUE),'nichenet_analysis_status.tsv');put(rbindlist(activities,fill=TRUE),'nichenet_ligand_activities.tsv.gz');put(rbindlist(bgout,fill=TRUE),'nichenet_backgrounds.tsv.gz');put(rbindlist(links,fill=TRUE),'nichenet_target_links.tsv.gz');cat('NICHENET COMPLETE\n')
