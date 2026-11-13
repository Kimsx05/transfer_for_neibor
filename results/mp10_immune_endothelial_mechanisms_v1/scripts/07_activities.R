source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R');source(file.path(S,'model_helpers.R'))
suppressPackageStartupMessages({library(edgeR);library(limma);library(decoupleR)})
b<-qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1);pr<-b$profiles;plan<-rd('model_plan.tsv');m20<-scan(file.path(P,'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'),what='character',quiet=TRUE)
nets<-list(CollecTRI=fread(file.path(REF,'collectri.tsv')),PROGENy=fread(file.path(REF,'progeny_top500.tsv')),DoRothEA_AB=fread(file.path(REF,'dorothea_AB.tsv')));coverage<-list();activities<-list();results<-list();stat<-list();idx<-0
for(family in c('Epi','Treg','Tfh/CXCL13-like','Arterial EC','Capillary EC')){
 pp<-pr[n_cells>=30 & if(family=='Epi')profile_type %in% c('L50','H50','H10','H05')else profile_type==family];if(nrow(pp)<3)next
 y<-DGEList(as.matrix(b$counts[,pp$profile_id,drop=FALSE]));keep<-filterByExpr(y,group=if(family=='Epi')factor(pp$profile_type)else rep(1,nrow(pp)));y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);mat<-cpm(y,log=TRUE,prior.count=2);mat<-mat-rowMeans(mat)
 for(resource in names(nets))for(mask in if(family=='Epi')c('Full','M20-masked')else'Full'){
 net<-copy(nets[[resource]]);mm<-mat;if(mask=='M20-masked'){net<-net[!target %in% m20];mm<-mm[!rownames(mm)%in%m20,,drop=FALSE]}
 cov<-net[,.(n_prior_targets=uniqueN(target),n_covered_targets=uniqueN(target[target %in% rownames(mm)]),n_M20_overlap=uniqueN(target[target %in% m20])),by=source];cov[,`:=`(family=family,resource=resource,mask=mask,min_required=if(resource=='PROGENy')10 else 5)];cov[,status:=ifelse(n_covered_targets>=min_required,'computed','insufficient_target_coverage')];coverage[[paste(family,resource,mask)]]<-cov
 cat('Activity',family,resource,mask,nrow(mm),ncol(mm),'\n')
 ac<-as.data.table(run_ulm(mm,net,.source='source',.target='target',.mor='mor',minsize=if(resource=='PROGENy')10 else 5));setnames(ac,c('condition','p_value'),c('profile_id','ulm_enrichment_p'));ac[,ulm_enrichment_q:=p.adjust(ulm_enrichment_p,'BH'),by=profile_id];ac[,`:=`(family=family,resource=resource,mask=mask)];activities[[paste(family,resource,mask)]]<-ac
 aw<-dcast(ac,source~profile_id,value.var='score');am<-as.matrix(aw[,-1]);rownames(am)<-aw$source
 subplan<-if(family=='Epi')plan[kind=='Epi']else plan[kind=='receiver'&target==family]
 for(i in seq_len(nrow(subplan))){p<-subplan[i];s<-select_model(p,pr);if(is.null(s))next
 am2<-am[,s$profiles$profile_id,drop=FALSE];valid<-apply(am2,1,function(x)all(is.finite(x))&&sd(x)>1e-10);if(!any(valid))next
 ft<-eBayes(lmFit(am2[valid,,drop=FALSE],s$design),trend=FALSE,robust=TRUE);tt<-as.data.table(topTable(ft,coef=ncol(s$design),number=Inf,sort.by='none'),keep.rownames='source');setnames(tt,c('logFC','P.Value','adj.P.Val'),c('activity_effect','sample_model_p','sample_model_q'));tt[,`:=`(model_id=p$model_id,family=family,resource=resource,mask=mask,n_specimens=uniqueN(s$profiles$sample_uid),n_recorded_patient_units=uniqueN(s$profiles$patient_uid_recorded),n_datasets=uniqueN(s$profiles$dataset))];results[[paste(p$model_id,resource,mask)]]<-tt
 }
 put(rbindlist(results,fill=TRUE),'activity_associations.tsv.gz');put(rbindlist(coverage,fill=TRUE),'activity_network_coverage.tsv')
 }
}
a<-rbindlist(activities,fill=TRUE);put(a,'profile_regulatory_activities.tsv.gz');r<-rbindlist(results,fill=TRUE);r[,joint_main_q:=NA_real_];for(resource in unique(r$resource))for(k in c('Epi','receiver')){ids<-plan[primary&kind==k,model_id];ix<-which(r$resource==resource&r$mask=='Full'&r$model_id %in% ids);r[ix,joint_main_q:=p.adjust(sample_model_p,'BH')]};put(r,'activity_associations.tsv.gz');cat('ACTIVITIES COMPLETE\n')
