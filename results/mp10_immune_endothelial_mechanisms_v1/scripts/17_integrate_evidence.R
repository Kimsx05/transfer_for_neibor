source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R');plan<-rd('model_plan.tsv');status<-rd('model_status.tsv');de<-rd('candidate_de_effects.tsv');b<-qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1);pr<-b$profiles;lr<-rd('nichenet_lr_prior.tsv.gz');ligands<-rd('nichenet_ligand_domain.tsv')$ligand;nn<-if(file.exists(tb('nichenet_ligand_activities.tsv.gz')))rd('nichenet_ligand_activities.tsv.gz')else data.table();nns<-if(file.exists(tb('nichenet_analysis_status.tsv')))rd('nichenet_analysis_status.tsv')else data.table();ac<-rd('activity_associations.tsv.gz');all<-list()
# Preserve all LR candidates; key-axis table is a view, not a candidate filter.
for(r in c('Treg','Tfh/CXCL13-like','Arterial EC','Capillary EC')){
 id<-plan[primary&kind=='receiver'&target==r,model_id];st<-status[model_id==id];pp<-pr[main_patient&profile_type==r&n_cells>=50];e<-unique(lr[,.(ligand=from,receptor=to)]);e[,`:=`(sender='Epi',receiver=r,prior_source='NicheNet v2 human / Zenodo7074291',in_ligand_target_prior=ligand%in%ligands,receiver_model_id=id,receiver_model_status=st$status,n_receiver_specimens=nrow(pp),n_receiver_recorded_units=uniqueN(pp$patient_uid_recorded),n_receiver_datasets=uniqueN(pp$dataset),HD_support='pending',functional_support='not_tested')]
 for(recp in unique(e$receptor)){
  if(!recp%in%rownames(b$counts)){e[receptor==recp,`:=`(receptor_status='feature_absent',receptor_median_detection=NA_real_,receptor_mean_RNA=NA_real_)];next}
  if(!nrow(pp)){e[receptor==recp,`:=`(receptor_status='insufficient_sample_coverage',receptor_median_detection=NA_real_,receptor_mean_RNA=NA_real_)];next}
  det<-as.numeric(b$detected_cells[recp,pp$profile_id])/pp$n_cells;expr<-as.numeric(b$normalized_sum[recp,pp$profile_id])/pp$n_cells;e[receptor==recp,`:=`(receptor_status=if(all(det==0))'observed_zero'else if(median(det)<.1)'below_primary_detection_rule'else'detected_primary',receptor_median_detection=median(det),receptor_mean_RNA=mean(expr))]
 }
 for(h in c('H50','H10','H05')){
  did<-plan[primary&kind=='Epi'&target==h,model_id];dt<-rd(paste0('de/',did,'.tsv.gz'))[,.(ligand=gene,logFC,FDR,joint_main_q)];setnames(dt,c('logFC','FDR','joint_main_q'),paste0(h,c('_logFC','_q','_joint_q')));e<-merge(e,dt,by='ligand',all.x=TRUE)
 }
 if(st$status=='completed'){
  dt<-rd(paste0('de/',id,'.tsv.gz'));e[,`:=`(receiver_tested_genes=nrow(dt),receiver_positive_DE=sum(dt$FDR<.05&dt$logFC>0),receiver_negative_DE=sum(dt$FDR<.05&dt$logFC<0))]
 }else e[,`:=`(receiver_tested_genes=NA_integer_,receiver_positive_DE=NA_integer_,receiver_negative_DE=NA_integer_)]
 for(d in c('positive','negative','union')){
  z<-nn[model_id==id&mask=='Full'&target_direction==d&detection_threshold==.1,.(ligand=test_ligand,aupr_corrected,rank_sender_agnostic,rank_focused,low_information)];if(nrow(z)){setnames(z,setdiff(names(z),'ligand'),paste0('NN_',d,'_',setdiff(names(z),'ligand')));e<-merge(e,z,by='ligand',all.x=TRUE)}
 }
 aa<-ac[model_id==id&resource=='CollecTRI'&mask=='Full'];pg<-ac[model_id==id&resource=='PROGENy'&mask=='Full'];e[,`:=`(receiver_CollecTRI_q05=if(nrow(aa))sum(aa$sample_model_q<.05)else NA_integer_,receiver_PROGENy_q05=if(nrow(pg))sum(pg$sample_model_q<.05)else NA_integer_,alternative_explanations='Epi composition/quality; correlated specimen programs; treatment/unmeasured covariates; sparse receiver coverage; prior branch ambiguity',interpretation_limit='RNA and predicted shared targets do not establish receptor occupancy, ligand-to-TF causality, immunosuppression or angiogenesis')]
 all[[r]]<-e
}
all<-rbindlist(all,fill=TRUE);put(all,'complete_candidate_evidence.tsv.gz');key<-all[(ligand%in%c('NECTIN2','PVR')&receptor%in%c('TIGIT','PVRIG','CD226')&receiver%in%c('Treg','Tfh/CXCL13-like'))|(ligand%in%c('VEGFA','VEGFB','VEGFC','PGF','ANGPT1','ANGPT2')&receptor%in%c('FLT1','KDR','FLT4','NRP1','NRP2','TEK')&receiver%in%c('Arterial EC','Capillary EC'))];put(key,'key_axis_evidence.tsv')
# Candidate effect stability retains every prespecified setting, not a best-P score.
z<-de[kind=='Epi'&only_dataset=='ALL'];main<-z[primary==TRUE,.(target,gene,main_logFC=logFC,main_q=FDR)];z<-merge(z,main,by=c('target','gene'),all.x=TRUE);z[,same_direction:=ifelse(is.na(logFC)|is.na(main_logFC),NA,sign(logFC)==sign(main_logFC))];put(z,'candidate_sensitivity_effects.tsv')
# Existing CellChat DB cofactor/complex definitions are context, not new inference.
cc<-fread(file.path(BASE,'tables/cellchat_LR_definitions.tsv'));put(cc[ligand%in%c('VEGFA','VEGFB','VEGFC','PGF','NECTIN2','PVR','TIGIT','CD226')|receptor%in%c('NECTIN2','PVR','TIGIT','PVRIG','CD226','FLT1','KDR','FLT1_KDR')],'candidate_CellChat_DB_context.tsv')
put(data.table(input_symbol=c('NECTIN2','PVRL2'),canonical_symbol='NECTIN2',mapping_type=c('approved_symbol','previous_symbol'),feature_present=c('NECTIN2','PVRL2')%in%rownames(b$counts),prior_canonical_present='NECTIN2'%in%ligands,action=c('use existing canonical feature','do not add zero row or duplicate canonical counts')),'candidate_symbol_mapping.tsv')
cat('EVIDENCE INTEGRATION COMPLETE\n')
