source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R')
a<-rd('analysis_cell_metadata.tsv.gz');q<-rd('specimen_exposure_eligibility.tsv');pr<-rd('pseudobulk_profiles.tsv');m<-rd('pseudobulk_membership.tsv.gz');ct<-rd('full_celltype_counts.tsv')
receivers<-c('Treg','Tfh/CXCL13-like','Arterial EC','Capillary EC');all<-list();res<-list()
for(r in receivers)for(denom in c('lineage','all_TME')){
 major<-if(r %in% c('Treg','Tfh/CXCL13-like'))'T/NK'else'Endothelial'
 d<-copy(q);n<-ct[analysis_celltype_c2l_final==r,.(sample_uid,numerator=n_cells)];dd<-ct[if(denom=='lineage')analysis_major_c2l_final==major else analysis_major_c2l_final!='Epithelial',.(denominator=sum(n_cells)),by=sample_uid]
 d<-merge(d,n,by='sample_uid',all.x=TRUE);d<-merge(d,dd,by='sample_uid',all.x=TRUE);d[is.na(numerator),numerator:=0];d[is.na(denominator),denominator:=0];d[,`:=`(receiver=r,denominator_type=denom,status=ifelse(denominator>0,'computed','no_denominator_cells'))];d[,proportion:=ifelse(denominator>0,numerator/denominator,NA_real_)];all[[paste(r,denom)]]<-d
 for(min_cells in c(50,30))for(patient in c('main','replacement'))for(exp in c('Q50','Q95'))for(omit in c('none',sort(unique(d$dataset)))){
  s<-d[denominator>=min_cells & (if(patient=='main')main_patient else replacement_patient)];if(omit!='none')s<-s[dataset!=omit];s[,x:=get(paste0(exp,'_z'))];s[,ds:=factor(dataset)]
  z<-data.table(receiver=r,denominator_type=denom,min_cells=min_cells,patient=patient,exposure=exp,omit_dataset=omit,n_specimens=nrow(s),n_recorded_patient_units=uniqueN(s$patient_uid_recorded),n_datasets=uniqueN(s$dataset),status='not_evaluable',reason='',estimate=NA_real_,se=NA_real_,p=NA_real_,dispersion=NA_real_)
  des<-if(uniqueN(s$dataset)>1)~ds+x else ~x
  if(nrow(s)<4){z[,reason:='insufficient_coverage']}else {mm<-model.matrix(des,s);if(qr(mm)$rank<ncol(mm)||nrow(mm)-ncol(mm)<3)z[,reason:='insufficient_residual_df_or_rank']else tryCatch({f<-glm(update(des,cbind(numerator,denominator-numerator)~.),data=s,family=quasibinomial());co<-summary(f)$coefficients['x',];z[,`:=`(status='completed',estimate=co[1],se=co[2],p=co[4],dispersion=summary(f)$dispersion)]},error=function(e)z[,reason:=conditionMessage(e)])}
  res[[length(res)+1]]<-z
 }
}
put(rbindlist(all,fill=TRUE),'abundance_counts_proportions.tsv');rr<-rbindlist(res,fill=TRUE);rr[,q:=p.adjust(p,'BH'),by=.(denominator_type,min_cells,patient,exposure,omit_dataset)];put(rr,'abundance_associations.tsv')
# Existing Epi clusters are descriptive annotations, not new malignant labels.
eold<-fread(file.path(P,'mp10_subtype_screen_round1/tables/epithelial_existing_metadata.tsv.gz'),select=c('cell_id','epi_allgene_cluster','annotation'));setnames(eold,'cell_id','barcode');em<-merge(m[profile_type %in% c('H50','H10','H05','L50')],eold,by='barcode',all.x=TRUE)
put(em[,.(n_cells=.N),by=.(sample_uid,profile_type,epi_allgene_cluster,annotation)],'epi_existing_cluster_composition.tsv')
qc<-rd('profile_quality.tsv');z<-merge(pr,qc,by='profile_id');qs<-list()
for(h in c('H50','H10','H05'))for(v in c('median_nCount','median_nFeature','median_mt')){hi<-z[main_patient & profile_type==h & n_cells>=50,.(sample_uid,high=get(v))];lo<-z[main_patient & profile_type=='L50' & n_cells>=50,.(sample_uid,low=get(v))];j<-merge(hi,lo,by='sample_uid');j[,`:=`(high_definition=h,measure=v,difference=high-low)];qs[[paste(h,v)]]<-j};put(rbindlist(qs),'epi_quality_paired.tsv')
cat('ABUNDANCE AND CONFOUNDING COMPLETE\n')
