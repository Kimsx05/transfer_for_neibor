.libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(limma)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1'
source(file.path(root,'scripts/speckle_getTransformedProps_official.R'))
s<-fread(file.path(root,'tables/samples.tsv'));m<-fread(file.path(root,'tables/cells.tsv.gz'));p<-fread(file.path(root,'tables/proportions.tsv.gz'));s[,MP10_eligible:=tolower(as.character(MP10_eligible))=='true']
# Official speckle transform is applied independently for each actual denominator.
families<-list(list(level='major',den='all_cells',parent='all',cells=m,type='major'),list(level='major',den='TME',parent='all',cells=m[major!='Epithelial'],type='major'),list(level='subtype',den='TME',parent='all',cells=m[major!='Epithelial'],type='subtype'))
for(pa in sort(unique(m[major!='Epithelial',major])))families[[length(families)+1]]<-list(level='subtype',den='parent',parent=pa,cells=m[major==pa],type='subtype')
comparisons<-c('Normal_Tumor','T1_T2','MP10_low_high','MP10_continuous','MP10_stage_complete_unadjusted','MP10_stage_adjusted','T_overall','Normal_Tumor_paired')
results<-list(); kk<-0L; designrows<-list()
for(f in families){
 tr<-getTransformedProps(clusters=f$cells[[f$type]],sample=f$cells$sample_uid,transform='logit')
 for(ds in c('ALL',sort(unique(s$dataset))))for(comp in comparisons){
  d<-copy(s);if(ds!='ALL')d<-d[dataset==ds]
  if(comp%in%c('Normal_Tumor','Normal_Tumor_paired')){d<-d[tissue_status%in%c('Normal','Tumor')];d[,predictor:=factor(tissue_status,levels=c('Normal','Tumor'))]}
  if(comp=='T1_T2'){d<-d[tissue_status=='Tumor' & T_category%in%c('T1','T2')];d[,predictor:=factor(T_category,levels=c('T1','T2'))]}
  if(comp=='MP10_low_high'){d<-d[MP10_eligible & MP10_group%in%c('low','high')];d[,predictor:=factor(MP10_group,levels=c('low','high'))]}
  if(grepl('^MP10_(continuous|stage)',comp)){d<-d[MP10_eligible & is.finite(MP10_sample)];if(comp!='MP10_continuous')d<-d[!is.na(T_category)];d[,predictor:=MP10_sample]}
  if(comp=='T_overall'){d<-d[tissue_status=='Tumor' & !is.na(T_category)];d[,predictor:=factor(T_category)]}
  d<-d[sample_uid%in%colnames(tr$Proportions)]
  if(comp=='Normal_Tumor_paired'){ids<-d[,.(n=uniqueN(tissue_status)),by=patient_uid_recorded][n==2,patient_uid_recorded];d<-d[patient_uid_recorded%in%ids]}
  if(nrow(d))d[,dataset:=droplevels(factor(dataset))]
  types<-rownames(tr$Proportions);y<-tr$TransformedProps[,d$sample_uid,drop=FALSE];raw<-tr$Proportions[,d$sample_uid,drop=FALSE]
  reason<-NA_character_; formula_text<-NA_character_;correlation<-NA_real_;blockmethod<-'independent_recorded_patient_ids';cf<-NULL
  if(nrow(d)<4)reason<-'fewer_than_4_specimens'
  if(!is.na(reason)){}else if(uniqueN(d$predictor)<2)reason<-'one_condition_or_constant_MP10'
  if(is.na(reason)){
   rhs<-c(if(uniqueN(d$dataset)>1)'dataset',if(comp=='MP10_stage_adjusted'&&uniqueN(d$T_category)>1)'factor(T_category)','predictor');formula_text<-paste('~',paste(rhs,collapse='+'));design<-model.matrix(as.formula(formula_text),data=d)
   if(qr(design)$rank<ncol(design))reason<-'dataset_stage_condition_confounded_rank_deficient'
   else if(nrow(design)-ncol(design)<2)reason<-'fewer_than_2_residual_df'
  }
  testable<-if(ncol(raw)>1)apply(raw,1,function(x)all(is.finite(x))&&diff(range(x))>1e-12) else rep(FALSE,nrow(raw))
  if(is.na(reason)&&any(testable)){
   fit<-tryCatch({
    yy<-y[testable,,drop=FALSE]
    if(anyDuplicated(d$patient_uid_recorded)){
     blockmethod<-'duplicateCorrelation_recorded_unverified_patient';dc<-duplicateCorrelation(yy,design,block=d$patient_uid_recorded);correlation<-dc$consensus.correlation
     if(!is.finite(correlation)){
      # Borrow the nuisance within-patient correlation from all atlas specimens
      # in the SAME denominator family, with dataset and tissue adjustment.
      da<-copy(s[sample_uid%in%colnames(tr$Proportions)])
      xa<-model.matrix(~factor(dataset)+factor(tissue_status),data=da)
      ya<-tr$TransformedProps[testable,da$sample_uid,drop=FALSE]
      correlation<-duplicateCorrelation(ya,xa,block=da$patient_uid_recorded)$consensus.correlation
      blockmethod<-'duplicateCorrelation_borrowed_full_atlas_same_denominator'
      if(!is.finite(correlation))stop('patient_correlation_not_estimable')
     }
     fit<-lmFit(yy,design,block=d$patient_uid_recorded,correlation=correlation)
    }else fit<-lmFit(yy,design)
    eBayes(fit,trend=FALSE,robust=nrow(yy)>=3)
   },error=function(e)e)
   if(inherits(fit,'error'))reason<-paste('fit_error',conditionMessage(fit)) else {
    idx<-grep('^predictor',colnames(design));cf<-topTable(fit,coef=idx,number=Inf,sort.by='none',confint=TRUE);cf<-as.data.table(cf,keep.rownames='subtype')
   }
  }
  z<-data.table(subtype=types,comparison=comp,scope=ds,level=f$level,denominator_type=f$den,parent_family=f$parent,n_specimens=nrow(d),n_patient_ids_recorded=uniqueN(d$patient_uid_recorded),n_patients_verified=0,n_datasets=uniqueN(d$dataset),model=formula_text,patient_method=blockmethod,patient_correlation=correlation,status=ifelse(testable,ifelse(is.na(reason),'estimated',reason),'constant_or_unavailable_composition'))
  z[,`:=`(coefficient=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,group0_mean=NA_real_,group1_mean=NA_real_,raw_difference_pp=NA_real_,raw_slope=NA_real_)]
  if(!is.null(cf)){
   ii<-match(cf$subtype,z$subtype);z[ii,P:=cf$P.Value]
   if('logFC'%in%names(cf))z[ii,`:=`(coefficient=cf$logFC,CI_low=cf$CI.L,CI_high=cf$CI.R)]
  }
  if(is.factor(d$predictor)&&nlevels(d$predictor)==2){
   a<-rowMeans(raw[,d$predictor==levels(d$predictor)[1],drop=FALSE]);b<-rowMeans(raw[,d$predictor==levels(d$predictor)[2],drop=FALSE]);z[,`:=`(group0_mean=a,group1_mean=b,raw_difference_pp=100*(b-a))]
  }else if(is.numeric(d$predictor)&&uniqueN(d$predictor)>1){z[,raw_slope:=apply(raw,1,function(v)coef(lm(v~d$predictor))[2])]}
  z[,effect_scale:=ifelse(comp=='T_overall','moderated_F_omnibus',ifelse(grepl('^MP10_(continuous|stage)',comp),'natural_logit_per_1_original_MP10','natural_logit_group1_minus_group0'))]
  kk<-kk+1L;results[[kk]]<-z
 }
}
r<-rbindlist(results,fill=TRUE);r[,q:=p.adjust(P,method='BH'),by=.(comparison,scope,level,denominator_type)]
map<-unique(p[,.(subtype,major)]);r<-merge(r,map,by='subtype',all.x=TRUE,sort=FALSE)
fwrite(r,file.path(root,'tables/proportion_models.tsv'),sep='\t',na='NA')
cat('Model statuses\n');print(r[,.N,by=.(comparison,scope,status)])
cat('PRIMARY MP10\n');print(r[scope=='ALL' & comparison=='MP10_continuous' & ((level=='major'&denominator_type=='TME')|(level=='subtype'&denominator_type=='parent'))][order(P),.(subtype,coefficient,CI_low,CI_high,P,q,n_specimens,status)])
