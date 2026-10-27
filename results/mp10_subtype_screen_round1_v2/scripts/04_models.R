.libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(limma)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1_v2';tt<-file.path(root,'tables');s<-fread(file.path(tt,'samples_v2.tsv'));p<-fread(file.path(tt,'proportions_v2.tsv.gz'))
s[,MP10_eligible:=tolower(as.character(MP10_eligible))=='true'];source(file.path(root,'scripts/speckle_getTransformedProps_official.R'))
families<-unique(p[,.(level,denominator_type,parent_family=ifelse(denominator_type=='parent',major,'all'))]);res<-list();membership<-list()
comps<-c('Normal_Tumor','T1_T2','sample_MP10_group_allEpi','MP10_continuous','MP10_stage_complete_unadjusted','MP10_stage_adjusted','T_overall','Normal_Tumor_paired')
makeD<-function(ds,co,eligible){d<-copy(s[sample_uid%in%eligible]);if(ds!='ALL')d<-d[dataset==ds]
 if(co%in%c('Normal_Tumor','Normal_Tumor_paired')){d<-d[tissue_status%in%c('Normal','Tumor')];d[,predictor:=factor(tissue_status,levels=c('Normal','Tumor'))]}
 if(co=='T1_T2'){d<-d[tissue_status=='Tumor' & T_category%in%c('T1','T2')];d[,predictor:=factor(T_category,levels=c('T1','T2'))]}
 if(co=='sample_MP10_group_allEpi'){d<-d[MP10_eligible==TRUE];d[,predictor:=factor(sample_MP10_group_allEpi,levels=c('low','high'))]}
 if(grepl('^MP10_',co)){d<-d[MP10_eligible & is.finite(MP10_sample_allEpi)];if(grepl('stage',co))d<-d[!is.na(T_category)];d[,predictor:=MP10_sample_allEpi]}
 if(co=='T_overall'){d<-d[tissue_status=='Tumor' & !is.na(T_category)];d[,predictor:=factor(T_category)]}
 if(co=='Normal_Tumor_paired'){ids<-d[,.(n=uniqueN(tissue_status)),by=patient_uid_recorded][n==2,patient_uid_recorded];d<-d[patient_uid_recorded%in%ids]}
 setorder(d,sample_uid);d}
mkdesign<-function(d,co){rhs<-c(if(uniqueN(d$dataset)>1 && co!='Normal_Tumor_paired')'factor(dataset)',if(co=='MP10_stage_adjusted'&&uniqueN(d$T_category)>1)'factor(T_category)',if(co=='Normal_Tumor_paired'&&uniqueN(d$patient_uid_recorded)>1)'factor(patient_uid_recorded)','predictor');model.matrix(as.formula(paste('~',paste(rhs,collapse='+'))),data=d)}
reasonD<-function(d,co,design){if(nrow(d)<4)return('fewer_than_4_specimens');if(uniqueN(d$predictor)<2)return('one_condition_or_constant_predictor');if(is.factor(d$predictor)&&any(vapply(levels(d$predictor),function(g)uniqueN(d[predictor==g,patient_uid_recorded])<2,logical(1))))return('fewer_than_2_recorded_patients_per_group');if(is.null(design)||qr(design)$rank<ncol(design))return('rank_deficient');if(nrow(design)-ncol(design)<2)return('fewer_than_2_residual_df');NA_character_}
for(fi in seq_len(nrow(families))){f<-families[fi];z<-p[level==f$level&denominator_type==f$denominator_type];if(f$parent_family!='all')z<-z[major==f$parent_family];types<-sort(unique(z$subtype));wide<-dcast(z,subtype~sample_uid,value.var='proportion');raw<-as.matrix(wide[,-1]);rownames(raw)<-wide$subtype
 # Exact official logit counts: add .5 to every category within each denominator family.
 counts<-dcast(z,subtype~sample_uid,value.var='numerator');cm<-as.matrix(counts[,-1]);rownames(cm)<-counts$subtype;pp<-sweep(cm+.5,2,colSums(cm+.5),'/');logit<-log(pp/(1-pp));asin<-asin(sqrt(raw))
 vars<-if(f$denominator_type=='parent')c('primary_asin','parent30_asin','parent30_logit','one_specimen_asin')else c('primary_asin','one_specimen_asin')
 # A shared fallback decision is made across transformations and both stage models BEFORE fitting.
 for(ds in c('ALL',sort(unique(s$dataset))))for(co in comps){mp<-co%in%c('sample_MP10_group_allEpi','MP10_continuous','MP10_stage_complete_unadjusted','MP10_stage_adjusted');for(variant in vars){if(variant!='primary_asin'&&!mp)next
  eligible<-unique(z[denominator>ifelse(grepl('parent30',variant),29,0),sample_uid]);d<-makeD(ds,co,eligible);method<-'independent_recorded_patient_ids';correl<-NA_real_
  transforms<-if(grepl('parent30',variant))list(asin,logit)else list(asin);pairedco<-if(grepl('stage',co))c('MP10_stage_complete_unadjusted','MP10_stage_adjusted')else co
  forceone<-variant=='one_specimen_asin'
  if(anyDuplicated(d$patient_uid_recorded)&&co!='Normal_Tumor_paired'&&!forceone){
   # Tiny repeat blocks cannot borrow information from the atlas. Require at least two repeated recorded IDs.
   nrep<-sum(table(d$patient_uid_recorded)>1);forceone<-nrep<2
   if(!forceone)for(pc in pairedco)for(ym in transforms){des<-tryCatch(mkdesign(d,pc),error=function(e)NULL);yy<-ym[,d$sample_uid,drop=FALSE];ok<-apply(yy,1,function(v)all(is.finite(v))&&diff(range(v))>1e-12);dc<-if(!is.null(des)&&is.na(reasonD(d,pc,des))&&any(ok))tryCatch(suppressWarnings(duplicateCorrelation(yy[ok,,drop=FALSE],des,block=d$patient_uid_recorded)$consensus.correlation),error=function(e)NA_real_)else NA_real_;if(!is.finite(dc))forceone<-TRUE}
  }
  if(forceone){d<-d[!duplicated(patient_uid_recorded)];method<-if(variant=='one_specimen_asin')'one_specimen_sensitivity' else 'one_specimen_fallback'}
  if(co=='Normal_Tumor_paired')method<-'paired_recorded_patient_fixed_effect'
  design<-tryCatch(mkdesign(d,co),error=function(e)NULL);reason<-reasonD(d,co,design);yy<-if(variant=='parent30_logit')logit[,d$sample_uid,drop=FALSE]else asin[,d$sample_uid,drop=FALSE];rr<-raw[,d$sample_uid,drop=FALSE];ok<-if(ncol(yy)>1)apply(yy,1,function(v)all(is.finite(v))&&diff(range(v))>1e-12)else rep(FALSE,nrow(yy));cf<-NULL
  if(is.na(reason)&&any(ok)){fit<-tryCatch({if(anyDuplicated(d$patient_uid_recorded)&&co!='Normal_Tumor_paired'){correl<-duplicateCorrelation(yy[ok,,drop=FALSE],design,block=d$patient_uid_recorded)$consensus.correlation;method<-'duplicateCorrelation_within_analysis';ff<-lmFit(yy[ok,,drop=FALSE],design,block=d$patient_uid_recorded,correlation=correl)}else ff<-lmFit(yy[ok,,drop=FALSE],design);ff<-eBayes(ff,trend=FALSE,robust=sum(ok)>=3);as.data.table(topTable(ff,coef=grep('^predictor',colnames(design)),number=Inf,sort.by='none',confint=TRUE),keep.rownames='subtype')},error=function(e)e);if(inherits(fit,'error'))reason<-paste('fit_error',conditionMessage(fit))else cf<-fit}
  out<-data.table(subtype=rownames(raw),analysis_variant=variant,comparison=co,scope=ds,annotation_level=f$level,denominator_type=f$denominator_type,parent_family=f$parent_family,n_specimens=nrow(d),n_patient_ids_recorded=uniqueN(d$patient_uid_recorded),n_datasets=uniqueN(d$dataset),patient_method=method,patient_correlation=correl,model=if(is.null(design))NA_character_ else paste(colnames(design),collapse=' + '),status=ifelse(ok,ifelse(is.na(reason),'estimated',reason),'constant_or_unavailable'),coefficient=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,raw_difference_pp=NA_real_)
  if(!is.null(cf)){ii<-match(cf$subtype,out$subtype);out[ii,P:=cf$P.Value];if('logFC'%in%names(cf))out[ii,`:=`(coefficient=cf$logFC,CI_low=cf$CI.L,CI_high=cf$CI.R)]}
  if(is.factor(d$predictor)&&nlevels(d$predictor)==2&&nrow(d))out[,raw_difference_pp:=100*(rowMeans(rr[,d$predictor==levels(d$predictor)[2],drop=FALSE])-rowMeans(rr[,d$predictor==levels(d$predictor)[1],drop=FALSE]))]
  out[,effect_scale:=ifelse(co=='T_overall','moderated_F_no_direction',paste0(ifelse(variant=='parent30_logit','natural_logit','asin_sqrt'),ifelse(is.numeric(d$predictor),'_per_1_original_MP10','_group1_minus_group0')))];res[[length(res)+1]]<-out
  if(nrow(d))membership[[length(membership)+1]]<-data.table(analysis_variant=variant,comparison=co,scope=ds,annotation_level=f$level,denominator_type=f$denominator_type,parent_family=f$parent_family,sample_uid=d$sample_uid,patient_uid_recorded=d$patient_uid_recorded,patient_method=method,fit_status=ifelse(is.na(reason),'design_eligible',reason))
 }}}
r<-rbindlist(res);r[,q:=p.adjust(P,'BH'),by=.(analysis_variant,comparison,scope,annotation_level,denominator_type)];r<-merge(r,unique(p[,.(subtype,major)]),by='subtype',all.x=TRUE,sort=FALSE);fwrite(r,file.path(tt,'proportion_models_v2.tsv'),sep='\t',na='NA');fwrite(rbindlist(membership),file.path(tt,'model_sample_membership.tsv'),sep='\t');cat('DONE',nrow(r),'results\n')
