source('/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_resistance_validation_v1/scripts/common.R')
source(file.path(B,'scripts/firth.R'))
# Cell-line backgrounds are the independent units.
a=rd('02_cell_lines/line_state_scores.tsv');w=dcast(a,cell_line~state,value.var='MP10_mean');d=melt(w,id.vars=c('cell_line','parental'),measure.vars=c('CisR','GemCisR','GemR'),variable.name='comparison',value.name='resistant_mean');d[,delta:=resistant_mean-parental];wt(d,'02_cell_lines/line_deltas.tsv');sumline=d[,{
 v=delta;boot=replicate(10000,mean(sample(v,length(v),replace=TRUE)));nz=v[v!=0];list(n_lines=.N,n_positive=sum(v>0),n_negative=sum(v<0),n_zero=sum(v==0),mean_delta=mean(v),median_delta=median(v),CI_low=quantile(boot,.025),CI_high=quantile(boot,.975),P=binom.test(sum(nz>0),length(nz),p=.5,alternative='two.sided')$p.value)
},by=comparison];sumline[,FDR_3_comparisons:=p.adjust(P,'BH')];wt(sumline,'02_cell_lines/overall_effects.tsv')
# Clinical metadata/endpoint construction occurs before models and keeps source labels.
g=rd('03_clinical/GSE169455_scored_metadata.tsv');g[,`:=`(cohort='GSE169455',patient_id=title,regimen='unknown_individual_regimen',pretreatment=TRUE,therapy_intent=get('neoadjuvant vs induction chemotherapy'),pathology_raw=get('pathologic tn stage rc specimen'),subtype=get('consensus classifier subtype'),cT_raw=get('clinical tnm-staging'))]
g[,cT:=factor(sub('.*[cC][tT]([0-9]).*','\\1',cT_raw))];g[,ypT:=sub('^[P]?T(.*)N.*$','\\1',toupper(pathology_raw))];g[ypT=='CIS',ypT:='IS'];g[,ypN:=toupper(sub('.*[nN]([0-9xX]+).*','\\1',pathology_raw))];g[,surgery:=grepl('^[pP][tT]',pathology_raw)];en=endpoint(g$ypT,g$ypN,g$surgery);g[,`:=`(good=en$good,pCR=en$pCR,strict_inclusion=FALSE,fallback_inclusion=therapy_intent=='neoadjuvant',author_good=NA_integer_,cycles=NA_real_,inclusion_reason=ifelse(therapy_intent=='induction','excluded_induction','strict_BLOCKED_individual_regimen_unavailable'),endpoint_evidence='GEO postoperative pTN after neoadjuvant/induction; treated as ypTN context')]
h=rd('03_clinical/GSE212810_scored_metadata.tsv');stopifnot(nrow(h)==39,!anyNA(h$patient_author),!anyDuplicated(h$patient_author));h[,`:=`(cohort='GSE212810',patient_id=patient_author,therapy_intent='neoadjuvant',cT=factor(cT),ypT=ifelse(ypT_code==777,'0',ifelse(ypT_code==0,'NONINVASIVE_OR_CIS',as.character(ypT_code))),ypN=ifelse(ypN_code==999,'NX',as.character(ypN_code)),surgery=procedure%in%c(3,4))];en=endpoint(h$ypT,h$ypN,h$surgery);h[,`:=`(good=en$good,pCR=en$pCR,strict_inclusion=regimen%in%c('GP_cisplatin','ddMVAC_cisplatin'),fallback_inclusion=TRUE,inclusion_reason=ifelse(regimen=='GCb_carboplatin','excluded_pure_carboplatin','included_strict_cisplatin'),endpoint_evidence='PMC10591055 Table S6A: 777=no residual tumor; 0=noninvasive/CIS; N999=unknown; no conversion author R to pCR')]
# GP verified as gemcitabine/cisplatin in article; GCb excluded.
all=rbindlist(list(g,h),fill=TRUE);all[,`:=`(nonresponse=1L-good,non_pCR=1L-pCR,author_nonresponse=1L-author_good)];all[,endpoint_reason:=fifelse(is.na(good),'insufficient_ypT_or_ypN_or_no_evaluable_surgery','known_pathology_sufficient')];all[,label_conflict:=!is.na(author_good)&!is.na(good)&author_good!=good]
# One reference scale per cohort, fixed before fitting and unchanged for sensitivity models.
scales=all[,{i=if(cohort[1]=='GSE169455')fallback_inclusion else strict_inclusion;list(n_scale=sum(i),mean=mean(MP10[i]),SD=sd(MP10[i]),scale_population=if(cohort[1]=='GSE169455')'NAC_individual_regimen_unknown_fallback'else 'strict_cisplatin_NAC')},by=cohort];all=merge(all,scales,by='cohort');all[,MP10_z:=(MP10-mean)/SD];stopifnot(all(is.finite(all$MP10_z)),all(scales$SD>0));wt(scales,'03_clinical/score_standardization.tsv');wt(all,'03_clinical/patient_inclusion_endpoint_audit.tsv')
wt(all[,.(n=.N,n_good=sum(good==1,na.rm=TRUE),n_nonresponse=sum(good==0,na.rm=TRUE),n_good_unknown=sum(is.na(good)),n_pCR=sum(pCR==1,na.rm=TRUE),n_non_pCR=sum(pCR==0,na.rm=TRUE),n_pCR_unknown=sum(is.na(pCR))),by=.(cohort,therapy_intent,regimen,strict_inclusion)],'03_clinical/enrollment_counts.tsv')
results=list();cells=list()
fitmodel=function(d,endpoint,scope,adjust='none',primary=FALSE){
 form=as.formula(paste(endpoint,'~MP10_z',if(adjust!='none')paste('+',adjust)else ''));cc=complete.cases(d[,all.vars(form),with=FALSE]);e=d[cc];rec=data.table(cohort=d$cohort[1],scope=scope,endpoint=endpoint,adjustment=adjust,primary=primary,n_enrolled=nrow(d),n=nrow(e),events=sum(e[[endpoint]]),missing=nrow(d)-nrow(e),OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,method='glm_logit_Wald',converged=FALSE,reason='not_estimable')
 if(nrow(e)<5||uniqueN(e[[endpoint]])<2){results[[length(results)+1]]<<-rec;return(invisible(NULL))}
 f=glm(form,data=e,family=binomial());co=summary(f)$coef;se=co['MP10_z','Std. Error'];bad=!f$converged||any(!is.finite(co))||any(abs(co[,'Std. Error'])>10)||any(f$fitted.values<1e-8|f$fitted.values>1-1e-8)
 if(bad){ff=firth_profile(model.matrix(f),e[[endpoint]],match('MP10_z',colnames(model.matrix(f))));rec[,`:=`(OR=exp(ff$beta),CI_low=exp(ff$low),CI_high=exp(ff$high),P=ff$P,converged=TRUE,reason='separation_or_instability_triggered_prespecified_penalty',method='Firth_Jeffreys_profile_likelihood')];f$fitted.values=ff$fitted}else{b=co['MP10_z','Estimate'];rec[,`:=`(OR=exp(b),CI_low=exp(b-1.95996398454005*se),CI_high=exp(b+1.95996398454005*se),P=co['MP10_z','Pr(>|z|)'],converged=TRUE,reason='OK')]}
 results[[length(results)+1]]<<-rec;cells[[length(cells)+1]]<<-data.table(cohort=rec$cohort,scope=scope,endpoint=endpoint,adjustment=adjust,patient_id=e$patient_id,GSM=e$GSM,response=e[[endpoint]],MP10_z=e$MP10_z,fitted_probability=f$fitted.values)
}
for(co in unique(all$cohort)){
 z=all[cohort==co];main=z[if(co=='GSE169455')fallback_inclusion else strict_inclusion];scope=if(co=='GSE169455')'NAC_regimen_unverified_sensitivity'else 'strict_cisplatin_NAC'
 for(ep in c('nonresponse','non_pCR')){
  fitmodel(main,ep,scope,primary=(co=='GSE212810'&&ep=='nonresponse'));fitmodel(main,ep,scope,'cT')
  if(co=='GSE169455')fitmodel(main,ep,scope,'subtype')
  if(co=='GSE212810'){fitmodel(z,ep,'all_NAC_including_GCb_sensitivity');fitmodel(main[cycles>=3],ep,'strict_cisplatin_cycles_ge3_sensitivity')}
 }
 if(co=='GSE212810')fitmodel(z,'author_nonresponse','author_NAC_labels_supplementary')
 if(co=='GSE169455')for(ep in c('nonresponse','non_pCR'))fitmodel(z[therapy_intent=='induction'],ep,'induction_separate_descriptive')
}
r=rbindlist(results,fill=TRUE);r[,FDR:=NA_real_];r[primary==TRUE&is.finite(P),FDR:=p.adjust(P,'BH',n=2)];r[cohort=='GSE212810'&scope=='strict_cisplatin_NAC'&endpoint=='non_pCR'&adjustment=='none'&is.finite(P),FDR:=p.adjust(P,'BH',n=2)]
# Explicit missing primary GSE169455 row, never silently substitute fallback.
r=rbind(r,data.table(cohort='GSE169455',scope='strict_cisplatin_NAC',endpoint='nonresponse',adjustment='none',primary=TRUE,n_enrolled=0,n=0,events=NA_integer_,missing=125,method='BLOCKED',converged=FALSE,reason='individual regimen unavailable; NAC unknown-regimen sensitivity reported separately'),fill=TRUE)
wt(r,'03_clinical/model_results.tsv');wt(rbindlist(cells,fill=TRUE),'03_clinical/model_patient_usage.tsv');print(sumline);print(r[,.(cohort,scope,endpoint,adjustment,n,events,OR,CI_low,CI_high,P,reason)]);cat('MODELS COMPLETE\n')
