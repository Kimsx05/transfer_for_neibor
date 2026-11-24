source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
core=gettab("fixed_core_members.tsv");ext=gettab("external_frozen_tails.tsv.gz")
ext[,program:=contrast];s=rbind(core[,.(program,gene,arm="UP",family="core")],ext[,.(program,gene,arm,family="external")])
s[program%in%c("E2F_TARGETS","G2M_CHECKPOINT","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE"),family:="comparators"]
share=intersect(core[program=="ATF4_growth",gene],core[program=="ATF4_ISR",gene]);cyc=core[program%in%c("E2F_TARGETS","G2M_CHECKPOINT"),gene]
q=gettab("project_sample_map.tsv");mats=readRDS(file.path(O,"tables/rank_backgrounds.rds"));scores=list();coverage=list()
for(key in names(mats)){
 mat=mats[[key]];for(mask in c("Full","M20","noShared","noCycle")){
  drop=switch(mask,Full=character(),M20=M20,noShared=share,noCycle=cyc)
  sig=if(mask%in%c("noShared","noCycle"))s[program%in%c("ATF4_growth","ATF4_ISR")] else copy(s)
  sig=sig[!gene%in%drop]; bg=setdiff(rownames(mat),drop);m=mat[bg,,drop=FALSE]
  if(mask!="Full")m=apply(m,2,function(x)(rank(x)-.5)/length(x))
  for(prog in unique(sig$program)){
   sub=sig[program==prog];value=rep(0,ncol(m));pass=TRUE
   for(ar in unique(sub$arm)){members=sub[arm==ar,gene];mapped=intersect(members,bg);id=paste(key,mask,prog,ar)
    coverage[[id]]=data.table(context=key,mask,program=prog,arm=ar,n_total=length(members),n_scored=length(mapped),bg_n=nrow(m),members=paste(mapped,collapse=";"))
    if(length(mapped)<ceiling(length(members)/2)||!length(mapped)){pass=FALSE;next}
    value=value+if(ar=="UP")colMeans(m[mapped,,drop=FALSE]) else -colMeans(m[mapped,,drop=FALSE])
   }
   scores[[paste(key,mask,prog)]]=data.table(context=key,mask,program=prog,family=sub$family[1],sample_uid=colnames(m),value=if(pass)value else NA_real_)
  }
 }
}
sc=rbindlist(scores);cv=rbindlist(coverage);put(cv,"score_coverage.tsv.gz")
sc[,type:=sub("_[0-9]+$","",context)];sc[,min_cells:=as.integer(sub(".*_","",context))]
sc=merge(sc,q[,.(type,sample_uid,n_cells,dataset,Q50_z,Q95_z,main_patient,replacement_patient,patient_uid_recorded,cap_fraction)],by=c("type","sample_uid"))
put(sc,"project_scores.tsv.gz")
source(file.path(O,"code/fit.R"))
model=list();leave=list();i=0L;full=sc[mask=="Full"];epi=tab(V,"sample_level_scores.tsv.gz")[family=="Epi_all"&mask=="M20"&variant=="standard"&signature%in%c("CHOL_SYNTH","HALLMARK_MTORC1_SIGNALING","FA_ELONG")]
for(z in split(sc,by=c("context","mask","program"),keep.by=TRUE)){
 settings=c("Q50",if(z$mask[1]=="Full")c("Q95","swap"));family=z$family[1]
 if(family=="core"&&z$mask[1]=="Full"&&z$min_cells[1]==50)settings=c(settings,"depth","cycle",if(z$type[1]=="Total EC")"composition",unique(epi$signature))
 for(setting in settings){
  d=if(setting=="swap")z[replacement_patient==TRUE] else z[main_patient==TRUE];stopifnot(!anyDuplicated(d$patient_uid_recorded))
  x=if(setting=="Q95")d$Q95_z else d$Q50_z
  if(setting%in%epi$signature){e=epi[signature==setting & n_cells>=50];x=e$score[match(d$sample_uid,e$sample_uid)]}
  cov=NULL
  if(setting=="depth"){
   # Actual library size is summed across available EC subtypes per specimen.
   pr=tab(B,"pseudobulk_profiles.tsv")
   pr=if(d$type[1]=="Total EC")pr[grepl(" EC$",profile_type),.(total_counts=sum(total_counts)),by=sample_uid] else pr[profile_type==d$type[1]]
   cov=log1p(pr$total_counts[match(d$sample_uid,pr$sample_uid)]/d$n_cells)
  }
  if(setting=="composition")cov=d$cap_fraction
  if(setting=="cycle"){cc=full[context==z$context[1]&program%in%c("E2F_TARGETS","G2M_CHECKPOINT"),.(cc=if(all(is.finite(value)))mean(value) else NA_real_),by=sample_uid];cov=cc$cc[match(d$sample_uid,cc$sample_uid)]}
  rr=fitone(x,d$value,d$dataset,cov)
  i=i+1L;model[[i]]=cbind(d[1,.(context,mask,program,family)],data.table(setting),rr)
  if(setting=="Q50"&&z$mask[1]=="Full"&&rr$status=="completed"){
   leave[[i]]=rbindlist(lapply(seq_len(nrow(d)),function(k)cbind(model[[i]][,.(context,mask,program,family,setting)],data.table(omitted=d$sample_uid[k]),fitone(x[-k],d$value[-k],d$dataset[-k]))))
  }
 }
}
res=rbindlist(model);res[,bh_family:=ifelse(setting%in%epi$signature,"Epi_M20_paired",setting)];res[,q:=p.adjust(P,"BH",n=.N),by=.(context,mask,family,bh_family)]
put(res,"patient_effects.tsv.gz");put(rbindlist(leave),"patient_leave_one.tsv.gz")
co=gettab("subtype_composition.tsv")[main_patient==TRUE &total>=50]
comp=co[,fitone(Q50_z,fraction,dataset),by=type];comp[,q:=p.adjust(P,"BH")];put(comp,"composition_effects.tsv")

