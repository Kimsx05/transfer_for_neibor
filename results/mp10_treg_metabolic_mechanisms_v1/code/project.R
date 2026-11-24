source(file.path(Sys.getenv('B_OUT'),'code/common.R'))
profiles<-readt(out,'profiles.tsv')
mat<-readRDS(file.path(out,'tables/treg_logCPM.rds'))
stopifnot(identical(colnames(mat),profiles$profile_id))
primary<-profiles$main_patient&profiles$n_cells>=50
ext<-readt(out,'perturbation_DE.tsv.gz')
sets<-readt(out,'uptake_gene_coverage.tsv')
hg2m<-readt(v2,'hallmark_gene_sets.tsv.gz')[signature=='HALLMARK_G2M_CHECKPOINT',unique(gene)]
members<-readt(out,'original_receiver_membership.tsv.gz')
masks<-list(Full=character(),M20=top,Downstream=unique(c(top,'SCAP','SREBF1','SREBF2','FASN','CD36')))
keys<-c('SCAP','FASN','TCR','CD36','SREBF2')
scores<-list();cover<-list();enrich<-list();modelrows<-list();effects<-list()
std<-function(v)as.numeric((v-mean(v[primary]))/sd(v[primary]))
percentile<-function(x)(rank(x,ties.method='average')-.5)/length(x)
for(ma in names(masks)){
 genes<-setdiff(rownames(mat),masks[[ma]])
 rnk<-apply(mat[genes,,drop=FALSE],2,percentile);rownames(rnk)<-genes
 for(k in c(keys,'G2M',unique(sets$signature))){
   if(k%in%keys){
     a<-ext[contrast==k];u<-a[human%in%genes];w<-percentile(u$oriented_t)-.5
     v<-as.numeric(crossprod(w,rnk[u$human,,drop=FALSE]))/sum(abs(w))
     n_input<-nrow(a);n_mapped<-sum(!is.na(a$human));n_det<-sum(a$human%in%rownames(mat));n_used<-nrow(u)
   }else{
     expected<-if(k=='G2M')hg2m else sets[signature==k,unique(gene)]
     g<-intersect(genes,expected);n_input<-length(expected);n_mapped<-n_input;n_det<-sum(expected%in%rownames(mat));n_used<-length(g)
     v<-if(length(g)>=5 && length(g)>=.5*length(setdiff(expected,masks[[ma]])))colMeans(rnk[g,,drop=FALSE]) else rep(NA_real_,ncol(rnk))
   }
   cover[[length(cover)+1]]<-data.table(mask=ma,program=k,input=n_input,mapped=n_mapped,detected_in_CPM_background=n_det,used=n_used,rank_universe=length(genes),scale_mean=mean(v[primary]),scale_SD=sd(v[primary]))
   scores[[length(scores)+1]]<-cbind(profiles[,.(profile_id,sample_uid,dataset,patient_uid_recorded,n_cells,main_patient,replacement_patient,Q50_z,Q95_z)],data.table(mask=ma,program=k,raw=as.numeric(v),score=std(v)))
 }
 for(mid in unique(members$model_id)){
   de<-readt(old,paste0('de/',mid,'.tsv.gz'))[!gene%in%masks[[ma]]]
   for(k in keys){
     a<-ext[contrast==k];u<-merge(a[!is.na(human)],de[,.(human=gene,signed_stat)],by='human')
     for(set in c('full_rank','strict_positive','strict_negative')){
        sel<-switch(set,full_rank=rep(TRUE,nrow(u)),strict_positive=u$adj.P.Val<.05 & u$oriented_t>0,strict_negative=u$adj.P.Val<.05 & u$oriented_t<0)
        stat<-if(set=='full_rank')cor(u$oriented_t,u$signed_stat,method='spearman') else if(sum(sel)>=1 && sum(!sel)>=1)(mean(rank(u$signed_stat)[sel])-(sum(sel)+1)/2)/sum(!sel)-.5 else NA_real_
        external<-if(set=='full_rank')rep(TRUE,nrow(a)) else a$adj.P.Val<.05 & (if(set=='strict_positive')a$oriented_t>0 else a$oriented_t<0)
        enrich[[length(enrich)+1]]<-data.table(model_id=mid,mask=ma,program=k,set=set,external_input=sum(external),mapped=sum(external&!is.na(a$human)),shared_background=nrow(u),used=sum(sel),stat=stat,stat_type=if(set=='full_rank')'Spearman' else 'AUC_minus_0.5',recorded_units=uniqueN(members[model_id==mid,patient_uid_recorded]))
     }
   }
 }
}
sco<-rbindlist(scores);saveT(sco,'specimen_scores.tsv.gz');saveT(rbindlist(cover),'program_coverage.tsv');saveT(rbindlist(enrich),'directional_enrichment.tsv')
fit<-function(d,adjust=FALSE,hc3=FALSE){
 d<-d[is.finite(x)&is.finite(y)];if(adjust)d<-d[is.finite(z)]
 stopifnot(!anyDuplicated(d$patient_uid_recorded));n<-nrow(d)
 ans<-data.table(n_specimens=n,n_recorded_units=uniqueN(d$patient_uid_recorded),n_datasets=uniqueN(d$dataset),df=NA_integer_,beta=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,reason='insufficient_units_or_df')
 if(n<4)return(ans)
 m<-lm(reformulate(c(if(uniqueN(d$dataset)>1)'dataset','x',if(adjust)'z'),'y'),d)
 ans[,df:=df.residual(m)];if(m$rank<ncol(model.matrix(m))||ans$df<3)return(ans)
 se<-summary(m)$coef['x',2]
 if(hc3){h<-hatvalues(m);if(max(h)>=1-1e-10){ans[,reason:='HC3_undefined_unit_leverage'];return(ans)}
 X<-model.matrix(m);V<-solve(crossprod(X));se<-sqrt((V%*%crossprod(X,X*(resid(m)/(1-h))^2)%*%V)['x','x'])}
 bt<-unname(coef(m)['x']);df0<-df.residual(m);ans[,c('beta','lo','hi','p','reason'):=list(bt,bt-qt(.975,df0)*se,bt+qt(.975,df0)*se,2*pt(-abs(bt/se),df0),'estimated')];ans
}
for(ma in names(masks)){
 w<-dcast(sco[mask==ma],profile_id~program,value.var='score')
 for(setting in c(unique(members$model_id),paste0('within_',unique(profiles$dataset)))){
   mid<-if(startsWith(setting,'within_'))'M0046' else setting
   d<-merge(members[model_id==mid,.(profile_id,sample_uid,patient_uid_recorded,dataset,x=exposure_z)],w,by='profile_id')
   if(startsWith(setting,'within_'))d<-d[dataset==sub('within_','',setting)]
   for(k in c(keys,'CHOL_IMPORT','LCFA_IMPORT','CHOL_SYNTH','G2M')){
     fam<-if(k%in%keys[1:3])'core' else if(k%in%keys[4:5])'supplement' else 'process'
     d[,y:=get(k)]
     for(adj in c('total',if(k!='TCR' && k%in%keys)'TCR',if(k!='G2M'&& k%in%keys)'G2M')){
       if(adj!='total')d[,z:=get(adj)]
       modelrows[[length(modelrows)+1]]<-cbind(data.table(setting=setting,mask=ma,program=k,adjustment=adj,family=fam),fit(d,adj!='total'))
     }
   }
 }
}
mod<-rbindlist(modelrows);mod[,q:=p.adjust(p,method='BH',n=if(family[1]=='core')3L else if(family[1]=='supplement')2L else 4L),by=.(setting,mask,adjustment,family)]
mod[,status:=fcase(!is.finite(p),'不可评估',q<.05 & beta>0,'完成且支持',default='完成但不支持')];saveT(mod,'model_status.tsv')
# Compare like-for-like H05 sender with repaired receiver. HC3 primary and registered multi-specimen sensitivity.
sender<-readt(v2,'sample_level_scores.tsv.gz')[family=='Epi_all'&profile_type=='H05'&variant=='standard'&mask=='Full'&main_patient==TRUE&n_cells>=50]
mm<-list();sm<-list()
for(k in c('CHOL_IMPORT','LCFA_IMPORT')){
 sx<-if(k=='CHOL_IMPORT')'CHOL_EFFLUX' else 'FA_ACYL_COA'
 d<-merge(sender[signature==sx,.(sample_uid,x=score,n_sender=n_cells)],sco[program==k&mask=='Full'&main_patient==TRUE&n_cells>=50,.(sample_uid,patient_uid_recorded,dataset,y=score,n_receiver=n_cells)],by='sample_uid')
 d[,c('sender','receiver'):=list(sx,k)];sm[[k]]<-d
 for(st in c('main','multi_specimen_datasets')){
   z<-copy(d);if(st!='main')z<-z[dataset%in%z[,.N,by=dataset][N>=2,dataset]]
   mm[[length(mm)+1]]<-cbind(data.table(setting=st,sender=sx,receiver=k),fit(z,hc3=TRUE))
 }
}
routes<-rbindlist(mm);routes[,q:=p.adjust(p,'BH',n=2),by=setting];saveT(routes,'recovered_sender_receiver.tsv');saveT(rbindlist(sm),'sender_receiver_plot_data.tsv')
writeLines(capture.output(sessionInfo()),file.path(out,'logs/projection_session.txt'))
print(mod[setting=='M0046'&mask=='Full'&adjustment=='total']);print(routes)
