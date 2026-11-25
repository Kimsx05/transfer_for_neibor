source('results/mp10_nectin_external_validation_v1/code/common.R')
famil<-unlist(lapply(names(sets)[1:4],function(a)lapply(familygenes,function(b)c(a,b))),recursive=FALSE)
aux<-list(c('NECTIN2','TIGIT'),c('NECTIN2','Treg'),c('MP10','Treg'))
old<-c(list(c('SREBF1','SREBF2')),lapply(c(core,list(c('SREBF1','SREBF2'))),function(v)ifelse(v%in%c('SREBF1','SREBF2'),paste0(v,'_RNA'),v)))
old<-Filter(function(v)any(grepl('RNA',v))||identical(v,c('SREBF1','SREBF2')),old)
pairs<-list(core=core,family=famil,Treg=aux,legacy_RNA=old)
stat<-function(d,ax,by,cov=NULL,spear=FALSE){
 x<-d[[ax]];y<-d[[by]];if(is.null(x)||is.null(y))stop('missing variable ',ax,' ',by)
 keep<-is.finite(x)&is.finite(y);if(!is.null(cov))keep<-keep & !is.na(d[[cov]]) & nzchar(as.character(d[[cov]]))
 x<-x[keep];y<-y[keep];n<-length(x);res<-data.table(n=n,r=NA_real_,partial_r=NA_real_,rho=NA_real_,beta=NA_real_,lo=NA_real_,hi=NA_real_,p=NA_real_,max_loo_delta=NA_real_,most_influential=NA_character_,status='NE: insufficient values or variance')
 if(n<4||sd(x)==0||sd(y)==0)return(res);df<-n-2
 if(!is.null(cov)){z<-d[[cov]][keep];if(uniqueN(z)<2)return(res);mm<-model.matrix(~z);fit<-lm.fit(mm,cbind(x,y));if(n-fit$rank<4)return(res);df<-n-fit$rank-1;x<-fit$residuals[,1];y<-fit$residuals[,2];if(sd(x)<1e-12||sd(y)<1e-12)return(res)}
 if(spear){x<-rank(x);y<-rank(y)}
 rr<-cor(x,y);ci<-tanh(atanh(rr)+c(-1,1)*qnorm(.975)/sqrt(df-1));res[,c('lo','hi','p','status'):=list(ci[1],ci[2],2*pt(-abs(rr)*sqrt(df/(1-rr^2)),df),'estimated')]
 if(spear)res[,rho:=rr] else if(!is.null(cov))res[,partial_r:=rr] else res[,r:=rr]
 if(!spear)res[,beta:=rr*sd(y)/sd(x)]
 x<-x-mean(x);y<-y-mean(y);k<-n/(n-1);loo<-(sum(x*y)-k*x*y)/sqrt((sum(x*x)-k*x*x)*(sum(y*y)-k*y*y));j<-which.max(abs(loo-rr));if(length(j))res[,c('max_loo_delta','most_influential'):=list(abs(loo[j]-rr),d$GSM[keep][j])];res
}
out<-list();sns<-list()
for(a in c('TCGA',cohorts)){
 td<-file.path(TB,'bulk',a);d<-fread(file.path(td,'scores.tsv'));if(a=='GSE31684')d[,subtype_used:=as.character(cluster)]
 for(f in names(pairs))for(pair in pairs[[f]]){
  x<-pair[1];y<-pair[2];base<-data.table(cohort=a,family_id=f,x=x,y=y,analysis_id=paste(a,'Full',x,y,sep='|'))
  est<-stat(d,x,y);out[[length(out)+1]]<-cbind(base,est)
  add<-function(typ,dt=d,ax=x,by=y,cov=NULL,spear=FALSE){v<-stat(dt,ax,by,cov,spear);sns[[length(sns)+1]]<<-cbind(base,model=typ,v)}
  add('raw');add('Spearman',spear=TRUE);if(y!='Treg'){add('T_lineage',cov='Tlineage');if(all(!is.finite(d$Tlineage)))add('CD3D_composition_proxy',cov='CD3D')}
  if(a=='TCGA')add('DNA_purity',cov='purity') else{
   for(g in c('NMIBC','MIBC')){add(g,d[invasion==g]);if(a=='GSE13507')add(paste0(g,'_T_category'),d[T_category_invasion==g])};dt<-d[invasion!='unknown'];dt[,Tstage:=sub('.*?(T[AaSsIi0-4]+).*','\\1',toupper(stage_used))];add('T_stage',dt,cov='Tstage');if(any(!is.na(d$subtype_used)&nzchar(as.character(d$subtype_used))))add('subtype_or_study_cluster',cov='subtype_used')
   if(a=='GSE31684')add('without_preRC_chemo',d[tolower(treatment_used)=='no'])
  }
  mname<-function(n)if(paste0(n,'_M20')%in%names(d))paste0(n,'_M20') else n
  add('M20',ax=mname(x),by=mname(y))
  ap<-paste0(x,'__',y,'_A');bp<-paste0(x,'__',y,'_B');add('remove_shared_members',ax=if(ap%in%names(d))ap else x,by=if(bp%in%names(d))bp else y)
 }
}
z<-rbindlist(out);z[,q_cohort:=p.adjust(p,'BH',n=.N),by=.(cohort,family_id)];z[,reference:=cohort=='TCGA'];z[,q:=p.adjust(p,'BH',n=.N),by=.(reference,family_id)]
ss<-rbindlist(sns);ss[,q:=p.adjust(p,'BH',n=.N),by=.(cohort,family_id,model)]
for(a in unique(z$cohort)){td<-file.path(TB,'bulk',a);for(f in names(pairs))put(z[cohort==a&family_id==f],file.path(td,paste0(f,'.tsv')));put(ss[cohort==a],file.path(td,'sensitivities.tsv'))}
put(z,file.path(TB,'all_cohorts.tsv'));put(ss,file.path(TB,'all_sensitivities.tsv.gz'))
meta<-function(d){
 u<-d[is.finite(r)&n>3];k<-nrow(u);res<-data.table(k=k,included=paste(u$cohort,collapse=';'),excluded=paste(setdiff(cohorts,u$cohort),collapse=';'),r=NA_real_,lo=NA_real_,hi=NA_real_,Q=NA_real_,Q_p=NA_real_,I2_pct=NA_real_,tau2=NA_real_,status='NE: fewer than two estimable cohorts')
 if(k<2)return(res);zz<-atanh(u$r);v<-1/(u$n-3);fn<-function(t){w<-1/(v+t);m<-weighted.mean(zz,w);sum(log(v+t))+log(sum(w))+sum(w*(zz-m)^2)}
 op<-optimize(fn,c(0,max(10,var(zz)*10)));t<-if(fn(0)<=op$objective)0 else op$minimum;w<-1/(v+t);m<-weighted.mean(zz,w);ci<-tanh(m+c(-1,1)*qnorm(.975)/sqrt(sum(w)));qstat<-sum((zz-weighted.mean(zz,1/v))^2/v)
 res[,c('r','lo','hi','Q','Q_p','I2_pct','tau2','status'):=list(tanh(m),ci[1],ci[2],qstat,pchisq(qstat,k-1,lower.tail=FALSE),if(qstat==0)0 else 100*max(0,(qstat-k+1)/qstat),t,'external random effects')];res
}
mt<-z[reference==FALSE,meta(.SD),by=.(family_id,x,y)];put(mt,file.path(TB,'meta_external.tsv'))

