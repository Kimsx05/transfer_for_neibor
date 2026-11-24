source("mp10_immune_endothelial_mechanisms_v1/scripts/common.R")
suppressPackageStartupMessages({library(Matrix);library(data.table);library(edgeR);library(qs)})
setDTthreads(1)
O=commandArgs(TRUE)[1]; B='results/mp10_immune_endothelial_mechanisms_v1'
rd=function(n) fread(file.path(B,'tables',n))
saveT=function(x,n) fwrite(x,file.path(O,n),sep='\t',na='NA')
modules=fread(file.path(O,'scripts/modules.tsv'));sets=setNames(strsplit(modules$genes,' '),modules$module);top=sets$MP10
obj=qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1)
pr=as.data.table(obj$profiles)[profile_type %in% c('H05','H10','H50','L50','Epi_all')]
saveT(merge(pr,rd('patient_relationship_audit.tsv'),by=c('sample_uid','dataset','patient_uid_recorded')),'tables/project_mapping.tsv')
ac=rd('profile_regulatory_activities.tsv.gz')[family=='Epi' & resource=='PROGENy' & source %in% c('MAPK','PI3K')]
cov=rd('activity_network_coverage.tsv')[family=='Epi' & resource=='PROGENy' & source %in% c('MAPK','PI3K')]
saveT(cov,'tables/PROGENy_coverage.tsv')
w=fread('resources/mp10_regulatory_priors/progeny_top500.tsv')[source %in% c('MAPK','PI3K')]
targets=unique(c(unlist(sets),'SREBF2','SCAP','CHMP4B','CHMP2A'))
selected=intersect(targets,rownames(obj$counts))
node_counts=as.data.table(as.table(as.matrix(obj$detected_cells[selected,pr$profile_id,drop=FALSE])))
setnames(node_counts,c('gene','profile_id','n_detected')); node_counts=merge(node_counts,pr,by='profile_id')
node_counts[,det_fraction:=n_detected/n_cells];saveT(node_counts,'tables/candidate_profile_detection.tsv.gz')
scores=list(); genes=list()
for(state in c('H05','H10','H50','continuous')){
 if(state=='continuous'){p=pr[main_patient==TRUE & eligible_50==TRUE & profile_type=='Epi_all']}else{
 p=pr[main_patient==TRUE & eligible_50==TRUE & profile_type %in% c(state,'L50')];p=p[sample_uid %in% p[,.N,by=sample_uid][N==2,sample_uid]]}
 y=DGEList(as.matrix(obj$counts[,p$profile_id]));y=calcNormFactors(y);mat=cpm(y,log=TRUE,prior.count=2)
 det=as.matrix(obj$detected_cells[,p$profile_id]); gf=filterByExpr(y,group=p$profile_type)
 for(mask in c('Full','M20-masked')){
  for(m in names(sets)){
   goal=sets[[m]];if(mask!='Full')goal=setdiff(goal,top)
   g=intersect(goal,rownames(mat)[rowSums(y$counts)>0])
   sc=if(length(g))colMeans(mat[g,,drop=FALSE])else rep(NA_real_,ncol(mat))
   scores[[length(scores)+1]]=data.table(p,module=m,mask=mask,score=sc,contrast=state,covered=length(g),expected=length(goal))
  }
  mask_it=mask
  for(m in c('MAPK','PI3K')){
   a=ac[source==m & get('mask')==mask_it] # replaced below for data.table scope
   s=a$score[match(p$profile_id,a$profile_id)]
   if(state=='continuous'){
    # Existing Epi_all profiles lack stored ULM. Run same signed univariate model in this cohort.
    bg=which(gf & (mask=='Full' | !rownames(mat)%in%top));mm=mat[bg,,drop=FALSE];mm=mm-rowMeans(mm)
    ww=w[source==m & target %in% rownames(mm)];act=decoupleR::run_ulm(mm,ww,.source='source',.target='target',.mor='mor',minsize=10);s=act$score[match(p$profile_id,act$condition)]
   }
   cv=cov[source==m & get('mask')==mask_it]
   scores[[length(scores)+1]]=data.table(p,module=m,mask=mask,score=s,contrast=state,covered=if(state=='continuous')nrow(ww)else cv$n_covered_targets,expected=cv$n_prior_targets)
  }
 }
 if(state!='continuous'){
  h=p[profile_type==state] # align by sample
  l=p[profile_type=='L50'][match(h$sample_uid,sample_uid)]
  sel=intersect(unique(c(unlist(sets),w$target,'CHMP2A')),rownames(mat))
  dif=mat[sel,h$profile_id,drop=FALSE]-mat[sel,l$profile_id,drop=FALSE]
  info=data.table(gene=sel,contrast=state,delta=rowMeans(dif),det_High=rowMeans(sweep(det[sel,h$profile_id,drop=FALSE],2,h$n_cells,'/')),det_Low=rowMeans(sweep(det[sel,l$profile_id,drop=FALSE],2,l$n_cells,'/')),n=nrow(h))
  genes[[state]]=info
 }
}
ss=rbindlist(scores);saveT(ss,'plot_data/project_profiles.tsv.gz')
gs=rbindlist(genes);saveT(gs,'tables/project_gene_drivers.tsv.gz')
saveT(merge(w,gs,by.x='target',by.y='gene')[,signed_delta:=mor*delta][order(contrast,source,-abs(signed_delta))],'tables/PROGENy_weighted_member_differences.tsv.gz')
# Signed weight times expression contrast describes contributors, not an exact decomposition of ULM t.
estimate=function(f,term,n){
 t=coef(summary(f));if(!term %in% rownames(t))return(NULL)
 ci=confint(f,term);data.table(n=n,effect=t[term,1],lo=ci[1],hi=ci[2],p=t[term,4])
}
res=list();i=0
for(k in unique(ss$contrast))for(v in unique(ss$mask))for(m in unique(ss$module)){
 a=ss[contrast==k & mask==v & module==m & is.finite(score)]
 if(!nrow(a))next
 fit=function(f,term,label,n){
  est=estimate(f,term,n);if(is.null(est))return();i<<-i+1;res[[i]]<<-data.table(contrast=k,module=m,mask=v,model=label,est)}
 if(k=='continuous'){fit(lm(score~Q95_z+factor(dataset),a),'Q95_z','dataset_adjusted',nrow(a));next}
 h=a[profile_type==k];l=a[profile_type=='L50'][match(h$sample_uid,sample_uid)]
 d=data.frame(y=h$score-l$score,ds=h$dataset,depth=log2(h$counts_per_cell/l$counts_per_cell))
 getdelta=function(mod){
  z=ss[contrast==k & mask=='Full' & module==mod]
  z$score[match(h$profile_id,z$profile_id)]-z$score[match(l$profile_id,z$profile_id)]}
 d$cyc=getdelta('cycling');d$stress=getdelta('stress')
 fit(lm(y~1,d),'(Intercept)','paired',nrow(d))
 if(k=='H05'){
  if(!m%in%c('cycling','stress'))fit(lm(y~cyc+stress+depth,d),'(Intercept)','at_zero_confounder_delta',nrow(d))
  for(ds in unique(d$ds)){dd=d[d$ds!=ds,];fit(lm(y~1,dd),'(Intercept)',paste0('omit_',ds),nrow(dd))}
 }
}
rr=rbindlist(res);rr[,q:=p.adjust(p,'BH'),by=.(contrast,mask,model)]
saveT(rr,'tables/project_statistics.tsv')
plan=rd('model_plan.tsv')[kind=='Epi' & primary==TRUE]
dl=rbindlist(lapply(plan$model_id,function(id){a=rd(paste0('de/',id,'.tsv.gz'));a[gene %in% unique(c(unlist(sets),'CHMP2A'))][,model_id:=id]}))
saveT(merge(dl,plan[,.(model_id,target)],by='model_id'),'tables/reused_edgeR_genes.tsv')
capture.output(sessionInfo(),file=file.path(O,'logs/R_session.txt'))

