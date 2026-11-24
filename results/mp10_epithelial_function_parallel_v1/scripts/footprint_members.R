source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R')
suppressPackageStartupMessages(library(edgeR))
O=commandArgs(TRUE)[1]
b=qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1)
p=b$profiles[n_cells>=30 & profile_type %in% c('L50','H50','H10','H05')]
y=DGEList(as.matrix(b$counts[,p$profile_id]));keep=filterByExpr(y,group=factor(p$profile_type))
y=calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);mat=cpm(y,log=TRUE,prior.count=2);mat=mat-rowMeans(mat)
w=fread(file.path(REF,'progeny_top500.tsv'))
top=strsplit(fread(file.path(O,'scripts/modules.tsv'))[module=='MP10',genes],' ')[[1]]
ac=fread(file.path(B,'tables/profile_regulatory_activities.tsv.gz'))[family=='Epi' & resource=='PROGENy']
h=p[profile_type=='H05' & main_patient==TRUE & n_cells>=50]
l=p[profile_type=='L50' & n_cells>=50][match(h$sample_uid,sample_uid)]
stopifnot(nrow(h)==24, !anyNA(l$profile_id))
ans=list()
for(msk in c('Full','M20-masked')){
 mm=mat[msk=='Full' | !rownames(mat)%in%top,,drop=FALSE]
 mm=sweep(mm,2,colMeans(mm),'-')
 for(path in c('MAPK','PI3K')){
  prior=w[source==path & target%in%rownames(mm)]
  wt=numeric(nrow(mm));names(wt)=rownames(mm);wt[prior$target]=prior$mor
  numerator=mm[prior$target,,drop=FALSE]*prior$mor
  denominator=sqrt(sum((wt-mean(wt))^2)*colSums(mm^2))
  corr=colSums(numerator)/denominator
  contrib=sweep(numerator,2,sqrt((nrow(mm)-2)/(1-corr^2))/denominator,'*')
  ref=ac[source==path & mask==msk];error=max(abs(colSums(contrib)-ref$score[match(colnames(contrib),ref$profile_id)]))
  stopifnot(error<1e-6)
  ans[[length(ans)+1]]=data.table(module=path,mask=msk,gene=prior$target,weight=prior$mor,n=nrow(h),mean_contribution=rowMeans(contrib[,h$profile_id]-contrib[,l$profile_id]),max_score_reconstruction_error=error)
 }
}
fwrite(rbindlist(ans)[order(module,mask,-abs(mean_contribution))],file.path(O,'tables/PROGENy_H05_exact_contributions.tsv'),sep='\t')

