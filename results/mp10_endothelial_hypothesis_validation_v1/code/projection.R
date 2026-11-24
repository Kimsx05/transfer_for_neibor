source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
library(fgsea)
de=gettab("external_gene_statistics.tsv.gz")[!is.na(gene)&is.finite(stat)];de[,oriented:=stat*orientation]
ref=gettab("external_frozen_tails.tsv.gz");c1=tab(B,"de/M0073.tsv.gz");c2=gettab("Total EC_50_Q50_DE.tsv.gz")
dom=Reduce(intersect,lapply(split(de$gene,de$contrast),unique));sets=split(ref$gene,paste(ref$contrast,ref$arm,sep="::"))
ranks=list(Capillary=c1,Total=c2);rankout=list();direction=list();bglog=list()
for(ty in names(ranks))for(mask in c("Full","M20")){
 x=ranks[[ty]][is.finite(signed_stat)&(mask=="Full"|!gene%in%M20)];bg=intersect(dom,x$gene);x=x[gene%in%bg];setorder(x,-signed_stat,gene)
 bglog[[paste(ty,mask)]]=x[,.(type=ty,mask,gene,signed_stat)]
 set.seed(123);z=as.data.table(fgseaMultilevel(sets,setNames(x$signed_stat,x$gene),minSize=10,maxSize=500,eps=1e-10,nPermSimple=10000,nproc=1))
 z[,leadingEdge:=vapply(leadingEdge,paste,collapse=";",FUN.VALUE="")];z=merge(data.table(pathway=names(sets)),z,by="pathway",all.x=TRUE)
 z[,c("contrast","arm"):=tstrsplit(pathway,"::",fixed=TRUE)];z[,q:=p.adjust(pval,"BH",n=length(sets))];z[, `:=`(type=ty,mask=mask,n_background=nrow(x),aligned_NES=ifelse(arm=="UP",NES,-NES),status=ifelse(is.finite(NES),"tested","NE size or numerical"))]
 rankout[[paste(ty,mask)]]=z
 direction[[paste(ty,mask)]]=merge(de,x[,.(gene,signed_stat)],by="gene")[,.(n=.N,rho=cor(oriented,signed_stat,method="spearman"),sign_agreement=mean(sign(oriented)==sign(signed_stat)),type=ty,mask=mask),by=contrast]
}
put(rbindlist(rankout),"perturbation_enrichment.tsv");put(rbindlist(direction),"perturbation_direction.tsv");put(rbindlist(bglog),"projection_universe.tsv.gz")
# Human literature sets on the full human DE domain; retain large ISR set.
core=gettab("fixed_core_members.tsv");core=core[program%in%c("ATF4_growth","ATF4_ISR","Fan_Red_ATF4")]
shared=intersect(core[program=="ATF4_growth",gene],core[program=="ATF4_ISR",gene])
cyc=gettab("fixed_core_members.tsv")[program%in%c("E2F_TARGETS","G2M_CHECKPOINT"),gene]
out=list()
for(ty in names(ranks))for(mask in c("Full","M20","noShared","noCycle")){
 drop=switch(mask,Full=character(),M20=M20,noShared=shared,noCycle=cyc)
 d=ranks[[ty]][is.finite(signed_stat)&!gene%in%drop];se=split(core$gene,core$program);set.seed(123)
 z=as.data.table(fgseaMultilevel(se,setNames(d$signed_stat,d$gene),minSize=5,maxSize=1000,nPermSimple=10000,nproc=1));z[,leadingEdge:=vapply(leadingEdge,paste,collapse=";",FUN.VALUE="")]
 z[, `:=`(type=ty,mask=mask,n_background=nrow(d))];out[[paste(ty,mask)]]=z
}
put(rbindlist(out),"literature_rank_enrichment.tsv")

