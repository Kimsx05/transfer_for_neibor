source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
source(file.path(O,"code/fit.R"))
growth=gettab("fixed_core_members.tsv")[program=="ATF4_growth",gene]
db=org.Hs.eg.db::org.Hs.eg.db
go=as.data.table(AnnotationDbi::select(db,keys=c("GO:0006418","GO:0003333","GO:0008652","GO:0002474"),keytype="GOALL",columns="SYMBOL"))[,.(GOALL,SYMBOL)];go=unique(go[!is.na(SYMBOL)])
sets=list(Anabolic_functional=intersect(growth,go[GOALL!="GO:0002474",SYMBOL]),MHC_I=unique(go[GOALL=="GO:0002474",SYMBOL]),JLR_selected=c("IL6","CXCL1","CXCL8"),PNAS_CXCL_locus=c("CXCL1","CXCL2","CXCL6","CXCL8"))
members=rbindlist(lapply(names(sets),function(p)data.table(program=p,gene=sets[[p]])));put(members,"secondary_prediction_members.tsv");put(go[SYMBOL%in%members$gene],"secondary_GO_membership.tsv")
q=gettab("project_sample_map.tsv");allsc=gettab("project_scores.tsv.gz");mat=readRDS(file.path(O,"tables/rank_backgrounds.rds"));res=list();covs=list();samples=list()
for(ty in c("Capillary EC","Total EC")){
 key=paste0(ty,"_50");m=mat[[key]];meta=q[type==ty & n_cells>=50 & main_patient==TRUE]
 for(n in names(sets)){
  for(mask in c("Full","M20",if(n=="Anabolic_functional")c("noCycle","noShared"))){
  core=gettab("fixed_core_members.tsv");drop=switch(mask,Full=character(),M20=M20,noCycle=core[program%in%c("G2M_CHECKPOINT","E2F_TARGETS"),gene],noShared=intersect(core[program=="ATF4_growth",gene],core[program=="ATF4_ISR",gene]))
  genes=setdiff(sets[[n]],drop);bg=setdiff(rownames(m),drop);ms=m[bg,meta$sample_uid,drop=FALSE];if(mask!="Full")ms=apply(ms,2,function(v)(rank(v)-.5)/length(v))
  hit=intersect(genes,bg);k=paste(ty,n,mask)
  covs[[k]]=data.table(type=ty,program=n,mask,expected=length(genes),covered=length(hit),members=paste(hit,collapse=";"))
  z=if(length(hit)&&length(hit)>=length(genes)/2)colMeans(ms[hit,,drop=FALSE]) else rep(NA_real_,ncol(ms))
  samples[[k]]=data.table(type=ty,program=n,mask,sample_uid=meta$sample_uid,value=z)
  for(x in c("Q50","Q95",if(n=="Anabolic_functional")"cycle")){
   co=NULL;if(x=="cycle"){ss=allsc[context==key&mask=="Full"&program%in%c("E2F_TARGETS","G2M_CHECKPOINT"),.(val=mean(value)),by=sample_uid];co=ss$val[match(meta$sample_uid,ss$sample_uid)]}
   res[[paste(k,x)]]=cbind(covs[[k]],data.table(setting=x),fitone(if(x=="Q95")meta$Q95_z else meta$Q50_z,z,meta$dataset,co))
  }
 }}
}
st=rbindlist(res);st[,q:=p.adjust(P,"BH",n=.N),by=.(type,mask,setting)];put(st,"secondary_prediction_effects.tsv");put(rbindlist(samples),"secondary_scores.tsv");put(rbindlist(covs),"secondary_coverage.tsv")
# Audit same-target knockdown and external reported compensation pathways.
ext=gettab("external_gene_statistics.tsv.gz");genes=unique(c(unlist(sets),gettab("candidate_M0073.tsv")$gene,"UQCRQ"))
put(ext[gene%in%genes],"external_prediction_genes.tsv.gz")
x=ext[contrast=="GSE121770__QPC_KOminusWT"&!is.na(gene)&is.finite(stat)]
hm=fread(file.path(R,"hallmark.tsv.gz"));set.seed(123);z=as.data.table(fgsea::fgseaMultilevel(split(hm$gene_symbol,hm$gs_name),setNames(x$stat,x$gene),minSize=15,maxSize=500,nPermSimple=10000,nproc=1));z[,leadingEdge:=vapply(leadingEdge,paste,collapse=";",FUN.VALUE="")]
put(z,"QPC_loss_Hallmark.tsv")

