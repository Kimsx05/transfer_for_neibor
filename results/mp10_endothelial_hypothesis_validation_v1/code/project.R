source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
suppressPackageStartupMessages({library(qs);library(fgsea)})
p=qread(file.path(B,"objects/pseudobulk.qs"),nthreads=1);pr=as.data.table(p$profiles)
q=tab(B,"specimen_exposure_eligibility.tsv")
cl=tab(B,"full_celltype_counts.tsv")[analysis_major_c2l_final=="Endothelial",.(sample_uid,type=analysis_celltype_c2l_final,n_cells)]
types=sort(unique(cl$type));grid=merge(CJ(sample_uid=q$sample_uid,type=types),cl,all.x=TRUE);grid[is.na(n_cells),n_cells:=0L];grid[,total:=sum(n_cells),by=sample_uid];grid[,fraction:=n_cells/total]
composition=merge(grid,q[,.(sample_uid,dataset,patient_uid_recorded,main_patient,replacement_patient,Q50_z,Q95_z)],by="sample_uid");put(composition,"subtype_composition.tsv")
counts=rbind(grid[,.(sample_uid,type,n_cells)],unique(grid[,.(sample_uid,type="Total EC",n_cells=total)])); meta=merge(counts,q,by="sample_uid")
meta[,cap_fraction:=grid[type=="Capillary EC",fraction][match(sample_uid,grid[type=="Capillary EC",sample_uid])]];put(meta,"project_sample_map.tsv")
ecpr=pr[profile_type%in%types];ec=p$counts[,ecpr$profile_id,drop=FALSE];tot=ec%*%sparseMatrix(i=seq_len(nrow(ecpr)),j=match(ecpr$sample_uid,q$sample_uid),x=1,dims=c(nrow(ecpr),nrow(q)),dimnames=list(NULL,q$sample_uid))
genes=strsplit("MYC ATF3 ATF4 HSPA5 HSPA9 DDIT3 HERPUD1 XBP1 ASNS PSAT1 SLC7A5 MTHFD2 EIF2AK3 PDK4 PDHA1 ACLY LONP1 STARD10 SREBF2 KLF6 CXCL1 CXCL2 CXCL3 CXCL5 CXCL6 CXCL8 CCL2 CCL5 CXCL10 CXCL11 IL6 ICAM1 VCAM1 HLA-A HLA-B HLA-C B2M TAP1 ISG15 IFIT1 MKI67 FOXO1"," ")[[1]]
hm=fread(file.path(R,"hallmark.tsv.gz"));hsets=split(hm$gene_symbol,hm$gs_name)
fgr=function(stats,sets){set.seed(123);z=as.data.table(fgseaMultilevel(sets,stats,minSize=15,maxSize=500,nPermSimple=10000,nproc=1,eps=1e-10));z[,leadingEdge:=vapply(leadingEdge,paste,collapse=";",FUN.VALUE="")];z}
legacy=tab(B,"gsea_results.tsv.gz")[model_id=="M0073"&library=="Hallmark"];m20=tab(V,"EC_M20_Hallmark_complete.tsv.gz");path=rbind(legacy[,.(pathway,mask,NES,q=padj,pval,size,leadingEdge)],m20[,.(pathway,mask,NES,q,pval,size,leadingEdge)])
put(path,"legacy_full_m20.tsv");path1=legacy[,setNames(strsplit(leadingEdge,";",fixed=TRUE),pathway)]
over=as.data.table(t(combn(names(path1),2)));over[,members:=mapply(function(x,y)paste(intersect(path1[[x]],path1[[y]]),collapse=";"),V1,V2)];over[,n:=ifelse(members=="",0L,lengths(strsplit(members,";",fixed=TRUE)))];put(over,"leading_edge_overlap.tsv")
sens=tab(B,"gsea_sensitivity_results.tsv.gz")[target=="Capillary EC"&library=="Hallmark"];put(sens,"legacy_sensitivity.tsv.gz")
candidate=list();status=list();backgrounds=list(); matrices=list()
for(ty in c("Capillary EC","Total EC","Venous EC","Arterial EC","Lymphatic EC")){
 if(ty=="Total EC")mat=tot else {pp=ecpr[profile_type==ty];mat=ec[,pp$profile_id,drop=FALSE];colnames(mat)=pp$sample_uid}
 for(nmin in c(50,30)){
  main=meta[type==ty & main_patient==TRUE & n_cells>=nmin];ix=main$sample_uid
  for(xx in c("Q50","Q95")){
   ds=factor(main$dataset);x=main[[paste0(xx,"_z")]];D=if(length(x)<2)matrix(1,length(x),1) else if(nlevels(ds)==1)model.matrix(~x) else model.matrix(~ds+x)
   df=length(ix)-qr(D)$rank;id=paste(ty,nmin,xx,sep="_")
   ok=length(ix)>=4&&df>=3&&qr(D)$rank==ncol(D);status[[id]]=data.table(type=ty,min_cells=nmin,exposure=xx,n=length(ix),datasets=uniqueN(ds),df,status=if(ok)"completed" else "NE: insufficient units/rank/df")
   if(!ok)next
   d=DGEList(as.matrix(mat[,ix,drop=FALSE]));keep=filterByExpr(d,design=D);d=calcNormFactors(d[keep,,keep.lib.sizes=FALSE]);key=paste(ty,nmin,sep="_")
   if(xx=="Q50"){
    allix=meta[type==ty & n_cells>=nmin,sample_uid];mm=as.matrix(mat[rownames(d),allix,drop=FALSE]);matrices[[key]]=apply(mm,2,function(v)(rank(v,ties.method="average")-.5)/length(v));rownames(matrices[[key]])=rownames(d)
    backgrounds[[key]]=data.table(context=key,gene=rownames(d));cand=intersect(genes,rownames(mat));cstat=data.table(gene=genes,passed=genes%in%rownames(d),UMI=NA_real_)
    cstat[gene%in%cand,UMI:=as.numeric(Matrix::rowSums(mat[gene,ix,drop=FALSE]))];cstat[,status:=fcase(passed,"tested",is.na(UMI),"absent_feature",UMI==0,"observed_zero",default="nonzero_filtered")];cstat[,context:=key];candidate[[key]]=cstat
   }
   if(ty=="Capillary EC")next
   d=estimateDisp(d,D,robust=TRUE);fit=glmQLFit(d,D,robust=TRUE);tt=as.data.table(topTags(glmQLFTest(fit,coef=ncol(D)),n=Inf)$table,keep.rownames="gene");tt[,signed_stat:=sign(logFC)*sqrt(pmax(F,0))];put(tt,paste0(id,"_DE.tsv.gz"));srt=setNames(tt$signed_stat,tt$gene);g=fgr(srt,hsets);put(g,paste0(id,"_Hallmark.tsv"))
  }
 }
}
put(rbindlist(status),"coverage_model_status.tsv");put(rbindlist(backgrounds),"score_background.tsv.gz");put(rbindlist(candidate),"candidate_detection.tsv")
saveRDS(matrices,file.path(O,"tables/rank_backgrounds.rds"))
md=tab(B,"de/M0073.tsv.gz");put(merge(candidate[["Capillary EC_50"]],md,all.x=TRUE,by="gene"),"candidate_M0073.tsv")
writeLines(capture.output(sessionInfo()),file.path(O,"logs/project_session.txt"))

