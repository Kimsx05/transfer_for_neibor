source(file.path(Sys.getenv("OUT"),"code/transfer.R"))
states<-which(meta$state_ref & !is.na(meta$H05_truth));parts<-list()
for(s in c("A1R2","A1R4","A1R5","A1R6","A2R2")){
 d<-file.path(out,"xenium",s)
 f<-Sys.glob(file.path(root,"highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5",paste0("*_",s,"_*h5")))
 x<-Read10X_h5(f,use.names=TRUE,unique.features=FALSE)[["Gene Expression"]]
 z<-fread(file.path(d,"state_transfer.tsv.gz")); u<-fread(file.path(d,"ucell.tsv.gz")); stopifnot(identical(z$id,u$id))
 w<-transfer(ref[,states],x[,z$id,drop=FALSE],meta$dataset[states])
 stopifnot(identical(rownames(w),z$id))
 z$chol<-u$Chol_biosynthesis_8;rt<-rank(z$H05_transfer)/nrow(z);rc<-rank(z$chol)/nrow(z)
 z$group<-ifelse(rt>=.8 & rc<=.2,"high_transfer_low_measured",ifelse(rt<=.2 & rc>=.8,"low_transfer_high_measured","remainder"))
 z$leading_ref_source<-colnames(w)[max.col(w,ties.method="first")];z$leading_ref_support<-apply(w,1,max)
 fwrite(data.frame(z,w,check.names=FALSE),file.path(d,"reference_source_diagnostic.tsv.gz"),sep="\t")
 for(ds in colnames(w))for(gr in unique(z$group)){
  j<-z$group==gr;parts[[length(parts)+1]]<-data.frame(sample=s,discordance=gr,reference_source=ds,n=sum(j),mean_source_weight=mean(w[j,ds]),dominant_source_fraction=mean(z$leading_ref_source[j]==ds),H05_support_spearman=if(sd(w[,ds])>0)cor(z$H05_transfer,w[,ds],method="spearman") else NA_real_)
 }
}
fwrite(do.call(rbind,parts),file.path(out,"reports/reference_source_review.tsv"),sep="\t",na="NA")

