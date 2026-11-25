source(file.path(Sys.getenv("OUT"),"code/transfer.R"))
s<-commandArgs(TRUE)[1];d<-file.path(out,"xenium",s)
f<-Sys.glob(file.path(root,"highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5",paste0("*_",s,"_*h5")))
x<-Read10X_h5(f,use.names=TRUE,unique.features=FALSE)[["Gene Expression"]]
x<-x[readLines(file.path(out,"shared/common_panel.txt")),,drop=FALSE]
flags<-fread(file.path(d,"flags.tsv.gz"));stopifnot(identical(colnames(x),flags$cell_id))
ri<-which(meta$identity_ref);valid<-colSums(x)>0
mw<-transfer(ref[,ri],x[,valid,drop=FALSE],meta$analysis_major_c2l_final[ri])
saveRDS(mw,file.path(d,"major_weights.rds"),compress=TRUE)
z<-data.frame(id=colnames(x),major_candidate=NA_character_,major_support=NA_real_,major_second=NA_character_,major_runner=NA_real_,major="Unknown",candidate=NA_character_,support=NA_real_,second=NA_character_,runner=NA_real_,supported="Unknown",flag="uncertain",marker_detected=NA_integer_,anchor_detected=NA_integer_)
z[valid,2:5]<-top2(mw)
for(k in unique(z$major_candidate[valid])){
 j<-which(z$major_candidate==k);raw<-marker_count(x[,j,drop=FALSE],major_markers[[k]])
 alt<-integer(length(j)); for(h in unique(z$major_second[j])) {a<-which(z$major_second[j]==h);alt[a]<-marker_count(x[,j[a],drop=FALSE],major_markers[[h]])}
 mix<-z$major_runner[j]>=.25 & alt>=1 & raw>=1 & z$major_support[j]-z$major_runner[j]<.25
 good<-raw>=1 & z$major_support[j]>=.5 & z$major_support[j]-z$major_runner[j]>=.12 & !mix
 z$major[j[good]]<-k; z$supported[j[good]]<-k; z$flag[j[good]]<-"major_only"
 z$major[j[mix]]<-"Mixed";z$supported[j[mix]]<-"Mixed";z$flag[j[mix]]<-"mixed"
}
for(k in c("T/NK","Myeloid","Endothelial","Fibroblast")){
 j<-which(z$major==k);ri<-which(meta$identity_ref & meta$analysis_major_c2l_final==k)
 if(length(j)<8)next
 w<-transfer(ref[,ri],x[,j,drop=FALSE],meta$analysis_celltype_c2l_final[ri]); z[j,7:10]<-top2(w)
 for(lab in unique(z$candidate[j])){
  a<-j[z$candidate[j]==lab];row<-tax[tax$label==lab,]
  n<-marker_count(x[,a,drop=FALSE],strsplit(row$markers," ",fixed=TRUE)[[1]]); anc<-marker_count(x[,a,drop=FALSE],strsplit(row$anchors," ",fixed=TRUE)[[1]])
  z$marker_detected[a]<-n;z$anchor_detected[a]<-anc
  ok<-n>=1 & anc>=1 & z$support[a]>=.50 & z$support[a]-z$runner[a]>=.15
  if(lab=="Treg") ok<-ok & (n>=2 | as.numeric(x["FOXP3",a])>=1)
  z$supported[a[ok]]<-lab;z$flag[a[ok]]<-"subtype_supported"
 }
}
# Treat degenerate boundaries as uncertain. Area outliers remain visible and omit from state analysis.
bad<-flags$boundary_degenerate; z$flag[bad]<-"uncertain_boundary";z$supported[bad]<-"Unknown";z$major[bad]<-"Unknown"
z$epi_eligible<-z$major=="Epithelial" & !flags$area_outlier & !bad
fwrite(z,file.path(d,"identity.tsv.gz"),sep="\t",na="NA")
j<-which(z$epi_eligible);ri<-which(meta$state_ref & !is.na(meta$H05_truth))
w<-transfer(ref[,ri],x[,j],meta$H05_truth[ri]);fwrite(data.frame(id=z$id[j],H05_transfer=w[,"H05"]),file.path(d,"state_transfer.tsv.gz"),sep="\t")
sets<-jsonlite::read_json(file.path(out,"shared/programs.json"),simplifyVector=TRUE)
scores<-UCell::ScoreSignatures_UCell(x[,j],features=sets,maxRank=min(1500,nrow(x)),ties.method="average",missing_genes="skip",ncores=1,chunk.size=400,name="")
fwrite(data.frame(id=rownames(scores),scores),file.path(d,"ucell.tsv.gz"),sep="\t")
# Independent expression clusters on the real panel for marker/label review.
obj<-FindVariableFeatures(mk(x[,valid]),nfeatures=900,verbose=FALSE);obj<-ScaleData(obj,verbose=FALSE);obj<-RunPCA(obj,npcs=15,verbose=FALSE,seed.use=123)
obj<-FindNeighbors(obj,dims=1:15,verbose=FALSE);obj<-FindClusters(obj,resolution=.45,random.seed=123,verbose=FALSE)
fwrite(data.frame(id=colnames(obj),cluster=as.character(Idents(obj))),file.path(d,"expression_clusters.tsv.gz"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"logs",paste0(s,"_R_session.txt")))

