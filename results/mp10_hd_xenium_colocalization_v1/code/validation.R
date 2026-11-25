source(file.path(Sys.getenv("OUT"),"code/transfer.R"))
d<-file.path(out,"shared/validation");thin<-sparse_input(d);sets<-jsonlite::read_json(file.path(out,"shared/programs.json"),simplifyVector=TRUE)[c("Chol_biosynthesis_8","MP10_covered")]
stopifnot(all(unlist(sets)%in%rownames(ref)))
for(ds in unique(meta$dataset[meta$validate])){
 train<-which(meta$dataset!=ds & meta$state_ref & !is.na(meta$H05_truth)); test<-which(meta$dataset==ds & meta$validate)
 stopifnot(!length(intersect(meta$sample_uid[train],meta$sample_uid[test])),!length(intersect(na.omit(meta$patient_uid_recorded[train]),na.omit(meta$patient_uid_recorded[test]))))
 for(mode in c("panel","thinned")){
  q<-if(mode=="panel") ref[,test] else thin[,meta$barcode[test]]
  a<-transfer(ref[,train],q,meta$H05_truth[train])
  v<-UCell::ScoreSignatures_UCell(q,features=sets,maxRank=min(1500,nrow(q)),ties.method="average",missing_genes="skip",ncores=1,chunk.size=400,name="")
  stopifnot(identical(rownames(v),meta$barcode[test]))
  tab<-data.frame(meta[test,c("barcode","dataset","sample_uid","patient_uid_recorded","MP10_cell","H05_truth")],mode=mode,v,H05_transfer=a[,"H05"],panel_depth=colSums(q))
  fwrite(tab,file.path(d,paste0("test_",ds,"_",mode,".tsv.gz")),sep="\t",na="NA")
 }
}

