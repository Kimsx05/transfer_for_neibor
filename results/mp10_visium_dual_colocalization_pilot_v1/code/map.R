source("code/init.R")
tmp <- ".work"
seed <- 20261126L
meta <- fread("reference_manifest.tsv.gz")
sections <- fread("tables/join_audit.tsv")$section_id
top <- scan("../mp10_visium_h05_mapping_pilot_v1/features/top20.txt", what="",quiet=TRUE)
sct <- function(id) {
 x <- read_matrix(id)
 x <- x[Matrix::rowSums(x)>0,,drop=FALSE]
 obj <- CreateSeuratObject(x)
 message("SCT ",id," ",nrow(x)," genes x ",ncol(x)," observations ",Sys.time())
 SCTransform(obj,assay="RNA",vst.flavor="v2",ncells=min(1000,ncol(x)),variable.features.n=1500,
             conserve.memory=TRUE,do.correct.umi=FALSE,return.only.var.genes=TRUE,seed.use=seed,verbose=FALSE)
}
for (mode in c("Epi","TME")) {
 ix <- which(meta$model==mode)
 stopifnot(identical(meta$barcode[ix],readLines(file.path(tmp,mode,"ids"))))
}
if (Sys.getenv("DUAL_RESUME")=="1") {
 complete <- TRUE
 for (s in sections) for (mode in c("Epi","TME")) {
  dest <- file.path("spots",paste0(s,".",mode,".tsv.gz"))
  fp <- file.path("features",paste0(s,".",mode,".tsv"))
  if(!file.exists(dest)||!file.exists(fp)){complete <- FALSE;next}
  d <- fread(dest);fx <- fread(fp); ax <- readLines(file.path(tmp,s,"ids"))
  sc <- as.matrix(d[,-c("section_id","barcode")])
  stopifnot(all(d$section_id==s),identical(d$barcode,ax),ncol(sc)==if(mode=="Epi")2L else 29L,
            all(is.finite(sc)),min(sc)>=0,max(abs(rowSums(sc)-1))<1e-7,nrow(fx)>30,!anyDuplicated(fx$gene))
  message("verified ",s," ",mode," scores and feature axes")
 }
 if(complete){writeLines(capture.output(sessionInfo()),"logs/cache_verification_R.txt");quit(save="no",status=0)}
}
reference <- list()
for(model in c("Epi","TME")) reference[[model]] <- sct(model)
for(s in sections) {
 if(Sys.getenv("DUAL_RESUME")=="1" && all(file.exists(file.path("spots",paste0(s,".",c("Epi","TME"),".tsv.gz"))))) next
 query <- sct(s)
 for(model in names(reference)){
  dest <- file.path("spots",paste0(s,".",model,".tsv.gz"))
  if(Sys.getenv("DUAL_RESUME")=="1" && file.exists(dest) && file.exists(file.path("features",paste0(s,".",model,".tsv"))))next
  ref <- reference[[model]]
  features <- intersect(VariableFeatures(ref),rownames(query[["SCT"]]))
  features <- features[!startsWith(features,"PilotUnusedFeature")]
  ref <- RunPCA(ref,features=features,npcs=30,seed.use=seed,verbose=FALSE)
  actual <- rownames(Loadings(ref[["pca"]]))
  set.seed(seed)
  anc <- FindTransferAnchors(reference=ref,query=query,normalization.method="SCT",reference.assay="SCT",query.assay="SCT",
            reference.reduction="pca",reduction="pcaproject",features=actual,dims=1:30,npcs=30,k.anchor=5,k.filter=NA,
            k.score=30,n.trees=50,recompute.residuals=TRUE,verbose=FALSE)
  stopifnot(setequal(actual,anc@anchor.features))
  keep_cells <- which(meta$model==model); m <- meta[keep_cells];stopifnot(identical(m$barcode,colnames(ref)))
  types <- sort(unique(m$label));id <- sprintf("L%02d",seq_along(types))
  pred <- TransferData(anc,refdata=setNames(id[match(m$label,types)],m$barcode),weight.reduction="pcaproject",dims=1:30,k.weight=50,verbose=FALSE)
  stopifnot(setequal(rownames(pred),colnames(query)))
  mat <- as.matrix(pred[colnames(query),paste0("prediction.score.",id)])
  stopifnot(all(is.finite(mat)),min(mat)>=0,max(abs(rowSums(mat)-1))<1e-6)
  colnames(mat) <- if(model=="TME") paste0("TME_score__",types) else ifelse(types=="Epi_H05","H05_separate_transfer_score","nonH05_separate_transfer_score")
  fwrite(data.table(gene=actual,Top20=actual%in%top),file.path("features",paste0(s,".",model,".tsv")),sep="\t")
  fwrite(data.table(section_id=s,barcode=colnames(query),mat),dest,sep="\t")
  message(s," ",model," features=",length(actual)," Top20=",sum(actual%in%top)," anchors=",nrow(anc@anchors)," ",Sys.time())
  rm(anc,ref,pred,mat);gc()
 }
 rm(query);gc()
}
writeLines(capture.output(warnings()),"logs/map_warnings.txt")
writeLines(capture.output(sessionInfo()),"versions_R.txt")

