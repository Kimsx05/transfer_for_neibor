source(file.path(Sys.getenv("OUT"),"code/io.R"))
suppressPackageStartupMessages(library(Seurat))
options(future.globals.maxSize=8*1024^3);future::plan("sequential")
rd <- file.path(out,"shared/reference"); ref <- sparse_input(rd); meta <- fread(file.path(rd,"membership.tsv.gz"),na.strings="NA")
stopifnot(identical(meta$barcode,colnames(ref)))
mk <- function(x) NormalizeData(CreateSeuratObject(x,min.cells=0,min.features=0),verbose=FALSE)
# Transfer on a reference PCA newly fitted after panel restriction. Neutral codes preserve original subtype punctuation.
transfer <- function(rx,qx,y,do_cluster=FALSE) {
  genes <- intersect(rownames(rx),rownames(qx));rx<-rx[genes,,drop=FALSE];qx<-qx[genes,,drop=FALSE]
  ids<-colnames(qx); if(length(unique(y))==1) return(matrix(1,nrow=length(ids),dimnames=list(ids,unique(y))))
  stopifnot(length(y)==ncol(rx))
  rr<-FindVariableFeatures(mk(rx),nfeatures=min(1500,length(genes)),verbose=FALSE)
  feats <- VariableFeatures(rr); feats<-feats[rowSums(qx[feats,,drop=FALSE])>0]
  dims <- seq_len(min(20,length(feats)-1,ncol(rx)-2,ncol(qx)-2))
  if(length(dims)<2)stop("Too few cells or shared features for anchors")
  rr<-ScaleData(rr,features=feats,verbose=FALSE);rr<-RunPCA(rr,features=feats,npcs=max(dims),verbose=FALSE,seed.use=123)
  qq<-mk(qx)
  aa<-FindTransferAnchors(reference=rr,query=qq,reference.reduction="pca",reduction="pcaproject",features=feats,dims=dims,k.anchor=min(5,max(1,floor(ncol(qx)/4))),k.score=min(30,max(2,floor(min(ncol(qx),ncol(rx))/3))),verbose=FALSE)
  labels<-sort(unique(y)); enc<-paste0("type",match(y,labels)); names(enc)<-colnames(rx)
  res<-TransferData(anchorset=aa,refdata=enc,dims=dims,k.weight=min(40,nrow(aa@anchors)-1L),verbose=FALSE)
  col<-paste0("prediction.score.type",seq_along(labels)); mat<-as.matrix(res[ids,col,drop=FALSE]);colnames(mat)<-labels
  mat
}
top2 <- function(w) {
  first <- max.col(w,ties.method="first");v<-w[cbind(seq_len(nrow(w)),first)]
  w[cbind(seq_len(nrow(w)),first)] <- -Inf
  second<-max.col(w,ties.method="first")
  data.frame(first=colnames(w)[first],support=v,second=colnames(w)[second],runner=w[cbind(seq_len(nrow(w)),second)])
}
tax <- unique(fread(file.path(out,"shared/taxonomy.tsv"))[,.(label,major,markers,anchors)])
major_markers<-jsonlite::read_json(file.path(out,"shared/major_markers.json"),simplifyVector=TRUE)
marker_count<-function(x,markers) as.numeric(colSums(x[intersect(rownames(x),markers),,drop=FALSE]>0))

