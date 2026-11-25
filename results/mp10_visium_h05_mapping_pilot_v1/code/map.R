.libPaths(unique(c(path.expand("~/R/x86_64-pc-linux-gnu-library/4.5"),path.expand("~/R/x86_64-pc-linux-gnu-library/4.4"),"/refdir/Rlib",.libPaths())))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,RCPP_PARALLEL_NUM_THREADS=1)
if(requireNamespace("RhpcBLASctl",quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
options(Ncpus=1,mc.cores=1)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(Seurat)})
setDTthreads(1)
options(future.globals.maxSize=8*1024^3)
future::plan("sequential")
seed <- 20261125L
set.seed(seed)
tmp <- Sys.getenv("PILOT_TMP")
top <- scan("features/top20.txt",what="",quiet=TRUE)
meta <- fread("reference_sampling_ids.tsv.gz")
sections <- fread("barcode_audit.tsv")$section_id
label_key <- c(EpiHigh="Epi_H05", EpiOther="Epi_nonH05", TNK="T/NK")
label <- meta$label
for (k in names(label_key)) label[label==label_key[k]] <- k
names(label) <- meta$barcode
read_matrix <- function(name) {
 p <- file.path(tmp,name)
 load_int <- function(k) readBin(file.path(p,k),integer(),n=file.size(file.path(p,k))/4,endian="little")
 genes <- readLines(file.path(p,"genes")); ids <- readLines(file.path(p,"ids"))
 # Preserve total UMI without allowing Seurat symbol rewriting to merge identifiers.
 bad <- grepl("[_|]",genes)
 genes[bad] <- paste0("PilotUnusedFeature",which(bad))
 stopifnot(!anyDuplicated(genes),!anyDuplicated(ids))
 new("dgCMatrix",i=load_int("i"),p=load_int("p"),x=as.numeric(load_int("x")),Dim=as.integer(c(length(genes),length(ids))),Dimnames=list(genes,ids))
}
reference <- read_matrix("reference")
stopifnot(identical(colnames(reference),names(label)))
shared <- Reduce(intersect,lapply(sections,function(s)readLines(file.path(tmp,s,"genes"))))
universe <- intersect(shared,rownames(reference)[Matrix::rowSums(reference>0)>=5])
universe <- universe[!grepl("[_|]",universe)]
candidate <- CreateSeuratObject(reference[universe,])
candidate <- FindVariableFeatures(candidate,selection.method="vst",nfeatures=1200,verbose=FALSE)
selected <- unique(c(VariableFeatures(candidate),intersect(top,universe)))
rm(candidate);gc()
build_object <- function(m,exclude) {
 if(length(exclude)) m <- m[!rownames(m)%in%exclude,,drop=FALSE]
 NormalizeData(CreateSeuratObject(m),normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
}
map <- function(s,m,branch,output) {
 set.seed(seed)
 exclude <- if(branch=="noTop20") top else character()
 keep <- setdiff(intersect(selected,rownames(m)[Matrix::rowSums(m)>0]),exclude)
 ref <- build_object(reference,exclude); query <- build_object(m,exclude)
 ref <- ScaleData(ref,features=keep,verbose=FALSE)
 ref <- RunPCA(ref,features=keep,npcs=25,seed.use=seed,verbose=FALSE)
 actual <- rownames(Loadings(ref[["pca"]]))
 anc <- FindTransferAnchors(reference=ref,query=query,normalization.method="LogNormalize",reduction="pcaproject",reference.reduction="pca",features=actual,dims=1:25,npcs=25,k.anchor=5,k.score=30,k.filter=NA,nn.method="annoy",n.trees=50,verbose=FALSE)
 stopifnot(setequal(anc@anchor.features,actual),!any(actual %in% exclude),nrow(anc@anchors)>=50)
 pred <- TransferData(anchorset=anc,refdata=label,weight.reduction="pcaproject",dims=1:25,k.weight=50,verbose=FALSE)
 stopifnot(setequal(rownames(pred),colnames(m)))
 pred <- pred[colnames(m),,drop=FALSE]
 class_names <- sort(unique(label))
 score <- as.matrix(pred[,paste0("prediction.score.",class_names),drop=FALSE]); stopifnot(all(is.finite(score)),all(score>=0),max(abs(rowSums(score)-1))<1e-6)
 decode <- function(x) {for(k in names(label_key)) x[x==k] <- label_key[[k]];x}
 colnames(score) <- paste0("support_",decode(class_names))
 df <- data.table(section_id=s,barcode=colnames(m),predicted_label=decode(as.character(pred$predicted.id)),score)
 df[,H05_support:=support_Epi_H05];df[,Epi_support:=support_Epi_H05+support_Epi_nonH05]
 df[,H05_conditional_support:=fifelse(Epi_support>0,H05_support/Epi_support,NA_real_)];df[,low_denominator:=Epi_support<0.10]
 fwrite(df,output,sep="\t",na="NA")
 fwrite(data.table(gene=actual,MP10_member=(actual %in% top),in_PCA=TRUE,in_anchors=TRUE),file.path("features",paste0(s,"_",branch,".tsv")),sep="\t")
 message(s," ",branch," query=",ncol(m)," features=",length(actual)," anchors=",nrow(anc@anchors)," top20=",sum(actual%in%top))
}
writeLines(capture.output(sessionInfo()),"versions_R.txt")
if (!("--pseudo-only" %in% commandArgs(TRUE))) for (s in sections) {
 m <- read_matrix(s)
 for (branch in c("main","noTop20")) {
  message("Starting ",s," ",branch)
  map(s,m,branch,file.path("spots",paste0(s,".",branch,".tsv.gz")));gc()
 }
 rm(m);gc()
}
# Hold the reference feature set fixed across paired depth realizations.
if (all(dir.exists(file.path(tmp,c("pseudo_medium","pseudo_low"))))) {
 for (depth in c("medium","low")) {
  pm <- read_matrix(paste0("pseudo_",depth))
  selected <- intersect(selected,rownames(pm)[Matrix::rowSums(pm)>0])
 }
 rm(pm)
}
for(depth in c("medium","low")) {
 s <- paste0("pseudo_",depth)
 if(dir.exists(file.path(tmp,s))) tryCatch(map(s,read_matrix(s),"main",paste0(s,".tsv.gz")),error=function(e) writeLines(conditionMessage(e),file.path("logs",paste0(s,"_failure.txt"))))
}

