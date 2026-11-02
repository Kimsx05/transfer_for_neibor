# Run from the project root. Use the installed UCell library and local MSigDB cache.
suppressPackageStartupMessages(library(msigdbr))
sets <- rbind(msigdbr(species="Homo sapiens",collection="H"),
   subset(msigdbr(collection="C2",subcollection="CP:REACTOME"), gs_name %in% paste0("REACTOME_",c("CHOLESTEROL_BIOSYNTHESIS","EXTRACELLULAR_MATRIX_ORGANIZATION","CO_INHIBITION_BY_CTLA4"))),
   subset(msigdbr(collection="C7",subcollection="IMMUNESIGDB"), gs_name %in% paste0("GSE15659_RESTING_VS_ACTIVATED_TREG_",c("DN","UP"))))
.libPaths(c(file.path(Sys.getenv("HOME"),"R/x86_64-pc-linux-gnu-library/4.5"),.libPaths()))
suppressPackageStartupMessages(library(Matrix))
 data.table::setDTthreads(2)
out <- "mp10_paper_distance_pilot_v1"
members <- unique(sets[,c("gs_name","gene_symbol","db_version","gs_collection","gs_subcollection","gs_id")])
members_out <- gzfile(file.path(out,"tables/gene_sets.tsv.gz"),"wt")
write.table(members,members_out,sep="\t",row.names=FALSE,quote=FALSE)
close(members_out)
gs <- lapply(split(members$gene_symbol,members$gs_name),unique)
stopifnot(length(gs)==55)
manifest <- read.delim(file.path(out,"tables/input_manifest.tsv"))
for (i in seq_len(nrow(manifest))) {
 r <- manifest[i,]; sec <- r$section_id
 path <- file.path(out,"spots",paste0(sec,"_programs.tsv.gz"))
 if(file.exists(path)) next
 message(sec)
 con <- if(grepl("gz$",r$counts_path)) gzfile(r$counts_path,"rt") else file(r$counts_path,"rt")
 x <- as(readMM(con),"CsparseMatrix"); close(con)
 genes <- readLines(r$genes_path); cells <- readLines(r$barcodes_path)
 stopifnot(nrow(x)==length(genes),ncol(x)==length(cells),!anyDuplicated(cells))
 if(anyDuplicated(genes)){
  u <- unique(genes)
  x <- sparseMatrix(i=match(genes,u),j=seq_along(genes),x=1,dims=c(length(u),length(genes))) %*% x
  genes <- u
 }
 dimnames(x) <- list(genes,cells)
 avail <- lapply(gs,intersect,genes)
 matched <- lengths(avail)
 coverage <- data.frame(section_id=sec,gs_name=names(gs),n_original=lengths(gs),n_matched=matched,n_universe=nrow(x),db_version=unique(sets$db_version),status=ifelse(matched>0,"scored","empty"))
 write.table(coverage,file.path(out,"tables",paste0(sec,"_match.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 uc <- UCell::ScoreSignatures_UCell(matrix=x,features=avail[matched>0],name=NULL,ncores=1,chunk.size=50)
 scores <- data.frame(section_id=sec,barcode=cells,uc[cells,,drop=FALSE],check.names=FALSE)
 for(k in names(gs)[matched==0])scores[[k]] <- NA_real_
 if(i<=2){
  obj <- Seurat::NormalizeData(Seurat::CreateSeuratObject(x,min.cells=0,min.features=0),verbose=FALSE)
  cmp <- c("REACTOME_CHOLESTEROL_BIOSYNTHESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_HYPOXIA")
  obj <- Seurat::AddModuleScore(obj,features=avail[cmp],name="ref",seed=1)
  scores[paste0("AMS_",cmp)] <- obj[[]][cells,paste0("ref",1:3)]
  rm(obj)
 }
 con <- gzfile(paste0(path,".tmp"),"wt")
 write.table(scores,con,sep="\t",quote=FALSE,row.names=FALSE);close(con);file.rename(paste0(path,".tmp"),path)
 rm(x,uc,scores);gc()
}
