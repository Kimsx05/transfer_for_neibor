.libPaths(c(file.path(Sys.getenv("HOME"),"R/x86_64-pc-linux-gnu-library/4.5"),.libPaths()))
suppressPackageStartupMessages(library(Matrix))
 data.table::setDTthreads(2)
args <- commandArgs(TRUE)
out <- args[1]
members <- read.delim(file.path(out,"tables/new_gene_sets.tsv"))
stopifnot(as.character(packageVersion("UCell"))=="2.14.0")
gs <- lapply(split(members$gene_symbol,members$gs_name),unique)

manifest <- read.delim(file.path(out,"tables/input_manifest.tsv"))
score_section <- function(i) {
 r <- manifest[i,]; sec <- r$section_id
 path <- file.path(out,"scores",paste0(sec,"_programs.tsv.gz"))
 if(file.exists(path)) return(invisible(NULL))
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
 detected <- vapply(avail,function(z)sum(Matrix::rowSums(x[z,,drop=FALSE])>0),integer(1))
 det <- do.call(rbind,lapply(names(avail),function(n){z=avail[[n]]; data.frame(section_id=sec,gs_name=n,gene_symbol=z,n_detected_spots=Matrix::rowSums(x[z,,drop=FALSE]>0),mean_count=Matrix::rowMeans(x[z,,drop=FALSE]))}))
 write.table(det,file.path(out,"scores",paste0(sec,"_gene_detection.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 coverage <- data.frame(section_id=sec,gs_name=names(gs),n_original=lengths(gs),n_matched=matched,n_detected=detected,n_universe=nrow(x),db_version="literature + MSigDB 2026.1.Hs",status=ifelse(matched>0,"scored","empty"))
 write.table(coverage,file.path(out,"scores",paste0(sec,"_match.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
 uc <- UCell::ScoreSignatures_UCell(matrix=x,features=avail[matched>0],name=NULL,maxRank=1500,ties.method="average",ncores=1,chunk.size=50)
 scores <- data.frame(section_id=sec,barcode=cells,uc[cells,,drop=FALSE],check.names=FALSE)
 for(k in names(gs)[matched==0])scores[[k]] <- NA_real_
 con <- gzfile(paste0(path,".tmp"),"wt")
 write.table(scores,con,sep="\t",quote=FALSE,row.names=FALSE);close(con);file.rename(paste0(path,".tmp"),path)
 rm(x,uc,scores);gc()
}

for (i in seq(as.integer(args[2]),nrow(manifest),by=as.integer(args[3]))) {
  tryCatch(score_section(i),error=function(e){
    message(manifest$section_id[i]," score unavailable: ",conditionMessage(e))
    writeLines(conditionMessage(e),file.path(out,"logs",paste0(manifest$section_id[i],"_score_failure.txt")))
  })
}
