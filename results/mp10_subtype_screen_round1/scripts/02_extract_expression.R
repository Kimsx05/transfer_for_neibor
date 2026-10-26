.libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(qs);library(Seurat);library(CellChat);library(Matrix)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1'
p<-file.path(dirname(root),'subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs')
cat('Loading one final annotated atlas',format(Sys.time()),'\n');obj<-qread(p,nthreads=2)
cat('assays/layers\n');print(Assays(obj));print(Layers(obj[['RNA']]));print(dim(obj));print(names(obj@meta.data))
DB<-CellChatDB.human
# Only database signalling genes are materialized; retain sparse normalized expression.
genes<-extractGene(DB); lay<-Layers(obj[['RNA']]); stopifnot('data'%in%lay)
x<-LayerData(obj,assay='RNA',layer='data',features=intersect(genes,rownames(obj[['RNA']])) )
stopifnot(inherits(x,'sparseMatrix'),min(x@x)>=0)
cat('RNA data',dim(x),'range',range(x@x),'\n')
qs::qsave(x,file.path(root,'objects/signaling_RNA_data.qs'),preset='fast',nthreads=2)
fwrite(data.table(gene=genes,present_in_RNA=genes%in%rownames(x)),file.path(root,'tables/CellChat_gene_coverage.tsv'),sep='\t')
fwrite(data.table(package=c('R','Seurat','CellChat','limma','qs'),version=c(as.character(getRversion()),sapply(c('Seurat','CellChat','limma','qs'),function(z)as.character(packageVersion(z))))),file.path(root,'tables/package_versions.tsv'),sep='\t')
writeLines(c(p,'RNA:data; existing nonnegative LogNormalize expression; no scaled residuals',capture.output(sessionInfo())),file.path(root,'logs/expression_source_session.txt'))
cat('Done',format(Sys.time()),'\n')
