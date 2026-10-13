suppressPackageStartupMessages({library(qs);library(Seurat);library(Matrix)})
OUT <- '/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_mp10_calibration_v5'
stopifnot(basename(OUT)=='neighborhood_mp10_calibration_v5')
p <- '/home/data/t070721/Metabolistic_scRNA/Multi2/3.4 malignant cell process/PCA_ALLGENE1/epi_cycling_allgene_harmony_clustered_final_after_small_removal.qs'
s <- qread(p); cat('object',dim(s),'assays',names(s@assays),'\n');print(colnames(s@meta.data));print(Layers(s[['RNA']]))
source('/home/data/t070721/codex_workspace/Bladder Metabolism/ST/mapping/12_epimp_fulltranscriptome/00_utils.R')
x <- bind_assay_layers(s,'RNA','counts');stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)));cat('counts',dim(x),'nnz',length(x@x),'\n')
d <- file.path(OUT,'01_score_fidelity','counts_cache');dir.create(d,showWarnings=FALSE)
writeBin(x@i,file.path(d,'indices.bin'),size=4,endian='little');writeBin(x@p,file.path(d,'indptr.bin'),size=4,endian='little');writeBin(x@x,file.path(d,'data.bin'),size=8,endian='little');writeLines(rownames(x),file.path(d,'genes.txt'));writeLines(colnames(x),file.path(d,'cells.txt'))
md <- s@meta.data[colnames(x),,drop=FALSE];md$cell_id <- rownames(md);write.table(md,gzfile(file.path(d,'metadata.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE)
u <- qread('/home/data/t070721/codex_workspace/Bladder Metabolism/ST/mapping/08_cell_EpiMP_scores.qs');write.table(data.frame(cell_id=colnames(u$usage_raw),MP10_usage=as.numeric(u$usage_raw['EpiMP10',])),file.path(OUT,'01_score_fidelity','original_usage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c(p,paste('n_genes',nrow(x)),paste('n_cells',ncol(x)),paste('nnz',length(x@x)),paste('usage_basis',u$basis_name),paste('usage_projection',u$projection)),file.path(d,'EXPORT_COMPLETE.txt'))
cat('AUCell available',requireNamespace('AUCell',quietly=TRUE),'\n')
