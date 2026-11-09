# No changes to original cell objects. Paths resolve from the result directory.
load_bin8 <- function(section_dir, attach_roi_rctd=TRUE) {
  section_dir <- normalizePath(section_dir)
  x <- readRDS(file.path(section_dir,'seurat.rds'))
  SeuratObject::DefaultAssay(x) <- 'RNA'
  ids <- x@misc$analysis_bin_ids
  stopifnot(identical(colnames(x[['RNA']]),ids))
  if(attach_roi_rctd && file.exists(file.path(section_dir,'rctd_ROI_results.tsv.gz'))) {
    d <- read.delim(gzfile(file.path(section_dir,'rctd_ROI_results.tsv.gz')),check.names=FALSE)
    j <- match(ids,d$bin_id);m<-data.frame(row.names=ids,rctd_status=ifelse(is.na(j),'not_tested',d$rctd_status[j]))
    for(k in intersect(c('first_type','second_type','weight_1','weight_2'),names(d)))m[[paste0('rctd_',k)]]<-d[[k]][j]
    x<-SeuratObject::AddMetaData(x,m)
  }
  x@misc$resolved_spatial_assets <- normalizePath(file.path(section_dir,x@misc$external_assets_relative),mustWork=FALSE)
  x
}
# source('code/load_result.R'); x <- load_bin8('Cervilla')
# SeuratObject::LayerData(x, assay='RNA', layer='counts')
# SeuratObject::Embeddings(x, reduction='umap.full')
# x$singler_blca_pruned_labels are diagnostic candidates, not accepted cell identities.
