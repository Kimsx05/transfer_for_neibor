# No expression copies. All additions use pc_ prefix and exact original cell IDs.
load_precommunication <- function(root, section='Cervilla', include_all=TRUE) {
 root<-normalizePath(root);if(section!='Cervilla')stop('D1 is light evaluability only; no formal result object was generated')
 py<-file.path(dirname(root),'hd_epi_niche_multimethod_v2/.venv/bin/python')
 tmp<-tempfile(fileext='.tsv');on.exit(unlink(tmp),add=TRUE)
 status<-system2(py,c(shQuote(file.path(root,'code/12_loaders.py')),shQuote(tmp)));if(status!=0)stop('Metadata conversion failed')
 d<-data.table::fread(tmp,data.table=FALSE);stopifnot(!anyDuplicated(d$cell_id));rownames(d)<-d$cell_id;d$cell_id<-NULL;names(d)<-paste0('pc_',names(d))
 source(file.path(dirname(root),'segmentation_assignment_v2/code/load_result.R'),local=TRUE)
 x<-load_assignment(file.path(dirname(root),'segmentation_assignment_v2',section),qc=if(include_all)'all' else 'analysis')
 stopifnot(all(rownames(d)%in%colnames(x)));x<-Seurat::AddMetaData(x,d);x$pc_in_analysis<-colnames(x)%in%rownames(d)
 x@misc$precommunication<-list(root=root,MP10_unit='fixed connected epithelial-only 64um region; see regional table, no replicated single-cell association score',new_counts=FALSE)
 x
}
load_precommunication_regions <- function(root,radius=50) {
 data.table::fread(file.path(root,'tables',paste0('region_result_r',radius,'.tsv.gz')),data.table=FALSE)
}
