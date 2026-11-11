# Original counts are read once from a native H5 or the sole MPBC sparse 8um aggregate.
load_extension <- function(root,section,attach_results=TRUE) {
 root<-normalizePath(root);o<-file.path(root,section)
 .libPaths(c(file.path(normalizePath(file.path(root,'../../..')),'mp10_pathology_regions_v2/Rlib'),file.path(dirname(root),'hd_annotation_diagnosis_singler_v1/Rlib'),.libPaths()))
 if(section%in%c('Cervilla','D1_Johnson')){
  e<-new.env(parent=globalenv());sys.source(file.path(dirname(root),'hd_bin8_banksy_v1/code/load_result.R'),envir=e);x<-e$load_banksy(file.path(dirname(root),'hd_bin8_banksy_v1',section))
  f<-file.path(o,'supplement_metadata.tsv.gz');if(file.exists(f)){d<-data.table::fread(f,data.table=FALSE);stopifnot(!anyDuplicated(d$bin_id));d<-d[match(colnames(x[['RNA']]),d$bin_id),];stopifnot(!anyNA(d$bin_id));rownames(d)<-d$bin_id;d$bin_id<-NULL;stopifnot(!any(names(d)%in%colnames(x[[]])));x<-SeuratObject::AddMetaData(x,d)};return(x)
 }
 info<-jsonlite::read_json(file.path(o,'input_summary.json'),simplifyVector=TRUE);d<-data.table::fread(file.path(o,'analysis_metadata.tsv.gz'),data.table=FALSE);ids<-d$bin_id;stopifnot(!anyDuplicated(ids));c<-Seurat::Read10X_h5(file.path(o,info$counts_relative),use.names=FALSE);stopifnot(all(ids%in%colnames(c)));c<-c[,ids,drop=FALSE];rownames(d)<-ids;x<-SeuratObject::CreateSeuratObject(c,assay='RNA',meta.data=d,min.cells=0,min.features=0)
 g<-data.table::fread(file.path(o,'gene_mapping.tsv'),data.table=FALSE);stopifnot(identical(rownames(x),g$feature_id));rownames(g)<-g$feature_id;x[['RNA']]<-SeuratObject::AddMetaData(x[['RNA']],g)
 qfile<-file.path(o,'normalized.rds');if(file.exists(qfile)){q<-readRDS(qfile);stopifnot(identical(dimnames(q),dimnames(c)));SeuratObject::LayerData(x,assay='RNA',layer='data')<-q;rm(q)}
 hv<-file.path(o,'frozen_HVG.tsv');if(file.exists(hv))SeuratObject::VariableFeatures(x[['RNA']])<-data.table::fread(hv)$feature_id
 x@misc$analysis_bin_ids<-ids;x@misc$source_counts<-info$counts_relative;x@misc$coordinate_note<-jsonlite::read_json(file.path(o,'transform.json'));x@misc$in_tissue_known<-info$in_tissue_known
 if(attach_results){f<-file.path(o,'metadata.tsv.gz');if(file.exists(f)){m<-data.table::fread(f,data.table=FALSE);stopifnot(identical(m$bin_id,ids));m<-m[,setdiff(names(m),names(x[[]])),drop=FALSE];rownames(m)<-ids;x<-SeuratObject::AddMetaData(x,m)}
  for(tag in c('l00','l02','l08')){f<-file.path(o,paste0(tag,'_embeddings.rds'));if(!file.exists(f))next;z<-readRDS(f);for(kind in c('pca','umap')){em<-z[[kind]];stopifnot(identical(rownames(em),ids),all(is.finite(em)));key<-paste0('B',sub('l','',tag),toupper(kind),'_');colnames(em)<-paste0(key,seq_len(ncol(em)));x[[paste0('banksy_',tag,'_',kind)]]<-SeuratObject::CreateDimReducObject(embeddings=em,key=key,assay='RNA')}}
 }
 f<-file.path(o,'supplement_metadata.tsv.gz');if(attach_results&&file.exists(f)){d<-data.table::fread(f,data.table=FALSE);stopifnot(identical(d$bin_id,ids));rownames(d)<-d$bin_id;d$bin_id<-NULL;stopifnot(!any(names(d)%in%colnames(x[[]])));x<-SeuratObject::AddMetaData(x,d)}
 SeuratObject::DefaultAssay(x)<-'RNA';x
}
