# Incremental attachment. Original RNA counts/data and original labels are retained.
load_banksy <- function(section_dir, branches=c('l00','l02','l08'), attach_roi_rctd=TRUE) {
 section_dir<-normalizePath(section_dir)
 source_dir<-file.path(dirname(dirname(section_dir)),'hd_bin8_comparison_v1',basename(section_dir))
 e<-new.env(parent=globalenv());sys.source(file.path(dirname(source_dir),'code/load_result.R'),envir=e)
 x<-e$load_bin8(source_dir,attach_roi_rctd=attach_roi_rctd)
 ids<-colnames(x[['RNA']]);d<-read.delim(gzfile(file.path(section_dir,'metadata.tsv.gz')),check.names=FALSE,stringsAsFactors=FALSE)
 stopifnot(!anyDuplicated(d$bin_id),setequal(ids,d$bin_id));d<-d[match(ids,d$bin_id),,drop=FALSE]
 stopifnot(identical(d$bin_id,ids),max(abs(d$x_um-x$x_um))<1e-7,max(abs(d$y_um-x$y_um))<1e-7)
 d$bin_id<-NULL;d$x_um<-d$y_um<-NULL
 # All added fields carry banksy_ or coarse_ prefixes; original metadata never overwritten.
 names(d)<-ifelse(grepl('^(banksy_|coarse_)',names(d)),names(d),paste0('banksy_',names(d)))
 stopifnot(!any(names(d)%in%names(x[[]])));rownames(d)<-ids;x<-SeuratObject::AddMetaData(x,d)
 for(tag in branches){z<-readRDS(file.path(section_dir,paste0(tag,'_embeddings.rds')));for(kind in c('pca','umap')){v<-z[[kind]];stopifnot(identical(rownames(v),ids),all(is.finite(v)));nm<-paste0('banksy_',tag,'_',kind);key<-paste0('B',sub('l','',tag),toupper(kind),'_');colnames(v)<-paste0(key,seq_len(ncol(v)));x[[nm]]<-SeuratObject::CreateDimReducObject(embeddings=v,key=key,assay='RNA')}}
 SeuratObject::DefaultAssay(x)<-'RNA';x@misc$banksy_incremental_root<-section_dir;x@misc$banksy_original_source<-normalizePath(source_dir);x@misc$banksy_note<-'Full valid bins; cluster IDs are branch-specific. RNA counts/data remain original. Coarse labels are reference candidates, not consensus or calibrated confidence.';x
}
# source('analysis/hd_bin8_banksy_v1/code/load_result.R')
# x <- load_banksy('analysis/hd_bin8_banksy_v1/Cervilla')
# x$banksy_l08_cluster; x$coarse_blca; x$coarse_hpca
# SeuratObject::Embeddings(x,'banksy_l08_umap')
