# No input labels or expression layers are changed; no duplicate Seurat object is saved.
attach_diagnosis <- function(object, section_dir) {
  d <- data.table::fread(file.path(section_dir, 'diagnosis.tsv.gz'), data.table=FALSE,na.strings=c('','NA'))
  stopifnot(!anyDuplicated(d$cell_id), length(unique(d$section))==1)
  ids <- colnames(object); ix <- match(ids,d$cell_id)
  if(anyNA(ix))stop('Object contains cells outside the fixed analysis denominator; subset explicitly first.')
  if('section_id'%in%colnames(object[[]]))stopifnot(all(as.character(object$section_id)==unique(d$section)))
  keep <- grep('^(diag_|singler_)',names(d),value=TRUE)
  md <- d[ix,keep,drop=FALSE];rownames(md)<-ids
  existing <- intersect(keep,colnames(object[[]]))
  if(length(existing))stop('Diagnostic fields already exist: load an original object to avoid overwriting.')
  SeuratObject::AddMetaData(object,metadata=md)
}
load_diagnosis <- function(section_dir) data.table::fread(file.path(section_dir,'diagnosis.tsv.gz'),na.strings=c('','NA'))
# Example:
# source('../hd_major_annotation_v1/code/load_result.R')
# x <- load_major('../hd_major_annotation_v1/Cervilla')
# x <- attach_diagnosis(x, 'Cervilla')
# Original x$final_major is unchanged. Candidate: x$singler_blca_pruned_labels.
