attach_hd_niches_v2 <- function(object, root, section) {
 root<-normalizePath(root);f<-tempfile(fileext='.tsv');on.exit(unlink(f))
 st<-system2(file.path(root,'.venv/bin/python'),c(shQuote(file.path(root,'code/load_result.py')),shQuote(root),shQuote(section),shQuote(f)));stopifnot(st==0)
 d<-data.table::fread(f,data.table=FALSE);stopifnot(!anyDuplicated(d$bin_id),all(d$bin_id%in%colnames(object[['RNA']])));rownames(d)<-d$bin_id;d$bin_id<-NULL;stopifnot(!any(names(d)%in%names(object[[]])))
 d<-d[match(colnames(object),rownames(d)),,drop=FALSE];rownames(d)<-colnames(object);d$niche_v2_anchor[is.na(d$niche_v2_anchor)]<-FALSE
 SeuratObject::AddMetaData(object,d)
}
