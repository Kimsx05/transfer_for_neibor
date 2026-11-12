attach_hd_niches <- function(object, root, section) {
 root <- normalizePath(root)
 temp <- tempfile(fileext='.tsv');on.exit(unlink(temp))
 status <- system2(file.path(root,'.venv/bin/python'),c(shQuote(file.path(root,'code/load_result.py')),shQuote(root),shQuote(section),shQuote(temp)))
 stopifnot(status==0)
 d<-data.table::fread(temp,data.table=FALSE);stopifnot(!anyDuplicated(d$bin_id),all(d$bin_id%in%colnames(object[['RNA']])))
 rownames(d)<-d$bin_id;d$bin_id<-NULL;stopifnot(!any(names(d)%in%names(object[[]])))
 d<-d[match(colnames(object),rownames(d)),,drop=FALSE];rownames(d)<-colnames(object);d$epi_anchor[is.na(d$epi_anchor)]<-FALSE
 SeuratObject::AddMetaData(object,d)
}
