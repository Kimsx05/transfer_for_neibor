load_treg_mp10 <- function(root, section, kind=c('T32','T16','Epi64','membership')) {
 kind<-match.arg(kind);root<-normalizePath(root)
 py<-file.path(dirname(root),'hd_epi_niche_multimethod_v2/.venv/bin/python')
 tmp<-tempfile(fileext='.tsv');on.exit(unlink(tmp),add=TRUE)
 status<-system2(py,c(shQuote(file.path(root,'code/load_result.py')),shQuote(section),kind,shQuote(tmp)))
 if(status!=0)stop('Result loader failed')
 data.table::fread(tmp,data.table=FALSE)
}
attach_treg_mp10_membership <- function(object, root, section, assay='RNA') {
 d<-load_treg_mp10(root,section,'membership');ids<-colnames(object[[assay]])
 stopifnot(!anyDuplicated(ids),!anyDuplicated(d$bin_id),all(ids%in%d$bin_id));d<-d[match(ids,d$bin_id),];rownames(d)<-d$bin_id
 d<-d[,grep('^trmp_',names(d)),drop=FALSE];stopifnot(!any(names(d)%in%names(object[[]])))
 SeuratObject::AddMetaData(object,d)
}
# Scores belong to regions. Membership attachment does not replicate regional scores as independent bins.
