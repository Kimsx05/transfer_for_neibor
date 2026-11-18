load_followup <- function(root, section, kind=c('regions','paired','membership','neighbors')) {
 kind <- match.arg(kind); root <- normalizePath(root)
 py <- file.path(dirname(root),'hd_epi_niche_multimethod_v2/.venv/bin/python')
 tmp <- tempfile(fileext='.tsv'); on.exit(unlink(tmp),add=TRUE)
 status <- system2(py,c(shQuote(file.path(root,'code/load_result.py')),shQuote(section),kind,shQuote(tmp)))
 if(status!=0)stop('Followup loader failed')
 read.delim(tmp,check.names=FALSE,stringsAsFactors=FALSE,na.strings='NA')
}
attach_followup_region_evidence <- function(table, root, section, id_column='region_id') {
 d<-load_followup(root,section);ids<-table[[id_column]]
 stopifnot(!anyDuplicated(ids),!anyDuplicated(d$region_id),all(ids%in%d$region_id))
 d<-d[match(ids,d$region_id),];fields<-c('old_evidence','followup_evidence','calibration_status','label_change','score_zero_reason','formal_communication_ready')
 for(f in fields){key<-paste0('fu2_',f);if(key%in%names(table))stop('Refusing overwrite');table[[key]]<-d[[f]]}
 table
}
# Membership attaches only the regional ID to bins. Regional scores are not independent bin observations.
