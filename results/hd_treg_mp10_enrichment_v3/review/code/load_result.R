# Region-level data, not per-bin cell identities. Does not load/copy expression matrices.
.v3_loader_root <- dirname(dirname(normalizePath(sys.frame(1)$ofile)))
load_enrichment <- function(section, root=.v3_loader_root) {
  stopifnot(section %in% c('Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'))
  reader <- function(kind, name) {
    p <- file.path(root, section, name)
    if (requireNamespace('arrow',quietly=TRUE)) return(as.data.frame(arrow::read_parquet(p)))
    py <- file.path(dirname(root),'hd_epi_niche_multimethod_v2','.venv','bin','python')
    out <- tempfile(fileext='.tsv');on.exit(unlink(out))
    status <- system2(py,c(shQuote(file.path(root,'code','load_result.py')),'--section',shQuote(section),'--kind',kind,'--tsv',shQuote(out)))
    if(status != 0L)stop('Python Parquet reader failed')
    d <- read.delim(out,check.names=FALSE,na.strings='NA')
    d[] <- lapply(d,function(x) if(is.character(x) && all(is.na(x) | x %in% c('True','False'))) ifelse(is.na(x),NA,x=='True') else x)
    d
  }
  z <- list(T=reader('T','T_regions.parquet'),Epi=reader('Epi','Epi_regions.parquet'),neighbors=reader('neighbors','neighborhoods.parquet'))
  stopifnot(!anyDuplicated(z$T$region_id),!anyDuplicated(z$Epi$region_id))
  z$effects <- read.delim(file.path(root,'tables','enrichment_effects.tsv'))
  z$effects <- z$effects[z$effects$section %in% section,,drop=FALSE]
  z
}
attach_enrichment_regions <- function(metadata,section,kind='T',id_col='region_id',root=.v3_loader_root) {
  d<-load_enrichment(section,root)[[kind]];ids<-if(id_col %in% names(metadata))metadata[[id_col]]else rownames(metadata)
  stopifnot(!anyDuplicated(ids));ii<-match(ids,d$region_id);cols<-setdiff(names(d),'region_id');dest<-paste0('v3_',cols)
  stopifnot(!any(dest %in% names(metadata)));metadata[dest]<-d[ii,cols,drop=FALSE];metadata
}
