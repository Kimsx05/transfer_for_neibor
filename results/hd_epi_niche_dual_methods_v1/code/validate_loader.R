args<-commandArgs(TRUE);root<-normalizePath(args[1]);ext<-file.path(dirname(root),'hd_bin8_banksy_extension_v1')
source(file.path(ext,'code/load_result.R'));source(file.path(root,'code/load_result.R'))
rows<-list()
for(s in c('Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2')) {
 x<-load_extension(ext,s);oldids<-colnames(x[['RNA']]);cts<-SeuratObject::LayerData(x,assay='RNA',layer='counts');total<-sum(cts);qc<-x[[]][,intersect(c('qc_keep','qc_keep_expression','in_tissue_original'),names(x[[]])),drop=FALSE]
 x<-attach_hd_niches(x,root,s)
 stopifnot(identical(oldids,colnames(x[['RNA']])),identical(cts,SeuratObject::LayerData(x,assay='RNA',layer='counts')),identical(qc,x[[]][,names(qc),drop=FALSE]))
 a<-x[[]]$epi_anchor;stopifnot(!anyNA(a));niche_cols<-grep('^niche_.*_niche$',names(x[[]]),value=TRUE);stopifnot(length(niche_cols)>=8,all(is.na(as.matrix(x[[]][!a,niche_cols,drop=FALSE]))))
 rows[[s]]<-data.frame(section=s,n_original=length(oldids),anchors=sum(a),UMI=total,IDs_aligned=TRUE,original_counts_unchanged=TRUE,original_QC_unchanged=TRUE,nonanchors_not_labeled=TRUE,loader_executed=TRUE)
 print(rows[[s]]);rm(x,cts);gc()
}
data.table::fwrite(do.call(rbind,rows),file.path(root,'tables/loader_validation.tsv'),sep='\t')
