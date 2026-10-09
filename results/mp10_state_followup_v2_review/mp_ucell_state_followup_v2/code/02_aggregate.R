source('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_followup_v2/code/00_common.R')
suppressPackageStartupMessages(library(Seurat))
md<-load_meta();ids<-md$cell_id
if(file.exists(out('02_pseudobulk_markers/aggregate_COMPLETE'))) {logmsg('Aggregate checkpoint reused');quit(save='no')}
logmsg('Loading frozen original Seurat once')
s<-qread(input,nthreads=2);stopifnot(setequal(colnames(s),ids));layers<-Layers(s[['RNA']]);cl<-grep('^counts($|\\.)',layers,value=TRUE);stopifnot(length(cl)==1)
x<-LayerData(s[['RNA']],layer=cl)[,ids,drop=FALSE];stopifnot(inherits(x,'sparseMatrix'),!anyDuplicated(rownames(x)),all(is.finite(x@x)),all(x@x>=0))
nonint<-sum(x@x!=floor(x@x));wt(data.table(n_features=nrow(x),n_cells=ncol(x),nnz=length(x@x),noninteger_entries=nonint,counts_layer=cl),'00_audit/counts_check.tsv');if(nonint>0)stop('Counts are noninteger; do not round')
wt(data.table(gene_id=rownames(x)),'00_audit/RNA_features.tsv.gz')
rednames<-Reductions(s);writeLines(rednames,out('00_audit/original_reductions.txt'))
if('umap_allgene'%in%rednames) {emb<-Embeddings(s[['umap_allgene']]);stopifnot(setequal(rownames(emb),ids),ncol(emb)==2);qsafe(emb[ids,],'00_audit/original_RNA_umap.qs');writeLines('Frozen original Seurat reduction umap_allgene, read only, cell IDs exact.',out('00_audit/RNA_umap_source.txt'))}
capture.output(s@commands,file=out('00_audit/original_Seurat_commands.txt'))
lib<-Matrix::colSums(x);stopifnot(all(lib>0));qsafe(setNames(lib,ids),'02_pseudobulk_markers/cell_library_sizes.qs')
# Verify existing data against ordinary LogNormalize over all features of 25 fixed cells.
dl<-grep('^data$',layers,value=TRUE);use_existing<-FALSE
if(length(dl)==1) {d<-LayerData(s[['RNA']],layer=dl);if(identical(rownames(d),rownames(x))&&setequal(colnames(d),ids)) {ix<-unique(round(seq(1,length(ids),length.out=25)));test<-x[,ix,drop=FALSE]%*%Diagonal(x=1e4/lib[ix]);test@x<-log1p(test@x);maxdiff<-max(abs(d[,ids[ix],drop=FALSE]-test));use_existing<-is.finite(maxdiff)&&maxdiff<1e-6} else maxdiff<-NA_real_} else maxdiff<-NA_real_
if(use_existing) {y<-d[,ids,drop=FALSE];method<-'Verified existing RNA data: log1p(counts/library*10000), all features of 25 deterministic cells; tolerance1e-6'} else {y<-x%*%Diagonal(x=1e4/lib);dimnames(y)<-dimnames(x);y@x<-log1p(y@x);method<-'Computed ordinary LogNormalize scale.factor=10000 in memory from counts; no original mutation'}
wt(data.table(method=method,max_abs_difference=maxdiff,existing_data_used=use_existing),'00_audit/expression_display_input.tsv')
rm(s);if(exists('d'))rm(d);gc()
detected<-x;detected@x[]<-1
for(r in c('r04','r06')) {p<-paste0('02_pseudobulk_markers/aggregate_',r,'.qs');if(file.exists(out(p)))next
 b<-unique(md[,.(sample,dataset,source,cluster=get(r))]);setorder(b,sample,cluster);b[,bin:=.I];key<-paste(b$sample,b$cluster,sep='\r');j<-match(paste(md$sample,md[[r]],sep='\r'),key);stopifnot(!anyNA(j));ind<-sparseMatrix(i=seq_along(ids),j=j,x=1,dims=c(length(ids),nrow(b)));b[,n_cells:=as.integer(Matrix::colSums(ind))];b[,library_size:=as.numeric(crossprod(ind,lib))]
 logmsg('Aggregating',r,nrow(b),'sample-cluster bins');obj<-list(meta=b,counts=x%*%ind,expression_sum=y%*%ind,detected_sum=detected%*%ind);qsafe(obj,p);wt(b,paste0('02_pseudobulk_markers/aggregate_',r,'_metadata.tsv'));rm(obj);gc()
}
writeLines('complete',out('02_pseudobulk_markers/aggregate_COMPLETE'));capture.output(sessionInfo(),file=out('logs/sessionInfo_aggregate.txt'));logmsg('Sparse aggregation complete')
