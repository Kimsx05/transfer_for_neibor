source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R')
suppressPackageStartupMessages(library(Seurat));dir.create(file.path(B,'objects/scenic'),showWarnings=FALSE)
a<-rd('analysis_cell_metadata.tsv.gz');o<-qread(objpath,nthreads=1);x<-LayerData(o,assay='RNA',layer='data');c<-LayerData(o,assay='RNA',layer='counts');out<-list()
for(lin in c('Endothelial','T/NK','Epithelial')){
 name<-c('Endothelial'='Endothelial','T/NK'='T_lineage','Epithelial'='Epi')[[lin]];z<-a[analysis_major_c2l_final==lin & (lin!='T/NK' | analysis_celltype_c2l_final!='NK')];setorder(z,barcode);set.seed(123);z[,priority:=runif(.N)];setorder(z,sample_uid,priority,barcode);z[,train:=seq_len(.N)<=100,by=sample_uid]
 tr<-z[train==TRUE];expr<-x[,tr$barcode,drop=FALSE];keep<-Matrix::rowSums(expr>0)>=max(3,ceiling(.01*ncol(expr)));expr<-expr[keep,,drop=FALSE];out[[name]]<-data.table(lineage=name,n_full=nrow(z),n_training=nrow(tr),n_genes=nrow(expr),n_specimens=uniqueN(z$sample_uid),gene_hash=digest(rownames(expr),algo='sha256'),training_hash=digest(tr$barcode,algo='sha256'))
 put(z,paste0('scenic_membership_',name,'.tsv.gz'));put(data.table(gene=rownames(expr)),paste0('scenic_genes_',name,'.tsv'))
 fwrite(data.table(gene=rownames(expr),as.data.table(as.matrix(expr))),file.path(B,'objects/scenic',paste0(name,'_training.tsv.gz')),sep='\t',compress='gzip');cat(name,nrow(expr),ncol(expr),'\n')
}
put(rbindlist(out),'scenic_training_audit.tsv');cat('SCENIC TRAINING INPUTS COMPLETE\n')
