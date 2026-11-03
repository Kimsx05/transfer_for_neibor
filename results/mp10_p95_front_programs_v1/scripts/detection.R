.libPaths(c(file.path(Sys.getenv('HOME'),'R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
suppressPackageStartupMessages({library(Matrix);library(data.table)});setDTthreads(1)
a=commandArgs(TRUE);o='results/mp10_p95_front_programs_v1';m=fread(file.path(o,'tables/input_manifest.tsv'));g=fread(file.path(o,'tables/new_gene_sets.tsv'));g=g[gs_name %in% c('PURAM_PEMT_TOP30','FLORES_CRISPR_PRO_INVASION_10')]
for(i in seq(as.integer(a[1]),nrow(m),by=6)){
 r=m[i,];p=file.path(o,'scores',paste0(r$section_id,'_signature_detection.tsv.gz'));if(file.exists(p))next
 con=gzfile(r$counts_path,'rt');x=as(readMM(con),'CsparseMatrix');close(con);genes=readLines(r$genes_path);cells=readLines(r$barcodes_path)
 if(anyDuplicated(genes)){u=unique(genes);x=sparseMatrix(i=match(genes,u),j=seq_along(genes),x=1,dims=c(length(u),length(genes)))%*%x;genes=u};rownames(x)=genes
 rows=list()
 for(n in unique(g$gs_name)){idx=intersect(g[gs_name==n,gene_symbol],genes);xx=x[idx,,drop=FALSE];total=Matrix::colSums(xx);top=apply(as.matrix(xx),2,max);rows[[n]]=data.table(section_id=r$section_id,barcode=cells,program=n,n_detected_genes=Matrix::colSums(xx>0),signature_counts=total,top1_count_share=ifelse(total>0,top/total,NA_real_))}
 fwrite(rbindlist(rows),p,sep='\t');rm(x);gc();message(r$section_id)
}
