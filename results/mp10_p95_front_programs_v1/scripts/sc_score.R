.libPaths(c(file.path(Sys.getenv('HOME'),'R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
suppressPackageStartupMessages({library(Matrix);library(data.table)})
setDTthreads(1)
o='results/mp10_p95_front_programs_v1'
x=as(readMM(file.path(o,'scores/sc_counts_temporary.mtx')),'CsparseMatrix');genes=readLines(file.path(o,'scores/sc_genes.txt'))
if(anyDuplicated(genes)){u=unique(genes);x=sparseMatrix(i=match(genes,u),j=seq_along(genes),x=1,dims=c(length(u),length(genes)))%*%x;genes=u}
rownames(x)=genes;m=fread(file.path(o,'tables/sc_cells.tsv'));colnames(x)=m$cell
gs=fread(file.path(o,'tables/new_gene_sets.tsv'));gs=gs[gs_name %in% c('PURAM_PEMT_TOP30','FLORES_CRISPR_PRO_INVASION_10')]
sets=lapply(split(gs$gene_symbol,gs$gs_name),function(z)intersect(unique(z),genes))
u=UCell::ScoreSignatures_UCell(x,features=sets,name=NULL,ncores=1,chunk.size=50,maxRank=1500,ties.method='average');z=cbind(m,as.data.table(u[m$cell,]));fwrite(z,file.path(o,'tables/sc_cell_scores.tsv.gz'),sep='\t')
fwrite(z[,lapply(.SD,median),by=.(sample,type,lineage,dataset),.SDcols=names(sets)],file.path(o,'tables/sc_sample_scores.tsv'),sep='\t')
a=list();j=1
for(n in names(sets))for(g in sets[[n]]){v=as.numeric(x[g,]);a[[j]]=cbind(m[,.(sample,type,lineage,dataset)],gs_name=n,gene=g,detected=v>0,logCP10K=log1p(v/Matrix::colSums(x)*10000));j=j+1}
a=rbindlist(a);fwrite(a[,.(fraction_detected=mean(detected),mean_logCP10K=mean(logCP10K),n_cells=.N),by=.(sample,type,lineage,dataset,gs_name,gene)],file.path(o,'tables/sc_gene_source.tsv'),sep='\t')
