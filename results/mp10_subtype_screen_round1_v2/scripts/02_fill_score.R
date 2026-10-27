.libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(qs);library(Seurat);library(Matrix)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1_v2';p<-dirname(root);tt<-file.path(root,'tables')
a<-fread(file.path(tt,'score_fill_cells.tsv'));v<-fread(file.path(tt,'score_validation_cells.tsv'));ids<-c(a$barcode,v$barcode)
genes<-readLines(file.path(p,'neighborhood_mp10_calibration_v5/01_score_fidelity/official_genes.txt'));g2<-readLines(file.path(p,'neighborhood_mp10_calibration_v5/01_score_fidelity/counts_cache/genes.txt'));stopifnot(identical(genes,g2),length(genes)==90011)
top<-fread(file.path(tt,'frozen_MP10_Top20.tsv'))$gene
cat('Loading final atlas',format(Sys.time()),'\n');obj<-qread(file.path(p,'subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs'),nthreads=2)
stopifnot(all(genes%in%rownames(obj[['RNA']])),all(ids%in%colnames(obj)))
x<-LayerData(obj,assay='RNA',layer='counts',features=genes,cells=ids);x<-x[genes,ids];stopifnot(inherits(x,'dgCMatrix'))
cat('Score input',dim(x),'\n')
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector score_sparse(IntegerVector ip,IntegerVector ri,NumericVector xx,IntegerVector top,int G,int K,double den){
 int nc=ip.size()-1;NumericVector ans(nc);
 for(int c=0;c<nc;c++){std::vector<double> vals;std::vector<double> tv(top.size(),0);double total=0;
 for(int z=ip[c];z<ip[c+1];z++){double val=xx[z];total+=val;if(val!=0)vals.push_back(val);for(int j=0;j<top.size();j++)if(ri[z]==top[j])tv[j]=val;}
 if(total==0){ans[c]=NA_REAL;continue;}std::sort(vals.begin(),vals.end());double score=0;
 for(double val:tv){double left=std::lower_bound(vals.begin(),vals.end(),val)-vals.begin();double right=std::upper_bound(vals.begin(),vals.end(),val)-vals.begin();double gt=vals.size()-right;double eq=val==0?G-vals.size():right-left;double rank=1+gt+(eq-1)/2;score+=std::max(0.0,K-rank)/den;}ans[c]=score;
 }return ans;}',cacheDir=file.path(root,'objects/rcpp_score'))
z<-score_sparse(x@p,x@i,x@x,match(top,genes)-1L,length(genes),4501L,89810)
fwrite(data.table(barcode=ids,filled_score=z,total_counts=Matrix::colSums(x),n_detected=Matrix::colSums(x>0)),file.path(tt,'score_recalculation.tsv'),sep='\t',na='NA')
writeMM(x[,v$barcode],file.path(root,'objects/validation_counts.mtx'));writeLines(genes,file.path(root,'objects/validation_genes.txt'));writeLines(v$barcode,file.path(root,'objects/validation_cells.txt'))
# Frozen B identity, Naive/Memory, GC/IFN and contamination markers; use all B cells, no relabeling.
markers<-c('CD79A','MS4A1','CD19','EPCAM','KRT8','KRT19','IGHD','TCL1A','FCER2','CD27','TNFRSF13B','AIM2','BCL6','AICDA','RGS13','ISG15','STAT1','MX1','MKI67','TOP2A','UBE2C')
m<-fread(file.path(tt,'cells_pre_score.tsv.gz'))[major=='B'];present<-markers%in%rownames(obj[['RNA']]);fwrite(data.table(gene=markers,measured=present),file.path(tt,'B_marker_gene_coverage.tsv'),sep='\t')
b<-LayerData(obj,assay='RNA',layer='data',features=markers[present],cells=m$barcode);b<-b[markers[present],m$barcode]
m[,aggregate_id:=paste(sample_uid,subtype,sep='|||')];gr<-unique(m$aggregate_id);gi<-match(m$aggregate_id,gr);cnt<-tabulate(gi,nbins=length(gr));gm<-sparseMatrix(i=seq_len(nrow(m)),j=gi,x=1)
me<-as.matrix(b%*%gm%*%Diagonal(x=1/cnt));de<-as.matrix((b>0)%*%gm%*%Diagonal(x=1/cnt));meta<-m[,.(dataset=dataset[1],sample_uid=sample_uid[1],subtype=subtype[1],n_cells=.N),by=aggregate_id];rr<-rbindlist(lapply(seq_along(gr),function(j)data.table(meta[match(gr[j],aggregate_id)],gene=rownames(b),mean=me[,j],detection=de[,j])))
fwrite(rr,file.path(tt,'B_marker_sample_support.tsv'),sep='\t');writeLines(c('Frozen full background: 90011 genes; K4501; denominator89810',capture.output(sessionInfo())),file.path(root,'logs/score_session.txt'));cat('DONE',format(Sys.time()),'\n')
