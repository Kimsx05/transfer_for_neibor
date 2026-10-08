.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/refdir/Rlib','/usr/local/lib/R/library'))
Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1)
options(stringsAsFactors=FALSE, future.globals.maxSize=20*1024^3)
suppressPackageStartupMessages({library(qs);library(data.table);library(Matrix);library(Seurat);library(digest)})
data.table::setDTthreads(2)
root <- '/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1'
frozen <- '/home/data/t070721/Metabolistic_scRNA/Multi2/4.3 NMF allgene/output_full_cc_seed5_v2'
input <- fread(file.path(frozen,'01_input_summary.csv'))$input_qs
mp <- sprintf('EpiMP%02d',1:14)
out <- function(...) file.path(root,...)
write_tsv <- function(x,path) {dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);fwrite(x,path,sep='\t',na='NA')}
qsave_atomic <- function(x,path) {qsave(x,paste0(path,'.tmp'),nthreads=1);stopifnot(file.rename(paste0(path,'.tmp'),path))}
logmsg <- function(...) {cat(format(Sys.time()),..., '\n');flush.console()}
config <- list(version='mp_ucell_state_v1',frozen=frozen,input=input,basis='sample_balanced_mean',mp=mp,top_n=50L,tie_order='C locale gene name ascending',assay='RNA',layer='counts',background='all RNA features',maxRank=1500L,ties.method='average',missing_genes='impute',chunk.size=500L,workers=2L,smoothed=FALSE,standardization='FULL column mean and sample SD, no clipping',sample_cap=1000L,sampling_seeds=42:44,analysis_seed=42L,nn_method='annoy',n_trees=50L,k=20L,metric='euclidean',l2.norm=FALSE,prune.SNN=1/15,algorithm=4L,leiden_method='leidenbase',objective='RBConfigurationVertexPartition',n.iter=10L,group.singletons=TRUE,resolutions_full_z=c(.4,.2,.6),resolution_other=.4,umap='uwot',umap_n_neighbors=30L,umap_min_dist=.3,umap_n_components=2L,umap_threads=1L,pca=FALSE,harmony=FALSE)
config_hash <- digest(config,algo='sha256')
if(file.exists(out('config.qs'))) stopifnot(identical(qread(out('config.qs')),config)) else {qsave_atomic(config,out('config.qs'));dput(config,file=out('config.R'));writeLines(config_hash,out('config.sha256'))}
