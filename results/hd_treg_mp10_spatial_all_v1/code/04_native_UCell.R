source('highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1/code/common.R')
root<-file.path(P,'highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1');library(UCell);if(Sys.getenv('HD_ALLOW_DAYTIME')!='1')gate()
cfg<-yaml::read_yaml(file.path(root,'frozen_config.yaml'));sig<-cfg$signatures;M<-cfg$UCell_maxRank
x<-as(readMM(file.path(root,'temporary/native_UCell_check.mtx')),'dgCMatrix');rownames(x)<-readLines(file.path(root,'temporary/native_UCell_genes.txt'));colnames(x)<-readLines(file.path(root,'temporary/native_UCell_ids.txt'))
a<-ScoreSignatures_UCell(x,features=sig,maxRank=M,name='',ncores=1,chunk.size=100,ties.method='average',missing_genes='skip')
a<-as.data.frame(a);a$region_id<-rownames(a);fwrite(a,file.path(root,'tables/UCell_native_raw_check.tsv'),sep='\t')
# Staged native UCell2.14 core. Public precalc.ranks infers maxRank=max(rank)+1,
# which is ambiguous for restricted gene caches and half-integer tied ranks.
# Call its exact native scoring core with explicit frozen cap, on full-background ranks.
files<-list.files(file.path(root,'temporary'),pattern='_ucell_ranks.mtx$',full.names=TRUE)
for(f in files){
 q<-as(readMM(f),'dgCMatrix');stem<-sub('_ucell_ranks.mtx$','',f);rownames(q)<-readLines(paste0(stem,'_ucell_genes.txt'));colnames(q)<-readLines(paste0(stem,'_ucell_ids.txt'))
 zz<-list()
 for(ii in split(seq_len(ncol(q)),ceiling(seq_len(ncol(q))/500))){
  dense<-as.data.table(as.matrix(q[,ii,drop=FALSE]),keep.rownames=TRUE);setkey(dense,'rn',physical=FALSE)
  z<-UCell:::u_stat_signature_list(sig,dense,maxRank=M,sparse=TRUE,missing_genes='skip',w_neg=1)
  colnames(z)<-names(sig);z<-as.data.frame(z);z$region_id<-colnames(q)[ii];zz[[length(zz)+1]]<-z
 }
 fwrite(rbindlist(zz),file.path(root,'tables',paste0(basename(stem),'_native_UCell.tsv.gz')),sep='\t');cat('SCORED',basename(stem),ncol(q),'\n');flush.console()
}
writeLines(capture.output(sessionInfo()),file.path(root,'logs/UCell_sessionInfo.txt'))
