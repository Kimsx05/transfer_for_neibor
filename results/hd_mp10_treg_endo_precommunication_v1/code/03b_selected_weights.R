source('highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1/code/common.R');library(spacexr)
for(unit in c('bin8','bin16')) {
 d<-rd(paste0(unit,'_rctd_results.tsv.gz'));w<-rd(paste0(unit,'_rctd_weights.tsv.gz'));rownames(d)<-d$bin_id;rownames(w)<-w$bin_id;w$bin_id<-NULL;w<-as.matrix(w[d$bin_id,,drop=FALSE]);types<-sort(unique(c(d$first_type,d$second_type)));b<-spacexr:::get_beta_doublet(d$bin_id,types,d,w)
 put(cbind(data.frame(bin_id=rownames(b)),as.data.frame(b,check.names=FALSE)),paste0(unit,'_native_selected_weights.tsv.gz'))
}
cat('NATIVE_SELECTED_WEIGHT_DECODE_COMPLETE\n')
