Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_THREAD_LIMIT=1)
.libPaths(c('/refdir/Rlib','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
library(data.table);setDTthreads(1)
root<-normalizePath('highres_VisiumHD/analysis/hd_treg_mp10_enrichment_v3');skill<-'/home/data/t070721/.agents/skills/jins-routine-plotting'
for(f in c('common/theme.R','common/palettes.R','common/export.R','plot_types/Violin_distribution/plot.R'))source(file.path(skill,f))
jin_setup();d<-fread(file.path(root,'tables/distribution_plot_data.tsv.gz'));sections<-c('Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2');d[,display:=paste(section,ifelse(group=='FOXP3-supported','FOXP3+','Other T'),sep='\n')];ord<-as.vector(rbind(paste(sections,'FOXP3+',sep='\n'),paste(sections,'Other T',sep='\n')));colors<-setNames(rep(c('#9c3985','#aeb8bf'),5),ord)
ix<-list()
for(v in c('Treg_main','log10_UMI','log1p_FOXP3','log1p_IL2RA','log1p_CTLA4')){
 p<-plot_violin_distribution(d,'display',v,method=v,order=ord,colors=colors,ylab=v)+theme(axis.text.x=element_text(angle=45,hjust=1))
 name<-paste0('distribution_',v);jin_export(p,file.path(root,'figures',name),11,4.8);ix[[length(ix)+1]]<-data.frame(figure=name,description='All T-eligible; FOXP3-supported versus remaining evaluable T. Frozen values; no new scores or tests.',PNG=paste0(name,'.png'),source='distribution_plot_data.tsv.gz; jins-routine-plotting Violin_distribution')
}
fwrite(rbindlist(ix),file.path(root,'tables/figure_index_R.tsv'),sep='\t');writeLines(c(paste('R',getRversion()),paste('ggplot2',packageVersion('ggplot2')),paste('ragg',packageVersion('ragg'))),file.path(root,'logs/R_plot_versions.txt'))
