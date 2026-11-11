source(file.path(commandArgs(TRUE)[1],'code/plot_common.R'));z<-readRDS(file.path(out,'tmp/supplement_plot.rds'));d<-z$metadata;v<-z$log;ro<-z$rois
make<-function(genes,ii=seq_len(nrow(d))){rbindlist(lapply(genes,function(g)data.table(x_um=d$x_um[ii],y_um=d$y_um[ii],value=v[g,ii],gene=g)))}
if(s=='D1_Johnson'){
 for(k in names(z$panels)){p<-spfinish(ggplot(make(z$panels[[k]]),aes(x_um,-y_um,color=value))+geom_point_rast(size=.09,raster.dpi=300)+facet_wrap(~gene,ncol=2)+scale_color_viridis_c(name='RNA LogNormalize'));export(p,paste0('D1_whole_',k),11,if(length(z$panels[[k]])>4)14 else 10,9)}
 selected<-c('conflict_focus','reference_agree_control')
}else selected<-if(s=='Cervilla')c('far_UMI_hotspot','stroma') else ro$roi
for(nm in selected){b<-ro[roi==nm];ii<-which(d$x_um>=b$x&d$x_um<b$X&d$y_um>=b$y&d$y_um<b$Y);p<-spfinish(ggplot(make(rownames(v),ii),aes(x_um,-y_um,color=value))+geom_point_rast(size=.65,raster.dpi=300)+facet_wrap(~gene,ncol=4)+scale_color_viridis_c(name='RNA LogNormalize'));p<-roi_fixed(p,b);export(p,paste0('ROI_single_genes_',nm),14,13,12)}
m<-fread(file.path(out,'supplement_marker_statistics.tsv'));m<-m[context=='l08'];m[,group:=factor(group,levels=as.character(sort(as.integer(unique(group)))))];m[,symbol:=factor(symbol,levels=rownames(v))];p<-ggplot(m,aes(symbol,group,size=pct_positive*100,color=mean_log))+geom_point()+scale_size_area(max_size=4,limits=c(0,100),name='Detected (%)')+scale_color_viridis_c(name='Mean RNA\nlog expression')+jin_theme()+theme(axis.text.x=element_text(angle=60,hjust=1))+labs(x=NULL,y='lambda 0.8 cluster');export(p,'expanded_marker_clusters_l08',12,max(7,length(unique(m$group))*.22),10)
# Whole-section fixed expanded-panel maps for cross-section reading, unsmoothed measured expression.
w<-rbindlist(lapply(names(z$panels),function(k)data.table(x_um=d$x_um,y_um=d$y_um,value=colMeans(v[z$panels[[k]],,drop=FALSE]),panel=k)));p<-spfinish(ggplot(w,aes(x_um,-y_um,color=value))+geom_point_rast(size=.08,raster.dpi=300)+facet_wrap(~panel,nrow=1)+scale_color_viridis_c(name='Mean RNA\nlog expression'));export(p,'spatial_expanded_original_panels',18,6,14)
cat('SUPPLEMENT_PLOTS_COMPLETE',s,'\n')
