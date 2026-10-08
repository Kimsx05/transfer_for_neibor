source('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1/code/00_common.R')
suppressPackageStartupMessages({library(ggplot2);library(ggrastr)})
M<-qread(out('02_ucell_scores/M_RAW.qs'));md<-qread(out('00_input_audit/cell_metadata.qs'));md[,MP10:=M[cell_id,'EpiMP10']];u<-qread(out('02_ucell_scores/original_usage_raw.qs'));md[,usage_MP10:=u[cell_id,'EpiMP10']]
th<-theme_classic(base_size=10)+theme(plot.background=element_rect(fill='transparent',colour=NA),panel.background=element_rect(fill='transparent',colour=NA),legend.background=element_rect(fill='transparent',colour=NA))
pal<-setNames(hcl.colors(uniqueN(md$dataset),'Dark 3'),sort(unique(md$dataset)))
savep<-function(p,name,w=8,h=6) {base<-out('07_mp10_review',name);ggsave(paste0(base,'.pdf'),p,width=w,height=h,bg='transparent',useDingbats=FALSE);ggsave(paste0(base,'.png'),p,width=w,height=h,dpi=300,bg='transparent')}
for(v in c('nFeature_counts','nCount_counts','percent.mt','usage_MP10')) {
 p<-ggplot(md,aes(.data[[v]],MP10,colour=dataset))+geom_point_rast(size=.15,alpha=.3,raster.dpi=300)+scale_colour_manual(values=pal)+labs(x=v,y='MP10 raw UCell')+th
 savep(p,paste0('FULL_MP10_vs_',v))
 savep(p+facet_wrap(~dataset,ncol=3),paste0('dataset_MP10_vs_',v),13,8)
}
for(v in c('dataset','source','annotation','original_cluster')) {
 p<-ggplot(md,aes(.data[[v]],MP10))+geom_violin(fill='#B8CBDD',linewidth=.2,scale='width')+geom_boxplot(width=.1,outlier.shape=NA,fill='white')+labs(x=NULL,y='MP10 raw UCell')+th+theme(axis.text.x=element_text(angle=30,hjust=1))
 savep(p,paste0('FULL_MP10_by_',v),max(8,uniqueN(md[[v]])*.35),6)
}
cm<-cor(M,method='spearman');write_tsv(as.data.table(cm,keep.rownames='MP'),out('07_mp10_review/score_spearman_matrix.tsv'))
z<-as.data.table(as.table(cm));setnames(z,c('MP1','MP2','rho'));p<-ggplot(z,aes(MP1,MP2,fill=rho))+geom_tile()+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),name='Spearman')+labs(x=NULL,y=NULL)+th+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5));savep(p,'FULL_14MP_score_correlation',9,8)
writeLines('complete',out('logs/COVARIATE_FIGURES_COMPLETE'))
