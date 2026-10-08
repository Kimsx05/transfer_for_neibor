source('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1/code/00_common.R')
suppressPackageStartupMessages({library(ggplot2);library(viridis)})
md<-qread(out('00_input_audit/cell_metadata.qs'))
th<-theme_classic(base_size=10)+theme(plot.background=element_rect(fill='transparent',colour=NA),panel.background=element_rect(fill='transparent',colour=NA),legend.background=element_rect(fill='transparent',colour=NA),axis.text.x=element_text(angle=30,hjust=1))
savep<-function(p,name,w=9,h=6) {base<-out('00_input_audit',name);ggsave(paste0(base,'.pdf'),p,width=w,height=h,bg='transparent',useDingbats=FALSE);ggsave(paste0(base,'.png'),p,width=w,height=h,bg='transparent',dpi=300)}
p<-ggplot(md,aes(dataset,nFeature_counts,fill=dataset))+geom_violin(scale='width',linewidth=.2)+geom_boxplot(width=.1,outlier.shape=NA,fill='white')+geom_hline(yintercept=1500,linetype=2)+guides(fill='none')+labs(x=NULL,y='Detected RNA features (counts)')+th
savep(p,'nFeature_by_dataset')
p<-ggplot(md,aes(sample,nFeature_counts,fill=dataset))+geom_boxplot(outlier.shape=NA,linewidth=.2)+coord_flip()+geom_hline(yintercept=1500,linetype=2)+labs(x=NULL,y='Detected RNA features (counts)')+th
savep(p,'nFeature_by_sample',12,13)
ov<-fread(out('01_signatures/overlap_jaccard_long.tsv'));p<-ggplot(ov,aes(MP1,MP2,fill=Jaccard))+geom_tile()+geom_text(aes(label=overlap),size=2.5)+scale_fill_viridis_c(option='plasma',limits=c(0,1))+labs(x=NULL,y=NULL)+th+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
savep(p,'signature_Jaccard_overlap',9,8)
det<-qread(out('02_ucell_scores/signature_detected_counts.qs'));d<-as.data.table(det);d[,dataset:=md$dataset];long<-melt(d,id.vars='dataset',variable.name='MP',value.name='detected');p<-ggplot(long,aes(dataset,detected,fill=dataset))+geom_boxplot(outlier.shape=NA,linewidth=.15)+facet_wrap(~MP,ncol=4)+guides(fill='none')+labs(x=NULL,y='Detected signature genes / 50')+th+theme(axis.text.x=element_text(angle=60,hjust=1,size=7))
savep(p,'signature_detection_by_dataset',14,12)
