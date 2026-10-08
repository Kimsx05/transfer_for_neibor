# Finalize branch-specific categorical colors without changing any scores, labels or embeddings.
s <- paste(readLines('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1/code/04_figures.R'),collapse='\n')
prefix <- strsplit(s,'for(i in seq_len(nrow(inst))) {',fixed=TRUE)[[1]][1]
eval(parse(text=prefix))
stopifnot(file.exists(out('logs/FIGURES_COMPLETE')))
saveforce<-function(p,base,w=7,h=5) {ggsave(paste0(base,'.pdf'),p,width=w,height=h,bg='transparent',useDingbats=FALSE,limitsize=FALSE);ggsave(paste0(base,'.png'),p,width=w,height=h,dpi=300,bg='transparent',limitsize=FALSE)}
for(i in seq_len(nrow(inst))) {
 nm<-inst$instance[i];dir<-out(inst$directory[i]);logmsg('Finalize independent palette',nm)
 emb<-qread(file.path(dir,'umap.qs'));cl<-fread(file.path(dir,'clusters_res0.4.tsv.gz'))
 d<-merge(as.data.table(emb,keep.rownames='cell_id'),cl,by='cell_id',sort=FALSE);d<-merge(d,md[,.(cell_id,dataset,sample)],by='cell_id',sort=FALSE);xy<-d
 levels_cl<-paste0('C',sort(as.integer(sub('C','',unique(cl$cluster)))))
 pal<-setNames(grDevices::hcl(h=(seq(0,360,length.out=length(levels_cl)+1)[seq_along(levels_cl)]+37*i)%%360,c=60,l=65),levels_cl)
 write_tsv(data.table(cluster=names(pal),color=pal),file.path(dir,'cluster_colors.tsv'))
 extensive<-nm %in% c('FULL_Z','FULL_RAW','BAL_Z_seed42','BAL_RAW_seed42')
 for(ds in if(extensive)c('ALL',sort(unique(d$dataset))) else 'ALL') {
  dd<-if(ds=='ALL')copy(d) else d[dataset==ds];base<-if(ds=='ALL')file.path(dir,'figures')else out('08_figures_by_dataset',ds,nm);prefix<-file.path(base,paste0(nm,'__',ds,'__res0.4__'))
  saveforce(catplot(dd,'cluster',pal,xy),paste0(prefix,'cluster'))
  if(nm=='FULL_RAW' && ds=='ALL') {saveforce(continuous(dd,'EpiMP10',xy),paste0(prefix,'MP10'),8,6);saveforce(wrap_plots(lapply(mp,function(m) continuous(dd,m,xy)),ncol=4),paste0(prefix,'all14_raw_UCell'),20,16)}
  if(extensive) {
   dd[,MP10:=M[cell_id,'EpiMP10']];dd[,cluster:=factor(cluster,levels=levels_cl)]
   p<-ggplot(dd,aes(cluster,MP10,fill=cluster))+geom_violin(scale='width',linewidth=.15)+geom_boxplot(width=.12,outlier.shape=NA,fill='white',linewidth=.2)+scale_fill_manual(values=pal)+guides(fill='none')+labs(x=NULL,y='MP10 raw UCell')+th
   saveforce(p,paste0(prefix,'MP10_distribution'),max(8,length(levels_cl)*.45),5)
   ct<-dd[,.(n=.N),by=.(sample,cluster)];ct[,proportion:=n/sum(n),by=sample]
   p<-ggplot(ct,aes(sample,proportion,fill=cluster))+geom_col()+scale_fill_manual(values=pal)+coord_flip()+th+labs(x=NULL,y='Proportion within sampled set / FULL branch')
   saveforce(p,paste0(prefix,'sample_cluster_proportions'),10,max(5,uniqueN(ct$sample)*.20))
  }
 }
}
writeLines('complete',out('logs/PALETTE_FINALIZATION_COMPLETE'))
