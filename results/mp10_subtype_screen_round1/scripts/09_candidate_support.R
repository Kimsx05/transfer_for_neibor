.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1';tb<-function(n)file.path(root,'tables',n)
a<-fread(tb('CellChat_epithelial_comparisons.tsv.gz'));db<-fread(tb('CellChat_DB_interactions.tsv'));cx<-fread(tb('CellChat_DB_complexes.tsv'));s<-fread(tb('samples.tsv'));cache<-qread(file.path(root,'objects/sample_expression_support.qs'));meta<-cache$metadata
# Candidate discovery keeps both directions and the union across all datasets.
u<-unique(a[comparison=='MP10_low_high' & (strength0>0|strength1>0),.(source,target,LR_ID,pathway,ligand,receptor,category,LR_display)])
u[,candidate_subtype:=ifelse(source=='Malignant_like_epithelial',target,source)]
# Exclude epithelial/technical labels from PRIORITY selection only; retained in full network outputs.
u[,priority_interpretable:=!grepl('epithelial|dataset_bias|excluded_|unresolved|boundary',candidate_subtype)]
# A small, recorded axis shortlist per subtype; all networks and candidate directions remain in the full tables.
ranked<-a[comparison=='MP10_low_high' & scope!='ALL',.(n_dataset_evaluable=sum(is.finite(difference)),n_dataset_up=sum(difference>0,na.rm=TRUE),n_dataset_down=sum(difference<0,na.rm=TRUE),median_abs_network_difference=median(abs(difference),na.rm=TRUE)),by=.(source,target,LR_ID)]
u<-merge(u,ranked,by=c('source','target','LR_ID'),all.x=TRUE)
u[,max_direction_agreement:=pmax(n_dataset_up,n_dataset_down)]
setorderv(u,c('candidate_subtype','max_direction_agreement','n_dataset_evaluable','median_abs_network_difference','LR_ID'),c(1,-1,-1,-1,1))
fwrite(u,tb('candidate_discovery_all_axes.tsv'),sep='\t')
u<-u[,head(.SD,3),by=c('candidate_subtype')]
expand_gene<-function(g){if(g%in%cx$complex){v<-unlist(cx[complex==g,-1,with=FALSE],use.names=FALSE);v[!is.na(v)&nzchar(v)]}else g}
rows<-list();gene_rows<-list();k<-0L
for(i in seq_len(nrow(u))){
 z<-u[i];lg<-expand_gene(z$ligand);rg<-expand_gene(z$receptor)
 for(j in seq_len(nrow(s))){
  sm<-s[j];ls<-paste(sm$sample_uid,z$source,sep='|||');rs<-paste(sm$sample_uid,z$target,sep='|||');li<-match(ls,meta$aggregate_id);ri<-match(rs,meta$aggregate_id)
  nl<-if(is.na(li))0L else meta$n_cells[li];nr<-if(is.na(ri))0L else meta$n_cells[ri]
  evaluable<-nl>=30&&nr>=30&&all(c(lg,rg)%in%rownames(cache$mean))
  lm<-ld<-rm<-rd<-NA_real_
  if(evaluable){
   lm<-exp(mean(log(as.numeric(cache$mean[lg,ls]))));rm<-exp(mean(log(as.numeric(cache$mean[rg,rs]))));ld<-min(as.numeric(cache$detection[lg,ls]));rd<-min(as.numeric(cache$detection[rg,rs]))
  }
  k<-k+1L;rows[[k]]<-data.table(LR_ID=z$LR_ID,source=z$source,target=z$target,candidate_subtype=z$candidate_subtype,sample_uid=sm$sample_uid,dataset=sm$dataset,patient_uid=sm$patient_uid_recorded,tissue_status=sm$tissue_status,T_category=sm$T_category,MP10_group=sm$MP10_group,MP10_sample=sm$MP10_sample,n_sender=nl,n_receiver=nr,evaluable=evaluable,ligand_mean=lm,receptor_mean=rm,ligand_min_detection=ld,receptor_min_detection=rd,expression_supported=evaluable&&ld>=.1&&rd>=.1)
 }
}
r<-rbindlist(rows);fwrite(r,tb('candidate_union_sample_expression_support.tsv.gz'),sep='\t',na='NA')
fwrite(u,tb('candidate_union_axes.tsv'),sep='\t')
cat('Sample checks for',nrow(u),'axes and',nrow(r),'axis-specimens\n')
