.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
library(data.table);library(qs)
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1';tb<-function(n)file.path(root,'tables',n);a<-fread(tb('selected_axes.tsv'));cx<-fread(tb('CellChat_DB_complexes.tsv'));s<-fread(tb('samples.tsv'));c<-qread(file.path(root,'objects/sample_expression_support.qs'));out<-list();k<-0L
expand<-function(x){if(x%in%cx$complex){v<-unlist(cx[complex==x,-1,with=FALSE],use.names=FALSE);v[!is.na(v)&nzchar(v)]}else x}
for(i in seq_len(nrow(a)))for(role in c('ligand','receptor')){
 genes<-expand(a[[role]][i]);ct<-if(role=='ligand')a$source[i] else a$target[i]
 for(j in seq_len(nrow(s)))for(g in genes){
  col<-paste(s$sample_uid[j],ct,sep='|||');id<-match(col,c$metadata$aggregate_id);nc<-if(is.na(id))0 else c$metadata$n_cells[id];present<-g%in%rownames(c$mean)&&col%in%colnames(c$mean)
  k<-k+1L;out[[k]]<-data.table(candidate_subtype=a$candidate_subtype[i],LR_ID=a$LR_ID[i],role=role,cell_type=ct,gene=g,sample_uid=s$sample_uid[j],patient_uid=s$patient_uid_recorded[j],dataset=s$dataset[j],MP10_group=s$MP10_group[j],n_cells=nc,mean_expression=if(present)as.numeric(c$mean[g,col]) else NA_real_,detection=if(present)as.numeric(c$detection[g,col]) else NA_real_,status=if(nc==0)'population_absent' else if(nc<30)'low_cell_coverage_descriptive_only' else if(!present)'gene_unavailable' else 'evaluable')
 }
}
fwrite(rbindlist(out),tb('selected_LR_gene_sample_support.tsv'),sep='\t',na='NA')
