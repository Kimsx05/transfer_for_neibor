.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1';tab<-function(n)file.path(root,'tables',n)
files<-list.files(file.path(root,'objects/networks'),pattern='qs$',full.names=TRUE);nets<-setNames(lapply(files,qread),sub('.qs$','',basename(files)))
db<-fread(tab('CellChat_DB_interactions.tsv'));map<-fread(tab('subtype_parent.tsv'));major_map<-setNames(map$major,map$subtype);major_map[c('Malignant_like_epithelial','Normal_epithelial','Other_epithelial_unresolved')]<-'Epithelial'
edge_rows<-list();major_rows<-list();coverage<-list()
for(key in names(nets)){
 n<-nets[[key]];idx<-which(n$net$prob>0 & n$net$pval<.05,arr.ind=TRUE);dn<-dimnames(n$net$prob)
 if(nrow(idx)){
  e<-data.table(source=dn[[1]][idx[,1]],target=dn[[2]][idx[,2]],LR_ID=dn[[3]][idx[,3]],probability=n$net$prob[idx],internal_P=n$net$pval[idx],scope=n$scope,comparison=n$comparison,condition=n$condition)
  edge_rows[[key]]<-merge(e,db[,.(LR_ID=interaction_name,pathway=pathway_name,ligand,receptor,category=annotation,LR_display=interaction_name_2)],by='LR_ID',all.x=TRUE)
  e[,`:=`(source_major=unname(major_map[source]),target_major=unname(major_map[target]))]
  major_rows[[key]]<-e[,.(interaction_count=.N,total_strength=sum(probability)),by=.(scope,comparison,condition,source_major,target_major)]
 }
 cov<-copy(n$counts);cov[,`:=`(scope=n$scope,comparison=n$comparison,condition=n$condition,evaluable=N>=30,data_max=n$params$data_max)];coverage[[key]]<-cov
}
edges<-rbindlist(edge_rows);fwrite(edges,tab('CellChat_detected_edges.tsv.gz'),sep='\t');fwrite(rbindlist(major_rows),tab('CellChat_major_aggregated.tsv'),sep='\t');fwrite(rbindlist(coverage),tab('CellChat_group_coverage.tsv'),sep='\t')
# Per-dataset edge comparisons use the UNION of LR/source/target discoveries.
# An LR removed by the gene-overexpression screen is distinct from an absent population.
read_edge<-function(net,src,tgt,lr){
 if(is.null(net))return(list(probability=NA_real_,internal_P=NA_real_,supported_strength=NA_real_,status='missing_condition'))
 dn<-dimnames(net$net$prob)
 if(!src%in%dn[[1]]||!tgt%in%dn[[2]])return(list(probability=NA_real_,internal_P=NA_real_,supported_strength=NA_real_,status='source_or_target_below_30'))
 if(!lr%in%dn[[3]])return(list(probability=NA_real_,internal_P=NA_real_,supported_strength=0,status='not_detected_LR_screen'))
 pr<-net$net$prob[src,tgt,lr];pv<-net$net$pval[src,tgt,lr]
 list(probability=pr,internal_P=pv,supported_strength=if(pv<.05)pr else 0,status=if(pr>0&&pv<.05)'detected' else 'not_detected_probability')
}
co_groups<-list(Normal_Tumor=c('Normal','Tumor'),T1_T2=c('T1','T2'),MP10_low_high=c('low','high'))
allcmp<-list();k<-0L
for(co in names(co_groups)){
 u<-unique(edges[comparison==co & (source=='Malignant_like_epithelial'|target=='Malignant_like_epithelial'),.(source,target,LR_ID,pathway,ligand,receptor,category,LR_display)])
 for(ds in c('ALL',sort(unique(edges[scope!='ALL',scope])))){
  n0<-nets[[paste(ds,co,co_groups[[co]][1],sep='__')]];n1<-nets[[paste(ds,co,co_groups[[co]][2],sep='__')]]
  a<-rbindlist(lapply(seq_len(nrow(u)),function(i){z<-u[i];l<-read_edge(n0,z$source,z$target,z$LR_ID);h<-read_edge(n1,z$source,z$target,z$LR_ID);data.table(probability0=l$probability,probability1=h$probability,P0=l$internal_P,P1=h$internal_P,strength0=l$supported_strength,strength1=h$supported_strength,status0=l$status,status1=h$status)}))
  if(nrow(a)){k<-k+1L;allcmp[[k]]<-cbind(u,a)[,`:=`(scope=ds,comparison=co,condition0=co_groups[[co]][1],condition1=co_groups[[co]][2],difference=a$strength1-a$strength0)]}
 }
}
fwrite(rbindlist(allcmp,fill=TRUE),tab('CellChat_epithelial_comparisons.tsv.gz'),sep='\t',na='NA')
cat('Exported',length(nets),'networks',nrow(edges),'detected edges\n')
# Read-back of saved numerical networks supersedes post-save logging errors from
# the initial status-column name collision (data.table reserves `key`).
status<-rbindlist(lapply(names(nets),function(k){n<-nets[[k]];stopifnot(identical(dim(n$net$prob),dim(n$net$pval)),all(is.finite(n$net$prob)),all(n$net$prob>=0),all(n$net$pval>=0&n$net$pval<=1),n$params$min.cells==30,n$params$nboot==100,!n$params$population.size);data.table(network_id=k,scope=n$scope,comparison=n$comparison,condition=n$condition,status='completed_readback_verified',n_cells=sum(n$contributions$n_cells),n_types=dim(n$net$prob)[1],n_LR=dim(n$net$prob)[3],data_max=n$params$data_max)}))
fwrite(status,tab('CellChat_network_status.tsv'),sep='\t')
