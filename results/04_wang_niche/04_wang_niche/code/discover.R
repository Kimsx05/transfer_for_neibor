#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
           VECLIB_MAXIMUM_THREADS='1', NUMEXPR_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
suppressPackageStartupMessages({library(scran);library(igraph);library(cluster);library(data.table);library(jsonlite)})
data.table::setDTthreads(1L)
args<-commandArgs(trailingOnly=TRUE);folder<-args[[1]];root<-dirname(dirname(folder))
params<-fromJSON(file.path(root,'parameters.json'));sig<-fromJSON(file.path(folder,'input_signature.json'))
seed<-as.integer(params$seed)
stamp<-function(...) cat(format(Sys.time(),'%FT%T%z'),paste(...,collapse=' '),'\n')
rss<-function(){z<-readLines('/proc/self/status');as.numeric(gsub('[^0-9.]','',z[startsWith(z,'VmHWM:')]))/1024^2}
write_json(list(R=R.version.string,packages=setNames(lapply(c('scran','bluster','igraph','cluster','data.table',
    'BiocNeighbors','BiocParallel','ggplot2','jsonlite','digest'),function(p)as.character(packageVersion(p))),
    c('scran','bluster','igraph','cluster','data.table','BiocNeighbors','BiocParallel','ggplot2','jsonlite','digest')),
    library_paths=.libPaths(),threads=1L),file.path(folder,'software_versions_R.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(getMethod('buildSNNGraph','ANY'),getFromNamespace('.setup_knn_data','scran'),
    bluster::makeSNNGraph,igraph::cluster_louvain,file=file.path(folder,'installed_interfaces.txt'))
d<-fread(file.path(folder,'composition_features.csv.gz'),check.names=FALSE)
ids<-d[[1]];stopifnot(identical(names(d)[-1],params$factor_order))
x<-as.matrix(d[,-1]);storage.mode(x)<-'double';rownames(x)<-ids
stopifnot(ncol(x)==30L,!anyNA(x),all(is.finite(x)),all(x>=0),max(abs(rowSums(x)-1))<1e-12)
evaluation<-fread(file.path(folder,'silhouette_evaluation_spots.tsv'))
ev<-evaluation$input_row_1based;stopifnot(identical(ids[ev],evaluation$spot_key),length(ev)<=5000L)
stamp('START',sig$dataset,'n=',nrow(x),'p=',ncol(x),'evaluation=',length(ev))
started<-proc.time()[3]
# Only this <=5000-spot Euclidean distance is ever formed, not all dataset pairs.
dist_ev<-dist(x[ev,,drop=FALSE],method='euclidean')
grid<-list();sil_rows<-list();size_rows<-list();graph_rows<-list();members<-list();i<-0L
worker_peaks<-numeric()
for(k in as.integer(params$k)){
  if(nrow(x)<=k){stamp('SKIP insufficient N for k=',k);next}
  set.seed(seed)
  stamp('GRAPH k=',k)
  graph<-scran::buildSNNGraph(x,k=k,d=NA_integer_,transposed=TRUE,type='rank',
      BNPARAM=BiocNeighbors::KmknnParam(distance='Euclidean'),num.threads=1L,
      BPPARAM=BiocParallel::SerialParam())
  w<-igraph::E(graph)$weight
  stopifnot(!igraph::is_directed(graph),igraph::vcount(graph)==nrow(x),length(w)>0,
      all(is.finite(w)),all(w>0))
  graph_rows[[length(graph_rows)+1L]]<-data.table(dataset=sig$dataset,k=k,vertices=vcount(graph),
      edges=ecount(graph),connected_components=components(graph)$no,
      weight_min=min(w),weight_max=max(w),weight_mean=mean(w),directed=is_directed(graph),
      SNN_weighting='rank',input_orientation='spots x factors; transposed=TRUE',d='NA_integer_',
      features=30L,metric='Euclidean',NN_engine='KmknnParam exact',threads=1L,
      graph_sha256=digest::digest(list(edges=as_edgelist(graph,names=FALSE),weights=w),algo='sha256'))
  fit_candidate<-function(res){
    key<-sprintf('k%d_r%.1f',k,res);set.seed(seed)
    fit<-igraph::cluster_louvain(graph,weights=w,resolution=res)
    lab<-as.integer(membership(fit))
    sizes<-table(lab);eval_sizes<-table(factor(lab[ev],levels=names(sizes)))
    reason<-'';value<-NA_real_;s<-NULL;ss<-NULL
    if(length(sizes)<2L)reason<-'ALL_ONE_CLASS'
    else if(all(sizes==1L))reason<-'ALL_SINGLETONS'
    else if(length(unique(lab[ev]))<2L)reason<-'EVALUATION_ONE_CLASS'
    else if(length(unique(lab[ev]))==length(ev))reason<-'EVALUATION_ALL_SINGLETONS'
    else {
      s<-tryCatch(cluster::silhouette(lab[ev],dist_ev),error=function(e)NULL)
      if(is.null(s)||any(!is.finite(s[,'sil_width'])))reason<-'INVALID_SILHOUETTE'
      else value<-mean(s[,'sil_width'])
    }
    if(!is.null(s)){
      sw<-data.table(raw_cluster=lab[ev],silhouette=s[,'sil_width'])
      ss<-sw[,.(n_evaluation=.N,mean_silhouette=mean(silhouette),median_silhouette=median(silhouette),
         min_silhouette=min(silhouette),n_negative=sum(silhouette<0)),by=raw_cluster]
      ss<-merge(data.table(raw_cluster=as.integer(names(sizes)),n_full=as.integer(sizes)),ss,
         by='raw_cluster',all.x=TRUE,sort=TRUE)
      ss[,`:=`(dataset=sig$dataset,candidate=key,k=k,resolution=res)]
    }
    size_row<-data.table(dataset=sig$dataset,candidate=key,k=k,resolution=res,
        raw_cluster=as.integer(names(sizes)),n_spots=as.integer(sizes),n_evaluation=as.integer(eval_sizes))
    grid_row<-data.table(dataset=sig$dataset,candidate=key,k=k,resolution=res,n_valid=nrow(x),
        n_evaluation=length(ev),n_clusters=length(sizes),n_clusters_evaluated=sum(eval_sizes>0),
        n_singletons=sum(sizes==1),min_cluster_size=min(sizes),median_cluster_size=median(sizes),
        max_cluster_size=max(sizes),mean_silhouette=value,modularity=modularity(fit),
        eligible=is.finite(value),exclusion_reason=reason,evaluation_mode=sig$evaluation_mode,
        seed=seed,status=if(is.finite(value))'PASS' else 'SKIPPED_NOT_ESTIMABLE')
    stamp('CANDIDATE',key,'clusters=',length(sizes),'silhouette=',round(value,6),'RSS_GiB=',round(rss(),3))
    list(key=key,lab=lab,silhouette=ss,size=size_row,grid=grid_row,worker_peak=rss())
  }
  # One dataset at a time; four candidate workers, each one BLAS/OpenMP thread.
  # Every worker explicitly resets the SAME seed, preserving serial candidate results.
  results<-parallel::mclapply(params$resolutions,fit_candidate,mc.cores=as.integer(params$candidate_workers),
      mc.preschedule=TRUE,mc.set.seed=FALSE,mc.allow.recursive=FALSE)
  for(result in results){
    if(inherits(result,'try-error'))stop('Candidate worker failed: ',result)
    i<-i+1L;members[[result$key]]<-result$lab;grid[[i]]<-result$grid;size_rows[[i]]<-result$size
    if(!is.null(result$silhouette))sil_rows[[length(sil_rows)+1L]]<-result$silhouette
    worker_peaks<-c(worker_peaks,result$worker_peak)
  }
  rm(graph,w,results,fit_candidate);gc()
}
g<-rbindlist(grid,fill=TRUE)
fwrite(g,file.path(folder,'silhouette_grid.tsv'),sep='\t',na='NA')
fwrite(rbindlist(sil_rows,fill=TRUE),file.path(folder,'per_cluster_silhouette_grid.tsv'),sep='\t',na='NA')
fwrite(rbindlist(size_rows),file.path(folder,'cluster_size_grid.tsv'),sep='\t')
fwrite(rbindlist(graph_rows),file.path(folder,'graph_diagnostics.tsv'),sep='\t')
fwrite(as.data.table(c(list(spot_key=ids),members)),file.path(folder,'all_candidate_memberships.csv.gz'),compress='gzip')
eligible_grid<-g[eligible==TRUE]
if(!nrow(eligible_grid)){
  write_json(list(dataset=sig$dataset,status='SKIPPED_NOT_ESTIMABLE',reason='No valid candidate',
      input_signature=sig$signature),file.path(folder,'discovery_status.json'),pretty=TRUE,auto_unbox=TRUE)
  quit(status=0L)
}
best<-max(eligible_grid$mean_silhouette)
ties<-eligible_grid[mean_silhouette>=best-params$tie_tolerance]
setorder(ties,resolution,k);chosen<-ties[1]
chosen[,`:=`(maximum_candidate_silhouette=best,tie_tolerance=params$tie_tolerance,
    n_tied_candidates=nrow(ties),tie_rule='lower resolution then lower k',
    features='all 30 untransformed subtype mean proportions',SNN_weighting='rank',
    algorithm='igraph::cluster_louvain',distance='Euclidean',input_signature=sig$signature)]
lab<-members[[chosen$candidate]]
raw<-sort(unique(lab));mapping<-data.table(raw_cluster=raw,n_spots=as.integer(table(lab)[as.character(raw)]),
    first_spot_key=vapply(raw,function(l)min(ids[lab==l]),character(1)))
setorder(mapping,-n_spots,first_spot_key);mapping[,niche_id:=sprintf('%s_N%02d',sig$dataset,seq_len(.N))]
labels<-mapping$niche_id[match(lab,mapping$raw_cluster)]
fwrite(data.table(spot_key=ids,raw_louvain_cluster=lab,niche_id=labels),file.path(folder,'selected_assignments.tsv'),sep='\t')
fwrite(mapping,file.path(folder,'niche_id_map.tsv'),sep='\t')
fwrite(chosen,file.path(folder,'selected_parameters.tsv'),sep='\t')
selected_sil<-rbindlist(sil_rows,fill=TRUE)[candidate==chosen$candidate]
selected_sil<-merge(selected_sil,mapping[,.(raw_cluster,niche_id)],by='raw_cluster',all.x=TRUE)
fwrite(selected_sil,file.path(folder,'per_cluster_silhouette.tsv'),sep='\t',na='NA')
write_json(list(dataset=sig$dataset,status='PASS',candidate=chosen$candidate,input_signature=sig$signature,
    n_valid=nrow(x),n_niches=length(raw),n_evaluation=length(ev),elapsed_seconds=unname(proc.time()[3]-started),
    peak_RSS_GiB=rss(),max_candidate_worker_RSS_GiB=max(worker_peaks),
    conservative_sum_peak_RSS_GiB=rss()+params$candidate_workers*max(worker_peaks),
    memory_measurement='VmHWM per process; sum is conservative and double-counts fork-shared pages',
    candidate_workers=params$candidate_workers,dataset_jobs=1L,BLAS_OpenMP_threads=1L,
    evaluation_mode=sig$evaluation_mode),file.path(folder,'discovery_status.json'),
    pretty=TRUE,auto_unbox=TRUE)
stamp('COMPLETE',sig$dataset,chosen$candidate,'niches=',length(raw),'seconds=',round(proc.time()[3]-started,1))
