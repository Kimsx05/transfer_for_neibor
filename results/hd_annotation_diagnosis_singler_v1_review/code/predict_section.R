source(file.path(commandArgs(TRUE)[1],'code/common.R'));s<-args[2];out<-file.path(root,s);dfile<-file.path(out,'diagnosis.tsv.gz');diag<-fread(dfile,data.table=FALSE);x<-load_major(file.path(old,s));ids<-x@misc$analysis_cell_ids;stopifnot(identical(ids,diag$cell_id),nrow(x[['RNA']])==18085,all(x$qc_keep_analysis),sum(x$nCount_RNA)==sum(diag$nCount_RNA),all(x$final_major==diag$final_major))
q<-LayerData(x,assay='RNA',layer='data');counts<-LayerData(x,assay='RNA',layer='counts');sym<-x[['RNA']][[]]$gene_symbol;stopifnot(identical(colnames(q),ids));input<-list(section=s,n=length(ids),genes=nrow(q),UMI=sum(counts),RNA_layers=Layers(x[['RNA']]),input_unchanged=TRUE)
# Only in-memory symbol view; if duplicates exist, sum counts before normalization.
sym[is.na(sym)|sym=='']<-rownames(q)[is.na(sym)|sym==''];dup<-sum(duplicated(sym));levels<-unique(sym)
if(dup>0){G<-sparseMatrix(i=match(sym,levels),j=seq_along(sym),x=1,dims=c(length(levels),length(sym)));qs<-G%*%counts;rownames(qs)<-levels;lib<-colSums(qs);qs@x<-log1p(qs@x*rep(10000/lib,diff(qs@p)))}else{qs<-q;rownames(qs)<-sym}
input$duplicate_symbol_rows_collapsed_in_temporary_view<-dup;write_json(input,file.path(out,'input_check.json'),auto_unbox=TRUE,pretty=TRUE)
panels<-readRDS(file.path(old,'reference/marker_panels.rds'));coverage<-list()
for(route in if(length(args)>=3)args[3] else c('blca','hpca')) {
 prefix<-paste0('singler_',route,'_');scorefile<-file.path(out,paste0(prefix,'scores.rds'))
 trained<-rebuildIndex(readRDS(file.path(root,'reference',paste0(route,'_trained.rds'))),num.threads=4)
 query<-if(route=='blca')q else qs;common<-intersect(rownames(query),rownames(trained$ref));query<-query[trained$options$test.genes,,drop=FALSE];used<-trained$markers$unique;stopifnot(all(used%in%rownames(query)))
 syms<-if(route=='blca')sym[match(common,rownames(q))] else common
 for(k in names(panels))coverage[[length(coverage)+1]]<-data.table(route=route,major=k,requested=length(panels[[k]]),shared_markers=sum(panels[[k]]%in%syms),present=paste(intersect(panels[[k]],syms),collapse=','),missing=paste(setdiff(panels[[k]],syms),collapse=','),all_shared_genes=length(common),model_selected_genes=length(used))
 if(paste0(prefix,'status')%in%colnames(diag)&&file.exists(scorefile)) {stopifnot(identical(rownames(readRDS(scorefile)$scores),ids));next}
 work<-query[used,,drop=FALSE];ss<-colSums(work);ss2<-work;ss2@x<-ss2@x^2;variance<-colSums(ss2)-ss^2/nrow(work);valid<-is.finite(variance)&variance>1e-12;rm(work,ss2);chunks<-split(which(valid),ceiling(seq_len(sum(valid))/5000));results<-vector('list',length(chunks));
 for(j in seq_along(chunks)){ci<-chunks[[j]];cat('CLASSIFY',s,route,j,length(chunks),'\n');flush.console();results[[j]]<-classifySingleR(query[,ci,drop=FALSE],trained,prune=FALSE,num.threads=4,BPPARAM=SerialParam());gc(FALSE)}
 z<-do.call(rbind,results);stopifnot(identical(rownames(z),ids[valid]));bad<-pruneScores(z);z$pruned.labels<-z$labels;z$pruned.labels[bad]<-NA_character_;dm<-getDeltaFromMedian(z);thresholds<-pruneScores(z,get.thresholds=TRUE)
 score<-matrix(NA_real_,length(ids),ncol(z$scores),dimnames=list(ids,colnames(z$scores)));score[valid,]<-z$scores
 saveRDS(list(scores=score,meaning='SingleR pre-fine-tuning correlation scores, not probabilities; not necessarily maximized by final labels'),scorefile,compress='gzip')
 put<-function(name,value){v<-rep(NA,length(ids));v[valid]<-value;diag[[paste0(prefix,name)]]<<-v}
 put('labels',z$labels);put('pruned_labels',z$pruned.labels);put('delta_next',z$delta.next);put('delta_median',dm);put('pruned',bad);put('initial_top_label',colnames(z$scores)[max.col(z$scores,ties.method='first')]);put('selected_initial_score',z$scores[cbind(seq_len(nrow(z)),match(z$labels,colnames(z$scores)))])
 diag[[paste0(prefix,'status')]]<-ifelse(valid,'evaluated','not_assessable_constant_expression_on_model_genes');diag[[paste0(prefix,'detected_model_genes')]]<-colSums(query[used,,drop=FALSE]>0)
 fwrite(data.table(label=names(thresholds),delta_median_threshold=as.numeric(thresholds)),file.path(out,paste0(route,'_pruning_thresholds.tsv')),sep='\t')
 write_json(list(route=route,section=s,all_objects=length(ids),evaluated=sum(valid),not_assessable=sum(!valid),pruned=sum(bad),nonpruned=sum(!bad),shared_genes=length(common),selected_model_genes=length(used),num_threads=4,chunk_size=5000,prune_once_after_merging=TRUE,nmads=3,min_diff_med='-Inf',min_diff_next=0,score_top_not_final=sum(diag[[paste0(prefix,'labels')]][valid]!=diag[[paste0(prefix,'initial_top_label')]][valid])),file.path(out,paste0(route,'_prediction_manifest.json')),auto_unbox=TRUE,pretty=TRUE)
 fwrite(diag,dfile,sep='\t');rm(z,results,score,trained);gc();cat('ROUTE_DONE',s,route,'\n')
}
fwrite(rbindlist(coverage),file.path(out,'gene_marker_coverage.tsv'),sep='\t');writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'));cat('SECTION_PREDICTION_COMPLETE',s,'\n')
