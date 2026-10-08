source('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1/code/00_common.R')
checks<-list();check<-function(name,value,detail='') {checks[[length(checks)+1]]<<-data.table(check=name,pass=isTRUE(value),detail=as.character(detail));if(!isTRUE(value)) warning(name)}
M<-qread(out('02_ucell_scores/M_RAW.qs'));Z<-qread(out('02_ucell_scores/M_Z.qs'));md<-qread(out('00_input_audit/cell_metadata.qs'));st<-qread(file.path(frozen,'01_input_state.qs'))
check('FULL_exact_frozen_ID_order',identical(rownames(M),st$cells));check('14_columns_in_fixed_order',identical(colnames(M),mp));check('finite_RAW_Z',all(is.finite(M))&&all(is.finite(Z)));check('exact_metadata_ID_order',identical(md$cell_id,rownames(M)))
check('Z_FULL_means_zero',max(abs(colMeans(Z)))<1e-12);check('Z_FULL_sample_SD_one',max(abs(apply(Z,2,sd)-1))<1e-12)
sc<-fread(out('02_ucell_scores/global_scaling.tsv'));expected<-M;for(j in 1:14) expected[,j]<-(M[,j]-sc$mean_FULL[j])/sc$sd_FULL[j]
check('Z_reconstructs_from_saved_scaling',max(abs(expected-Z))<1e-12)
inst<-fread(out('instance_manifest.tsv'));nlabels<-0L;runstats<-list()
for(i in seq_len(nrow(inst))) {
 nm<-inst$instance[i];dir<-out(inst$directory[i]);ids<-if(startsWith(nm,'FULL')) rownames(M) else qread(out('05_balanced',paste0('seed',sub('.*seed','',nm),'_cell_ids.qs')))
 u<-qread(file.path(dir,'umap.qs'));check(paste(nm,'UMAP IDs'),identical(rownames(u),ids));check(paste(nm,'UMAP finite 2D'),ncol(u)==2&&all(is.finite(u)))
 identity<-qread(file.path(dir,'instance_identity.qs'));X<-if(grepl('_RAW',nm)) M[ids,mp,drop=FALSE] else Z[ids,mp,drop=FALSE]
 check(paste(nm,'direct unchanged 14D matrix'),identical(identity$matrix_hash,digest(X,'sha256')))
 rr<-if(nm=='FULL_Z') c(.4,.2,.6) else .4
 for(res in rr) {cl<-fread(file.path(dir,paste0('clusters_res',res,'.tsv.gz')));check(paste(nm,res,'labels IDs'),identical(cl$cell_id,ids)&&!anyDuplicated(cl$cell_id));nlabels<-nlabels+1L;runstats[[length(runstats)+1]]<-data.table(instance=nm,resolution=res,n_cells=nrow(cl),n_clusters=uniqueN(cl$cluster),smallest_cluster=min(table(cl$cluster)),largest_cluster=max(table(cl$cluster)))}
 if(startsWith(nm,'BAL')) {ct<-table(md$sample[match(ids,md$cell_id)]);full<-table(md$sample);check(paste(nm,'specimen cap'),all(ct[names(full)]==pmin(full,1000)));check(paste(nm,'FULL original order'),identical(ids,rownames(M)[rownames(M)%in%ids]))}
}
check('8_UMAP_instances_10_labels',nrow(inst)==8&&nlabels==10)
f<-fread(out('00_input_audit/input_manifest.tsv'),colClasses=c(mtime='character'));stat<-file.info(f$path);check('original_inputs_size_mtime_unchanged',all(f$bytes==stat$size)&all(f$mtime==as.character(stat$mtime)))
check('figures_finished',file.exists(out('logs/FIGURES_COMPLETE')))
write_tsv(rbindlist(checks),out('00_input_audit/final_validation.tsv'));write_tsv(rbindlist(runstats),out('06_comparisons/run_summary.tsv'))
stopifnot(all(rbindlist(checks)$pass))
files<-list.files(out('code'),full.names=TRUE);files<-files[!file.info(files)$isdir];write_tsv(data.table(path=files,md5=unname(tools::md5sum(files))),out('00_input_audit/code_checksums.tsv'))
writeLines('complete',out('logs/VALIDATION_COMPLETE'))
