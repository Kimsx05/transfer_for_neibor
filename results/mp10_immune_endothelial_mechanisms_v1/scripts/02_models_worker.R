source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R');source(file.path(S,'model_helpers.R'));suppressPackageStartupMessages({library(edgeR);library(limma)})
args<-commandArgs(TRUE);worker<-as.integer(args[1]);workers<-as.integer(args[2]);plan<-rd('model_plan.tsv');b<-qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1);pr<-b$profiles;dir.create(file.path(B,'objects/model_parts'),showWarnings=FALSE)
inputkey<-digest(list(digest(file.path(B,'objects/pseudobulk.qs'),file=TRUE,algo='sha256'),digest(file.path(B,'config/run_config.json'),file=TRUE,algo='sha256'),edgeR=as.character(packageVersion('edgeR')),method='TMM_filterByExpr_design_robustQL_v1'),algo='sha256')
legacy<-rd('model_status_before_parallel_resume.tsv')
for(ii in seq_len(nrow(plan))[(seq_len(nrow(plan))-1)%%workers==worker]){
 p<-plan[ii];id<-p$model_id;s<-select_model(p,pr);key<-digest(list(inputkey,p,s),algo='sha256');out<-file.path(B,'objects/model_parts',paste0(id,'.qs'));statusfile<-tb(paste0('de/',id,'_status.tsv'))
 if(file.exists(out)&&file.exists(statusfile)){oo<-qread(out,nthreads=1);if(identical(oo$input_hash,key)){cat(id,'verified cache\n');next}}
 status<-copy(p);status[,`:=`(status='not_evaluable',reason='insufficient_coverage_or_residual_df_or_rank')]
 if(is.null(s)){status[,`:=`(n_profiles=0L,n_specimens=0L,n_recorded_patient_units=0L,n_datasets=0L)];put(status,paste0('de/',id,'_status.tsv'));next}
 sel<-s$profiles;design<-s$design;status[,`:=`(n_profiles=nrow(sel),n_specimens=uniqueN(sel$sample_uid),n_recorded_patient_units=uniqueN(sel$patient_uid_recorded),n_datasets=uniqueN(sel$dataset),datasets=paste(sort(unique(sel$dataset)),collapse=';'),design_rank=qr(design)$rank,design_columns=ncol(design),residual_df=nrow(design)-qr(design)$rank)]
 cat(id,p$kind,p$target,p$min_cells,p$omit_dataset,'\n');flush.console()
 tryCatch({
 y<-DGEList(as.matrix(b$counts[,sel$profile_id,drop=FALSE]));keep<-filterByExpr(y,design=design);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);depath<-tb(paste0('de/',id,'.tsv.gz'))
 if(id%in%legacy[status=='completed',model_id]&&file.exists(depath)){
 tt<-fread(depath);stopifnot(setequal(tt$gene,rownames(y)),all(tt$model_id==id),legacy[model_id==id,n_profiles]==nrow(sel));tt[,signed_stat:=sign(logFC)*sqrt(pmax(F,0))];cat('Reusing validated completed serial inference',id,'\n')
 }else{
 y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);tt<-as.data.table(topTags(glmQLFTest(fit,coef=ncol(design)),n=Inf,sort.by='none')$table,keep.rownames='gene');tt[,`:=`(model_id=id,signed_stat=sign(logFC)*sqrt(pmax(F,0)),joint_main_q=NA_real_)];tt[,mean_detection_fraction:=Matrix::rowMeans(sweep(b$detected_cells[gene,sel$profile_id,drop=FALSE],2,sel$n_cells,'/'))]
 }
 put(tt,paste0('de/',id,'.tsv.gz'));qsave(list(model=p,profiles=sel,design=design,logCPM=cpm(y,log=TRUE,prior.count=2),norm_factors=y$samples$norm.factors,genes=rownames(y),input_hash=key),paste0(out,'.tmp'),preset='high',nthreads=1);file.rename(paste0(out,'.tmp'),out);status[,`:=`(status='completed',reason='',input_hash=key)]
 },error=function(e){status[,`:=`(status='technically_blocked',reason=conditionMessage(e))]})
 put(status,paste0('de/',id,'_status.tsv'))
}
cat('WORKER COMPLETE',worker,'\n')
