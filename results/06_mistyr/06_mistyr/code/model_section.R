#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
args<-commandArgs(TRUE);root<-normalizePath(args[1]);folder<-normalizePath(args[2])
.libPaths(c(file.path(root,'runtime/Rlib'),.libPaths()))
suppressPackageStartupMessages({library(mistyR);library(data.table);library(jsonlite);library(future)})
plan(sequential);setDTthreads(1L)
p<-fromJSON(file.path(root,'parameters.json'));d<-fread(file.path(folder,'model_input.tsv'),na.strings='NA')
sec<-d$section_id[1]; ds<-d$dataset[1]
write_t<-function(x,name)fwrite(x,file.path(folder,name),sep='\t',na='NA')
active_dir<-'';audit_rows<-list();internal_rows<-list()
metrics<-function(y,pred){c(r2=if(var(y)>0)1-sum((y-pred)^2)/sum((y-mean(y))^2) else NA_real_,rmse=sqrt(mean((y-pred)^2)),mae=mean(abs(y-pred)))}
usable<-function(x)all(is.finite(x)) && is.finite(var(x)) && var(x)>0
raw_meta<-function(y,oob,seed,k){
  folds<-withr::with_seed(seed,caret::createFolds(y,k=k))
  rbindlist(lapply(seq_along(folds),function(i){te<-folds[[i]];tr<-setdiff(seq_along(y),te);fit<-lm(y~oob,data=data.frame(y=y[tr],oob=oob[tr]));pred<-predict(fit,newdata=data.frame(oob=oob[te]));m<-metrics(y[te],pred);data.table(fold=i,n_train=length(tr),n_test=length(te),r2=m['r2'],rmse=m['rmse'],mae=m['mae'])}))
}
# Same ranger learner as mistyR 1.18.0; additionally retains models/OOB predictions and enforces target exclusion.
audit_rf<-function(view_data,target,seed,...){
  predictors<-setdiff(names(view_data),target)
  stopifnot(length(predictors)>0,!target %in% predictors)
  for(v in predictors){
    stopifnot(!identical(view_data[[v]],view_data[[target]]))
    # A perfect monotone copy is also rejected. Ordinary correlations are reported separately.
    rho<-suppressWarnings(cor(view_data[[v]],view_data[[target]],method='spearman'))
    stopifnot(!is.finite(rho)||abs(rho)<1-1e-12)
  }
  fit<-ranger::ranger(x=as.data.frame(view_data[,predictors,drop=FALSE]),y=view_data[[target]],num.trees=p$num_trees,importance=p$importance,num.threads=1,seed=seed,min.node.size=p$min_node_size)
  saveRDS(fit,file.path(active_dir,paste0('forest_',target,'.rds')),compress=FALSE)
  fwrite(data.table(index=seq_len(nrow(view_data)),observed=view_data[[target]],oob_prediction=fit$predictions),file.path(active_dir,paste0('oob_',target,'.tsv')),sep='\t')
  m<-metrics(view_data[[target]],fit$predictions)
  f<-raw_meta(view_data[[target]],fit$predictions,seed,p$internal_meta_folds);f[,target:=target]
  fwrite(f,file.path(active_dir,paste0('unclipped_internal_meta_folds_',target,'.tsv')),sep='\t')
  internal_rows[[length(internal_rows)+1L]]<<-data.table(model=basename(active_dir),target=target,n=nrow(view_data),target_variance=var(view_data[[target]]),n_predictors=length(predictors),oob_r2=m['r2'],oob_rmse=m['rmse'],oob_mae=m['mae'],internal_meta_r2=mean(f$r2),internal_meta_rmse=mean(f$rmse),evaluation='OOB + random meta folds; not spatial CV')
  audit_rows[[length(audit_rows)+1L]]<<-data.table(model=basename(active_dir),target=target,predictor=predictors,target_excluded=TRUE,deterministic_copy_excluded=TRUE,status='PASS')
  fwrite(internal_rows[[length(internal_rows)]],file.path(active_dir,paste0('unclipped_summary_',target,'.tsv')),sep='\t')
  fwrite(audit_rows[[length(audit_rows)]],file.path(active_dir,paste0('exclusion_audit_',target,'.tsv')),sep='\t')
  list(unbiased.predictions=tibble::tibble(index=seq_len(nrow(view_data)),prediction=fit$predictions),importances=fit$variable.importance)
}
run_native<-function(tab,tag,targets=NULL){
  active_dir<<-file.path(folder,tag);dir.create(active_dir,showWarnings=FALSE)
  views<-create_initial_view(tibble::as_tibble(tab))
  stopifnot(identical(setdiff(names(views),'misty.uniqueid'),'intraview'))
  run_misty(views,results.folder=active_dir,seed=p$seed,target.subset=targets,cv.folds=p$internal_meta_folds,cached=FALSE,append=FALSE,model.function=audit_rf)
  native<-collect_results(active_dir)
  # future_map serializes closures even with sequential plan: collect audit side effects from files.
  internal_rows<<-c(internal_rows,lapply(list.files(active_dir,pattern='^unclipped_summary_.*tsv$',full.names=TRUE),fread))
  audit_rows<<-c(audit_rows,lapply(list.files(active_dir,pattern='^exclusion_audit_.*tsv$',full.names=TRUE),fread))
  saveRDS(native,file.path(active_dir,'native_collected_results.rds'))
  fwrite(native$improvements,file.path(active_dir,'native_internal_metrics_clipped.tsv'),sep='\t',na='NA')
  imp<-as.data.table(native$importances)
  imp[,`:=`(standardized_importance=Importance,model=tag)]
  # Never allow package NaN-to-zero handling to define non-estimable importance.
  for(t in unique(imp$Target)){
    raw<-fread(file.path(active_dir,paste0('importances_',t,'_intra.txt')))
    rawmap<-setNames(raw$imp,raw$target)
    imp[Target==t,raw_importance:=rawmap[Predictor]]
    if(!is.finite(sd(raw$imp))||sd(raw$imp)==0)imp[Target==t,standardized_importance:=NA_real_]
    imp[Target==t & Predictor==t,`:=`(raw_importance=NA_real_,standardized_importance=NA_real_)]
  }
  imp[,c('sample','Importance'):=NULL]
  imp
}
status<-list(); association<-list();
assoc<-function(tab,targets,preds,scope){
  rbindlist(lapply(targets,function(t)rbindlist(lapply(setdiff(preds,t),function(v){
    z<-is.finite(tab[[t]])&is.finite(tab[[v]]); a<-tab[[t]][z];b<-tab[[v]][z]
    data.table(scope=scope,target=t,predictor=v,n=sum(z),spearman_rho=if(length(a)>2&&var(a)>0&&var(b)>0)cor(a,b,method='spearman') else NA_real_,meaning='marginal association; separate from importance; no causal sign')
  }))))
}
tryCatch({
  av<-p$major_variables;a<-d[complete.cases(d[,..av])]
  valid<-av[vapply(a[,..av],usable,logical(1))];bad<-setdiff(av,valid)
  write_t(data.table(variable=av,variance=vapply(a[,..av],var,numeric(1)),valid=av %in% valid),'A_variable_status.tsv')
  if(nrow(a)<p$min_spots||length(valid)<2)stop('NOT_ESTIMABLE_A: insufficient spots or nonconstant majors')
  imp<-run_native(as.data.frame(a[,..valid]),'wang_style_cellcell_intra')
  perf<-rbindlist(internal_rows)[model=='wang_style_cellcell_intra']
  imp<-merge(imp,perf[,.(Target=target,oob_r2,internal_meta_r2)],by='Target',all.x=TRUE)
  imp[,importance_interpretable:=is.finite(internal_meta_r2)&internal_meta_r2>0&is.finite(oob_r2)&oob_r2>0]
  imp[,status:=fifelse(importance_interpretable,'PASS','PASS_WITH_WARNING')]
  write_t(imp,'cellcell_importance.tsv')
  association[[1]]<-assoc(a,valid,c(valid,'log1p_nCount','epi_fraction_gate'),'A_all_tissue')
  status[[1]]<-data.table(branch='A',status=if(length(bad))'PASS_WITH_WARNING' else 'PASS',reason=paste('nonconstant major targets',length(valid),'excluded',paste(bad,collapse=',')),n=nrow(a))
},error=function(e){status[[1]]<<-data.table(branch='A',status=if(grepl('NOT_ESTIMABLE',conditionMessage(e)))'SKIPPED_NOT_ESTIMABLE' else 'BLOCKED',reason=conditionMessage(e),n=nrow(d))})
tryCatch({
  bv<-c('MP10',p$technical_variables,p$tme_variables)
  b<-d[analysis_eligible_legacy==TRUE & complete.cases(d[,..bv])]
  if(nrow(b)<p$min_spots || !usable(b$MP10))stop('NOT_ESTIMABLE_B: insufficient complete cases or zero target variance')
  pred<-c(p$technical_variables,p$tme_variables); valid<-pred[vapply(b[,..pred],usable,logical(1))]
  tech<-intersect(p$technical_variables,valid);full<-valid
  stopifnot(length(tech)>0,all(!grepl('MP|rank|display|high',full,ignore.case=TRUE)))
  impbase<-run_native(as.data.frame(b[,c('MP10',tech),with=FALSE]),'mp10_technical_baseline',targets='MP10')
  impfull<-run_native(as.data.frame(b[,c('MP10',full),with=FALSE]),'mp10_target_extension',targets='MP10')
  foldrows<-list();predrows<-list();foldaudit<-list();
  for(k in sort(unique(b$spatial_fold))){
    tr<-which(b$spatial_fold!=k);te<-which(b$spatial_fold==k)
    if(length(tr)<p$min_train_spots||length(te)<p$min_test_spots||var(b$MP10[tr])==0){
      foldrows[[length(foldrows)+1L]]<-data.table(fold=k,model='both',n_train=length(tr),n_test=length(te),r2=NA_real_,rmse=NA_real_,mae=NA_real_,status='SKIPPED_NOT_ESTIMABLE');next
    }
    for(tag in c('baseline','tme')){
      vv<-if(tag=='baseline')p$technical_variables else pred
      # Predictor variability is assessed in training fold only.
      vv<-vv[vapply(b[tr,..vv],usable,logical(1))]
      if(!length(vv))stop('NOT_ESTIMABLE_B: no nonconstant predictors in training fold')
      fit<-ranger::ranger(x=as.data.frame(b[tr,..vv]),y=b$MP10[tr],num.trees=p$num_trees,importance=p$importance,num.threads=1,seed=p$seed+k,min.node.size=p$min_node_size)
      cal<-lm(y~oob,data=data.frame(y=b$MP10[tr],oob=fit$predictions))
      rf_pred<-predict(fit,data=as.data.frame(b[te,..vv]),num.threads=1)$predictions
      predout<-as.numeric(predict(cal,newdata=data.frame(oob=rf_pred)))
      stopifnot(all(is.finite(predout)),length(intersect(tr,te))==0)
      m<-metrics(b$MP10[te],predout)
      foldrows[[length(foldrows)+1L]]<-data.table(fold=k,model=tag,n_train=length(tr),n_test=length(te),r2=m['r2'],rmse=m['rmse'],mae=m['mae'],status='PASS')
      pr<-b[te,.(section_id,barcode,spatial_fold,MP10)];pr[,`:=`(model=tag,prediction=predout,rf_prediction=rf_pred,train_mean=mean(b$MP10[tr]))];predrows[[length(predrows)+1L]]<-pr
      foldaudit[[length(foldaudit)+1L]]<-data.table(fold=k,model=tag,predictors=paste(vv,collapse=';'),n_train=length(tr),n_test=length(te),target_excluded=TRUE,fold_overlap=0,calibration='lm trained on training-fold OOB only',status='PASS')
      saveRDS(list(forest=fit,calibrator=cal,predictors=vv,train_barcodes=b$barcode[tr],test_barcodes=b$barcode[te]),file.path(folder,paste0('spatial_fold_',k,'_',tag,'.rds')),compress=FALSE)
    }
  }
  fr<-rbindlist(foldrows,fill=TRUE);pr<-rbindlist(predrows,fill=TRUE)
  if(!nrow(pr)){
    write_t(fr,'mp10_spatial_fold_metrics.tsv')
    stop('NOT_ESTIMABLE_B: no valid spatial folds; internal models retained')
  }
  write_t(fr,'mp10_spatial_fold_metrics.tsv');write_t(pr,'mp10_spatial_predictions.tsv');write_t(rbindlist(foldaudit),'spatial_fold_audit.tsv')
  write_t(b[,.(section_id,barcode,spatial_fold,analysis_eligible_legacy)],'mp10_complete_case_spots.tsv')
  complete<-nrow(pr)==2*nrow(b) && uniqueN(pr$spatial_fold)==p$spatial_folds
  perf<-rbindlist(lapply(c('baseline','tme'),function(tag){
    z<-pr[model==tag];m<-if(nrow(z)>0)metrics(z$MP10,z$prediction) else rep(NA_real_,3)
    data.table(model=tag,n_complete_case=nrow(b),n_predictions=nrow(z),target_variance=var(b$MP10),spatial_r2=if(complete)m[1] else NA_real_,spatial_rmse=if(complete)m[2] else NA_real_,spatial_mae=if(complete)m[3] else NA_real_,spatial_q2_train_mean=if(complete)1-sum((z$MP10-z$prediction)^2)/sum((z$MP10-z$train_mean)^2) else NA_real_,spatial_cv_complete=complete)
  }))
  delta<-perf[model=='tme',spatial_r2]-perf[model=='baseline',spatial_r2]
  interpretable<-isTRUE(complete && perf[model=='tme',spatial_r2]>0 && delta>0)
  perf[,`:=`(delta_r2_tme_minus_baseline=delta,delta_rmse_baseline_minus_tme=perf[model=='baseline',spatial_rmse]-perf[model=='tme',spatial_rmse],importance_interpretable=interpretable,status=if(interpretable)'PASS' else 'PASS_WITH_WARNING')]
  write_t(perf,'mp10_baseline_vs_tme.tsv')
  impfull[,`:=`(importance_interpretable=interpretable,spatial_full_r2=perf[model=='tme',spatial_r2],delta_spatial_r2=delta,status=if(interpretable)'PASS' else 'PASS_WITH_WARNING')]
  write_t(impfull,'mp10_importance.tsv');write_t(impbase,'mp10_baseline_importance.tsv')
  association[[2]]<-assoc(b,'MP10',pred,'B_legacy_eligible_complete_cases')
  status[[2]]<-data.table(branch='B',status=if(complete && interpretable)'PASS' else 'PASS_WITH_WARNING',reason=if(!complete)'INCOMPLETE_SPATIAL_FOLDS; importance not interpreted' else if(!interpretable)'NONPOSITIVE_PERFORMANCE_OR_INCREMENT; retained, not technical failure' else 'positive spatial performance and increment',n=nrow(b))
},error=function(e){status[[2]]<<-data.table(branch='B',status=if(grepl('NOT_ESTIMABLE',conditionMessage(e)))'SKIPPED_NOT_ESTIMABLE' else 'BLOCKED',reason=conditionMessage(e),n=nrow(d[analysis_eligible_legacy==TRUE]))})
st<-rbindlist(status);st[,`:=`(section_id=sec,dataset=ds)];write_t(st,'section_status.tsv')
write_t(rbindlist(internal_rows),'internal_performance_unclipped.tsv');write_t(rbindlist(audit_rows),'predictor_exclusion_audit.tsv')
write_t(rbindlist(association),'abundance_mp_technical_associations.tsv')
capture.output(sessionInfo(),file=file.path(folder,'R_sessionInfo.txt'))
print(st)
if(any(st$status=='BLOCKED'))quit(status=2L)
