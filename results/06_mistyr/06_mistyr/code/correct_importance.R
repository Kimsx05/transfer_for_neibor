# Correct an observed mistyR 1.18.0 coefficient serialization defect without modifying its package or models.
args<-commandArgs(TRUE);root<-normalizePath(args[1]);Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
suppressPackageStartupMessages(library(data.table));setDTthreads(1L)
man<-fread(file.path(root,'section_input_summary.tsv'));audits<-list()
for(i in seq_len(nrow(man))){
 folder<-file.path(root,'by_dataset',man$dataset[i],man$section_id[i])
 for(tag in c('wang_style_cellcell_intra','mp10_target_extension','mp10_technical_baseline')){
  output<-switch(tag,wang_style_cellcell_intra='cellcell_importance.tsv',mp10_target_extension='mp10_importance.tsv',mp10_technical_baseline='mp10_baseline_importance.tsv')
  f<-file.path(folder,output);if(!file.exists(f)||file.info(f)$size<10)next
  imp<-fread(f,na.strings='NA');native_dir<-file.path(folder,tag)
  if(!'native_standardized_importance' %in% names(imp))imp[,native_standardized_importance:=standardized_importance]
  rows<-strsplit(readLines(file.path(native_dir,'coefficients.txt'))[-1],'[[:space:]]+')
  for(t in unique(imp$Target)){
   row<-rows[[which(vapply(rows,function(x)x[1]==t,logical(1)))]]
   stopifnot(length(row)==6L,is.na(suppressWarnings(as.numeric(row[4]))))
   native_header_p<-as.numeric(row[5]);serialized_slope_p<-as.numeric(row[6])
   oob<-fread(file.path(native_dir,paste0('oob_',t,'.tsv')))
   fit<-lm(observed~oob_prediction,data=oob);co<-coef(summary(fit));slope_p<-co[2,4]
   stopifnot(isTRUE(all.equal(slope_p,serialized_slope_p,tolerance=1e-8)),isTRUE(all.equal(co[1,4],native_header_p,tolerance=1e-8)))
   raw<-fread(file.path(native_dir,paste0('importances_',t,'_intra.txt')))
   # Preserve exact source p after verifying it against independently reconstructed OOB regression.
   z<-if(is.finite(sd(raw$imp))&&sd(raw$imp)>0)(raw$imp-mean(raw$imp))/sd(raw$imp) else rep(NA_real_,nrow(raw))
   expected_native<-z*(1-native_header_p);names(expected_native)<-raw$target
   z<-z*(1-serialized_slope_p);names(z)<-raw$target
   old<-imp[Target==t & Predictor %in% raw$target]
   stopifnot(isTRUE(all.equal(old$native_standardized_importance,unname(expected_native[old$Predictor]),tolerance=1e-8)))
   imp[Target==t,`:=`(standardized_importance=unname(z[Predictor]),native_p_intra_field=native_header_p,calibration_slope_p=serialized_slope_p,standardization_status='PROJECT_CORRECTED_NATIVE_COEFFICIENT_SERIALIZATION')]
   audits[[length(audits)+1L]]<-data.table(section_id=man$section_id[i],dataset=man$dataset[i],model=tag,target=t,header_field_count=5L,data_field_count=6L,native_p_intra_field=native_header_p,serialized_true_slope_p=serialized_slope_p,recomputed_oob_slope_p=slope_p,coefficient_p_verified=TRUE,native_importance_reproduced=TRUE,status='PASS_WITH_WARNING')
  }
  fwrite(imp,f,sep='\t',na='NA')
 }
}
fwrite(rbindlist(audits),file.path(root,'native_coefficient_serialization_audit.tsv'),sep='\t',na='NA')
cat('Verified and corrected',length(audits),'model-target records; original API outputs preserved.\n')
