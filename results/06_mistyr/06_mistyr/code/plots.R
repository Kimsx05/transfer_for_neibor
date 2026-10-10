args<-commandArgs(TRUE);root<-normalizePath(args[1]);Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
options(warn=1)
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(1L)
# Reuse project save_both function by extracting its AST only; do not execute its parent QC pipeline.
old_code<-file.path(dirname(dirname(root)),'ST/code/03_spatial_QC_summary.R')
expr<-parse(old_code);selected<-Filter(function(z)is.call(z)&&identical(z[[1]],as.name('<-'))&&identical(z[[2]],as.name('save_both')),as.list(expr));stopifnot(length(selected)==1L);eval(selected[[1]])
replace_dpi<-function(z){if(is.call(z)||is.pairlist(z)){for(j in seq_along(z))z[[j]]<-replace_dpi(z[[j]])};if(is.numeric(z)&&length(z)==1L&&z==180)z<-600;z}
body(save_both)<-replace_dpi(body(save_both));capture.output(save_both,file=file.path(root,'adapted_R_save_both.txt'))
index<-list()
savefig<-function(g,folder,stem,w,h){
 dir.create(file.path(folder,'figures'),showWarnings=FALSE);dir.create(file.path(folder,'thumbnails'),showWarnings=FALSE)
 save_both(g,file.path(folder,'figures'),stem,w,h)
 ggsave(file.path(folder,'thumbnails',paste0(stem,'.png')),g,width=w,height=h,dpi=110,bg='white',limitsize=FALSE)
 index[[length(index)+1L]]<<-data.table(stem=stem,pdf=file.path(folder,'figures',paste0(stem,'.pdf')),png=file.path(folder,'figures',paste0(stem,'.png')),thumbnail=file.path(folder,'thumbnails',paste0(stem,'.png')),png_dpi=600)
}
a<-fread(file.path(root,'cellcell_dataset_median_standardized_importance.tsv'),na.strings='NA')
b<-fread(file.path(root,'mp10_importance.tsv'),na.strings='NA')
perf<-fread(file.path(root,'mp10_baseline_vs_tme.tsv'),na.strings='NA')
man<-fread(file.path(root,'section_input_summary.tsv'),na.strings='NA')
dict<-fread(file.path(root,'variable_dictionary.tsv'),na.strings='NA')
dict[variable=='log1p_nCount',source_label:='log1p(nCount) [technical]']
dict[variable=='epi_fraction_gate',source_label:='Original Epi fraction [technical]']
for(ds in unique(a$dataset)){
 folder<-file.path(root,'by_dataset',ds);z<-a[dataset==ds]
 majors<-unique(c(z$Target_label,z$Predictor_label));grid<-CJ(Target_label=majors,Predictor_label=majors)
 z<-merge(grid,z,by=c('Target_label','Predictor_label'),all.x=TRUE)
 g<-ggplot(z,aes(Target_label,Predictor_label,fill=median_standardized_importance_reportable))+geom_tile()+
  scale_fill_gradient2(low='#2166ac',mid='#fff7bc',high='#b2182b',midpoint=0,na.value='#bdbdbd',name='Median standardized\nimportance')+
  geom_tile(data=z[!is.na(median_standardized_importance_reportable)&median_standardized_importance_reportable==0],fill='black')+
  geom_text(data=z[!is.na(median_standardized_importance_reportable)],aes(label=paste0('n=',n_reportable)),size=2.2,color='#222222')+
  theme_bw(base_size=10)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
  labs(x='Target major cell type',y='Predictor major cell type',caption='Internally reportable models only; no spatial validation for A.\nNegative = below-average importance, not inhibition. Zero black; NA grey.')
 savefig(g,folder,'cellcell_median_standardized_importance',8.5,7)
 z<-b[dataset==ds]
 grid<-CJ(section_id=man[dataset==ds]$section_id,Predictor=dict[role %in% c('B_TME_predictor','B_technical_predictor')]$variable)
 z<-merge(grid,z,by=c('section_id','Predictor'),all.x=TRUE)
 z[,Predictor_label:=dict$source_label[match(Predictor,dict$variable)]]
 z[,shown:=fifelse(importance_interpretable,standardized_importance,NA_real_)]
 g<-ggplot(z,aes(section_id,Predictor_label,fill=shown))+geom_tile()+
  scale_fill_gradient2(low='#2166ac',mid='#fff7bc',high='#b2182b',midpoint=0,na.value='#bdbdbd',name='Standardized\nimportance')+
  geom_tile(data=z[!is.na(shown)&shown==0],fill='black')+theme_bw(base_size=9)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=65,hjust=1))+
  labs(x=NULL,y=NULL,caption=paste0('MP10 target extension: ',uniqueN(z[importance_interpretable==TRUE]$section_id),'/',uniqueN(z$section_id),' reportable sections (spatial full R2 > 0 and delta R2 > 0).\nNegative is relative importance, not biological direction. Zero black; NA grey.'))
 savefig(g,folder,'mp10_predictive_importance_by_section',max(8,uniqueN(z$section_id)*.35+5),10)
 z<-perf[dataset==ds]
 g<-ggplot(z,aes(spatial_r2,section_id,color=model,group=section_id))+geom_vline(xintercept=0,linetype=2,color='grey60')+geom_line(color='grey65',na.rm=TRUE)+geom_point(size=2,na.rm=TRUE)+
  geom_text(data=z[model=='tme' & spatial_cv_complete==FALSE],aes(x=0,y=section_id,label='CV unavailable'),inherit.aes=FALSE,hjust=0,nudge_x=.03,color='grey40',size=3)+
  scale_color_manual(values=c(baseline='#31688e',tme='#d1495b'),labels=c(baseline='Epi fraction + log1p counts',tme='Technical + all 29 TME factors'))+
  theme_bw(base_size=10)+theme(legend.position='bottom')+labs(x='Spatial block out-of-fold R2 (unclipped)',y=NULL,color=NULL,caption='Same complete-case spots and fixed contiguous bands for both models.\nWithin-section held-out blocks; no buffer and no between-patient validation.')
 savefig(g,folder,'mp10_baseline_vs_tme_spatial_R2',9,max(4,uniqueN(z$section_id)*.3+2))
}
pilot<-fread(file.path(root,'by_dataset/GSE319536/GSE319536_BC15/model_input.tsv'))
g<-ggplot(pilot,aes(array_col/2,array_row*sqrt(3)/2,color=factor(spatial_fold),shape=analysis_eligible_legacy))+geom_point(size=.65)+coord_equal()+scale_y_reverse()+theme_bw(base_size=10)+labs(x='True lattice x',y='True lattice y',color='Held-out fold',shape='Legacy eligible',caption='BC15 pilot: geometry-only contiguous bands. All spots shown; B uses eligible complete cases.')
savefig(g,root,'pilot_spatial_fold_audit',8,7)
fwrite(rbindlist(index),file.path(root,'figure_index.tsv'),sep='\t')
