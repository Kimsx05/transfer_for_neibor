source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R')
b<-qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1);pr<-b$profiles
# Expression and detection for every tested profile, focused candidates, plus complete ligand/receiver prior audit downstream.
focus<-c('NECTIN2','PVRL2','PVR','VEGFA','VEGFB','VEGFC','PGF','ANGPT1','ANGPT2','TGFB1','TGFB2','CXCL12','MDK','MIF','LGALS9','CDH1','TIGIT','PVRIG','CD226','FLT1','KDR','FLT4','NRP1','NRP2','TEK')
ex<-rbindlist(lapply(focus,function(g){if(!g %in% rownames(b$counts))return(data.table(gene=g,status='feature_absent'));z<-copy(pr);z[,`:=`(gene=g,total_UMI=as.numeric(b$counts[g,profile_id]),positive_cells=as.numeric(b$detected_cells[g,profile_id]),normalized_mean=as.numeric(b$normalized_sum[g,profile_id])/n_cells,status='computed')];z[,positive_fraction:=positive_cells/n_cells];z} ),fill=TRUE);put(ex,'candidate_expression_detection.tsv.gz')
# Exposure variation and dataset explanatory fraction are descriptive, not High/Low tests.
q<-rd('specimen_exposure_eligibility.tsv');v<-rbindlist(lapply(c('Q50','Q95'),function(e)q[,.(exposure=e,n=.N,mean=mean(get(e)),sd=sd(get(e)),min=min(get(e)),max=max(get(e))),by=dataset]));put(v,'within_dataset_exposure_variation.tsv')
cat('CANDIDATE TABLES COMPLETE\n')
