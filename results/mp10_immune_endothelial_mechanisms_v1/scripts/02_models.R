source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R')
suppressPackageStartupMessages({library(edgeR);library(limma)})
b<-qread(file.path(B,'objects/pseudobulk.qs'),nthreads=1);pr<-b$profiles;cfg<-fromJSON(file.path(B,'config/run_config.json'));dir.create(tb('de'),showWarnings=FALSE)
receivers<-c('Treg','Tfh/CXCL13-like','Arterial EC','Capillary EC');datasets<-sort(unique(pr$dataset));plan<-list();counter<-0
add<-function(kind,target,threshold=50,patient='main',exposure='Q50',subset='own',omit='none',only='ALL') {counter<<-counter+1;plan[[counter]]<<-data.table(model_id=sprintf('M%04d',counter),kind=kind,target=target,min_cells=threshold,patient=patient,exposure=exposure,subset=subset,omit_dataset=omit,only_dataset=only)}
for(h in c('H50','H10','H05')) {add('Epi',h);add('Epi',h,30);add('Epi',h,patient='replacement');for(n in c(50,30))add('Epi',h,n,subset='common');for(d in datasets){add('Epi',h,omit=d);add('Epi',h,only=d)}}
for(r in receivers){add('receiver',r);add('receiver',r,30);add('receiver',r,exposure='Q95');add('receiver',r,patient='replacement');for(d in datasets)add('receiver',r,omit=d)}
plan<-rbindlist(plan);plan[,primary:=min_cells==50 & patient=='main' & exposure=='Q50' & subset=='own' & omit_dataset=='none' & only_dataset=='ALL'];put(plan,'model_plan.tsv')
statuses<-list(); fits<-list(); memberships<-list(); checks<-list(); summaries<-list()
for(ii in seq_len(nrow(plan))){p<-plan[ii];id<-p$model_id;cat(id,p$kind,p$target,p$min_cells,p$patient,p$subset,p$omit_dataset,'\n')
 ix<-if(p$patient=='main')pr$main_patient else pr$replacement_patient
 sel<-pr[ix & n_cells>=p$min_cells]
 if(p$omit_dataset!='none')sel<-sel[dataset!=p$omit_dataset]
 if(p$only_dataset!='ALL')sel<-sel[dataset==p$only_dataset]
 if(p$kind=='Epi'){
  if(p$subset=='common'){ok<-sel[profile_type %in% c('H50','H10','H05','L50'),.N,by=sample_uid][N==4,sample_uid];sel<-sel[sample_uid %in% ok]}
  ok<-sel[profile_type %in% c(p$target,'L50'),.N,by=sample_uid][N==2,sample_uid];sel<-sel[sample_uid %in% ok & profile_type %in% c(p$target,'L50')];setorder(sel,sample_uid,profile_type)
  sel[,specimen:=factor(sample_uid)];sel[,state:=factor(ifelse(profile_type=='L50','Low','High'),levels=c('Low','High'))]
 }else{sel<-sel[profile_type==p$target];setorder(sel,sample_uid);sel[,exposure_z:=get(paste0(p$exposure,'_z'))];sel[,dataset_factor:=factor(dataset)]}
 np<-uniqueN(sel$patient_uid_recorded);ns<-uniqueN(sel$sample_uid);nd<-uniqueN(sel$dataset)
 status<-copy(p);status[,`:=`(n_profiles=nrow(sel),n_specimens=ns,n_recorded_patient_units=np,n_datasets=nd,datasets=paste(sort(unique(sel$dataset)),collapse=';'),status='not_evaluable',reason='')]
 if(np<4){status[,reason:='fewer_than_4_recorded_patient_units'];statuses[[id]]<-status;next}
 design<-if(p$kind=='Epi')model.matrix(~specimen+state,sel)else if(nd>1)model.matrix(~dataset_factor+exposure_z,sel)else model.matrix(~exposure_z,sel)
 rank<-qr(design)$rank;df<-nrow(design)-rank;status[,`:=`(design_rank=rank,design_columns=ncol(design),residual_df=df)]
 memberships[[id]]<-cbind(data.table(model_id=id),sel)
 checks[[id]]<-data.table(model_id=id,profile_id=sel$profile_id,design_row=apply(design,1,function(z)paste(z,collapse=';')),design_columns=paste(colnames(design),collapse=';'))
 if(rank<ncol(design)||df<3){status[,reason:=if(rank<ncol(design))'rank_deficient_dataset_exposure_confounding'else'insufficient_residual_df'];statuses[[id]]<-status;next}
 outcome<-tryCatch({
  yy<-as.matrix(b$counts[,sel$profile_id,drop=FALSE]);y<-DGEList(yy);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y,method='TMM');y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);tt<-as.data.table(topTags(glmQLFTest(fit,coef=ncol(design)),n=Inf,sort.by='none')$table,keep.rownames='gene')
  tt[,`:=`(model_id=id,signed_stat=sign(logFC)*sqrt(pmax(F,0)),joint_main_q=NA_real_)];tt[,mean_detection_fraction:=Matrix::rowMeans(sweep(b$detected_cells[gene,sel$profile_id,drop=FALSE],2,sel$n_cells,'/'))]
  put(tt,paste0('de/',id,'.tsv.gz'))
  # Compact sample-level normalized expression supports downstream activities; raw pseudobulk is reused.
  fits[[id]]<-list(model=p,profiles=sel,design=design,logCPM=cpm(y,log=TRUE,prior.count=2),norm_factors=y$samples$norm.factors,genes=rownames(y))
  summaries[[id]]<-tt[,.(n_tested=.N,n_positive=sum(FDR<.05&logFC>0),n_negative=sum(FDR<.05&logFC<0),min_p=min(PValue))][,model_id:=id]
  status[,`:=`(status='completed',reason='')];TRUE
 },error=function(e){status[,`:=`(status='technically_blocked',reason=conditionMessage(e))];FALSE})
 statuses[[id]]<-status
 put(rbindlist(statuses,fill=TRUE),'model_status.tsv')
}
put(rbindlist(statuses,fill=TRUE),'model_status.tsv');put(rbindlist(memberships,fill=TRUE),'model_membership.tsv.gz');put(rbindlist(checks,fill=TRUE),'model_design_matrices.tsv.gz');put(rbindlist(summaries,fill=TRUE),'de_summary.tsv')
for(k in c('Epi','receiver')){ids<-plan[primary==TRUE & kind==k,model_id];tt<-rbindlist(lapply(ids,function(id){f<-tb(paste0('de/',id,'.tsv.gz'));if(file.exists(f))fread(f)else NULL}));if(nrow(tt)){tt[,joint_main_q:=p.adjust(PValue,'BH')];for(id in unique(tt$model_id))put(tt[model_id==id],paste0('de/',id,'.tsv.gz'))}}
qsave(fits,file.path(B,'objects/model_inputs.qs'),preset='high',nthreads=1)
# Expression and detection for every tested profile, focused candidates, plus complete ligand/receiver prior audit downstream.
focus<-c('NECTIN2','PVRL2','PVR','VEGFA','VEGFB','VEGFC','PGF','ANGPT1','ANGPT2','TGFB1','TGFB2','CXCL12','MDK','MIF','LGALS9','CDH1','TIGIT','PVRIG','CD226','FLT1','KDR','FLT4','NRP1','NRP2','TEK')
ex<-rbindlist(lapply(focus,function(g){if(!g %in% rownames(b$counts))return(data.table(gene=g,status='feature_absent'));z<-copy(pr);z[,`:=`(gene=g,total_UMI=as.numeric(b$counts[g,profile_id]),positive_cells=as.numeric(b$detected_cells[g,profile_id]),normalized_mean=as.numeric(b$normalized_sum[g,profile_id])/n_cells,status='computed')];z[,positive_fraction:=positive_cells/n_cells];z} ),fill=TRUE);put(ex,'candidate_expression_detection.tsv.gz')
# Exposure variation and dataset explanatory fraction are descriptive, not High/Low tests.
q<-rd('specimen_exposure_eligibility.tsv');v<-rbindlist(lapply(c('Q50','Q95'),function(e)q[,.(exposure=e,n=.N,mean=mean(get(e)),sd=sd(get(e)),min=min(get(e)),max=max(get(e))),by=dataset]));put(v,'within_dataset_exposure_variation.tsv')
cat('MODELS COMPLETE\n')
