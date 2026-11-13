source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R')
suppressPackageStartupMessages(library(Seurat))
cfg<-list(task='mp10_immune_endothelial_mechanisms_v1',baseline='c2lfinal_20261108T123140Z_0b5fd64f',seed=123,seed_scenic_sensitivity=456,threads_total=4,threads_per_process=1,min_cells=c(50,30),min_independent_units=4,min_residual_df=3,quantile_type=7,main_patient_specimen='T4463209a',replacement='T4463209b',exposures=c('Q50','Q95'),exposure_scale='one SD across original main-a specimens, fixed for all models',epi_model='~ specimen + state',receiver_model='~ dataset + exposure',normalization='edgeR TMM',DE='edgeR robust quasi-likelihood',gene_filter='edgeR filterByExpr with actual design',gsea=list(primary='Hallmark',secondary='Reactome',minSize=15,maxSize=500),progeny=list(top=500,min_targets=10),collectri=list(method='decoupleR ULM',min_targets=5,split_complexes=FALSE),dorothea_confidence=c('A','B'),nichenet=list(target_fdr=.05,target_min_informative=10,detection=c(.10,.05),detection_rule='median of eligible specimen nonzero-cell fractions',target_sets=c('positive','negative','union')),scenic=list(cap_per_specimen=100,min_detection_fraction=.01,min_cells_gene=3,aucell_auc_threshold=.05),patient_identity_caveat='recorded IDs are derived/provisional except confirmed repeated a/b relationship; not verified independent patients',frozen_before_effects=as.character(Sys.time()))
write_json(cfg,file.path(B,'config/run_config.json'),pretty=TRUE,auto_unbox=TRUE)
a<-fread(file.path(BASE,'tables/analysis_cells.tsv.gz')); sm<-fread(file.path(BASE,'tables/sample_metadata.tsv')); q<-fread(file.path(HDEF,'tables/specimen_quantiles_coverage.tsv'))
stopifnot(nrow(q)==40,all(q$sample_uid %in% sm[epi_state_eligible==TRUE,sample_uid]))
a<-a[sample_uid %in% q$sample_uid];stopifnot(!anyDuplicated(a$barcode),all(!is.na(a$analysis_celltype_c2l_final)))
q[,main_patient:=!grepl('T4463209b',specimen_id)];q[,replacement_patient:=!grepl('T4463209a',specimen_id)]
for(v in c('Q50','Q95')){mu<-mean(q[main_patient==TRUE,get(v)]);sdv<-sd(q[main_patient==TRUE,get(v)]);q[,(paste0(v,'_z')):=(get(v)-mu)/sdv];q[,(paste0(v,'_center')):=mu];q[,(paste0(v,'_scale')):=sdv]}
put(q,'specimen_exposure_eligibility.tsv');put(sm[sample_uid %in% q$sample_uid],'sample_metadata.tsv')
cat('Reading original object\n');o<-qread(objpath,nthreads=1)
md<-as.data.table(o@meta.data,keep.rownames='barcode');put(md,'source_cell_metadata.tsv.gz')
layers<-Layers(o[['RNA']]);stopifnot(all(c('counts','data') %in% layers))
ct<-LayerData(o,assay='RNA',layer='counts');ex<-LayerData(o,assay='RNA',layer='data')
stopifnot(!anyDuplicated(rownames(ct)),!anyDuplicated(colnames(ct)),all(ct@x>=0),all(abs(ct@x-round(ct@x))<1e-8),all(a$barcode %in% colnames(ct)),all(a$barcode %in% colnames(ex)),identical(rownames(ct),rownames(ex)))
put(data.table(layer=c('counts','data'),n_genes=c(nrow(ct),nrow(ex)),n_cells=c(ncol(ct),ncol(ex)),min=c(min(ct),min(ex)),max=c(max(ct),max(ex)),integer_values=c(all(ct@x==round(ct@x)),all(ex@x==round(ex@x)))), 'source_layer_audit.tsv')
stopifnot(all(md[match(a$barcode,barcode),dataset]==a$dataset))
qc<-md[match(a$barcode,barcode),.(barcode,nCount_RNA,nFeature_RNA,percent.mt,subtype_v2,state_v2,source_lineage)]
a<-merge(a,qc,by='barcode',sort=FALSE);a<-merge(a,q[,.(sample_uid,Q50,Q90,Q95)],by='sample_uid',sort=FALSE)
e<-a[analysis_major_c2l_final=='Epithelial'];stopifnot(nrow(e)==116448,all(is.finite(e$MP10_cell)))
m<-rbindlist(list(e[,.(barcode,sample_uid,profile_type='Epi_all')],e[MP10_cell<=Q50,.(barcode,sample_uid,profile_type='L50')],e[MP10_cell>Q50,.(barcode,sample_uid,profile_type='H50')],e[MP10_cell>Q90,.(barcode,sample_uid,profile_type='H10')],e[MP10_cell>Q95,.(barcode,sample_uid,profile_type='H05')],a[analysis_major_c2l_final!='Epithelial',.(barcode,sample_uid,profile_type=analysis_celltype_c2l_final)]))
m[,profile_id:=paste(sample_uid,profile_type,sep='|')];lev<-sort(unique(m$profile_id));m[,profile_order:=match(profile_id,lev)];m[,cell_order:=match(barcode,colnames(ct))]
put(m,'pseudobulk_membership.tsv.gz')
map<-sparseMatrix(i=m$cell_order,j=m$profile_order,x=1,dims=c(ncol(ct),length(lev)),dimnames=list(colnames(ct),lev))
cat('Aggregating counts\n');pb<-ct%*%map
cat('Aggregating detections\n');bn<-ct;bn@x[]<-1;det<-bn%*%map;rm(bn);gc()
cat('Aggregating original normalized expression\n');nmap<-map[match(colnames(ex),rownames(map)),,drop=FALSE];norm<-ex%*%nmap
pr<-m[,.(n_cells=.N),by=.(profile_id,sample_uid,profile_type)][match(lev,profile_id)]
pr[,total_counts:=Matrix::colSums(pb)];pr[,detected_genes:=Matrix::colSums(pb>0)];pr[,counts_per_cell:=total_counts/n_cells];pr<-merge(pr,q[,.(sample_uid,dataset,specimen_id,patient_uid_recorded,main_patient,replacement_patient,Q50,Q95,Q50_z,Q95_z)],by='sample_uid',sort=FALSE);pr<-pr[match(lev,profile_id)]
pr[,`:=`(eligible_50=n_cells>=50,eligible_30=n_cells>=30)];put(pr,'pseudobulk_profiles.tsv')
qsave(list(counts=pb,detected_cells=det,normalized_sum=norm,profiles=pr,feature_universe=rownames(ct)),file.path(B,'objects/pseudobulk.qs'),preset='high',nthreads=1)
put(data.table(gene=rownames(ct),original_order=seq_len(nrow(ct)),total_counts=Matrix::rowSums(pb)),'feature_universe.tsv.gz')
put(a[,.(barcode,sample_uid,dataset,patient_uid_recorded,analysis_celltype_c2l_final,analysis_major_c2l_final,MP10_cell,nCount_RNA,nFeature_RNA,percent.mt,subtype_v2,state_v2,source_lineage)],'analysis_cell_metadata.tsv.gz')
qcsummary<-merge(m,a[,.(barcode,nCount_RNA,nFeature_RNA,percent.mt,subtype_v2)],by='barcode',allow.cartesian=TRUE)[,.(median_nCount=as.numeric(median(nCount_RNA)),median_nFeature=as.numeric(median(nFeature_RNA)),median_mt=as.numeric(median(percent.mt)),n_subtype_annotated=sum(!is.na(subtype_v2)),subtypes=paste(sort(unique(na.omit(subtype_v2))),collapse=';')),by=profile_id];put(qcsummary,'profile_quality.tsv')
put(a[,.(n_cells=.N),by=.(sample_uid,dataset,analysis_major_c2l_final,analysis_celltype_c2l_final)],'full_celltype_counts.tsv')
files<-c(objpath,file.path(BASE,'tables/analysis_cells.tsv.gz'),file.path(BASE,'tables/sample_metadata.tsv'),file.path(HDEF,'tables/specimen_quantiles_coverage.tsv'),file.path(P,'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'))
manifest<-rbindlist(lapply(files,function(f){fi<-file.info(f);data.table(path=normalizePath(f),bytes=fi$size,mtime=as.character(fi$mtime),sha256=digest(f,file=TRUE,algo='sha256'))}));put(manifest,'input_manifest.tsv')
put(data.table(item=c('counts_barcodes','genes','frozen_scores','pseudobulk_membership','config'),sha256=c(digest(colnames(ct),algo='sha256'),digest(rownames(ct),algo='sha256'),digest(e[,.(barcode,MP10_cell)],algo='sha256'),digest(m,algo='sha256'),digest(cfg,algo='sha256'))),'input_content_hashes.tsv')
writeLines(capture.output(sessionInfo()),file.path(B,'logs/sessionInfo_pseudobulk.txt'))
cat('PSEUDOBULK COMPLETE:',nrow(pr),'profiles',nrow(pb),'genes\n')
