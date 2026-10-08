source('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1/code/00_common.R')
Sys.setlocale('LC_COLLATE','C')
logmsg('Audit starting')
versions <- rbindlist(lapply(c('UCell','Seurat','SeuratObject','leidenbase','igraph','uwot','qs','Matrix','BiocParallel','ggplot2','ggrastr','SCP'),function(p) data.table(package=p,version=if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else NA_character_)))
versions<-rbind(data.table(package='R',version=R.version.string),versions)
write_tsv(versions,out('00_input_audit/software_versions.tsv'))
capture.output(sessionInfo(),file=out('logs/sessionInfo_audit.txt'))
capture.output(formals(UCell::ScoreSignatures_UCell),Seurat:::FindNeighbors.default,Seurat:::FindClusters.default,Seurat:::RunLeiden,Seurat:::RunUMAP.default,file=out('00_input_audit/installed_API.txt'))
files <- c(input,file.path(frozen,c('01_input_summary.csv','00_parameters.csv','00_run_signature.txt','01_input_state.qs','01_analysis_sample_id_map.csv','01_sample_metadata_map.csv','07_EpiMP_gene_basis.qs','07_EpiMP_basis_sample_balanced.csv','07_EpiMP_top50_genes.csv','08_cell_EpiMP_scores.qs')),'/home/data/t070721/codex_workspace/Bladder Metabolism/cellchat/specimen_level_CellChat_input_audit.csv')
fi <- file.info(files)
manifest <- data.table(path=files,bytes=fi$size,mtime=as.character(fi$mtime),version_evidence=c('Original input path recorded in frozen 01_input_summary.csv',rep('output_full_cc_seed5_v2 / run_signature=7806d59bc293c240',10),'Existing specimen audit; patient status must be respected'))
write_tsv(manifest,out('00_input_audit/input_manifest.tsv'))
b <- qread(file.path(frozen,'07_EpiMP_gene_basis.qs'))
basis <- b$sample_balanced_mean
stopifnot(identical(dim(basis),c(3000L,14L)),identical(colnames(basis),mp),!anyDuplicated(rownames(basis)),all(is.finite(basis)),b$primary_name=='sample_balanced_mean')
bcsv <- read.csv(file.path(frozen,'07_EpiMP_basis_sample_balanced.csv'),row.names=1,check.names=FALSE)
stopifnot(identical(dimnames(basis),dimnames(as.matrix(bcsv))),max(abs(basis-as.matrix(bcsv)))<1e-14)
writeLines(c(paste('basis dimensions',paste(dim(basis),collapse=' x ')),paste('columns',paste(colnames(basis),collapse=',')),paste('CSV maximum absolute difference',max(abs(basis-as.matrix(bcsv))))),out('00_input_audit/basis_check.txt'))
st <- qread(file.path(frozen,'01_input_state.qs'))
logmsg('Loading original Seurat object')
seu <- qread(input,nthreads=2)
ids <- st$cells
stopifnot(!anyDuplicated(ids),setequal(ids,colnames(seu)),identical(st$analysis_cell_hash,digest(sort(ids),algo='xxhash64')))
md <- as.data.table(seu@meta.data[ids,,drop=FALSE],keep.rownames='cell_id')
md[, sample := paste(dataset,orig.ident,sep='|')]
stopifnot(identical(md$sample,unname(st$analysis_sample_id_by_cell[ids])))
# Original specimen audit explicitly designates orig.ident as specimen_id. Do not infer patients from names.
sa <- fread(tail(files,1))
sa[, sample := paste(dataset,specimen_id,sep='|')]
stopifnot(!anyDuplicated(sa$sample))
mi <- match(md$sample,sa$sample)
stopifnot(!anyNA(mi))
md[, patient_status := sa$patient_id_status[mi]]
md[, patient_source := sa$patient_id_source[mi]]
md[, patient_recorded := sa$patient_uid[mi]]
# Only independently confirmed mappings may enter primary patient statistics. Derived/provisional remain audit-only.
md[, patient := fifelse(patient_status %in% c('confirmed','verified'),patient_recorded,NA_character_)]
md[, source := as.character(Type)]
md[, annotation := as.character(major_celltype)]
md[, original_cluster := as.character(epi_allgene_cluster)]
write_tsv(unique(md[,.(sample,dataset,orig.ident,sample_id,source,patient,patient_recorded,patient_status,patient_source)]),out('00_input_audit/sample_patient_mapping.tsv'))
write_tsv(md[,.(n_cells=.N),by=.(sample,dataset,source,patient,patient_status)],out('00_input_audit/sample_cell_counts.tsv'))
write_tsv(md[,.(n_cells=.N,n_samples=uniqueN(sample)),by=.(dataset,source)],out('00_input_audit/dataset_source_counts.tsv'))
write_tsv(md,out('00_input_audit/cell_manifest.tsv.gz'))
qsave_atomic(md,out('00_input_audit/cell_metadata.qs'))
writeLines(c('FULL is exactly the original object cells and frozen 01_input_state.qs$cells; no annotation/source/sample size filtering.',paste('n_cells',length(ids),'n_samples',uniqueN(md$sample),'n_datasets',uniqueN(md$dataset)),paste('frozen_cell_hash',st$analysis_cell_hash),'Patient map is derived/provisional in existing audit; no clinical confirmation assumed. Patient-level replication unavailable unless status explicitly verified/confirmed.','orig.ident is specimen_id in existing specimen audit; analysis_sample = dataset|orig.ident validated cell by cell.'),out('00_input_audit/FULL_identity.txt'))
assay <- seu[['RNA']]
layers <- Layers(assay)
la <- rbindlist(lapply(layers,function(l) {a <- LayerData(assay,layer=l);data.table(layer=l,n_features=nrow(a),n_cells=ncol(a),class=class(a)[1],features_hash=digest(rownames(a),'sha256'),cells_hash=digest(colnames(a),'sha256'))}))
write_tsv(la,out('00_input_audit/RNA_layers.tsv'))
cl <- grep('^counts($|\\.)',layers,value=TRUE)
if(length(cl)!=1L) stop('Counts layer topology requires explicit reconciliation; stopping before score computation')
x <- LayerData(assay,layer=cl)
stopifnot(inherits(x,'sparseMatrix'),setequal(colnames(x),ids),!anyDuplicated(rownames(x)),all(x@x>=0),all(is.finite(x@x)))
x <- x[,ids,drop=FALSE]
write_tsv(data.table(gene=rownames(x)),out('00_input_audit/RNA_feature_universe.tsv.gz'))
sig <- rbindlist(lapply(mp,function(m) {ii <- order(-basis[,m],rownames(basis),method='radix')[1:50];data.table(MP=m,gene=rownames(basis)[ii],rank=1:50,loading=basis[ii,m])}))
sig[, `:=`(ID_before=gene,ID_after=gene,present_in_RNA=gene %in% rownames(x),mapping='exact identifier; no renaming',source=file.path(frozen,'07_EpiMP_gene_basis.qs:sample_balanced_mean'))]
stopifnot(all(sig[,uniqueN(gene),by=MP]$V1==50L),all(nzchar(sig$gene)))
features <- setNames(lapply(mp,function(m) sig[MP==m]$gene),mp)
sig[,signature_sha256 := digest(paste(gene,collapse='\n'),algo='sha256',serialize=FALSE),by=MP]
write_tsv(sig,out('01_signatures/top50_signatures.tsv'))
qsave_atomic(features,out('01_signatures/signatures.qs'))
old <- fread(file.path(frozen,'07_EpiMP_top50_genes.csv'))
diff <- merge(sig[,.(MP,rank,canonical_gene=gene)],old[,.(MP=EpiMP,rank,frozen_gene=gene)],by=c('MP','rank'),all=TRUE)
write_tsv(diff,out('01_signatures/frozen_top50_comparison.tsv'))
write_tsv(sig[,.(n_signature=.N,n_present=sum(present_in_RNA),n_missing=sum(!present_in_RNA),signature_sha256=signature_sha256[1]),by=MP],out('01_signatures/coverage.tsv'))
ov <- CJ(MP1=mp,MP2=mp);ov[,overlap := mapply(function(a,b) length(intersect(features[[a]],features[[b]])),MP1,MP2)];ov[,Jaccard := overlap/(100-overlap)]
write_tsv(ov,out('01_signatures/overlap_jaccard_long.tsv'))
write_tsv(dcast(ov,MP1~MP2,value.var='overlap'),out('01_signatures/overlap_matrix.tsv'))
write_tsv(dcast(ov,MP1~MP2,value.var='Jaccard'),out('01_signatures/jaccard_matrix.tsv'))
md[,nFeature_counts := Matrix::colSums(x>0)];md[,nCount_counts := Matrix::colSums(x)]
qsave_atomic(md,out('00_input_audit/cell_metadata.qs'))
depth <- rbindlist(lapply(c('dataset','sample'),function(g) md[,.(n=.N,min=min(nFeature_counts),q05=as.numeric(quantile(nFeature_counts,.05)),median=as.numeric(median(nFeature_counts)),q95=as.numeric(quantile(nFeature_counts,.95)),max=max(nFeature_counts),fraction_below_1500=mean(nFeature_counts<1500)),by=c(g)][,unit:=g]),fill=TRUE)
write_tsv(depth,out('02_ucell_scores/depth_by_dataset_sample.tsv'))
detected <- sapply(features,function(g) Matrix::colSums(x[intersect(g,rownames(x)),,drop=FALSE]>0));rownames(detected)<-ids
qsave_atomic(detected,out('02_ucell_scores/signature_detected_counts.qs'))
qsave_atomic(t(qread(file.path(frozen,'08_cell_EpiMP_scores.qs'))$usage_raw[mp,ids,drop=FALSE]),out('02_ucell_scores/original_usage_raw.qs'))
rm(seu,assay,b);gc()
identity <- list(config=config_hash,signature=digest(features,'sha256'),cells=digest(ids,'sha256'),features=digest(rownames(x),'sha256'),input_manifest=digest(manifest,'sha256'),UCell=as.character(packageVersion('UCell')))
if(file.exists(out('02_ucell_scores/M_RAW.qs'))) {
 stopifnot(identical(qread(out('02_ucell_scores/score_identity.qs')),identity));M<-qread(out('02_ucell_scores/M_RAW.qs'));logmsg('Reusing verified M_RAW')
} else {
 qsave_atomic(identity,out('02_ucell_scores/score_identity.qs'))
 logmsg('UCell once on FULL: all RNA features, 14 top50 sets, chunks500 workers2')
 M <- UCell::ScoreSignatures_UCell(matrix=x,features=features,maxRank=1500,ties.method='average',missing_genes='impute',name=NULL,chunk.size=500,ncores=2,BPPARAM=BiocParallel::MulticoreParam(workers=2),force.gc=TRUE)
 M <- as.matrix(M);stopifnot(setequal(rownames(M),ids),setequal(colnames(M),mp));M<-M[ids,mp,drop=FALSE]
 stopifnot(all(is.finite(M)),!anyDuplicated(rownames(M)),nrow(M)==length(ids))
 qsave_atomic(M,out('02_ucell_scores/M_RAW.qs'))
}
rm(x);gc()
qc <- data.table(MP=mp,min=apply(M,2,min),max=apply(M,2,max),mean=colMeans(M),sd=apply(M,2,sd),zero_fraction=colMeans(M==0))
write_tsv(qc,out('02_ucell_scores/raw_score_QC.tsv'))
dup<-duplicated(M)|duplicated(M,fromLast=TRUE)
write_tsv(data.table(n=nrow(M),all_zero_fraction=mean(rowSums(M)==0),duplicate_excess_fraction=mean(duplicated(M)),any_duplicate_vector_fraction=mean(dup),n_unique_vectors=nrow(unique(M))),out('02_ucell_scores/vector_QC.tsv'))
stopifnot(all(qc$sd>0),all(is.finite(qc$sd)))
Z <- M
for(j in seq_along(mp)) Z[,j] <- (M[,j]-qc$mean[j])/qc$sd[j]
qsave_atomic(Z,out('02_ucell_scores/M_Z.qs'))
write_tsv(qc[,.(MP,mean_FULL=mean,sd_FULL=sd)],out('02_ucell_scores/global_scaling.tsv'))
writeLines(digest(qc[,.(MP,mean,sd)],'sha256'),out('02_ucell_scores/global_scaling.sha256'))
for(seed in 42:44) {
 set.seed(seed)
 sel <- unlist(lapply(split(seq_len(nrow(md)),md$sample),function(ii) if(length(ii)>1000) sample(ii,1000,replace=FALSE) else ii),use.names=FALSE)
 sel<-sort(sel);bal<-md[sel,.(cell_id,dataset,sample,patient,patient_status)]
 write_tsv(bal,out('05_balanced',paste0('seed',seed,'_cell_manifest.tsv.gz')))
 qsave_atomic(bal$cell_id,out('05_balanced',paste0('seed',seed,'_cell_ids.qs')))
 counts<-merge(md[,.(before=.N),by=sample],bal[,.(after=.N),by=sample],by='sample',all=TRUE)
 stopifnot(all(counts$after==pmin(counts$before,1000)),!anyDuplicated(bal$cell_id))
 write_tsv(counts,out('05_balanced',paste0('seed',seed,'_sample_counts.tsv')))
}
capture.output(sessionInfo(),file=out('logs/sessionInfo_scoring.txt'))
writeLines('complete',out('02_ucell_scores/COMPLETE'))
logmsg('Audit and scoring complete')
