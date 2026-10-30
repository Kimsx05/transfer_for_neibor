P<-'/home/data/t070721/codex_workspace/Bladder Metabolism';V1<-file.path(P,'mp10_spatial_communication_v1');BASE<-file.path(P,'mp10_spatial_communication_v2_all_sections')
.libPaths(unique(c(file.path(V1,'Rlib'),'/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(SpatialCellChat)})
setDTthreads(1);future::plan(future::sequential);source(file.path(BASE,'scripts/scc_compat.R'));out<-file.path(BASE,'logs/scc_compat_qa');dir.create(out,recursive=TRUE,showWarnings=FALSE)
fixed<-scc_make_computeCellDistance();checks<-list()
# Actual Sample2 coordinates can demonstrate the dimension truncation without RNA/network inference.
if(!file.exists(file.path(out,'distance_parity.tsv')))for(sec in c('GSE171351_A1','Inhousedata_T4469479','CanDis_2026_Sample2')){
 xy<-as.matrix(fread(file.path(BASE,'input_review/sections',sec,'coords.tsv'))[,.(x_um,y_um)])
 old<-SpatialCellChat::computeCellDistance(xy,interaction.range=200,contact.range=105,ratio=1,tol=1e-6)
 new<-fixed(xy,interaction.range=200,contact.range=105,ratio=1,tol=1e-6);n<-nrow(xy);m<-nrow(old$d.spatial)
 stopifnot(identical(dim(new$d.spatial),c(n,n)),identical(dim(new$adj.contact),c(n,n)))
 de<-max(abs(old$d.spatial-new$d.spatial[seq_len(m),seq_len(m)]));ce<-max(abs(old$adj.contact-new$adj.contact[seq_len(m),seq_len(m)]));stopifnot(de==0,ce==0)
 if(m<n){stopifnot(n==10062L,m==10061L,sum(new$d.spatial[n,])==0,new$adj.contact[n,n]==1);err<-tryCatch({Matrix::Diagonal(n)*old$d.spatial;'NO_ERROR'},error=function(e)conditionMessage(e));stopifnot(grepl('non-conformable',err));sccdummy<-Matrix::Diagonal(n)*new$d.spatial;stopifnot(all(dim(sccdummy)==c(n,n)));writeLines(err,file.path(out,'original_error_reproduced.txt'))}
 checks[[sec]]<-data.table(section=sec,full_spots=n,native_dimension=m,patched_dimension=n,distance_submatrix_error=de,contact_submatrix_error=ce,triggered=m<n)
 rm(old,new);gc()
}
if(length(checks))fwrite(rbindlist(checks),file.path(out,'distance_parity.tsv'),sep='\t')
# Native LR selection depends on feature membership and measured signaling genes;
# a sparse dummy with the exact cached row names verifies this without RNA copies.
sec<-'CanDis_2026_Sample2';IN<-file.path(BASE,'input_review/sections',sec);OUT<-file.path(BASE,'spatialcellchat',sec,'main');f<-fread(file.path(OUT,'spatial_variable_gene_test.tsv'))
chat<-methods::new('SpatialCellChat');chat@data.signaling<-sparseMatrix(i=integer(),j=integer(),x=numeric(),dims=c(nrow(f),2),dimnames=list(f$gene,c('QA_A','QA_B')))
chat@DB<-SpatialCellChat::subsetDB(CellChatDB.human,search=c('Secreted Signaling','ECM-Receptor','Cell-Cell Contact'),non_protein=FALSE)
genes<-fread(file.path(IN,'genes.tsv'))$gene;spot<-fread(file.path(IN,'meta.tsv'))$spot_id
r<-scc_restore_features(chat,IN,OUT,genes,spot);stopifnot(r$restored,length(r$chat@var.features$features)==604,nrow(r$chat@LR$LRsig)==1841)
r2<-scc_restore_features(chat,IN,OUT,genes,spot);stopifnot(r2$restored,identical(r$chat@var.features,r2$chat@var.features),identical(r$chat@LR,r2$chat@LR))
fwrite(data.table(section=sec,features=length(r$chat@var.features$features),LR=nrow(r$chat@LR$LRsig),checkpoint_bytes=file.info(file.path(OUT,'feature_selection_checkpoint.rds'))$size,LR_replay_exact=TRUE,checkpoint_restore_exact=TRUE),file.path(out,'checkpoint_validation.tsv'),sep='\t')
cat('COMPATIBILITY_DISTANCE_AND_CHECKPOINT_QA_PASS\n')
