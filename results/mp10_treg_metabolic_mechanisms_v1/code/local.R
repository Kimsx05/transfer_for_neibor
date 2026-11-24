source(file.path(Sys.getenv('B_OUT'),'code/common.R'))
# The object is already aggregated. run.sh additionally locks this read through R process exit.
b<-qread(file.path(old,'objects/pseudobulk.qs'),nthreads=1)
p<-as.data.table(b$profiles);tr<-p[profile_type=='Treg']; ids<-tr$profile_id
stopifnot(!anyDuplicated(tr$sample_uid))
counts<-as.matrix(b$counts[,ids,drop=FALSE]);ref<-tr$main_patient&tr$n_cells>=50
stopifnot(all(counts==floor(counts)),all(colSums(counts)==tr$total_counts))
keep0<-rowSums(counts)>0
y<-calcNormFactors(DGEList(counts[keep0,,drop=FALSE]))
cpm0<-cpm(y);keep<-rowMeans(cpm0[,ref]>=1)>=.25
logmat<-cpm(y,log=TRUE,prior.count=2)[keep,,drop=FALSE]
saveRDS(logmat,file.path(out,'tables/treg_logCPM.rds'))
saveT(data.table(gene=rownames(counts),any_count=keep0,primary_CPM1_fraction=rowMeans(cpm(DGEList(counts,lib.size=y$samples$lib.size*y$samples$norm.factors))[,ref]>=1),rank_background=rownames(counts)%in%rownames(logmat)),'background.tsv.gz')
aud<-readt(old,'patient_relationship_audit.tsv');cv<-merge(aud,tr[,.(sample_uid,n_cells,total_counts,detected_genes,Q50,Q95,Q50_z,Q95_z)],by='sample_uid',all.x=TRUE)
cv[is.na(n_cells),n_cells:=0L];cv[,reason:=fcase(n_cells==0,'no_Treg_profile',!main_specimen,'confirmed_repeat_b',n_cells<30,'below_30',n_cells<50,'30_only',default='main_50')];cv[,counts_available:=n_cells>0]
anno<-readt(old,'analysis_cell_metadata.tsv.gz')[analysis_celltype_c2l_final=='Treg',.(subtype_counts=paste(names(table(subtype_v2)),as.integer(table(subtype_v2)),collapse=';'),original_state_counts=paste(names(table(state_v2)),as.integer(table(state_v2)),collapse=';'),final_label=unique(analysis_celltype_c2l_final)),by=sample_uid]
saveT(merge(cv,anno,by='sample_uid',all.x=TRUE),'coverage.tsv')
markers<-c('CD3D','CD3E','CD4','FOXP3','IL2RA','CTLA4','IKZF2','TIGIT','MKI67','CCR8','ABCA1','ABCG1','LDLR','SCARB1','NPC1','APOA1','APOA2','APOB','APOC1','APOE','CLU','CD36','FASN','HMGCR','SREBF1','SREBF2','SCAP')
q<-p[!profile_type%in%c('H05','H10','H50','L50')]
present<-intersect(markers,rownames(counts))
dt<-rbindlist(lapply(present,function(g){d<-copy(q[,.(sample_uid,dataset,patient_uid_recorded,main_patient,profile_type,profile_id,n_cells,total_counts)]);d[,gene:=g];d[,count:=as.numeric(b$counts[g,profile_id])];d[,positive_cells:=as.numeric(b$detected_cells[g,profile_id])];d[,detection:=positive_cells/n_cells];d[,cpm:=1e6*count/total_counts];d[,within_sample_count_share:=if(sum(count)>0)count/sum(count) else NA_real_,by=sample_uid];d}))
saveT(dt,'marker_sources_by_specimen.tsv.gz')
saveT(dt[profile_type=='Treg'],'treg_identity_and_lipid_markers.tsv')
saveT(dt[main_patient==TRUE,.(specimens=.N,median_detection=median(detection),median_CPM=median(cpm),median_count_share=median(within_sample_count_share,na.rm=TRUE),total_counts=sum(count)),by=.(gene,profile_type)][order(gene,-total_counts)],'source_summary.tsv')
sets<-readt(v2,'official_gene_sets.tsv.gz')[signature%in%c('CHOL_IMPORT','LCFA_IMPORT','CHOL_SYNTH')]
sets<-unique(sets,by=c('signature','gene'))
cov<-readt(v2,'signature_coverage.tsv')[family=='Treg'&mask=='Full']
for(s in unique(sets$signature)){
 ix<-which(sets$signature==s);gs<-sets$gene[ix];og<-strsplit(cov[signature==s,covered_genes],';',fixed=TRUE)[[1]]
 sets[ix,old_used:=gs%in%og];sets[ix,count_feature:=gs%in%rownames(counts)];sets[ix,any_detected:=gs%in%rownames(cpm0)];sets[ix,new_used:=gs%in%rownames(logmat)]
}
saveT(sets,'uptake_gene_coverage.tsv')
saveT(sets[,.(expected=.N,old=sum(old_used),raw_feature=sum(count_feature),any_detect=sum(any_detected),new_CPM_background=sum(new_used),eligible=sum(new_used)>=5&mean(new_used)>=.5),by=signature],'uptake_coverage.tsv')
mem<-readt(old,'model_membership.tsv.gz')[model_id%in%sprintf('M%04d',46:54)]
saveT(mem,'original_receiver_membership.tsv.gz')
saveT(tr,'profiles.tsv')
saveT(data.table(gene=rownames(logmat),TMM_logCPM_background=TRUE),'rank_background.tsv.gz')
cat('BACKGROUND',nrow(logmat),'primary',sum(ref),'of',nrow(tr),'\n')
print(unique(sets$signature))
