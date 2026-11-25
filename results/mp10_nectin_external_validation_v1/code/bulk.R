source('results/mp10_nectin_external_validation_v1/code/common.R')
a<-commandArgs(TRUE)[1];td<-file.path(TB,'bulk',a);dir.create(td,recursive=TRUE,showWarnings=FALSE)
readSeries<-function(path){v<-fread(path,skip='"ID_REF"',na.strings=c('null','NA',''));v[!grepl('^!',v[[1]])]}
if(a=='TCGA'){
  tc<-function(n)if(file.exists(file.path('/tmp',n)))file.path('/tmp',n) else file.path(I,n)
  raw<-fread(tc('mp10_tcga_complete.txt'));clinical<-fread(tc('mp10_clinical.txt'),skip='PATIENT_ID\tSAMPLE_ID')
  md<-clinical[SAMPLE_TYPE=='Primary'&SAMPLE_ID%in%names(raw)&grepl('-01$',SAMPLE_ID)][order(SAMPLE_ID)][!duplicated(PATIENT_ID)]
  pur<-fread(tc('mp10_abs.tsv'))[is.finite(purity),.(purity=median(purity)),by=array];md<-merge(md,pur,by.x='SAMPLE_ID',by.y='array',all.x=TRUE)
  md[,GSM:=SAMPLE_ID];md[,patient_sample_id:=PATIENT_ID];ids<-md$GSM;g<-raw$Hugo_Symbol;g[g%in%'PVRL2']<-'NECTIN2';ok<-!is.na(g)&nzchar(g);x<-log2(as.matrix(raw[ok,..ids])+1);g<-g[ok];m<-libmatrix(rowsum(x,g,reorder=TRUE)/as.numeric(table(g)));me<-m;rownames(me)<-rename(rownames(me))
  put(md,file.path(td,'cohort_manifest/sample_inclusion.tsv'))
}else{
 md<-fread(file.path(td,'cohort_manifest/sample_inclusion.tsv'))[included==1];ids<-md$GSM
 v<-readSeries(file.path(I,paste0(a,'_series_matrix.txt.gz')))
 if(a=='GSE31684'){an<-fread('results/mp10_resistance_validation_v2_expanded/00_plan_audit/GPL570_annotation.tsv',select=c('ID','Gene Symbol'));setnames(an,c('id','raw_symbol'))}else{
 gpl<-md$platform[1];p<-file.path(I,paste0(gpl,'.soft')); if(gpl=='GPL6947' && file.exists(file.path(I,'GPL6947_family.soft.gz')))p<-file.path(I,'GPL6947_family.soft.gz');if(gpl=='GPL6947' && file.exists(file.path(I,'GPL6947_extracted.soft')))p<-file.path(I,'GPL6947_extracted.soft');an0<-fread(p,skip='ID\t',sep='\t',quote='',fill=TRUE);sym<-intersect(c('Symbol','Gene symbol','Gene Symbol','ILMN_Gene'),names(an0))[1];stopifnot(!is.na(sym));an<-an0[!startsWith(ID,'!'),.(id=ID,raw_symbol=get(sym))]
 }
 an[,gene:=vapply(strsplit(fifelse(is.na(raw_symbol),'',raw_symbol),' */// *| *// *| *; *'),function(z){s<-unique(rename(z));s<-s[!s%in%c('','---','NA','null')];if(length(s)==1)s else NA_character_},'')]
 pr<-an[match(v[[1]],id)];pr[,input_id:=v[[1]]];put(pr,file.path(td,'probe_mapping.tsv.gz'))
 x<-as.matrix(v[,..ids]);storage.mode(x)<-'double';me<-prof(x,pr);m<-libmatrix(me)
 if(a=='GSE32894'){
  # Tumour matrix and technical variables, without across-sample mean centering.
  raw<-fread(file.path(I,'GSE32894_non-normalized_308UCsamples.txt.gz'),skip='ID_REF\tUC_',na.strings=c('null','NA',''))
  names0<-sub('^Bladder Tumor ','',md$source);stopifnot(all(names0%in%names(raw)),!anyDuplicated(names0))
  x<-as.matrix(raw[,..names0]);colnames(x)<-ids;rownames(x)<-raw[[1]];ann<-an[match(raw[[1]],id)]
  proto<-fread(file.path(D,'config/extraction_32894.tsv'))[match(ids,GSM),protocol]
  group<-ifelse(grepl('Allprep',proto,ignore.case=TRUE),'AllPrep','RNeasy');plate<-md$labeling_plate
  stopifnot(length(unique(group))==2,'P4'%in%plate)
  aux<-fread(file.path(I,'GSE32894_reps_normals_preprocess_metadata.txt.gz'),skip='Sample name\ttitle',sep='\t',header=TRUE,fill=TRUE)
  aux<-aux[!is.na(get('Characteristics: Labeling_Plate')) & nzchar(get('Characteristics: Labeling_Plate'))]
  extra<-fread(file.path(I,'GSE32894_reps_normals_preprocess_matrix-non-normalized.txt.gz'),skip='ID_REF\tN_',na.strings=c('null','NA',''))
  put(data.table(metadata_id=aux[[1]],matched=aux[[1]]%in%names(extra)),file.path(td,'preprocess_aux_mapping.tsv'));aux<-aux[aux[[1]]%in%names(extra)];ax<-aux[[1]];stopifnot(all(ax%in%names(extra)),!anyDuplicated(ax));xx<-as.matrix(extra[match(raw[[1]],extra[[1]]),..ax]);x<-cbind(x,xx);group<-c(group,ifelse(grepl('Allprep',aux[['extract protocol']],ignore.case=TRUE),'AllPrep','RNeasy'));plate<-c(plate,aux[['Characteristics: Labeling_Plate']])
  keep<-rowSums(!is.finite(x))==0 ;x<-x[keep,,drop=FALSE];ann<-ann[keep]
  med<-function(x)matrixStats::rowMedians(x)
  correct<-function(x,group,ref){base<-med(x[,group==ref,drop=FALSE]);for(v0 in setdiff(unique(group),ref)){j<-group==v0;s<-med(x[,j,drop=FALSE]);good<-is.finite(s)&s>0 & base>0;x[good,j]<-x[good,j,drop=FALSE]*(base[good]/s[good]);x[!good,j]<-NA_real_};x}
  x<-correct(x,group,'RNeasy');x<-correct(x,plate,'P4');ok<-rowSums(!is.finite(x))==0;x<-x[ok,,drop=FALSE];ann<-ann[ok]
  x<-normalizeBetweenArrays(x,method='quantile');ok<-matrixStats::rowMins(x)>0 & !is.na(ann$gene);x<-log2(x[ok,,drop=FALSE]);ann<-ann[ok];intensity<-rowMeans(x);sel<-intensity>=median(intensity)
  m<-libmatrix(prof(x[sel,ids,drop=FALSE],ann[sel]));colnames(m)<-ids
  cg<-intersect(rownames(m),rownames(me));qc<-data.table(gene=cg,cor_with_published=as.numeric(vapply(cg,function(g)suppressWarnings(cor(m[g,],me[g,],use='pairwise.complete.obs')),0)));put(qc,file.path(td,'reconstruction_concordance.tsv.gz'))
  put(data.table(probe=ann$id,gene=ann$gene,precenter_mean=intensity,upper_half=sel),file.path(td,'reconstruction.tsv.gz'))
 }
}
stopifnot(ncol(m)==length(ids),!anyDuplicated(ids));out<-data.table(GSM=ids);coverage<-list()
for(g in unique(c(rna,'CD3D'))){name<-if(g%in%c('SREBF1','SREBF2'))paste0(g,'_RNA') else g;out[,(name):=if(g%in%rownames(me))as.numeric(me[g,]) else NA_real_]}
extras<-list()
for(pair in c(core,list(c('SREBF1','SREBF2')),unlist(lapply(names(sets)[1:4],function(n)lapply(familygenes,function(h)c(n,h))),recursive=FALSE))){
 n<-pair[1];b<-pair[2];s1<-sets[[n]];s2<-if(b%in%names(sets))sets[[b]] else b;shared<-intersect(s1,s2)
 if(length(shared)){key<-paste(n,b,sep='__');extras[[paste0(key,'_A')]]<-setdiff(s1,shared);if(b%in%names(sets))extras[[paste0(key,'_B')]]<-setdiff(s2,shared)}
}
for(ma in c('Full','M20')){
 s<-if(ma=='Full')c(sets,extras) else lapply(sets[-1],setdiff,top);mm<-if(ma=='Full')m else m[!rownames(m)%in%top,,drop=FALSE]
 z<-score(mm,s);z$cov[,mask:=ma];coverage[[ma]]<-z$cov
 for(n in rownames(z$mat)){name<-if(ma=='Full')n else paste0(n,'_M20');out[,(name):=as.numeric(z$mat[n,])]}
}
if(a=='TCGA'){
 for(ma in c('Full','M20')){
  orig<-fread(file.path(V,'plot_data',paste0('TCGA_primary_program_',ma,'_MP10_mTORC1.tsv')))[match(ids,sample_uid)]
  for(n in nodes){nm<-if(ma=='M20'&&n%in%c('SREBF1','SREBF2','mTORC1'))paste0(n,'_M20') else n
   delta<-max(abs(orig[[n]]-out[[nm]]));if(delta>1e-5)warning(nm,' differs ',delta);out[,(nm):=orig[[n]]]}
 }
}
cov<-rbindlist(coverage);cov[,c('cohort','background'):=list(a,ifelse(mask=='Full',nrow(m),sum(!rownames(m)%in%top)))];put(cov,file.path(td,'score_coverage.tsv'))
put(merge(md,out,by='GSM'),file.path(td,'scores.tsv'));cat(a,nrow(m),nrow(out),'\n')

