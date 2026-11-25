source('results/mp10_nectin_external_validation_v1/code/common.R')
suppressPackageStartupMessages(library(org.Hs.eg.db));tg<-AnnotationDbi::select(org.Hs.eg.db,keys=c('TIGIT','PVRIG','CD226'),keytype='SYMBOL',columns='ENTREZID')
# Resolve receptor symbols by unambiguous numeric Entrez IDs from the original platform.
for(a in cohorts){
 td<-file.path(TB,'bulk',a);d<-fread(file.path(td,'scores.tsv'));pr<-fread(file.path(td,'probe_mapping.tsv.gz'));md<-fread(file.path(td,'cohort_manifest/sample_inclusion.tsv'))
 if(a=='GSE13507'){d[,T_category_invasion:=ifelse(grepl('^T[234]',stage),'MIBC','NMIBC')];d[,invasion:=ifelse(invasiveness=='muscle invasive','MIBC','NMIBC')]}
 series<-fread(file.path(I,paste0(a,'_series_matrix.txt.gz')),skip='"ID_REF"');ids<-d$GSM
 if(a=='GSE31684'){pp<-fread('results/mp10_resistance_validation_v2_expanded/00_plan_audit/GPL570_annotation.tsv');idcol<-grep('ENTREZ|ENTREZ_GENE|Gene ID',names(pp),ignore.case=TRUE,value=TRUE)[1]}else{path<-file.path(I,if(a=='GSE32894')'GPL6947_extracted.soft' else 'GPL6102.soft');if(!file.exists(path))path<-file.path(I,'GPL6947.soft');pp<-fread(path,skip='ID\t',sep='\t',fill=TRUE,quote='');idcol<-'Entrez_Gene_ID'}
 map<-list()
 for(k in tg$SYMBOL){t<-tg[tg$SYMBOL==k,'ENTREZID'];probes<-if(!is.na(idcol))as.character(pp[[1]][trimws(as.character(pp[[idcol]]))%in%t]) else character()
  sel<-pr[(input_id%in%probes|gene==k)&!grepl('///',fifelse(is.na(raw_symbol),'',raw_symbol))];if(nrow(sel)){x<-as.matrix(series[match(sel$input_id,series[[1]]),..ids]);d[,(k):=apply(x,2,median,na.rm=TRUE)]}
  map[[k]]<-data.table(symbol=k,entrez=t,platform_probes=paste(probes,collapse=';'),used=paste(sel$input_id,collapse=';'),raw_symbol=paste(unique(sel$raw_symbol),collapse=';'),n_samples=sum(is.finite(d[[k]])))
 }
 put(rbindlist(map),file.path(td,'receptor_mapping.tsv'));put(d,file.path(td,'scores.tsv'))
 cv<-fread(file.path(td,'score_coverage.tsv'));univ<-unique(pr$gene)
 if(a=='GSE32894'){gg<-fread(file.path(I,'GPL6947_extracted.soft'),skip='ID\t',sep='\t',fill=TRUE,quote='');univ<-rename(gg$Symbol)}
 cv[,platform_annotated_n:=vapply(strsplit(members,';',fixed=TRUE),function(g)sum(g%in%univ),integer(1))];put(cv,file.path(td,'score_coverage.tsv'))
 amp<-d[,lapply(.SD,function(v)sum(is.finite(v))),.SDcols=c(rna[!rna%in%c('SREBF1','SREBF2')],'Tlineage','Treg')];put(melt(amp,measure.vars=names(amp),variable.name='feature',value.name='measured_samples'),file.path(td,'expression_coverage.tsv'))
 if(a=='GSE13507'){cc<-as.data.table(readxl::read_excel(file.path(I,'GSE13507_clinical_info.xls')));cc[,patient_sample_id:=get('Sample name')];ck<-merge(md[included==1,.(GSM,patient_sample_id,stage)],cc,by='patient_sample_id');stopifnot(nrow(ck)==165);put(ck,file.path(td,'clinical_crosscheck.tsv'))}
 if(a=='GSE31684'){cc<-fread(file.path(I,'GSE31684_table_of_clinical_details.txt.gz'));ck<-merge(md[,.(GSM,stage_used,treatment_used)],cc,by.x='GSM',by.y='GEO');stopifnot(nrow(ck)==93,all(tolower(ck$treatment_used)==tolower(ck$PreRC_Chemo)),all(ck$stage_used==ck$RC_Stage));put(ck,file.path(td,'clinical_crosscheck.tsv'))}
}
old<-rbindlist(lapply(list.files(file.path(V,'statistics'),pattern='^TCGA_primary',full.names=TRUE),function(f){v<-fread(f)[method=='pearson'];v[,source_file:=basename(f)];v}));put(old,file.path(TB,'TCGA_original_Pearson.tsv'))

