source('results/mp10_nectin_external_validation_v1/code/common.R')
suppressPackageStartupMessages(library(edgeR))
for(line in c('HCV29','KYSE150','TE5','HeLa')){
 a<-switch(line,HCV29='GSE236742',HeLa='GSE267018','GSE143803')
 td<-file.path(TB,'perturbation',a,line);md<-fread(file.path(dirname(td),'sample_inclusion.tsv'));dir.create(td,showWarnings=FALSE,recursive=TRUE)
 if(line=='HCV29'){
  v<-as.data.table(readxl::read_excel(file.path(I,'GSE236742_Rawcounts.xlsx'),sheet='HCV29Rawcounts'));id<-v[[1]];xs<-names(v)[-1];t0<-sub('^HCV29(C|Ra|RM)([0-9]+)$','\\1',xs);gr<-c(C='DMSO',Ra='Rapa',RM='RMC')[t0]
  anno<-as.data.table(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,keys=unique(id),keytype='ENSEMBL',columns='SYMBOL'))[,.(gene=if(uniqueN(na.omit(SYMBOL))==1)na.omit(SYMBOL)[1] else NA_character_),by=ENSEMBL];g<-anno[match(id,ENSEMBL),gene]
  titles<-paste0('HCV29, ',c(C='Ctrl',Ra='Rapamycin, 10nM',RM='RMC-6272, 3nM')[t0],', Rep',sub('.*([0-9])$','\\1',xs));cts<-c(Rapamycin='Rapa-DMSO',RMC6272='RMC-DMSO')
 }else if(line=='HeLa'){
  v<-fread(file.path(I,'GSE267018_hela_fpkm.csv.gz'));id<-v[[1]];g<-v[[2]];xs<-names(v)[-(1:2)];titles<-xs;gr<-paste0(ifelse(grepl('SREBP1',xs),'KO1',ifelse(grepl('SREBP2',xs),'KO2','WT')),ifelse(grepl('FBS',xs),'FBS','Dep'));cts<-c(KO1_FBS='KO1FBS-WTFBS',KO2_FBS='KO2FBS-WTFBS',KO1_depletion='KO1Dep-WTDep',KO2_depletion='KO2Dep-WTDep',KO1_interaction='(KO1Dep-WTDep)-(KO1FBS-WTFBS)',KO2_interaction='(KO2Dep-WTDep)-(KO2FBS-WTFBS)')
 }else{
  v<-fread(file.path(I,paste0(a,'_',line,'_siSREBF1_readCounts_genes.txt.gz')));id<-g<-v[[1]];xs<-names(v)[-1];gr<-ifelse(grepl('siNC',xs),'NC','KD');titles<-paste0(line,'_',sub('siNC_','siNC',xs),' RNA-Seq');cts<-c(siSREBF1='KD-NC')
 }
 gsm<-md[match(titles,title),GSM];stopifnot(!anyNA(gsm),!anyDuplicated(gsm));mp<-data.table(matrix_column=xs,title=titles,GSM=gsm,group=gr);put(mp,file.path(td,'sample_mapping.tsv'))
 g<-rename(g);ok<-!is.na(g)&!g%in%c('','---');put(data.table(input_id=id,symbol=g,annotation_valid=ok),file.path(td,'gene_mapping.tsv.gz'))
 m<-rowsum(as.matrix(v[ok,..xs]),g[ok],reorder=TRUE);stopifnot(all(m>=0),all(is.finite(m)));group<-factor(gr);design<-model.matrix(~0+group);colnames(design)<-levels(group);cm<-makeContrasts(contrasts=unname(cts),levels=design);colnames(cm)<-names(cts)
 if(line=='HeLa'){mat<-libmatrix(log2(m+1));fit<-lmFit(mat,design)}else{stopifnot(all(m==round(m)));y<-DGEList(m);y<-y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE];y<-calcNormFactors(y);vo<-voom(y,design,plot=FALSE);mat<-vo$E;fit<-lmFit(vo,design)}
 fit<-eBayes(contrasts.fit(fit,cm),trend=line=='HeLa')
 ss<-c(sets[1:4],list(G2M=hall[signature=='HALLMARK_G2M_CHECKPOINT',gene]))
 gm<-strsplit(readLines(file.path(I,'h.all.v2024.1.Hs.symbols.gmt')),'\t');names(gm)<-sapply(gm,'[',1)
 ss$Apoptosis<-gm[['HALLMARK_APOPTOSIS']][-(1:2)];ss$Chol_homeostasis<-gm[['HALLMARK_CHOLESTEROL_HOMEOSTASIS']][-(1:2)]
 ps<-score(mat,lapply(ss,setdiff,familygenes));put(ps$cov,file.path(td,'program_coverage.tsv'));put(cbind(mp,as.data.table(t(ps$mat))),file.path(td,'program_scores.tsv'))
 # Each pathway gets its own residual variance and ordinary t interval.
 pf<-contrasts.fit(lmFit(ps$mat,design),cm)
 all<-list();pr<-list()
 for(j in seq_along(cts)){
  k<-names(cts)[j];z<-as.data.table(topTable(fit,coef=j,number=Inf,sort.by='none',confint=TRUE),keep.rownames='feature');z[,contrast:=k];setnames(z,c('logFC','CI.L','CI.R','P.Value','adj.P.Val'),c('effect','low','high','p','q_genome'));all[[k]]<-z
  se<-pf$stdev.unscaled[,j]*pf$sigma;df<-pf$df.residual;ef<-pf$coefficients[,j];pr[[k]]<-data.table(feature=rownames(pf),contrast=k,effect=ef,low=ef-qt(.975,df)*se,high=ef+qt(.975,df)*se,p=2*pt(-abs(ef/se),df))
 }
 z<-rbindlist(all);p<-rbindlist(pr);p[,q:=p.adjust(p,'BH'),by=contrast];p[,scale:='ssGSEA units'];fg<-merge(CJ(feature=familygenes,contrast=names(cts)),z,all.x=TRUE,by=c('feature','contrast'));fg[,q:=p.adjust(p,'BH',n=5),by=contrast];fg[,scale:=if(line=='HeLa')'delta log2(FPKM+1)' else 'log2FC voom'];fg[,status:=ifelse(is.na(p),'NE: not measured or filtered','estimated')]
 for(k in names(cts)){take<-cm[,k]!=0;arm<-paste0(rownames(cm)[take],'=',as.integer(table(group)[rownames(cm)[take]]),collapse=';');fg[contrast==k,samples_by_group:=arm];p[contrast==k,samples_by_group:=arm]}
 put(z,file.path(td,'genome.tsv.gz'));put(z[,.(n_tested=.N,n_BH05=sum(q_genome<.05),median_absolute_effect=median(abs(effect))),by=contrast],file.path(td,'global.tsv'))
 put(fg,file.path(td,'family.tsv'));put(p,file.path(td,'programs.tsv'))
 if(line=='HeLa')for(cond in c('FBS','depletion','interaction')){for(tp in c('family','programs')){tab<-if(tp=='family')fg else p;put(tab[endsWith(contrast,cond)],file.path(td,cond,paste0(tp,'.tsv')))}}
}
