.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/refdir/Rlib',.libPaths()))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(GSVA);library(edgeR);library(org.Hs.eg.db);library(AnnotationDbi)})
setDTthreads(2);B='/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_resistance_validation_v1';P=dirname(dirname(B));S=file.path(B,'00_audit/sources');cfg=fromJSON(file.path(B,'config/frozen_v1.json'));set.seed(cfg$seed)
stopifnot(digest::digest(file=file.path(B,'config/frozen_v1.json'),algo='sha256')==readLines(file.path(B,'config/frozen_v1.sha256')))
wt=function(x,p)fwrite(x,file.path(B,p),sep='\t',na='NA');rd=function(p)fread(file.path(B,p))
map_score=function(x,original,symbol,cohort,counts=TRUE){
 stopifnot(nrow(x)==length(original),!anyDuplicated(colnames(x)),all(is.finite(x)))
 known=keys(org.Hs.eg.db,keytype='SYMBOL');mapped=symbol;method=ifelse(symbol%in%known,'exact_SYMBOL','unmapped_author_symbol')
 unknown=unique(symbol[!symbol%in%known & !is.na(symbol) & nzchar(symbol)]);valid=intersect(unknown,keys(org.Hs.eg.db,keytype='ALIAS'))
 if(length(valid)){al=as.data.table(AnnotationDbi::select(org.Hs.eg.db,valid,'SYMBOL','ALIAS'));al=unique(al[!is.na(SYMBOL)]);al=al[,if(uniqueN(SYMBOL)==1).(SYMBOL=SYMBOL[1])else .(SYMBOL=NA_character_),by=ALIAS];ii=match(symbol,al$ALIAS);ok=!is.na(ii)&!is.na(al$SYMBOL[ii])&!symbol%in%known;mapped[ok]=al$SYMBOL[ii[ok]];method[ok]='unique_ALIAS'}
 bad=is.na(mapped)|!nzchar(mapped);mapped[bad]=paste0('UNMAPPED_',original[bad]);method[bad]='original_ID_background'
 ma=data.table(original_id=original,author_symbol=symbol,mapped_symbol=mapped,mapping=method);wt(ma,paste0('00_audit/',cohort,'_gene_mapping.tsv'))
 if(counts){stopifnot(all(x>=0));xx=rowsum(x,mapped,reorder=FALSE);y=DGEList(xx);y=calcNormFactors(y,method='TMM');expr=cpm(y,log=TRUE,prior.count=1);wt(data.table(sample=colnames(x),library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors),paste0('00_audit/',cohort,'_normalization.tsv'))}else{u=unique(mapped);ix=split(seq_along(mapped),factor(mapped,levels=u));expr=t(vapply(ix,function(i)if(length(i)==1)x[i,] else apply(x[i,,drop=FALSE],2,median),numeric(ncol(x))));colnames(expr)=colnames(x);xx=expr}
 mp=cfg$genes$MP10;present=mp%in%rownames(expr);stopifnot(sum(present)>=16)
 wt(data.table(gene=mp,in_matrix=present,detection_rate=vapply(mp,function(g)if(g%in%rownames(xx))mean(xx[g,]>0)else NA_real_,numeric(1)),detection_definition=if(counts)'count>0'else 'positive normalized intensity, not biological detection'),paste0('00_audit/',cohort,'_MP10_coverage.tsv'))
 wt(data.table(cohort=cohort,n_original=nrow(x),n_exact_or_alias=sum(method%in%c('exact_SYMBOL','unique_ALIAS')),n_background=nrow(expr),n_signature=20,n_present=sum(present),missing=paste(mp[!present],collapse=';'),input_min=min(x),input_max=max(x),input_q50=median(x),preprocessing=if(counts)'TMM log2CPM prior.count=1'else 'author RMA quantile normalized ComBat; no re-log'),paste0('00_audit/',cohort,'_matrix_audit.tsv'))
 param=ssgseaParam(expr,list(MP10=mp),minSize=16,maxSize=100000,alpha=.25,normalize=TRUE,checkNA='yes');score=as.numeric(gsva(param,verbose=FALSE,BPPARAM=BiocParallel::SerialParam())[1,]);stopifnot(length(score)==ncol(expr),all(is.finite(score)),sd(score)>0)
 list(scores=data.table(matrix_id=colnames(expr),MP10=score),expr=expr)
}
endpoint=function(t,n,surgery=TRUE){
 t=toupper(as.character(t));n=toupper(as.character(n));good=pcr=rep(NA_integer_,length(t));npos=n%in%c('1','2','3','N1','N2','N3');nzero=n%in%c('0','N0');tzero=t%in%c('0','T0');tlow=t%in%c('A','IS','1','TA','TIS','T1','NONINVASIVE_OR_CIS');thigh=grepl('^(T)?[234]',t)
 good[npos|thigh]=0L;pcr[npos|thigh|tlow]=0L;good[(tzero|tlow)&nzero]=1L;pcr[tzero&nzero]=1L;good[!surgery]=NA_integer_;pcr[!surgery]=NA_integer_;data.table(good=good,pCR=pcr)
}
# Risk-specific endpoint unit checks, not simulated scientific results.
z=endpoint(c('0','TA','1','2','0','IS','TX','TX'),c('0','0','0','NX','NX','NX','0','1'));stopifnot(identical(z$good,c(1L,1L,1L,0L,NA_integer_,NA_integer_,NA_integer_,0L)),identical(z$pCR,c(1L,0L,0L,0L,NA_integer_,0L,NA_integer_,0L)))

z=endpoint(c('2B','3B','4B','IS','A'),c('0','0','0','0','X'));stopifnot(identical(z$good,c(0L,0L,0L,1L,NA_integer_)),identical(z$pCR,c(0L,0L,0L,0L,0L)))
