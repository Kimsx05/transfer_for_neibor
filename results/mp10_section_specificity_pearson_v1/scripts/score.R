.libPaths(c(file.path(Sys.getenv('HOME'),'R/x86_64-pc-linux-gnu-library/4.5'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(Seurat);library(Matrix)})
setDTthreads(1);future::plan('sequential');options(future.globals.maxSize=20*1024^3)
o='results/mp10_section_specificity_pearson_v1';a=commandArgs(TRUE);m=fread(file.path(o,'tables/input_manifest.tsv'));gs=fread(file.path(o,'tables/gene_sets.tsv'));sets=lapply(split(gs$gene_symbol,gs$gs_name),unique);newnames=unique(fread(file.path(o,'tables/new5_gene_sets.tsv'))$gs_name)
alias=c(CCN2='CTGF',CAVIN3='PRKCDBP',ACKR3='CXCR7',GSDME='DFNA5',P3H2='LEPREL1')
score_one=function(i){
 r=m[i,];sec=r$section_id;p=file.path(o,'scores',paste0(sec,'_SCT_AddModuleScore.tsv.gz'));if(file.exists(p)&&file.exists(file.path(o,"scores",paste0(sec,"_gene_coverage.tsv"))))return(NULL);message(sec)
 con=gzfile(r$counts_path,'rt');x=as(readMM(con),'CsparseMatrix');close(con);genes=readLines(r$genes_path);cells=readLines(r$barcodes_path)
 if(anyDuplicated(genes)){u=unique(genes);x=sparseMatrix(i=match(genes,u),j=seq_along(genes),x=1,dims=c(length(u),length(genes)))%*%x;genes=u};dimnames(x)=list(genes,cells)
 elig=fread(file.path(o,'scores',paste0(sec,'_eligibility.tsv')));valid=cells %in% elig[valid_tissue==TRUE,barcode];stopifnot(sum(valid)>1)
 mapset=function(v,pool){z=v;for(k in seq_along(z))if(!z[k]%in%pool && z[k]%in%names(alias) && alias[[z[k]]]%in%pool)z[k]=alias[[z[k]]];intersect(unique(z),pool)}
 usets=lapply(sets[newnames],mapset,pool=genes);uc=UCell::ScoreSignatures_UCell(x,features=usets,name=NULL,ncores=1,chunk.size=50,maxRank=1500,ties.method='average');fwrite(data.table(section_id=sec,barcode=cells,uc[cells,,drop=FALSE]),file.path(o,'scores',paste0(sec,'_new5_UCell.tsv.gz')),sep='\t')
 obj=CreateSeuratObject(x[,valid,drop=FALSE],min.cells=0,min.features=0,assay='RNA');rm(uc);gc()
 set.seed(20261104)
 obj=SCTransform(obj,assay='RNA',new.assay.name='SCT',vst.flavor='v2',ncells=min(5000,ncol(obj)),seed.use=20261104,return.only.var.genes=TRUE,conserve.memory=TRUE,verbose=FALSE)
 layer=GetAssayData(obj,assay='SCT',layer='data');pool=rownames(layer);aset=lapply(sets,mapset,pool=pool);good=lengths(aset)>0
 # data is log1p corrected counts; entire available SCT gene pool, not variable features.
 stopifnot(nrow(layer)>length(VariableFeatures(obj)),all(dim(layer)==dim(GetAssayData(obj,assay='SCT',layer='counts'))))
 co=GetAssayData(obj,assay='SCT',layer='counts');take=seq_len(min(20,nrow(co)));stopifnot(isTRUE(all.equal(as.matrix(layer[take,,drop=FALSE]),log1p(as.matrix(co[take,,drop=FALSE])),check.attributes=FALSE)))
 obj=AddModuleScore(obj,features=aset[good],pool=pool,nbin=24,ctrl=100,assay='SCT',slot='data',seed=20261104,search=FALSE,name='AMS85_')
 ans=as.data.table(obj[[]][,paste0('AMS85_',seq_len(sum(good))),drop=FALSE]);setnames(ans,paste0('AMS_SCT__',names(aset)[good]));ans=cbind(data.table(section_id=sec,barcode=colnames(obj)),ans)
 for(n in names(aset)[!good])ans[[paste0('AMS_SCT__',n)]]=NA_real_
 fwrite(ans,paste0(p,'.tmp'),sep='\t',compress='gzip');file.rename(paste0(p,'.tmp'),p)
 coverage=rbindlist(lapply(names(sets),function(n){matched=aset[[n]];raw=mapset(sets[[n]],genes);data.table(section_id=sec,program=n,n_original=length(sets[[n]]),n_RNA_available=length(raw),n_RNA_detected=sum(Matrix::rowSums(x[raw,valid,drop=FALSE])>0),n_SCT_available=length(matched),missing_SCT=paste(setdiff(sets[[n]],matched),collapse='|'),used_symbols=paste(matched,collapse='|'),n_SCT_pool=length(pool),n_scored_spots=sum(valid),assay='SCT',layer='data',nbin=24,ctrl=100,seed=20261104,status=ifelse(length(matched)>0,'scored','empty_signature'))}))
 fwrite(coverage,file.path(o,'scores',paste0(sec,'_gene_coverage.tsv')),sep='\t');fwrite(data.table(section_id=sec,Seurat=as.character(packageVersion('Seurat')),sctransform=as.character(packageVersion('sctransform')),UCell=as.character(packageVersion('UCell')),SCT_flavor='v2',fit='native R NB fallback; glmGamPoi unavailable',min_cells=5,ncells=min(5000,sum(valid)),full_SCT_data_genes=length(pool),n_tissue=sum(valid)),file.path(o,'scores',paste0(sec,'_software.tsv')),sep='\t')
 rm(obj,x,layer,co);gc();message(sec,' COMPLETE')
}
for(i in seq(as.integer(a[1]),nrow(m),by=as.integer(a[2])))tryCatch(score_one(i),error=function(e){message('FAIL ',m$section_id[i],': ',conditionMessage(e));writeLines(conditionMessage(e),file.path(o,'logs',paste0(m$section_id[i],'_score_failure.txt')))})
