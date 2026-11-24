suppressPackageStartupMessages({library(data.table);library(edgeR)})
setDTthreads(1);O=commandArgs(TRUE)[1];raw=fread(Sys.glob('cache/mp10_public_resources_v1/GSE182638/*raw_readcounts*'))
map=fread(file.path(O,'tables/GSE182638_samples.tsv')); stopifnot(nrow(map)==18)
z=fread(file.path(O,'scripts/modules.tsv'));sets=setNames(strsplit(z$genes,' '),z$module);top=sets$MP10
x=as.matrix(raw[,-1]);x=rowsum(x,raw[[1]]);stopifnot(all(x>=0),all(x==round(x)))
res=list();detail=list();profiles=list()
for(line in unique(map[['cell line']])){
 m=map[get('cell line')==line];y=calcNormFactors(DGEList(x[,m$sample]));mat=cpm(y,log=TRUE,prior.count=2)
 genes=intersect(unique(c(unlist(sets),'CHMP2A')),rownames(mat)[rowSums(y$counts)>0])
 for(mask in c('Full','M20-masked'))for(mod in names(sets)){
  members=sets[[mod]];if(mask!='Full')members=setdiff(members,top);g=intersect(members,genes);if(!length(g))next
  v=colMeans(mat[g,,drop=FALSE]);df=data.table(m,score=v);profiles[[length(profiles)+1]]=data.table(df,module=mod,mask=mask)
  w=dcast(df,`cell passage`~arm,value.var='score');stopifnot(nrow(w)==3)
  for(ab in list(c('RSL3','UT'),c('FR','RSL3'))){
   dif=w[[ab[1]]]-w[[ab[2]]];tt=t.test(dif)
   res[[length(res)+1]]=data.table(line,module=mod,mask,contrast=paste(ab,collapse='-'),n_passages=length(dif),effect=mean(dif),lo=tt$conf.int[1],hi=tt$conf.int[2],p=tt$p.value,members=paste(g,collapse=' '),closer_to_UT=abs(mean(w$FR-w$UT))<abs(mean(w$RSL3-w$UT)),closer_passages=sum(abs(w$FR-w$UT)<abs(w$RSL3-w$UT)))
  }
 }
 design=model.matrix(~factor(m[['cell passage']])+factor(m$arm,levels=c('UT','RSL3','FR')))
 keep=filterByExpr(y,design=design);yy=calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);fit=glmQLFit(estimateDisp(yy,design,robust=TRUE),design,robust=TRUE)
 for(j in 1:2){co=rep(0,ncol(design));co[4:5]=if(j==1)c(1,0)else c(-1,1)
 tab=as.data.table(topTags(glmQLFTest(fit,contrast=co),n=Inf)$table,keep.rownames='gene')
 detail[[length(detail)+1]]=tab[gene %in% genes][,list(line=line,contrast=c('RSL3-UT','FR-RSL3')[j],gene,logFC,PValue,FDR)]
 }
}
r=rbindlist(res);r[,q:=p.adjust(p,'BH'),by=.(line,mask,contrast)]
fwrite(r,file.path(O,'tables/GSE182638_modules.tsv'),sep='\t');fwrite(rbindlist(detail),file.path(O,'tables/GSE182638_genes.tsv'),sep='\t');fwrite(rbindlist(profiles),file.path(O,'plot_data/GSE182638_profiles.tsv'),sep='\t')

