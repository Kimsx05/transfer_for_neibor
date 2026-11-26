suppressPackageStartupMessages(library(Seurat))
library(data.table)
setDTthreads(1)
for (s in commandArgs(TRUE)){
 p=file.path("niche",s); z=fread(file.path(p,"niche_input.tsv.gz")); ids=z$id
 x=Matrix::sparseMatrix(i=rep(1L,length(ids)),j=seq_along(ids),x=1,dimnames=list("placeholder",ids));obj=CreateSeuratObject(x)
 obj$ct=z$ct;co=data.frame(x=z$x,y=z$y,row.names=ids)
 obj[["fov"]]=CreateFOV(coords=co,type="centroids",assay="RNA")
 yy=list()
 for(seed in c(261126,261127,261128)){
  set.seed(seed);m=BuildNicheAssay(obj,fov="fov",group.by="ct",neighbors.k=20,niches.k=4,iter.max=100)
  yy[[as.character(seed)]]=as.integer(m$niches)
 }
 fwrite(data.frame(id=ids,niche=yy[[1]]),file.path(p,"BuildNicheAssay.tsv.gz"),sep="\t")
 fwrite(data.frame(seed_pair=c("1-2","1-3","2-3"),ARI=combn(3,2,function(i)mclust::adjustedRandIndex(yy[[i[1]]],yy[[i[2]]]))),file.path(p,"BuildNicheAssay_stability.tsv"),sep="\t")
}
