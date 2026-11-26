suppressPackageStartupMessages({library(sf);library(BiocParallel);library(data.table)})
setDTthreads(1)
# Use the actual source at the recorded commit when local R package installation is unavailable.
source("code/CRAWDAD/R/processing_inputs.R");source("code/CRAWDAD/R/finding_trends.R")
for(s in commandArgs(TRUE)){
 hd=s %in% c("Cervilla","D1_Johnson","UC_SCC");d=file.path(if(hd)"hd" else "xenium",s);t=fread(file.path(d,"crawdad_input.tsv.gz"));keys=if(hd)"MP10" else c("MP10_covered","Chol_biosynthesis_8","H05_transfer");xy=data.frame(x=t$x_um,y=t$y_um)
 # Native calibration field has all units within a complete 400 x 400 um square.
 # Choose the Treg closest to tissue median to ensure that the focal subset exists.
 rr=which(t$Treg);dist=(xy$x[rr]-median(xy$x))^2+(xy$y[rr]-median(xy$y))^2;center=rr[which.min(dist)];keep=abs(xy$x-xy$x[center])<=200 & abs(xy$y-xy$y[center])<=200
 z=t[keep];coords=xy[keep,];pat=toSF(coords,factor(ifelse(z[[keys[1]]]==1,"high",ifelse(z[[keys[1]]]==0,"low","other")),levels=c("high","low","other")),verbose=FALSE)
 pat$celltypes=droplevels(pat$celltypes)
 folds=makeShuffledCells(pat,scales=c(100,200,400,800),perms=3,ncores=1,verbose=FALSE)
 raw=findTrends(pat,neighDist=50,shuffleList=folds,ncores=1,subsetList=list(Treg_component=which(z$Treg)),verbose=FALSE,returnMeans=TRUE)
 uni=as.data.frame(raw[[1]]);if(!"high" %in% names(uni))uni$high=NA_real_;uni$scale=c(100,200,400,800);uni$n=length(which(keep));uni$focal=sum(z$Treg);uni$cx=xy$x[center];uni$cy=xy$y[center];fwrite(uni,file.path(d,"CRAWDAD_native_calibration.tsv"),sep="\t")
}
