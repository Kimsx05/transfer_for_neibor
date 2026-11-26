library(sf);library(BiocParallel);library(data.table);setDTthreads(1);source("code/CRAWDAD/R/processing_inputs.R");source("code/CRAWDAD/R/finding_trends.R");dyn.load("code/crawdad_counts.so")
d="xenium/A1R2";z=fread(file.path(d,"crawdad_input.tsv.gz"));cc=fread(file.path(d,"CRAWDAD_native_calibration.tsv"));t=z[abs(x_um-cc$cx[1])<=200&abs(y_um-cc$cy[1])<=200];xy=data.frame(x=t$x_um,y=t$y_um);lab=factor(ifelse(t$MP10_covered==1,"high",ifelse(t$MP10_covered==0,"low","other")));a=toSF(xy,lab,verbose=FALSE);sh=makeShuffledCells(a,scales=c(100,200),perms=3,ncores=1,verbose=FALSE)
buf=st_buffer(a[t$Treg,],dist=50);q=st_intersection(a,buf$geometry);inc=lengths(st_intersects(a,buf))>0
print(c(nrow(q),length(intersect(rownames(q),rownames(a))),sum(inc)));print(bplapply(1,function(i)RNGkind(),BPPARAM=SnowParam(1)))
for(r in c(100,200)){
 print(sapply(sh[[as.character(r)]],function(v)sum(v[inc]=="high")))
 for(kind in c("Mersenne-Twister","L'Ecuyer-CMRG")){
 RNGkind(kind);mark=matrix(as.integer(lab=="high"));print(c(kind,.Call("tile_counts",xy$x,xy$y,mark,inc,as.double(r),3L)))
 }
}
