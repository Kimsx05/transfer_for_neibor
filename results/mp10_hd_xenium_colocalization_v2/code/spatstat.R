suppressPackageStartupMessages({library(spatstat.geom);library(spatstat.random);library(spatstat.explore);library(data.table)})
setDTthreads(1)
for(s in commandArgs(TRUE)){
 hd=s %in% c("Cervilla","UC_SCC","D1_Johnson");d=file.path(if(hd)"hd" else "xenium",s);t=fread(file.path(d,"spatstat_input.tsv.gz"));keys=if(hd)"MP10" else c("MP10_covered","Chol_biosynthesis_8","H05_transfer")
 win=owin(range(t$x_um)+c(-1,1),range(t$y_um)+c(-1,1));Y=ppp(t$x_um,t$y_um,window=win,marks=seq_len(nrow(t)),check=FALSE)
 g=factor(t$stratum);set.seed(261126);p=rlabel(Y,group=g);stopifnot(all(g==g[marks(p)]),all(sort(marks(p))==seq_len(nrow(t))))
 stats=function(X,r,...){
  idx=marks(X);h=t[[key]][idx]==1;dist=t$distance
  a=sapply(r,function(radius)mean(dist[t$cf200 & h]<=radius)-mean(dist[t$cf200 & !h]<=radius))
  fv(data.frame(r=r,obs=a),argu="r",valu="obs",ylab=quote(Delta*F(r)),yexp=quote(Delta*F(r)),fmla=.~r,alim=range(r),labl=c("r","obs"),desc=c("radius","high minus low"))
 }
 for(key in keys){
  set.seed(261126);rr=c(25,50,100,150,200);env=envelope(Y,fun=stats,r=rr,simulate=expression(rlabel(Y,group=g)),nsim=199,nrank=10,global=TRUE,verbose=FALSE)
  fwrite(as.data.frame(env),file.path(d,paste0("spatstat_envelope_",key,".tsv")),sep="\t")
 }
}
