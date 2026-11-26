library(spatstat.geom);library(spatstat.explore);library(data.table);setDTthreads(1)
for(s in commandArgs(TRUE)){
 d=file.path(if(s %in% c("Cervilla","D1_Johnson","UC_SCC"))"hd" else "xenium",s);xy=fread(file.path(d,"k_points.tsv.gz"))
 dn=scan(file.path(d,"k_shape.txt"),quiet=TRUE);m=matrix(as.integer(readBin(file.path(d,"k_mask.bin"),"raw",n=prod(dn)))==1,nrow=dn[1]);window=owin(xrange=c(-4,8*dn[2]-4),yrange=c(-4,8*dn[1]-4),mask=m);tt=xy[Treg==TRUE]
 result=list()
 for(k in intersect(names(xy),c("MP10","MP10_covered","Chol_biosynthesis_8","H05_transfer"))){
  ep=xy[get(k)==1 & Treg==FALSE];if(nrow(tt)<30 || nrow(ep)<20){
   result[[k]]=data.frame(state=k,r=NA,Kcross_inhom=NA,n_Treg=nrow(tt),n_high=nrow(ep),reason="sparse for intensity fitting");next
  }
  a=rbind(data.frame(x=tt$x,y=tt$y,type="Treg"),data.frame(x=ep$x,y=ep$y,type="high"))
  X=ppp(a$x,a$y,window=window,marks=factor(a$type),check=TRUE)
  rr=seq(0,200,10)
  K=Kcross.inhom(X,i="Treg",j="high",r=rr,correction="border",sigma=200)
  result[[k]]=data.frame(state=k,r=rr,Kcross_inhom=K$border,n_Treg=nrow(tt),n_high=nrow(ep),reason="native Euclidean raster window; exploratory")
 }
 fwrite(rbindlist(result),file.path(d,"spatstat_inhomogeneous_K.tsv"),sep="\t")
}
