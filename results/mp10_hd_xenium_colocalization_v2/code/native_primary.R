library(spatstat.geom);library(spatstat.random);library(data.table);setDTthreads(1)
for(s in commandArgs(TRUE)){
 d=file.path(if(s %in% c("Cervilla","D1_Johnson","UC_SCC"))"hd" else "xenium",s)
 t=fread(file.path(d,"spatstat_input.tsv.gz"));keys=intersect(names(t),c("MP10","MP10_covered","Chol_biosynthesis_8","H05_transfer"))
 win=owin(range(t$x_um)+c(-1,1),range(t$y_um)+c(-1,1));Y=ppp(t$x_um,t$y_um,window=win,marks=seq_len(nrow(t)),check=FALSE);g=factor(t$stratum)
 H=as.matrix(t[,..keys]);y=as.matrix(t[,.(density,Tfrac)]);tot=colSums(y,na.rm=TRUE);n=colSums(is.finite(y));finite=is.finite(y);y[!finite]=0;storage.mode(finite)="double"
 obsstat=function(idx){h=H[idx,,drop=FALSE];sm=crossprod(h,y);nh=crossprod(h,finite);lo=matrix(tot,nrow=ncol(h),ncol=2,byrow=TRUE)-sm;nl=matrix(n,nrow=ncol(h),ncol=2,byrow=TRUE)-nh
  as.vector(sm/nh-lo/nl)
 }
 set.seed(261126);ob=obsstat(seq_len(nrow(t)));null=matrix(NA_real_,4999,length(ob));last=0
 while(last<4999){b=min(32,4999-last);sim=rlabel(Y,group=g,nsim=b,drop=FALSE)
  for(j in seq_len(b)){ix=marks(sim[[j]]);if(j==1)stopifnot(identical(g,g[ix]));null[last+j,]=obsstat(ix)}
  last=last+b
 }
 center=colMeans(null);sd=apply(null,2,sd);p=(1+colSums(abs(sweep(null,2,center))>=matrix(abs(ob-center)-1e-12,nrow=4999,ncol=length(ob),byrow=TRUE)))/5000;p[sd<1e-10]=NA_real_
 fwrite(data.frame(sample=s,state=rep(keys,2),endpoint=rep(c("density","Treg/T"),each=length(keys)),observed=ob,null_mean=center,null_sd=sd,P=p,perms=4999),file.path(d,"spatstat_primary.tsv"),sep="\t")
}
