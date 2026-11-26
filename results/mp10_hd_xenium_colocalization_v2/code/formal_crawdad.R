library(sf);library(data.table);setDTthreads(1)
source("code/CRAWDAD/R/processing_inputs.R");source("code/CRAWDAD/R/finding_trends.R")
dyn.load("code/crawdad_counts.so")
kind=BiocParallel::bplapply(1,function(i)RNGkind(),BPPARAM=BiocParallel::SnowParam(workers=1))[[1]]
do.call(RNGkind,as.list(unname(kind)));writeLines(kind,"audit/CRAWDAD_rngkind.txt")
for(s in commandArgs(TRUE)){
 d=file.path(if(s %in% c("Cervilla","D1_Johnson","UC_SCC"))"hd" else "xenium",s);z=fread(file.path(d,"crawdad_input.tsv.gz"))
 keys=intersect(c("MP10","MP10_covered","Chol_biosynthesis_8","H05_transfer"),names(z))
 evaluate=function(t,B,pilot=FALSE){
  xy=data.frame(x=as.double(t$x_um),y=as.double(t$y_um));a=as.matrix(t[,..keys])==1
  if(pilot){a=cbind(high=a[,1],low=t[[keys[1]]]==0,other=t[[keys[1]]]<0)}
  usekeys=colnames(a);storage.mode(a)="integer";sf=toSF(xy,factor(a[,1]),verbose=FALSE)
  # The named subset is not an exclusive cell-type label. The upstream deduplicated sf intersection includes focal coordinates.
  # Keep this native zero-distance behavior separate from primary no-self statistics.
  include=lengths(st_intersects(sf,st_buffer(sf[t$Treg,],dist=50)))>0
  n=sum(include);obs=colSums(a[include,,drop=FALSE]);rows=list()
  for(r in c(100,200,400,800)){
   sim=.Call("tile_counts",xy$x,xy$y,a,include,as.double(r),as.integer(B))
   for(j in seq_along(usekeys)){nu=sim[,j];p0=(nu+obs[j])/(2*n);zz=(obs[j]-nu)/n/sqrt(p0*(1-p0)*2/n);valid=is.finite(zz)
    rows[[length(rows)+1]]=data.table(sample=s,state=usekeys[j],neighDist=50,tile_um=r,perms=B,buffer_nodes=n,focal=sum(t$Treg),observed_high=obs[j],expected_high=mean(nu),mean_proportion_Z=if(all(valid))mean(zz) else NA_real_,finite_perm_Z=sum(valid),mean_Z_conditional=if(any(valid)) mean(zz[valid]) else NA_real_,P=if(sd(nu)>0)(1+sum(abs(nu-mean(nu))>=abs(obs[j]-mean(nu))-1e-10))/(B+1) else NA_real_)
   }
  }
  rbindlist(rows)
 }
 calibration=fread(file.path(d,"CRAWDAD_native_calibration.tsv"));cx=calibration$cx[1];cy=calibration$cy[1];local=z[abs(x_um-cx)<=200 & abs(y_um-cy)<=200];check=evaluate(local,3,TRUE)
 cols=intersect(c("high","low","other"),names(calibration));nc=melt(calibration,id.vars="scale",measure.vars=cols,variable.name="state",value.name="native_Z")
 test=merge(check,nc,by.x=c("state","tile_um"),by.y=c("state","scale"));test[,absolute_error:=abs(mean_proportion_Z-native_Z)];test[,NA_agrees:=is.na(mean_proportion_Z)==is.na(native_Z)];fwrite(test,file.path(d,"CRAWDAD_native_check.tsv"),sep="\t")
 full=evaluate(z,999);full[,native_calibration_max_error:=if(any(is.finite(test$absolute_error))) max(test$absolute_error,na.rm=TRUE) else NA_real_];fwrite(full,file.path(d,"CRAWDAD_formal.tsv"),sep="\t")
}
