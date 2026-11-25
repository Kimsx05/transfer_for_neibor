suppressPackageStartupMessages({library(data.table);library(Matrix);library(scTenifoldKnk);library(scTenifoldNet)})
setDTthreads(1);RhpcBLASctl::blas_set_num_threads(as.integer(Sys.getenv("OPENBLAS_NUM_THREADS","3")))
options(mc.cores=1);args=commandArgs(TRUE);id=args[1];seed=if(id=="alt")456L else 123L
dir.create("cache/jobs",showWarnings=FALSE);dir.create("qc/jobs",showWarnings=FALSE)
out=function(name)paste0("qc/jobs/",id,"_",name,".tsv.gz")
put=function(x,name)fwrite(x,out(name),sep="\t",na="")
timing=function(label,expr) {
 t=proc.time();disk=function()sum(file.info(list.files(c("cache/jobs","qc/jobs"),full.names=TRUE))$size,na.rm=TRUE);before=disk()
 on.exit({z=proc.time()-t; rss=as.numeric(sub(".*?([0-9]+) kB.*","\\1",grep("^VmHWM",readLines("/proc/self/status"),value=TRUE)))
 write.table(data.frame(stage=label,wall_s=z[3],cpu_s=sum(z[1:2]),process_highwater_RSS_KiB=rss,output_byte_change=disk()-before),paste0("logs/time_",id,".tsv"),append=file.exists(paste0("logs/time_",id,".tsv")),col.names=!file.exists(paste0("logs/time_",id,".tsv")),row.names=FALSE,quote=FALSE,sep="\t")})
 force(expr)
}
genes=fread("tables/model_panel.tsv")$gene; m=fread("cache/metadata.tsv.gz"); flag=if(id=="pilot")"pilot" else if(id=="alt")"alt" else "main"
take=which(m[[flag]])
if(startsWith(id,"sub")) {
 set.seed(12000L+as.integer(sub("sub","",id)));ps=sort(unique(m$patient_id[take]));take=take[m$patient_id[take]%in%sample(ps,ceiling(.8*length(ps)))]
}
if(startsWith(id,"drop_"))take=take[m$dataset[take]!=sub("drop_","",id)]
if(startsWith(id,"state_"))take=take[m$state[take]==sub("state_","",id)]
m0=m[take];put(unique(m0[,.(dataset,patient_id,state,barcode)]),"membership")
K=c("SREBF1","SREBF2"); nPat=uniqueN(m0$patient_id)
if(nrow(m0)<80 || nPat<4){put(data.table(reason="state cells <80 or patients <4",n_cells=nrow(m0),n_patients=nPat),"status");quit(save="no")}
ct=matrix(readBin("cache/cells_by_genes_int32.bin","integer",n=as.double(nrow(m))*length(genes),size=4,endian="little"),nrow=length(genes))
ct=ct[,take,drop=FALSE];rownames(ct)=genes;colnames(ct)=m0$barcode
orient=fread("qc/orientation.tsv.gz")[take];stopifnot(identical(colnames(ct),orient$barcode),all(colSums(ct)==orient$panel_UMI),all(ct>=0))
pos=rowSums(ct>0);pp=Reduce("+",lapply(split(seq_len(ncol(ct)),m0$patient_id),function(j)as.integer(rowSums(ct[,j,drop=FALSE])>0)))
rawmean=rowMeans(ct);variance=rowMeans(ct^2)-rawmean^2;good=pos>=10 & pp>=3 & variance>0
cov=data.table(gene=genes,n_positive=pos,n_positive_patients=pp,raw_mean=rawmean,raw_variance=variance,detected_fraction=pos/ncol(ct),in_network=good)
ct=as(ct[good,,drop=FALSE],"dgCMatrix");genes=rownames(ct);stopifnot(all(colSums(ct)>0))
set.seed(seed)
wt=timing("reference_fit",{
 norm=scTenifoldNet::cpmNormalization(ct)
 nets=scTenifoldNet::makeNetworks(norm,nNet=10,nCells=min(ncol(ct),200L),nComp=3,q=.9,nCores=1)
 tn=scTenifoldNet::tensorDecomposition(nets,K=3,maxError=1e-5,maxIter=1000,nDecimal=3,seed=seed)$X
 w=as.matrix(scTenifoldKnk:::strictDirection(tn,lambda=0));diag(w)=0;t(w)
})
rm(tn,norm,nets,w);invisible(gc());stopifnot(identical(rownames(wt),genes),identical(colnames(wt),genes),all(is.finite(wt)))
cov[good,outdegree:=rowSums(wt!=0)];cov[good,indegree:=colSums(wt!=0)];put(cov,"coverage")
if(id=="main")saveRDS(as(wt,"dgCMatrix"),"cache/jobs/main_WT.rds")
ko=intersect(K,genes)
if(id=="main"){
 tf=scan("../../resources/mp10_regulatory_priors/human_tfs.txt",what="",quiet=TRUE)
 reserved=unique(c(fread("tables/program_members.tsv")$gene,"NECTIN1","NECTIN2","NECTIN3","NECTIN4","PVR",K))
 cand=cov[in_network & outdegree>0 & gene%in%tf & (!gene%in%reserved | gene%in%K)]
 v=as.matrix(cand[,.(mean=log1p(raw_mean),detected_fraction,degree=log1p(outdegree))]);v=scale(v);v[,!is.finite(attr(v,"scaled:scale"))|attr(v,"scaled:scale")==0]=0
 matches=rbindlist(lapply(intersect(K,cand$gene),function(k){
 ix=match(k,cand$gene);b=copy(cand);b[,match_distance:=sqrt(rowSums(sweep(v,2,v[ix,],"-")^2))];b=b[!gene%in%K][order(match_distance,gene)][1:min(5,.N)];b[,matched_to:=k];b
 }))
 put(matches,"controls");ko=unique(c(ko,matches$gene))
}
graph=igraph::graph_from_adjacency_matrix(wt!=0,mode="directed",diag=FALSE)
reach=lapply(ko,function(k)igraph::as_ids(igraph::subcomponent(graph,k,mode="out")));names(reach)=ko
H=timing("WT_heat_kernel",scTenifoldKnk::heatKernel(wt,t=10))
source("code/direction.R")
score_direction=function(depth) {
 if(startsWith(id,"sub")) gram_direction(ct,genes,ko,depth) else scTenifoldKnk::knockoutDirection(ct,as.list(ko),genes,t=5,regressLibSize=depth)
}
dirs=timing("direction_scores",score_direction(FALSE))
adj=timing("depth_sensitivity",score_direction(TRUE))
# Compare directions after removing PC1 from depth-residualized logCPM expression.
pcs=timing("PC1_sensitivity",{
 S=scTenifoldKnk:::.logCPM(ct,genes);L=scTenifoldKnk:::.regressLibSize(S)
 L=t(scale(t(L)));L[!is.finite(L)]=0;R=tcrossprod(L)/max(1,ncol(ct)-1)
 ev=RSpectra::eigs_sym(R,1,which="LA");R=R-ev$values[1]*tcrossprod(ev$vectors[,1])
 sg=sqrt(pmax(diag(R),0));R=R/outer(sg,sg);R[!is.finite(R)]=0
 HK=scTenifoldKnk::heatKernel(R,t=5);ks=as.list(ko);names(ks)=ko;scTenifoldKnk:::.kernelEffect(HK,S,ks)
})
rm(S,L,R,HK);gc()
all=rbindlist(lapply(ko,function(k)timing(paste0("KO_",k),{
 ds=scTenifoldKnk::hkManifoldAlignment(wt,ct,list(k),t=10,H=H)[1,]
 dr=as.data.table(scTenifoldKnk:::.drStatistics(ds,names(ds),names(ds)==k,noiseLevel=sqrt(.Machine$double.eps)*max(ds),direction=dirs[k,]))
 dr[,depth_direction_score:=adj[k,gene]];dr[,pc1_direction_score:=pcs[k,gene]];dr[,reachable:=gene%in%reach[[k]]];dr[,perturbed_gene:=k];dr[,rank:=NA_real_];dr[gene!=k,rank:=rank(-distance,ties.method="average")];dr
})))
timing("export",put(if(id%in%c("main","pilot")||startsWith(id,"state_"))all else all[gene=="NECTIN2"],"effects"))
if(id=="main"){
 # Patient-conditioned magnitudes on the fixed WT kernel, without reusing a pooled direction.
 records=list()
 for(p in unique(m0$patient_id))for(stratum in c("All","H05","H10","H50","L50")){
  ix=which(m0$patient_id==p & (if(stratum=="All")rep(TRUE,nrow(m0)) else m0[[stratum]]))
  if(length(ix)<10 || !"NECTIN2"%in%genes)next
  x=ct[,ix,drop=FALSE];z=scTenifoldKnk::hkManifoldAlignment(wt,x,as.list(intersect(K,ko)),H=H)
  for(k in intersect(K,ko)){
   valid="NECTIN2"%in%reach[[k]] && all(apply(as.matrix(x[c(k,"NECTIN2"),,drop=FALSE]),1,var)>0)
   records[[length(records)+1]]=data.table(patient_id=p,dataset=m0$dataset[ix[1]],stratum=stratum,n_cells=length(ix),perturbed_gene=k,effect_value=if(valid)z[k,"NECTIN2"]else NA_real_,reason=if(valid)"" else "constant gene or no directed path")
  }
 }
 put(rbindlist(records),"patient")
}
put(data.table(reason="fitted",n_cells=nrow(m0),n_patients=nPat,n_genes=length(genes)),"status")

