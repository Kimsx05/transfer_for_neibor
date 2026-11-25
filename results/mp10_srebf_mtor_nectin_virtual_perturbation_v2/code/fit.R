suppressPackageStartupMessages({library(Matrix);library(data.table);library(scTenifoldKnk);library(scTenifoldNet)})
setDTthreads(1);options(mc.cores=1);RhpcBLASctl::blas_set_num_threads(1)
cfg=jsonlite::fromJSON("config/main.json");K=cfg$source;arg=commandArgs(TRUE);id=arg[1];stage=arg[2]
d=file.path("cache/jobs",id);dir.create(d,recursive=TRUE,showWarnings=FALSE)
put=function(x,n)fwrite(x,file.path(d,n),sep="\t",na="")
P=fread("tables/panel_members.tsv");meta=fread("cache/metadata.tsv.gz");m=copy(meta)
seed=if(startsWith(id,"seed"))as.integer(sub("seed","",id)) else cfg$main_seed
if(startsWith(id,"sub")){set.seed(12000L+as.integer(sub("sub","",id)));u=sort(unique(m$patient_id));m=m[patient_id%in%sample(u,ceiling(length(u)*cfg$patient_fraction))]}
if(startsWith(id,"drop_"))m=m[dataset!=sub("drop_","",id)]
if(startsWith(id,"cluster_"))m=m[state==sub("cluster_","",id)]
put(m[,.(barcode,patient_id,specimen_id,dataset,state)],"membership.tsv.gz")
if(startsWith(id,"cluster_")&&(nrow(m)<cfg$cluster_min_cells||uniqueN(m$patient_id)<cfg$cluster_min_patients)){
 put(data.table(status="NE",reason="cluster <500 cells or <4 patients"),"status.tsv");quit(save="no")
}
ng=nrow(P);z=matrix(readBin("cache/panel.bin","integer",size=4,endian="little",n=ng*nrow(meta)),nrow=ng)
take=which(if(id=="full")P$P_full else P$P_clean)
ct=z[take,match(m$barcode,meta$barcode),drop=FALSE];rm(z);rownames(ct)=P$gene[take];colnames(ct)=m$barcode
dets=rowSums(ct>0);np=Reduce("+",lapply(split(seq_len(ncol(ct)),m$patient_id),function(ix)as.integer(rowSums(ct[,ix,drop=FALSE])>0)))
mu=rowMeans(ct);vr=rowMeans(ct^2)-mu^2
cov=data.table(gene=rownames(ct),cells=dets,patients=np,mean=mu,det=dets/nrow(m),variance=vr,ok=dets>=cfg$raw_gene_min_cells & np>=cfg$raw_gene_min_patients & vr>0)
ct=as(ct[cov$ok,,drop=FALSE],"dgCMatrix");g=rownames(ct)
# Complement is exclusively a denominator carrier. It never enters PCR, correlation, heat, or ranks.
remainder=as.numeric(m$nCount_RNA-colSums(ct));stopifnot(all(remainder>=0))
raw=rbind(ct,Matrix(matrix(remainder,nrow=1,dimnames=list("__excluded_UMI__",m$barcode)),sparse=TRUE))
stopifnot(all(colSums(raw)==m$nCount_RNA))
if(stage=="reference" || !id%in%c("full","clean")){
 set.seed(seed);N=scTenifoldNet::cpmNormalization(raw)[g,,drop=FALSE]
 subs=scTenifoldNet::makeNetworks(N,nNet=cfg$nNet,nCells=cfg$nCells,nComp=cfg$nComp,q=cfg$q,nCores=1)
 T=scTenifoldNet::tensorDecomposition(subs,K=cfg$td_K,maxError=cfg$td_maxError,maxIter=cfg$td_maxIter,nDecimal=cfg$td_nDecimal,seed=seed)$X
 W=as.matrix(scTenifoldKnk:::strictDirection(T,lambda=cfg$lambda));diag(W)=0;W=t(W);rm(subs,T,N);invisible(gc())
 if(id%in%c("full","clean"))saveRDS(as(W,"dgCMatrix"),file.path(d,"WT.rds"))
} else W=as.matrix(readRDS(file.path(d,"WT.rds")))
stopifnot(identical(rownames(W),g),identical(colnames(W),g))
cov[ok==TRUE,outdegree:=rowSums(W!=0)];cov[ok==TRUE,indegree:=colSums(W!=0)];put(cov,"coverage.tsv.gz")
if(stage=="reference")quit(save="no")
controls=fread("tables/control_matching.tsv");legacy=fread("../mp10_srebf_nectin_virtual_perturbation_v1/tables/matched_controls.tsv")$gene
lmap=fread("tables/gene_mapping.tsv",select=c("feature_id","gene"))
allKO=unique(c(K,controls$control,lmap$gene[match(legacy,lmap$feature_id)]));allKO=allKO[!is.na(allKO)]
ks=intersect(allKO,cov[ok==TRUE & outdegree>0,gene]);ko=as.list(ks);names(ko)=ks
put(data.table(ko=allKO,status=ifelse(allKO%in%ks,"EVALUABLE","NE"),reason=ifelse(allKO%in%ks,"",ifelse(allKO%in%g,"no effective outgoing edge","local expression/patient/variance gate or absent panel"))),"KO_status.tsv")
G=igraph::graph_from_adjacency_matrix(W!=0,mode="directed",diag=FALSE)
paths=lapply(ks,function(k)igraph::as_ids(igraph::subcomponent(G,k,mode="out")));names(paths)=ks;rm(G)
S=scTenifoldKnk:::.logCPM(raw,g)
stopifnot(identical(as.numeric(S$libSize),as.numeric(m$nCount_RNA)))
direction=function(depth=FALSE,pc=FALSE){
 L=if(depth)scTenifoldKnk:::.regressLibSize(S) else S$L
 L=L-rowMeans(L);den=sqrt(rowSums(L^2));den[den==0]=Inf
 R=tcrossprod(L/den)
 if(pc){ev=RSpectra::eigs_sym(R,1,which="LA");R=R-ev$values[1]*tcrossprod(ev$vectors[,1]);sig=sqrt(pmax(diag(R),0));R=R/outer(sig,sig);R[!is.finite(R)]=0}
 scTenifoldKnk:::.kernelEffect(heatKernel(R,t=cfg$directionT),S,ko)
}
dirs=direction();depth=direction(TRUE);pc=direction(TRUE,TRUE)
if(id=="clean"){
 # Validate numerical identity with official APIs on a test block using the complete denominator.
 small=g[unique(c(match(ks,g),seq_len(min(80,length(g)))))]
 aa=scTenifoldKnk::knockoutDirection(raw,ko,genes=small,t=cfg$directionT)
 SS=scTenifoldKnk:::.logCPM(raw,small);L=SS$L-SS$mu;bb=scTenifoldKnk:::.kernelEffect(heatKernel(tcrossprod(L/sqrt(rowSums(L^2))),t=cfg$directionT),SS,ko)
 stopifnot(max(abs(aa-bb))<1e-8);put(data.table(max_abs_diff=max(abs(aa-bb))),"direction_QC.tsv")
}
H=heatKernel(W,t=cfg$heatT);distances=hkManifoldAlignment(W,raw,ko,t=cfg$heatT,H=H)
formatDR=function(dt,k,method){
 dt=as.data.table(dt);dt[,is_KO:=gene==k];dt[,rank:=NA_real_];dt[is_KO==FALSE,rank:=frank(-distance,ties.method="average")];dt[,rank_denominator:=length(g)-1L]
 dt[,rank_fraction:=(rank-1)/max(1L,length(g)-2L)];dt[,reachable:=gene%in%paths[[k]]];dt[,depth_score:=depth[k,gene]];dt[,pc1_score:=pc[k,gene]];dt[is_KO==TRUE,c("depth_score","pc1_score"):=.(NA_real_,NA_real_)]
 dt[,ko:=k];dt[,analysis:=id];dt[,branch:=method];dt[,status:="EVALUABLE"];dt[,coverage_cells:=cov$cells[match(gene,cov$gene)]];dt[,coverage_patients:=cov$patients[match(gene,cov$gene)]]
 dt
}
effects=rbindlist(lapply(ks,function(k){ds=distances[k,];formatDR(scTenifoldKnk:::.drStatistics(ds,g,g==k,noiseLevel=sqrt(.Machine$double.eps)*max(ds),direction=dirs[k,]),k,"heat")}))
put(effects,"heat.tsv.gz")
if(id=="clean"){
 manifoldK=intersect(unique(c(K,controls[match_order==1,control])),ks)
 ma=rbindlist(lapply(manifoldK,function(k){Y=W;Y[k,]=0;set.seed(seed);A=scTenifoldNet::manifoldAlignment(W,Y,d=cfg$manifold_dim,nCores=1);formatDR(scTenifoldKnk::dRegulation(A,gKO=k,direction=dirs[k,]),k,"manifold")}))
 put(ma,"manifold.tsv.gz")
 focus=intersect(g,unique(c(cfg$endpoints,K,fread("tables/program_members.tsv")$gene)))
 rows=list()
 for(p in unique(m$patient_id))for(s in c("All","H05","H10","H50","L50")){
  j=which(m$patient_id==p & if(s=="All")rep(TRUE,nrow(m))else m[[s]])
  n=length(j);x=raw[,j,drop=FALSE]
  if(n>=10){sx=scTenifoldKnk:::.logCPM(x,g);vk=intersect(K,ks);zz=hkManifoldAlignment(W,x,as.list(vk),H=H)}
  for(k in K){
   valid=rep(FALSE,length(focus));v=rep(NA_real_,length(focus))
   if(n>=10 && k%in%ks){valid=is.finite(sx$sd[focus]) & is.finite(sx$sd[k]) & focus!=k;v[valid]=zz[k,focus[valid]]}
   rows[[length(rows)+1L]]=data.table(patient_id=p,stratum=s,ko=k,gene=focus,n_cells=n,distance=v,status=ifelse(valid,"EXPLORATORY","NE"),reason=ifelse(valid,"",if(n<10)"<10 cells" else "constant/self/unavailable source or target"))
  }
 }
 put(rbindlist(rows),"patient.tsv.gz")
}
put(data.table(status="complete",reason="",cells=nrow(m),patients=uniqueN(m$patient_id),genes=length(g),kos=length(ks)),"status.tsv")

