 .libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4','/refdir/Rlib',.libPaths())))
 options(mc.cores=1,Ncpus=1)
 suppressPackageStartupMessages({library(data.table);library(limma);library(GSVA)})
 setDTthreads(1);if(requireNamespace('RhpcBLASctl',quietly=TRUE)){RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)}
 D<-'results/mp10_nectin_external_validation_v1';I<-file.path(D,'inputs');TB<-file.path(D,'tables');F0<-file.path(D,'figures');PD<-file.path(D,'plot_data')
 cohorts<-c('GSE32894','GSE31684','GSE13507')
 V<-'results/mp10_epithelial_treg_tcga_visualization_v2';H<-'results/mp10_hypothesis_followup_v2/tables';B<-'results/mp10_treg_metabolic_mechanisms_v1/tables'
 put<-function(z,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);fwrite(z,path,sep='\t',na='NA',nThread=1)}
 top<-scan('neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt',what='',quiet=TRUE)
 stopifnot(setequal(top,strsplit('HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1',' ')[[1]]))
 hall<-fread(file.path(H,'hallmark_gene_sets.tsv.gz'));ed<-fread(file.path(H,'scenic_corrected_regulon_edges.tsv'))
 immun<-jsonlite::fromJSON('highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1/models/signatures.json')
 sets<-list(MP10=top,mTORC1=hall[signature=='HALLMARK_MTORC1_SIGNALING',unique(gene)],SREBF1=ed[TF=='SREBF1',unique(target)],SREBF2=ed[TF=='SREBF2',unique(target)],Treg=immun$Treg_main,Tlineage=immun$T_lineage)
 familygenes<-c(paste0('NECTIN',1:4),'PVR');rna<-c(familygenes,'SREBF1','SREBF2','TIGIT','PVRIG','CD226','MTOR');nodes<-c('MP10','mTORC1','SREBF1','SREBF2','NECTIN2')
 aliases<-setNames(paste0('NECTIN',1:4),paste0('PVRL',1:4))
 rename<-function(x){x<-trimws(x);idx<-which(x%in%names(aliases));x[idx]<-unname(aliases[x[idx]]);x}
 core<-list(c('MP10','mTORC1'),c('MP10','SREBF1'),c('MP10','SREBF2'),c('MP10','NECTIN2'),c('mTORC1','SREBF1'),c('mTORC1','SREBF2'),c('mTORC1','NECTIN2'),c('SREBF1','NECTIN2'),c('SREBF2','NECTIN2'))
 score<-function(mat,ss){
   cov<-rbindlist(lapply(names(ss),function(n){g<-unique(ss[[n]]);u<-intersect(g,rownames(mat));data.table(variable=n,expected=length(g),used=length(u),members=paste(g,collapse=';'),retained=paste(u,collapse=';'),missing=paste(setdiff(g,u),collapse=';'),valid=length(u)>=5&&length(u)>=.5*length(g))}))
   out<-matrix(NA_real_,length(ss),ncol(mat),dimnames=list(names(ss),colnames(mat)))
   ok<-cov[valid==TRUE,variable];if(length(ok)){v<-gsva(ssgseaParam(mat,ss[ok],minSize=5,alpha=.25,normalize=FALSE),BPPARAM=BiocParallel::SerialParam(),verbose=FALSE);out[rownames(v),]<-v}
   list(mat=out,cov=cov)
 }
 libmatrix<-function(x)x[apply(x,1,function(v)all(is.finite(v)) && max(v)>min(v)),,drop=FALSE]
 prof<-function(x,pr,method='median'){
  stopifnot(nrow(x)==nrow(pr));ix<-which(!is.na(pr$gene)&nzchar(pr$gene));a<-as.data.table(x[ix,,drop=FALSE]);a[,gene:=pr$gene[ix]];r<-a[,lapply(.SD,median,na.rm=TRUE),by=gene];m<-as.matrix(r[,-1]);rownames(m)<-r$gene;m
 }
