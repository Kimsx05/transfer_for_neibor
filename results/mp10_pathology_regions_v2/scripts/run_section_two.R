args=commandArgs(TRUE);sec=args[1];cores=if(length(args)>1)as.integer(args[2]) else 2L
cores=min(cores,6L)
ROOT=normalizePath('.');OUT=file.path(ROOT,'mp10_pathology_regions_v2');.libPaths(unique(c(file.path(OUT,'Rlib'),'/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4','/refdir/Rlib',.libPaths())))
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1);RhpcBLASctl::blas_set_num_threads(1);RhpcBLASctl::omp_set_num_threads(1)
suppressPackageStartupMessages({library(Matrix);library(data.table);library(SpaCET);library(Banksy);library(SingleCellExperiment)})
setDTthreads(1);set.seed(20261027)
man=fread(file.path(OUT,'tables/section_manifest.tsv'));r=man[section_id==sec];stopifnot(nrow(r)==1);dest=r$output_dir;meta=fread(file.path(dest,'input_spots.tsv'));markers=list(epithelial=c('EPCAM','KRT8','KRT18','KRT19','KRT7'),urothelial=c('UPK1A','UPK1B','UPK2','UPK3A','UPK3B','KRT20'),basal=c('KRT5','KRT14','TP63','KRT17'),smooth_muscle=c('MYH11','ACTG2','CNN1','TAGLN','DES','MYLK'),fibroblast=c('COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA'),immune=c('PTPRC','CD3D','CD3E','CD79A','LYZ','LST1'),endothelial=c('PECAM1','VWF','CDH5','KDR'),proliferation=c('MKI67','TOP2A','UBE2C','CENPF','TYMS'),pericyte=c('RGS5','CSPG4','PDGFRB','MCAM'))
cat(format(Sys.time()),'reading original counts',sec,'\n');con=if(grepl('gz$',r$counts_path))gzfile(r$counts_path) else file(r$counts_path);x=as(readMM(con),'CsparseMatrix');rownames(x)=readLines(r$genes_path);colnames(x)=readLines(r$barcodes_path);stopifnot(setequal(colnames(x),meta$barcode));x=x[,meta$barcode];stopifnot(all(x@x>=0),all(x@x==round(x@x)),!anyDuplicated(meta$SpaCET_id))
if(anyDuplicated(rownames(x)))stop('Duplicate gene symbols require explicit aggregation')
ng=colSums(x>0);umi=colSums(x);keep=ng>=100&umi>0;meta[,algorithm_qc:=ifelse(keep,'retained','less_than_100_detected_genes')];meta[,nCount:=umi];meta[,nFeature:=ng];fwrite(meta,file.path(dest,'algorithm_qc.tsv'),sep='\t');xm=x[,keep];mm=meta[keep];colnames(xm)=mm$SpaCET_id
# Explicitly keep both directions of software IDs; original barcodes never shortened.
fwrite(mm[,.(dataset_id,section_id,barcode,SpaCET_id,array_row,array_col,x,y)],file.path(dest,'barcode_crosswalk.tsv'),sep='\t')
if(!file.exists(file.path(dest,'spacet.done'))){
 tryCatch({
 co=data.frame(pixel_row=mm$y,pixel_col=mm$x,array_row=mm$array_row,array_col=mm$array_col,row.names=mm$SpaCET_id);md=data.frame(barcode=mm$barcode,row.names=mm$SpaCET_id)
 obj=create.SpaCET.object(xm,co,metaData=md,imagePath=r$image_path,platform='Visium',organism='human');obj=SpaCET.quality.control(obj,min.genes=100)
 obj=SpaCET.deconvolution(obj,cancerType='BLCA',signatureType=NULL,adjacentNormal=FALSE,coreNo=cores)
 obj=SpaCET.identify.interface(obj,MalignantCutoff=0.5)
 pr=obj@results$deconvolution$propMat;major=unique(c('Malignant',names(obj@results$deconvolution$Ref$lineageTree),'Unidentifiable'));major=intersect(major,rownames(pr));sig=obj@results$deconvolution$malRes$sig;lab=major[max.col(t(pr[major,,drop=FALSE]),ties.method='first')];res=data.frame(barcode=mm$barcode,malignant_fraction=as.numeric(pr['Malignant',mm$SpaCET_id]),SpaCET_signature_source=paste(sig,collapse=':'),SpaCET_major_label=lab,SpaCET_native_interface=as.character(obj@results$CCI$interface[1,mm$SpaCET_id]));fwrite(res,file.path(dest,'spacet_spots.tsv'),sep='\t');fwrite(data.frame(barcode=mm$barcode,t(pr[major,mm$SpaCET_id,drop=FALSE]),check.names=FALSE),file.path(dest,'spacet_major_fractions.tsv'),sep='\t');saveRDS(list(malRes=obj@results$deconvolution$malRes,lineageTree=obj@results$deconvolution$Ref$lineageTree),file.path(dest,'spacet_signature_results.rds'),compress=TRUE);writeLines(as.character(sig),file.path(dest,'spacet_signature_source.txt'));writeLines('complete',file.path(dest,'spacet.done'));rm(obj);gc()
 },error=function(e){writeLines(conditionMessage(e),file.path(dest,'spacet.error'));cat('SpaCET failed:',conditionMessage(e),'\n')})
}
# Marker expression is annotation evidence, never independent validation.
norm=xm%*%Diagonal(x=1e4/Matrix::colSums(xm));rownames(norm)=rownames(xm);colnames(norm)=colnames(xm);lognorm=norm;lognorm@x=log1p(lognorm@x);mk=data.frame(barcode=mm$barcode)
for(nm in names(markers)){gs=intersect(markers[[nm]],rownames(lognorm));mk[[paste0('marker_',nm)]]=if(length(gs))Matrix::colMeans(lognorm[gs,,drop=FALSE]) else NA_real_;mk[[paste0('detected_',nm)]]=if(length(gs))Matrix::colSums(xm[gs,,drop=FALSE]>0) else 0L}
for(gn in unique(unlist(markers)))mk[[paste0('gene_',gn)]]=if(gn%in%rownames(xm))as.numeric(lognorm[gn,]) else NA_real_
fwrite(mk,file.path(dest,'marker_evidence.tsv'),sep='\t');fwrite(data.frame(marker_set=rep(names(markers),lengths(markers)),gene=unlist(markers),available=unlist(markers)%in%rownames(xm)),file.path(dest,'marker_coverage.tsv'),sep='\t')
if(!file.exists(file.path(dest,'banksy.done'))){
 tryCatch({
 # Select 3000 HVGs from the full transcriptome; no metabolic/MP10 gene restriction.
 varfit=scran::modelGeneVar(lognorm);hvg=scran::getTopHVGs(varfit,n=min(3000,nrow(varfit)));hvg=hvg[is.finite(varfit[hvg,'bio'])];writeLines(hvg,file.path(dest,'banksy_HVGs.txt'))
 normbanksy=xm[hvg,,drop=FALSE]%*%Diagonal(x=mean(Matrix::colSums(xm))/Matrix::colSums(xm));colnames(normbanksy)=mm$SpaCET_id;se=SingleCellExperiment(assays=list(normcounts=normbanksy),colData=data.frame(gx=mm$array_col/2,gy=mm$array_row*sqrt(3)/2,row.names=mm$SpaCET_id))
 ed=fread(file.path(dest,'banksy_two_ring_edges.tsv'));map=rep(NA_integer_,nrow(meta));map[which(keep)]=seq_len(sum(keep));ed[,from:=map[from]];ed[,to:=map[to]];ed=ed[!is.na(from)&!is.na(to)];loc=cbind(mm$array_col/2,mm$array_row*sqrt(3)/2);ed[,distance:=sqrt(rowSums((loc[from,,drop=FALSE]-loc[to,,drop=FALSE])^2))];ed[,weight:=exp(-distance^2/median(distance)^2),by=from];ed[,weight:=weight/sum(weight),by=from];ed[,phi:=atan2(loc[to,2]-loc[from,2],loc[to,1]-loc[from,1])];ed[phi<0,phi:=phi+2*pi];ed=ed[,.(from,to,weight,phi)]
 # Project extension: exact <=2-hop retained-tissue graph replaces unrestricted kNN.
 # Official BANKSY harmonics, scaling, PCA and Leiden routines remain unchanged.
 orig=getFromNamespace('computeNeighbors','Banksy');assignInNamespace('computeNeighbors',function(...)copy(ed),ns='Banksy')
 se=computeBanksy(se,assay_name='normcounts',coord_names=c('gx','gy'),compute_agf=TRUE,k_geom=18,parallel=FALSE)
 assignInNamespace('computeNeighbors',orig,ns='Banksy')
 se=runBanksyPCA(se,use_agf=TRUE,lambda=0.2,npcs=30,seed=20261027,lazy=FALSE)
 res=if(r$pilot)c(0.4,0.6) else 0.6
 se=clusterBanksy(se,use_agf=TRUE,lambda=0.2,npcs=30,k_neighbors=30,resolution=res,seed=20261027)
 labs=as.data.frame(colData(se));cols=grep('^clust',names(labs),value=TRUE);fwrite(data.frame(barcode=mm$barcode,labs[,cols,drop=FALSE]),file.path(dest,'banksy_domains.tsv'),sep='\t');pcs=reducedDim(se,1);saveRDS(pcs,file.path(dest,'banksy_PCA.rds'),compress=TRUE);fwrite(data.frame(barcode=mm$barcode,spatial_neighbors=tabulate(ed$from,nbins=nrow(mm))),file.path(dest,'banksy_neighbor_coverage.tsv'),sep='\t');writeLines('complete',file.path(dest,'banksy.done'));rm(se);gc()
 },error=function(e){writeLines(conditionMessage(e),file.path(dest,'banksy.error'));cat('BANKSY failed:',conditionMessage(e),'\n')})
}
writeLines(capture.output(sessionInfo()),file.path(dest,'software_session.txt'));cat(format(Sys.time()),'section algorithms finished',sec,'\n')
