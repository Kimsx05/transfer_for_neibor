.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(CellChat)})
source('mp10_subtype_screen_round1/scripts/fast_trimean.R');root<-fastroot
set.seed(10);x<-matrix(sample(c(0,0,0,1:20)/3,10000,TRUE),nrow=1000);gr<-factor(sample(letters[1:8],1000,TRUE));a<-stats::aggregate(x,list(gr),CellChat::triMean);b<-fastenv$aggregate(x,list(gr),CellChat::triMean);err1<-max(abs(as.matrix(a[,-1])-as.matrix(b[,-1])));stopifnot(err1<1e-14)
x<-qread(file.path(root,'objects/signaling_RNA_data.qs'));m<-fread(file.path(root,'tables/cells.tsv.gz'));m<-m[communication_type%in%c('SPP1 TAM','Malignant_like_epithelial')];m<-m[,head(.SD,250),by=communication_type];x<-x[,m$barcode]
cc<-createCellChat(x,meta=data.frame(labels=factor(m$communication_type),row.names=m$barcode),group.by='labels');cc@DB<-CellChatDB.human;cc<-subsetData(cc);cc<-identifyOverExpressedGenes(cc);cc<-identifyOverExpressedInteractions(cc);cc@LR$LRsig<-head(cc@LR$LRsig,20)
a<-computeCommunProb(cc,type='triMean',nboot=10,seed.use=123,population.size=FALSE);b<-computeCommunProb_fast(cc,type='triMean',nboot=10,seed.use=123,population.size=FALSE)
err2<-max(abs(a@net$prob-b@net$prob));err3<-max(abs(a@net$pval-b@net$pval));stopifnot(err2<1e-12,err3==0)
fwrite(data.table(check=c('exact_type7_aggregate_max_abs_error','CellChat_probability_max_abs_error','CellChat_permutation_P_max_abs_error'),value=c(err1,err2,err3)),file.path(root,'tables/triMean_acceleration_validation.tsv'),sep='\t')
