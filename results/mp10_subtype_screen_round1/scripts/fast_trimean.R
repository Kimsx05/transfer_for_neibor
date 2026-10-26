# Only replace stats::aggregate's matrix/group iteration with exact compiled
# triMean aggregation. CellChat inference, permutation labels and parameters
# remain the installed original function. Used only after numerical comparison.
fastroot<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1'
Rcpp::sourceCpp(file.path(fastroot,'scripts/fast_trimean.cpp'),cacheDir=file.path(fastroot,'objects/rcpp_cache'),showOutput=FALSE)
computeCommunProb_fast<-CellChat::computeCommunProb
fastenv<-new.env(parent=environment(CellChat::computeCommunProb))
fastenv$aggregate<-function(x,by,FUN,...){
 if(!is.matrix(x)||length(by)!=1||!identical(FUN,CellChat::triMean))return(stats::aggregate(x,by,FUN,...))
 f<-droplevels(factor(by[[1]]));mat<-grouped_trimean_cpp(x,as.integer(f),nlevels(f));colnames(mat)<-colnames(x)
 cbind(data.frame(Group.1=factor(levels(f),levels=levels(f))),as.data.frame(mat,check.names=FALSE))
}
environment(computeCommunProb_fast)<-fastenv
