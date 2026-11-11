# Same Seurat 5.5.1 GroupSingletons rule; cache membership indices instead of rescanning all IDs per cluster.
# Sequential membership updates, original column order, numeric mean(), and seed-1 ties are retained.
fast_group_singletons <- function(ids,SNN,group.singletons=TRUE,verbose=TRUE){
 singles<-names(which(table(ids)==1));singles<-intersect(unique(ids),singles)
 if(!group.singletons){ids[ids%in%singles]<-'singleton';return(ids)}
 if(!length(singles))return(ids)
 cls<-setdiff(as.character(unique(ids)),singles)
 membership<-setNames(lapply(cls,function(j)which(ids==j)),cls)
 single_index<-setNames(match(singles,ids),singles)
 for(i in singles){ii<-single_index[[as.character(i)]];row<-as.numeric(SNN[ii,,drop=TRUE]);connectivity<-setNames(vapply(membership,function(j)mean(row[j]),numeric(1)),cls);set.seed(1);m<-max(connectivity,na.rm=TRUE);mi<-which(connectivity==m,arr.ind=TRUE);closest<-sample(names(connectivity[mi]),1);ids[ii]<-closest;membership[[closest]]<-sort(c(membership[[closest]],ii))}
 if(verbose)message(length(singles),' singletons identified. ',length(unique(ids)),' final clusters (cached identical rule).')
 ids
}
