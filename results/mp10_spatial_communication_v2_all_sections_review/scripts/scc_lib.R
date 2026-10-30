# Shared aggregation and IO; no mutation of v1 inputs or position matrices.
scc_atomic_lines <- function(lines,path) {tmp<-paste0(path,'.tmp.',Sys.getpid());writeLines(lines,tmp);stopifnot(file.rename(tmp,path))}
scc_put <- function(x,name,out) {tmp<-file.path(out,paste0(name,'.tsv.tmp.',Sys.getpid()));data.table::fwrite(x,tmp,sep='\t',na='NA');stopifnot(file.rename(tmp,file.path(out,paste0(name,'.tsv'))))}
scc_support <- function(base,inp,section,Qlist) {
 p<-file.path(inp,'type_support.tsv');if(!file.exists(p))p<-file.path(base,'input_review','type_support.tsv')
 if(!file.exists(p))stop('Required type_support.tsv is missing; Q>0 cannot substitute eligibility')
 z<-data.table::fread(p);if('section_id'%in%names(z))z<-z[section_id==section]
 if(!'level'%in%names(z)&&'network_level'%in%names(z))data.table::setnames(z,'network_level','level')
 stopifnot(all(c('level','node','type_evaluable')%in%names(z)),!anyDuplicated(z[,.(level,node)]))
 for(lev in names(Qlist))stopifnot(all(colnames(Qlist[[lev]])%in%z[level==lev,node]))
 z[,type_evaluable:=as.logical(type_evaluable)];z
}
scc_load_Q <- function(inp,spotid) {
 z<-lapply(c('major','fine'),function(lev){x<-data.table::fread(file.path(inp,paste0('Q_',lev,'.tsv')));stopifnot(identical(x$spot_id,spotid));q<-as.matrix(x[,-1]);storage.mode(q)<-'double';rownames(q)<-spotid;stopifnot(all(is.finite(q)),all(q>=0),max(abs(rowSums(q)-1))<1e-8);q});names(z)<-c('major','fine');z
}
scc_aggregate <- function(Plist,spotid,Qlist,support,out,section) {
 started<-Sys.time();lrstat<-vector('list',length(Plist));nativeflow<-crossflow<-diagtotal<-0
 if(length(Plist))for(i in seq_along(Plist)){
  p<-Matrix::drop0(Plist[[i]]);stopifnot(identical(dim(p),c(length(spotid),length(spotid))),all(is.finite(p@x)),all(p@x>=0));d<-diag(p);diag(p)<-0;p<-Matrix::drop0(p)
  lrstat[[i]]<-data.table::data.table(interaction=names(Plist)[i],cross_spot_links=length(p@x),cross_spot_flow=sum(p),same_spot_links=sum(d>0),same_spot_flow=sum(d),source_positions=sum(Matrix::rowSums(p)>0),target_positions=sum(Matrix::colSums(p)>0))
  crossflow<-crossflow+sum(p);diagtotal<-diagtotal+sum(d)
 }
 lrstat<-if(length(lrstat))data.table::rbindlist(lrstat) else data.table::data.table(interaction=character(),cross_spot_links=integer(),cross_spot_flow=numeric(),same_spot_links=integer(),same_spot_flow=numeric(),source_positions=integer(),target_positions=integer())
 scc_put(lrstat,'position_LR_support',out);checks<-list()
 for(lev in names(Qlist)) {
  Q<-Qlist[[lev]];H<-1*(Q>0);N<-colnames(Q);n<-length(N);ss<-support[level==lev][match(N,node)];stopifnot(identical(ss$node,N))
  good<-ss$type_evaluable & colSums(Q)>0;edgegood<-rep(good,n)&rep(good,each=n)
  total<-flow<-links<-diagflow<-matrix(0,n,n,dimnames=list(N,N));rows<-vector('list',length(Plist));diags<-vector('list',length(Plist));maxerr<-0
  if(length(Plist))for(i in seq_along(Plist)) {
   pp<-Matrix::drop0(Plist[[i]]);D<-Matrix::Diagonal(x=diag(pp));diag(pp)<-0;pp<-Matrix::drop0(pp)
   F<-as.matrix(crossprod(Q,pp%*%Q));K<-as.matrix(crossprod(H,(pp>0)%*%H));B<-SpatialCellChat::computeAvgCommunProb_LR_Visium_Avg(pp,Q,H)
   expected<-F/K;expected[is.nan(expected)]<-0;maxerr<-max(maxerr,max(abs(expected-B)));DF<-as.matrix(crossprod(Q,D%*%Q))
   total<-total+B;flow<-flow+F;links<-links+K;diagflow<-diagflow+DF
   b<-as.vector(B);f<-as.vector(F);k<-as.vector(K);df<-as.vector(DF);b[!edgegood]<-NA_real_;f[!edgegood]<-NA_real_;k[!edgegood]<-NA_real_;df[!edgegood]<-NA_real_
   rows[[i]]<-data.table::data.table(source=rep(N,n),target=rep(N,each=n),interaction=names(Plist)[i],mean_strength=b,flow=f,supported_position_pair_denominator=k,evaluable=edgegood)
   diags[[i]]<-data.table::data.table(source=rep(N,n),target=rep(N,each=n),interaction=names(Plist)[i],same_spot_flow=df,evaluable=edgegood)
  }
  checks[[lev]]<-data.table::data.table(level=lev,aggregation_helper_max_abs_error=maxerr,cross_flow_mass_error=sum(flow)-crossflow,same_spot_mass_error=sum(diagflow)-diagtotal,status='PASS')
  stopifnot(maxerr<1e-10,abs(sum(flow)-crossflow)<1e-7*max(1,crossflow),abs(sum(diagflow)-diagtotal)<1e-7*max(1,diagtotal))
  for(nm in c('total','flow','links','diagflow')){m<-get(nm);m[!good,]<-NA_real_;m[,!good]<-NA_real_;assign(nm,m)}
  reason<-ifelse(good,'computable','NE_insufficient_abundance_support')
  edges<-data.table::data.table(source=rep(N,n),target=rep(N,each=n),total_mean_strength=as.vector(total),total_flow=as.vector(flow),summed_LR_supported_pair_denominators=as.vector(links),same_spot_flow_excluded=as.vector(diagflow),source_effective_abundance=rep(colSums(Q),n),target_effective_abundance=rep(colSums(Q),each=n),source_support_spots=rep(colSums(H),n),target_support_spots=rep(colSums(H),each=n),source_eligible=rep(good,n),target_eligible=rep(good,each=n),source_eligibility=rep(reason,n),target_eligibility=rep(reason,each=n),evaluable=edgegood,n_LR=length(Plist),group_P_test='not performed; not used')
  for(col in intersect(c('n_abundance_ge_0_01','total_expected_abundance','identifiability_status'),names(ss))){edges[[paste0('source_',col)]]<-rep(ss[[col]],n);edges[[paste0('target_',col)]]<-rep(ss[[col]],each=n)}
  scc_put(edges,paste0(lev,'_edges'),out)
  empty<-data.table::data.table(source=character(),target=character(),interaction=character(),mean_strength=numeric(),flow=numeric(),supported_position_pair_denominator=numeric(),evaluable=logical())
  data.table::fwrite(if(length(rows))data.table::rbindlist(rows) else empty,file.path(out,paste0(lev,'_per_LR_aggregation.tsv.gz')),sep='\t',na='NA')
  data.table::fwrite(if(length(diags))data.table::rbindlist(diags) else data.table::data.table(source=character(),target=character(),interaction=character(),same_spot_flow=numeric(),evaluable=logical()),file.path(out,paste0(lev,'_same_spot_per_LR.tsv.gz')),sep='\t',na='NA')
  scc_put(data.table::data.table(node=N,total,evaluable=good),paste0(lev,'_total_mean_strength_matrix'),out);scc_put(data.table::data.table(node=N,flow),paste0(lev,'_total_flow_matrix'),out)
 }
 if(length(Plist))Ptotal<-Reduce('+',Plist) else Ptotal<-Matrix::sparseMatrix(i=integer(),j=integer(),dims=c(length(spotid),length(spotid)))
 sp<-data.table::data.table(spot_id=spotid,outgoing_including_same_spot=Matrix::rowSums(Ptotal),incoming_including_same_spot=Matrix::colSums(Ptotal),same_spot=diag(Ptotal));sp[,`:=`(outgoing_cross_spot=outgoing_including_same_spot-same_spot,incoming_cross_spot=incoming_including_same_spot-same_spot)]
 scc_put(sp,'spot_sender_receiver',out);scc_put(data.table::rbindlist(checks),'aggregation_checks',out)
 list(aggregation_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),n_LR=length(Plist),cross_spot_flow=crossflow,same_spot_flow=diagtotal)
}
