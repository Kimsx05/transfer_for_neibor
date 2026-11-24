# Additional descriptive overlap check, outside the registered primary families.
source(file.path(Sys.getenv('B_OUT'),'code/common.R'))
d<-readt(out,'perturbation_DE.tsv.gz')[!is.na(human)]
p<-readt(out,'profiles.tsv')[main_patient==TRUE&n_cells>=50]
mat<-readRDS(file.path(out,'tables/treg_logCPM.rds'))[,p$profile_id]
pct<-function(x)(rank(x)-.5)/length(x)
block<-unique(c(top,'SCAP','FASN','SREBF1','SREBF2','CD36'))
g2m<-readt(v2,'hallmark_gene_sets.tsv.gz')[signature=='HALLMARK_G2M_CHECKPOINT',unique(gene)]
tcr<-d[contrast=='TCR'&adj.P.Val<.05]
outcome<-list()
for(control in c('TCR','G2M')){
 genes<-intersect(setdiff(rownames(mat),block),if(control=='TCR')tcr$human else g2m)
 # Rank within disjoint gene backgrounds as well as disjoint program members.
 if(control=='TCR'){
  sign0<-tcr$oriented_t[match(genes,tcr$human)]
  z<-as.numeric(scale(as.numeric(crossprod(pct(sign0)-.5,apply(mat[genes,,drop=FALSE],2,pct)))))
 }else{
  # A fixed expression-level readout avoids borrowing any genes from the outcome.
  z<-as.numeric(scale(colMeans(mat[genes,,drop=FALSE])))
 }
 bg<-setdiff(rownames(mat),c(block,genes));r<-apply(mat[bg,,drop=FALSE],2,pct);rownames(r)<-bg
 for(k in c('SCAP','FASN','CD36','SREBF2')){
   u<-d[contrast==k&human%in%bg];w<-pct(u$oriented_t)-.5
   y<-as.numeric(scale(as.numeric(crossprod(w,r[u$human,,drop=FALSE]))))
   f<-lm(y~dataset+z+Q50_z,p);cf<-summary(f)$coef['Q50_z',];ci<-confint(f)['Q50_z',]
   outcome[[length(outcome)+1]]<-data.table(control=control,program=k,control_genes=length(genes),program_genes=nrow(u),overlap=length(intersect(genes,u$human)),patients_recorded=nrow(p),beta=cf[1],lo=ci[1],hi=ci[2],p=cf[4])
 }
}
o<-rbindlist(outcome);o[,q:=p.adjust(p,'BH'),by=control];saveT(o,'disjoint_control_sensitivity.tsv')
