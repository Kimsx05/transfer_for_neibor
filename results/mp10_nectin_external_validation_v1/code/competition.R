source('results/mp10_nectin_external_validation_v1/code/common.R')
td<-file.path(TB,'competition');p<-fread(file.path(B,'profiles.tsv'))[main_patient==TRUE&n_cells>=50]
og<-fread(file.path(H,'sample_level_scores.tsv.gz'))[profile_type=='Treg'&family=='Treg'&mask=='Full'&variant=='standard'&signature=='HALLMARK_OXIDATIVE_PHOSPHORYLATION'&main_patient==TRUE&n_cells>=50,.(sample_uid,y=score)]
p<-merge(p,og,by='sample_uid');stopifnot(!anyDuplicated(p$patient_uid_recorded))
m<-readRDS(file.path(B,'treg_logCPM.rds'))[,p$profile_id]
d<-fread(file.path(B,'perturbation_DE.tsv.gz'))[contrast=='TCR'&adj.P.Val<.05&!is.na(human)];pct<-function(x)(rank(x)-.5)/length(x)
off<-fread(file.path(H,'official_gene_sets.tsv.gz'));ox<-hall[signature=='HALLMARK_OXIDATIVE_PHOSPHORYLATION',gene]
ctrl<-list(TCR=d$human,G2M=hall[signature=='HALLMARK_G2M_CHECKPOINT',gene],CHOL_SYNTH=off[signature=='CHOL_SYNTH',gene]);mapping<-list()
for(k in names(ctrl)){
 g<-intersect(setdiff(ctrl[[k]],ox),rownames(m));if(k=='TCR'){w<-pct(d$oriented_t[match(g,d$human)])-.5;v<-as.numeric(crossprod(w,apply(m[g,,drop=FALSE],2,pct)))} else v<-colMeans(m[g,,drop=FALSE])
 p[,(k):=as.numeric(scale(v))];mapping[[k]]<-data.table(control=k,used=paste(g,collapse=';'),n=length(g),removed_overlap=paste(intersect(ctrl[[k]],ox),collapse=';'))
}
diag<-fread(file.path(V,'plot_data/diagnostic_specimen_profiles.tsv'))[mask=='Full'&profile_type=='H05'&eligible_50==TRUE&main_patient==TRUE&gene=='NECTIN2',.(sample_uid,H05_NECTIN2=TMM_logCPM)]
p<-merge(p,diag,by='sample_uid',all.x=TRUE);put(p,file.path(td,'paired_specimens.tsv'));put(rbindlist(mapping),file.path(td,'control_members.tsv'));out<-list()
for(ex in c('Q50_z','H05_NECTIN2'))for(pop in c('primary','replicated_datasets')){
 u<-copy(p[is.finite(get(ex))]);if(pop=='replicated_datasets')u<-u[,if(.N>1).SD,by=dataset];u[,exposure:=as.numeric(scale(get(ex)))];for(k in c('unadjusted','TCR','G2M','CHOL_SYNTH')){
 f<-lm(reformulate(c('dataset','exposure',if(k!='unadjusted')k),'y'),data=u);maxh<-max(hatvalues(f));df<-df.residual(f);r<-data.table(exposure=ex,population=pop,control=k,n=nrow(u),patients=uniqueN(u$patient_uid_recorded),datasets=uniqueN(u$dataset),rank=f$rank,df=df,max_leverage=maxh,beta=unname(coef(f)['exposure']),lo=NA_real_,hi=NA_real_,p=NA_real_,status='NE: unit leverage or rank/df')
 if(df>0&&maxh<1-1e-8&&f$rank==length(coef(f))){se<-sqrt(sandwich::vcovHC(f,type='HC3')['exposure','exposure']);r[,c('lo','hi','p','status'):=list(beta-qt(.975,df)*se,beta+qt(.975,df)*se,2*pt(-abs(beta/se),df),'estimated')]}
 out[[length(out)+1]]<-r
 }}
res<-rbindlist(out);res[,q:=p.adjust(p,'BH',n=.N),by=.(exposure,population)];put(res,file.path(td,'models.tsv'))

