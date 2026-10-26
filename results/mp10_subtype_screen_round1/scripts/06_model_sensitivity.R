# One deterministic specimen per recorded patient checks whether the sole repeated
# tumor-patient ID drives MP10 results, without changing the primary sample table.
.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
library(data.table);library(limma)
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1'
source(file.path(root,'scripts/speckle_getTransformedProps_official.R'))
s<-fread(file.path(root,'tables/samples.tsv'));s<-s[tolower(MP10_eligible)=='true'];setorder(s,patient_uid_recorded,sample_uid);s<-s[!duplicated(patient_uid_recorded)]
m<-fread(file.path(root,'tables/cells.tsv.gz'));out<-list()
fam<-list(list(level='major',den='TME',parent='all',cells=m[major!='Epithelial'],type='major'))
for(pa in unique(m[major!='Epithelial',major]))fam[[length(fam)+1]]<-list(level='subtype',den='parent',parent=pa,cells=m[major==pa],type='subtype')
for(f in fam){
 tr<-getTransformedProps(f$cells[[f$type]],f$cells$sample_uid,'logit');d<-s[sample_uid%in%colnames(tr$Proportions)];y<-tr$TransformedProps[,d$sample_uid,drop=FALSE];pr<-tr$Proportions[,d$sample_uid,drop=FALSE];ok<-apply(pr,1,function(z)diff(range(z))>1e-12)
 if(!any(ok))next
 X<-model.matrix(~factor(dataset)+MP10_sample,d);fit<-eBayes(lmFit(y[ok,,drop=FALSE],X),robust=sum(ok)>=3,trend=FALSE);tt<-as.data.table(topTable(fit,coef='MP10_sample',number=Inf,sort.by='none',confint=TRUE),keep.rownames='subtype');tt[,`:=`(level=f$level,denominator_type=f$den,n_specimens=nrow(d),n_patient_ids_recorded=uniqueN(d$patient_uid_recorded),method='lexicographically_first_eligible_specimen_per_recorded_patient')];out[[length(out)+1]]<-tt
}
r<-rbindlist(out);r[,q:=p.adjust(P.Value,'BH'),by=.(level,denominator_type)];fwrite(r,file.path(root,'tables/MP10_one_specimen_patient_sensitivity.tsv'),sep='\t')
