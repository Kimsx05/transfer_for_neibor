source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R');suppressPackageStartupMessages(library(fgsea))
plan<-rd('model_plan.tsv');m20<-scan(file.path(P,'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'),what='character',quiet=TRUE);all<-list();cov<-list();audit<-list();dir.create(tb('gsea'),showWarnings=FALSE)
for(lib in c('Hallmark','Reactome')){
 gs<-fread(file.path(REF,if(lib=='Hallmark')'hallmark.tsv.gz'else'reactome.tsv.gz'));pathways<-split(gs$gene_symbol,gs$gs_name);ids<-if(lib=='Reactome')plan[primary==TRUE,model_id]else plan$model_id
 for(id in ids){f<-tb(paste0('de/',id,'.tsv.gz'));if(!file.exists(f))next;p<-plan[model_id==id];de<-fread(f);de[,signed_stat:=sign(logFC)*sqrt(pmax(F,0))];audit[[id]]<-de[,.(model_id=id,n_negative_F=sum(F<0),min_F=min(F),negative_F_all_P1=all(PValue[F<0]==1))]
 for(mask in if(p$kind=='Epi')c('Full','M20-masked')else'Full'){
 out<-tb(paste0('gsea/',id,'_',lib,'_',mask,'.tsv.gz'));signature<-digest(list(de[,.(gene,signed_stat)],pathways,mask,'fgsea1.36.2',123),algo='sha256')

 dd<-de[is.finite(signed_stat)];if(mask=='M20-masked')dd<-dd[!gene%in%m20];setorder(dd,-signed_stat,gene);stats<-setNames(dd$signed_stat,dd$gene);pw<-if(mask=='M20-masked')lapply(pathways,setdiff,y=m20)else pathways
 cc<-rbindlist(lapply(names(pw),function(n)data.table(pathway=n,n_original=length(unique(pathways[[n]])),n_M20_overlap=length(intersect(pathways[[n]],m20)),n_mapped=length(intersect(pw[[n]],names(stats))))));cc[,`:=`(model_id=id,library=lib,mask=mask)];cov[[paste(id,lib,mask)]]<-cc
 if(file.exists(out)&&file.exists(paste0(out,'.sha256'))&&readLines(paste0(out,'.sha256'))==signature){all[[paste(id,lib,mask)]]<-fread(out);next}
 set.seed(123);cat(id,lib,mask,'\n');rr<-as.data.table(fgseaMultilevel(pw,stats,minSize=15,maxSize=500,nPermSimple=1000,eps=1e-10,nproc=1));rr[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))];rr[,`:=`(model_id=id,library=lib,mask=mask,db_version=paste(unique(gs$db_version),collapse=';'))];fwrite(rr,out,sep='\t',compress='gzip');writeLines(signature,paste0(out,'.sha256'));all[[paste(id,lib,mask)]]<-rr
 }
 }
}
z<-rbindlist(all,fill=TRUE);z[,joint_main_q:=NA_real_];for(lib in unique(z$library))for(k in c('Epi','receiver')){ii<-which(z$library==lib&z$mask=='Full'&z$model_id %in% plan[primary&kind==k,model_id]);z[ii,joint_main_q:=p.adjust(pval,'BH')]};put(z,'gsea_results.tsv.gz');put(rbindlist(cov,fill=TRUE),'gsea_coverage.tsv.gz');put(rbindlist(audit,fill=TRUE),'DE_numerical_boundary_audit.tsv');cat('GSEA CURRENT RESULTS COMPLETE\n')
