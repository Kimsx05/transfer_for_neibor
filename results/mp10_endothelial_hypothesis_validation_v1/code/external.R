source("results/mp10_endothelial_hypothesis_validation_v1/code/maps.R")
parse_soft=function(acc){
 f=file.path(C,acc,paste0(acc,"_family.soft.gz"));x=readLines(gzfile(f),warn=FALSE)
 p=grep("^\\^SAMPLE = ",x);rbindlist(lapply(seq_along(p),function(i){
  y=x[p[i]:(if(i<length(p))p[i+1]-1 else length(x))]
  field=function(k)paste(sub("^.*? = ","",grep(paste0("^!Sample_",k," = "),y,value=TRUE)),collapse=" | ")
  data.table(acc,gsm=sub("^.*= ","",y[1]),title=field("title"),description=field("description"),organism=field("organism_ch1"),attributes=field("characteristics_ch1"),source=field("source_name_ch1"),processing=field("data_processing"),donor_block="unreported; unpaired",supp=basename(field("supplementary_file_1")))
 }))
}
allmeta=list();allde=list();audits=list()
for(acc in c("GSE174785","GSE201466","GSE207787","GSE207919","GSE121770")){
 dir=file.path(C,acc);meta=parse_soft(acc);cat(acc,"\n")
 if(acc=="GSE174785"){
  tmp=tempfile();dir.create(tmp);files=untar(Sys.glob(file.path(dir,"*tar")),list=TRUE)
  stopifnot(all(files==basename(files)),setequal(files,meta$supp));untar(Sys.glob(file.path(dir,"*tar")),exdir=tmp)
  meta[,column:=supp];l=lapply(meta$column,function(n)fread(file.path(tmp,n)));ids=l[[1]]$target_id
  stopifnot(all(vapply(l,function(x)identical(ids,x$target_id),logical(1))),mean(startsWith(ids,"ENSMUST"))>.9)
  ann=fread(file.path(C,"GENCODE_mouse_M25/wgEncodeGencodeAttrsVM25.txt.gz"),header=FALSE,select=c(1,5))
  tx=unique(ann[,.(ENSEMBLTRANS=sub("[.].*","",V5),id=sub("[.].*","",V1))])
  tx=tx[,if(.N==1).(id=id[1]),by=ENSEMBLTRANS]
  ii=match(ids,tx$ENSEMBLTRANS);v=do.call(cbind,lapply(l,"[[","tpm"));x=rowsum(v[!is.na(ii),],tx$id[ii[!is.na(ii)]])
  audits[[acc]]=data.table(acc,type="kallisto TPM, summed unique GENCODE M25 transcript-to-gene",raw_rows=length(ids),mapped_rows=sum(!is.na(ii)),min_assigned_fraction=min(colSums(x)/colSums(v)),max_assigned_fraction=max(colSums(x)/colSums(v)))
  meta[,group:=paste0(ifelse(grepl("Oxidative",title),"Ox","Gly"),ifelse(grepl("WT",title),"WT","KO"))]
  ctr=c(ATF4_Ox_WTminusKO="OxWT-OxKO",ATF4_Gly_WTminusKO="GlyWT-GlyKO",ATF4_context_interaction="(OxWT-OxKO)-(GlyWT-GlyKO)",WT_OxminusGly="OxWT-GlyWT",KO_OxminusGly="OxKO-GlyKO")
  rm(l,v);unlink(tmp,recursive=TRUE)
 }else{
  files=Sys.glob(file.path(dir,"*txt.gz"));if(!length(files)){allmeta[[acc]]=meta;next}
  d=fread(files,sep="\t");nr=nrow(d)
  if(acc=="GSE121770"){
   meta[,column:=title];cols=meta$column;x=as.matrix(d[,..cols])
   rownames(x)=d[[1]];stopifnot(all(startsWith(rownames(x),"ENSMUSG")),all(x==round(x)),min(x)>=0)
   meta[,group:=sub("_.*","",title)];ctr=c(QPC_KOminusWT="KO-WT")
  }else{
   id=if(acc=="GSE207919")sub(".*\\[([^]]+)\\].*","\\1",meta$title) else meta$description
   cols=names(d);meta[,column:=vapply(id,function(z){k=cols[cols==z|startsWith(cols,paste0(z,"_"))];stopifnot(length(k)==1);k},character(1))]
   target=sub(".*si([A-Za-z0-9]+) .*","\\1",meta$title);target[target%in%c("Control","CTRL")]="CTRL"
   time=as.integer(sub(".*?([0-9]+)hr.*","\\1",meta$title));meta[,group:=paste0(target,time)]
   stopifnot(all(table(meta$group)==3),!anyDuplicated(meta$column))
   symbol_col=intersect(c("Gene Symbol","Gene.Symbol","gene_name","Gene_Symbol","Gene"),names(d))[1]
   sym=human(d[[symbol_col]]);cols=meta$column; xx=as.matrix(d[,..cols]);stopifnot(all(is.finite(xx)),min(xx)>=0)
   # Repeated multi-locus entries can have identical profiles. Average within mapped gene.
   x=rowsum(xx[!is.na(sym),],sym[!is.na(sym)])/as.numeric(table(factor(sym,levels=sort(unique(na.omit(sym))))));rm(xx)
   ctr=character();for(t in sort(unique(time[time>0]))){ctr[paste0("TNF_CTRL_",t)]=paste0("CTRL",t,"-CTRL0")}
   for(k in setdiff(unique(target),"CTRL")){
    for(t in sort(unique(time)))ctr[paste0(k,"_KD_",t)]=paste0(k,t,"-CTRL",t)
    for(t in sort(unique(time[time>0])))ctr[paste0(k,"_interaction_",t)]=paste0("(",k,t,"-",k,"0)-(CTRL",t,"-CTRL0)")
   }
  }
  audits[[acc]]=data.table(acc,type=if(acc=="GSE121770")"HTSeq counts" else "CPM",raw_rows=nr,mapped_rows=nrow(x))
 }
 colnames(x)=meta$column;stopifnot(!anyDuplicated(meta$column),!anyDuplicated(rownames(x)),ncol(x)==nrow(meta))
 design=model.matrix(~0+factor(meta$group));colnames(design)=levels(factor(meta$group));con=makeContrasts(contrasts=unname(ctr),levels=design);colnames(con)=names(ctr)
 if(acc=="GSE121770"){
  y=DGEList(x);keep=filterByExpr(y,design=design);y=calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);y=estimateDisp(y,design,robust=TRUE);fit=glmQLFit(y,design,robust=TRUE)
 }else{
  keep=rowSums(x>=1)>=min(table(meta$group));fit=eBayes(contrasts.fit(lmFit(log2(x[keep,,drop=FALSE]+0.5),design),con),trend=TRUE,robust=TRUE)
 }
 audits[[acc]][,tested:=sum(keep)];allmeta[[acc]]=meta
 for(j in seq_along(ctr)){
  if(acc=="GSE121770"){
   d=as.data.table(topTags(glmQLFTest(fit,contrast=con[,j]),n=Inf)$table,keep.rownames="id");d[,stat:=sign(logFC)*sqrt(pmax(F,0))];d[,gene:=m$gene[match(id,m$ensembl)]]
  }else{
   d=as.data.table(topTable(fit,coef=j,n=Inf,sort.by="none",confint=TRUE),keep.rownames="id");setnames(d,c("t","P.Value","adj.P.Val"),c("stat","PValue","FDR"))
   d[,gene:=if(acc=="GSE174785")m$gene[match(id,m$ensembl)] else id]
  }
  d[,contrast:=paste(acc,names(ctr)[j],sep="__")];d[,formula:=ctr[[j]]];d[,orientation:=ifelse(grepl("_KD_|_interaction_[0-9]",names(ctr)[j]),-1,1)]
  allde[[d$contrast[1]]]=d
 }
}
de=rbindlist(allde,fill=TRUE);put(de,"external_gene_statistics.tsv.gz");put(rbindlist(allmeta,fill=TRUE),"external_sample_map.tsv");put(rbindlist(audits,fill=TRUE),"external_matrix_audit.tsv")
de=de[!is.na(gene)];de[,oriented:=stat*orientation];setorder(de,contrast,-oriented,gene)
tails=de[,rbind(head(.SD[oriented>0],100)[,arm:="UP"],tail(.SD[oriented<0],100)[,arm:="DN"]),by=contrast]
put(tails,"external_frozen_tails.tsv.gz")

