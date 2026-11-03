library(data.table);setDTthreads(2);library(fgsea)
o="mp10_paper_distance_subtypes_v2"
members=fread(file.path(o,"tables/reactome.tsv"));stopifnot(unique(members$db_version)=="2026.1.Hs")
sets=lapply(split(members$gene_symbol,members$gs_name),unique)
m=fread("mp10_paper_distance_full_v1/tables/input_manifest.tsv")
all=list()
for(ds in unique(m$dataset_id)){
 rows=list()
 for(s in m[dataset_id==ds & source_context=="tumor_study_specimen_not_region_diagnosis",section_id]){
  f=file.path(o,"tables",paste0(s,"_gene_ranks.tsv.gz"))
  if(!file.exists(f))next
  z=fread(f);z[,section:=s];z[,value:=frank(stat,ties.method="average")/.N-.5];z[,unit:=m[section_id==s,analysis_unit_uid]]
  z[is.na(unit),unit:=paste0("section_unpaired_",s)];rows[[s]]=z
 }
 if(!length(rows))next
 u=rbindlist(rows)[,.(value=mean(value)),by=.(gene,unit)][,.(value=mean(value)),by=gene]
 u[,dataset:=ds];all[[ds]]=u
}
z=rbindlist(all);fwrite(z,file.path(o,"tables/reactome_dataset_rank.tsv.gz"),sep="\t")
all[["dataset_equal"]]=z[,.(value=mean(value)),by=gene]
res=list()
for(ds in names(all)){
 ranks=setNames(all[[ds]]$value,all[[ds]]$gene)
 covered=lapply(sets,intersect,names(ranks));ratio=lengths(covered)/lengths(sets)
 test=covered[lengths(covered)>=15 & lengths(covered)<=500 & ratio>=.5]
 set.seed(20261103)
 a=fgseaSimple(test,stats=sort(ranks,decreasing=TRUE),nperm=10000,nproc=1,minSize=15,maxSize=500)
 a[,dataset:=ds];a[,leadingEdge:=vapply(leadingEdge,paste,collapse="|",FUN.VALUE="")]
 res[[ds]]=a
}
fwrite(rbindlist(res),file.path(o,"tables/reactome_enrichment.tsv"),sep="\t")
fwrite(data.table(pathway=names(sets),n_original=lengths(sets)),file.path(o,"tables/reactome_test_universe.tsv"),sep="\t")
