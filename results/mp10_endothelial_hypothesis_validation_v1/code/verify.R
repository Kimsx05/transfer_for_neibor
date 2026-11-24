source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
# Cross-check independent files, including the published model membership.
g=gettab("external_sample_map.tsv");stopifnot(nrow(g)==76,!anyDuplicated(g$gsm),all(g[,uniqueN(column)==.N,by=acc]$V1))
stopifnot(g[acc=="GSE207787"&grepl("^PDHA1",group),.N]==6,g[acc=="GSE207919",.N]==12)
put(g[,.(libraries=.N),by=.(acc,organism,group)],"external_group_counts.tsv")
orig=tab(V,"mTORC1_M0073_membership.tsv");me=gettab("project_sample_map.tsv")[type=="Capillary EC"&main_patient==TRUE&n_cells>=50]
stopifnot(setequal(orig$sample_uid,me$sample_uid),uniqueN(me$patient_uid_recorded)==12)
put(tab(B,"patient_relationship_audit.tsv"),"patient_relationship_audit.tsv")
st=gettab("coverage_model_status.tsv");stopifnot(st[type=="Total EC"&min_cells==50&exposure=="Q50",n]==14)
a=gettab("legacy_full_m20.tsv")
coverage=tab(B,"gsea_coverage.tsv.gz")[model_id=="M0073"&library=="Hallmark"&mask=="Full"]
missing=coverage[!pathway%in%a[mask=="Full",pathway]]
stopifnot(all(missing$n_mapped<15 | missing$n_mapped>500))
if(nrow(missing))a=rbind(a,missing[,.(pathway,mask,size=n_mapped)],fill=TRUE)
a[,status:=ifelse(is.finite(NES),"computed","NE")]
a[,reason:=ifelse(is.finite(NES),NA_character_,"outside original GSEA size range; no inferential estimate")]
put(a,"legacy_full_m20.tsv")
stopifnot(all(a[,.N,by=mask]$N==50))
g1=function(k)strsplit(a[mask=="Full"&pathway==paste0("HALLMARK_",k),leadingEdge],";")[[1]]
stopifnot(identical(intersect(g1("MYC_TARGETS_V1"),g1("TNFA_SIGNALING_VIA_NFKB")),"MYC"),length(intersect(g1("MTORC1_SIGNALING"),g1("UNFOLDED_PROTEIN_RESPONSE")))==10)
s=gettab("external_frozen_tails.tsv.gz");stopifnot(uniqueN(s$contrast)==24,all(s[,.N,by=.(contrast,arm)]$N==100))
stopifnot(nrow(gettab("perturbation_enrichment.tsv"))==24*2*2*2)
d=gettab("candidate_M0073.tsv");stopifnot(d[gene=="MKI67",UMI]==28,d[gene=="MKI67",status]=="nonzero_filtered")
tfs=c("MYC","ATF3","ATF4","RELA","KLF6","SREBF2","FOXO1");put(tab(B,"activity_associations.tsv.gz")[model_id=="M0073"&source%in%tfs],"legacy_TF_activity.tsv")
sc=tab(B,"scenic_Endothelial_regulon_coverage.tsv");seed=tab(B,"scenic_second_seed_coverage.tsv")[lineage=="Endothelial"]
gate=data.table(TF=tfs);gate[,motif_targets:=sc$n_targets[match(TF,sc$TF)]];gate[,seed456_evaluable:=TF%in%seed[status=="completed",TF]];gate[,independent_downstream_endpoint:=FALSE];gate[,virtual_KO:="NE; no validated patient-balanced network stability or independent endpoint"];put(gate,"virtual_KO_status.tsv")
hi=fread(file.path(B,"plot_data/Hallmark_sensitivity_Capillary_EC.tsv.gz"),sep="\t");put(hi,"legacy_display_sensitivity.tsv.gz")
v=tab(V,"model_results.tsv.gz");put(v[grepl("Capillary EC",paste(family,receiver,y),fixed=TRUE)&grepl("MTOR|mTOR",paste(x,y,pathway))],"legacy_patient_mTOR_audit.tsv.gz")
ex=gettab("external_matrix_audit.tsv");stopifnot(ex[acc=="GSE174785",min_assigned_fraction]>.99)
put(gettab("patient_effects.tsv.gz")[context=="Capillary EC_50"&mask=="Full"&setting=="Q50"],"primary_patient_statistics.tsv")
f=list.files(file.path(O,"figures"),pattern="\\.png$",full.names=TRUE);stopifnot(length(f)==9)
validate_img=function(f){stem=sub("[.]png$","",f);tabfile=file.path(O,"plot_data",paste0(basename(stem),".tsv"));x=png::readPNG(f)
 pdf=pdftools::pdf_info(paste0(stem,".pdf"));stopifnot(pdf$pages==1,file.exists(tabfile),max(x)-min(x)>.6);data.table(stem=basename(stem),pdf_pages=pdf$pages,width=dim(x)[2],height=dim(x)[1],data_rows=nrow(fread(tabfile)))}
put(rbindlist(lapply(f,validate_img)),"figure_verification.tsv")
files=list.files(C,recursive=TRUE,full.names=TRUE)
pdfs=files[grepl("PMC(10515159|5380221|8432967|6521885)/.*[.]pdf$",files)]
if(length(pdfs))put(rbindlist(lapply(pdfs,function(f){pg=pdftools::pdf_text(f);data.table(file=basename(f),pages=length(pg),characters=sum(nchar(pg)),captions=paste(head(grep("figure|table",trimws(unlist(strsplit(pg,"\n"))),ignore.case=TRUE,value=TRUE),10),collapse=" | "))})),"supplement_pdf_readout.tsv")
packs=c("qs","data.table","edgeR","limma","fgsea","babelgene","org.Hs.eg.db","readxl","ggplot2","ragg","ROCR","pdftools")
put(data.table(package=packs,version=vapply(packs,function(n)as.character(packageVersion(n)),"")),"software_versions.tsv")
writeLines(capture.output(sessionInfo()),file.path(O,"logs/final_session.txt"))
cat("PASS: 76 libraries, 24 contrasts, M0073 membership, independent profile counts, marker filtering, full sets, 9 figure pairs.\n")

