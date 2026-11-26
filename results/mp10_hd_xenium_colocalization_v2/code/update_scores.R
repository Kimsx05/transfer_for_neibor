.libPaths(unique(c(file.path(path.expand("~"),"R/x86_64-pc-linux-gnu-library/4.5"),.libPaths())))
s=commandArgs(TRUE)[1];v=normalizePath(file.path(Sys.getenv("MP_V2_BASE"),"../mp10_hd_xenium_colocalization_v1"));Sys.setenv(OUT=v);source(file.path(v,"code/transfer.R"))
el=fread(file.path("review",s,"eligibility.tsv"));old=fread(file.path(v,"xenium",s,"state_transfer.tsv.gz"));added=setdiff(el[eligible_override==TRUE,id],old$id)
if(length(added)){
 h5=Sys.glob(file.path(root,"highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5",paste0("*_",s,"_*h5")))
 x=Read10X_h5(h5,use.names=TRUE)[["Gene Expression"]][readLines(file.path(v,"shared/common_panel.txt")),unique(c(old$id,added))]
 ri=which(meta$state_ref & !is.na(meta$H05_truth));w=transfer(ref[,ri],x,meta$H05_truth[ri]);sets=jsonlite::read_json(file.path(v,"shared/programs.json"),simplifyVector=TRUE)
 uc=UCell::ScoreSignatures_UCell(x[,added,drop=FALSE],features=sets,maxRank=1500,ties.method="average",missing_genes="skip",ncores=1,name="")
 tab=data.frame(id=added,H05_transfer=w[added,"H05"],uc[added,,drop=FALSE],check.names=FALSE)
 for(k in c("MP10_covered","Chol_biosynthesis_8")){a=as.matrix(x[sets[[k]],added,drop=FALSE]);tab[[paste0(k,"_detected")]]=colSums(a>0);tab[[paste0(k,"_max_gene_fraction")]]=apply(a,2,max)/pmax(1,colSums(a))}
 fwrite(tab,file.path("xenium",s,"additional_scores.tsv.gz"),sep="\t")
}
