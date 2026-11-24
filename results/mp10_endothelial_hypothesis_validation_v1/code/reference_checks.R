source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
logs=fread(file.path(O,"logs/resource_attempts.tsv"))
files=unique(logs[,.(resource,file)])
files[,path:=file.path(C,resource,file)]
books=files[file.exists(path)&grepl("[.]xlsx$",file)]
pdf=files[file.exists(path)&grepl("[.]pdf$",file)]
sheet_info=rbindlist(lapply(seq_len(nrow(books)),function(i){
  sheets=readxl::excel_sheets(books$path[i])
  data.table(accession=books$resource[i],file=books$file[i],type="xlsx",pages_or_sheets=length(sheets),readout=paste(sheets,collapse=" | "))
}))
pdf_info=rbindlist(lapply(seq_len(nrow(pdf)),function(i){
  pg=pdftools::pdf_text(pdf$path[i])
  data.table(accession=pdf$resource[i],file=pdf$file[i],type="pdf",pages_or_sheets=length(pg),readout=paste(sum(nchar(pg)),"text characters parsed"))
}))
put(rbind(sheet_info,pdf_info),"supplement_readout.tsv")
jlr=file.path(C,"PMC10407908","mmc1.xlsx")
peaks=as.data.table(readxl::read_excel(jlr,sheet=3,.name_repair="minimal"))
put(peaks[`Gene Name`=="KLF6"],"published_HUVEC_SREBP2_peaks_near_KLF6.tsv")
core=gettab("fixed_core_members.tsv")
signatures=core[,.(mapped_members=uniqueN(gene)),by=program]
signatures[,source:=fcase(program=="ATF4_growth","Torrence Fig 1 source data 1; sheets 2 and 3 intersection",program=="ATF4_ISR","Torrence Fig 1 source data 1; sheet 5",program=="Fan_Red_ATF4","Fan Data S1; sheet 2",default="MSigDB Hallmark 2026.1.Hs")]
signatures[,members:="fixed_core_members.tsv"]
ext=gettab("external_frozen_tails.tsv.gz")[,.(mapped_members=uniqueN(gene)),by=.(program=contrast,arm)]
ext[,source:=sub("__.*","",program)];ext[,members:="external_frozen_tails.tsv.gz"]
secondary=gettab("secondary_prediction_members.tsv")[,.(mapped_members=uniqueN(gene)),by=program]
secondary[,source:=fcase(program=="Anabolic_functional","ATF4_growth intersect GO:0006418/0003333/0008652",program=="MHC_I","GO:0002474",program=="JLR_selected","JLR 2023 text; selected cytokine readouts",default="PNAS 2023; selected CXCL locus")]
secondary[,members:="secondary_prediction_members.tsv"]
put(rbindlist(list(signatures,ext,secondary),fill=TRUE),"signature_registry.tsv")
