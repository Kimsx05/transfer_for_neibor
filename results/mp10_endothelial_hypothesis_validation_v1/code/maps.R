source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
h=fread(file.path(R,"hgnc_complete_set_20261118.txt"),select=c("symbol","prev_symbol","alias_symbol","status"),quote="")[status=="Approved"]
a=rbindlist(lapply(c("prev_symbol","alias_symbol"),function(k)h[nzchar(get(k)),.(old=unlist(strsplit(get(k),"|",fixed=TRUE))),by=symbol]))
a=unique(a)[,if(uniqueN(symbol)==1).(gene=symbol[1]),by=old]
human=function(x)ifelse(x%in%h$symbol,x,a$gene[match(x,a$old)])
m=as.data.table(get("orthologs_df",asNamespace("babelgene")))[taxon_id==10090 & support_n>=3]
m[,mouse_n:=uniqueN(human_entrez),by=entrez];m[,human_n:=uniqueN(entrez),by=human_entrez]
m[,gene:=human(human_symbol)];m=m[mouse_n==1 & human_n==1 & !is.na(gene)]
m=m[!duplicated(gene)&!duplicated(gene,fromLast=TRUE)]
put(m,"orthology_one_to_one.tsv.gz")
f=Sys.glob(file.path(C,"eLife*","*xlsx")); getset=function(i){x=as.data.frame(readxl::read_excel(f,i,.name_repair="minimal"));data.table(id=x[[2]],mouse=x[[3]],entrez=as.character(x[[4]]))}
ins=getset(3)[id%in%getset(2)$id];isr=getset(5); stopifnot(nrow(ins)==77,nrow(isr)==774,length(intersect(ins$id,isr$id))==61)
fan=as.character(readxl::read_excel(file.path(C,"PMC8432967/mmc2.xlsx"),2,skip=3,col_names=FALSE)[[1]])
s=rbindlist(list(ins[,program:="ATF4_growth"],isr[,program:="ATF4_ISR"],data.table(mouse=fan,program="Fan_Red_ATF4")),fill=TRUE)
s[,gene:=m$gene[match(entrez,m$entrez)]];s[program=="Fan_Red_ATF4",gene:=m$gene[match(mouse,m$symbol)]]
s[,shared_mouse:=mouse%in%intersect(ins$mouse,isr$mouse)]
put(s,"literature_original_sets.tsv")
hm=fread(file.path(R,"hallmark.tsv.gz"));keep=c("MYC_TARGETS_V1","MYC_TARGETS_V2","TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE","E2F_TARGETS","G2M_CHECKPOINT","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE")
sets=rbind(s[!is.na(gene),.(program,gene)],hm[sub("^HALLMARK_","",gs_name)%in%keep,.(program=sub("^HALLMARK_","",gs_name),gene=gene_symbol)])
put(unique(sets),"fixed_core_members.tsv")

