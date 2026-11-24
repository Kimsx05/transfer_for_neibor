suppressPackageStartupMessages(library(data.table));setDTthreads(1);o=commandArgs(TRUE)[1]
z=fread(file.path(o,'scripts/modules.tsv'));g=unique(c(unlist(strsplit(z$genes,' ')),'SREBF2','SCAP','CHMP2A'))
all=as.data.table(babelgene:::orthologs_df)[taxon_id==10090 & support_n>=3]
all[,n_mouse:=uniqueN(entrez),by=human_entrez];all[,n_human:=uniqueN(human_entrez),by=entrez]
out=merge(data.table(human_symbol=g),all,by='human_symbol',all.x=TRUE)
out[,eligible:=!is.na(entrez) & n_mouse==1 & n_human==1]
fwrite(out,file.path(o,'tables/mouse_orthology.tsv'),sep='\t',na='NA')

